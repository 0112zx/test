
library(tidyverse)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(circlize)
library(ComplexHeatmap)
library(ggh4x)
library(patchwork)

# primary glass WES sample data

glass_cnsig_data <- read_tsv('/boot3/bio_liaojl/GLASS_copy_number_signature_work/Result/CN_signature_identification/COSMIC/sigprofiler_output_wes_ft163/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Activities/CNV48_De-Novo_Activities_refit.txt') %>% 
  filter(str_detect(Samples, '-TP-')) %>% 
  mutate(patient = str_sub(Samples, 1, 12)) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  select(-patient) %>% 
  pivot_longer(-Samples, names_to = 'CNSig', values_to = 'weight') %>% 
  group_by(Samples) %>% 
  mutate(weight = weight/sum(weight)) %>% 
  ungroup() %>% 
  mutate(CNSig = str_replace(CNSig, 'CNV48', 'CNSig'))
# first time: 163
# primary: 160 samples, GLSS, 157; TCGA, 3



glass_cli_sam <- read_tsv('/boot3/bio_liaojl/GLASS_copy_number_signature_work/Data/GLASS/legacy/20190328_data_freeze/Tables/clinical_surgeries.tsv')
cli_case <- read_tsv('/boot3/bio_liaojl/GLASS_copy_number_signature_work/Data/GLASS/legacy/20190328_data_freeze/Tables/clinical_cases.tsv')
variants_titan_pp <- read_tsv('/boot3/bio_liaojl/GLASS_copy_number_signature_work/Data/GLASS/legacy/20190328_data_freeze/Tables/variants_titan_pp.tsv') %>% 
  mutate(sample_barcode = str_sub(pair_barcode, 1, 15)) %>% 
  select(sample_barcode, purity, ploidy) %>% 
  distinct(sample_barcode, .keep_all = TRUE)





# Similarity to COSMIC 21 -------------------------------------------------


load('/boot3/bio_liaojl/GLASS_copy_number_signature_work/Result/CN_signature_identification/COSMIC/all_sim_sorted.Rdata')

titan_wesft_cosmic21 <- all_sim_sorted %>% 
  filter(Method_A %in% 'Titan', Method_B %in% 'Ascat', Seq_A %in% 'WES', Cohort_A %in% 'Ft') %>% 
  select(SigName_A, SigName_B, cosine_similarity) %>% 
  rename(Signature = SigName_A, CN_signature = SigName_B) %>% 
  mutate(Signature = str_c('CNSig', Signature))

CN_sig_order <- titan_wesft_cosmic21 %>% distinct(CN_signature) %>% pull(CN_signature) %>% str_sort(numeric = TRUE)

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/cn_sig_similarity.pdf', width = 12, height = 3)

titan_wesft_cosmic21 %>% 
  mutate(cosine_similarity = round(cosine_similarity, 2), CN_signature = factor(CN_signature, levels = CN_sig_order)) %>% 
  ggplot(aes(CN_signature, Signature, fill = cosine_similarity)) +
  geom_tile() +
  geom_text(aes(label = cosine_similarity)) +
  labs(x = NULL, y = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()




# prevalence and weight across GLASS --------------------------------------


bar_p_data <- glass_cnsig_data %>% 
  group_by(CNSig) %>% 
  summarise(frac = mean(weight > 0), 
            exp1 = sum(weight > 0), 
            exp2 = n())

bar_p <- bar_p_data %>% 
  mutate(expression = map2(exp1, exp2, ~bquote(over(.(.x), .(.y))))) %>% 
  ggplot(aes(CNSig, frac, fill = CNSig)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = expression), parse = TRUE, vjust = -0.15, size = 3) +
  labs(x = NULL, y = 'Prevalence') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  ggsci::scale_fill_jco() +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())

box_w <- glass_cnsig_data %>% 
  filter(weight > 0) %>% 
  ggplot(aes(CNSig, weight, col = CNSig)) +
  geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  labs(x = NULL, y = 'Weight') +
  ggsci::scale_color_jco() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/CNSig_prevalence_comb.pdf', width = 5, height = 6)

box_w / bar_p

dev.off()


cnsig_stat <- glass_cnsig_data %>% 
  filter(weight > 0) %>% 
  group_by(CNSig) %>% 
  summarise(min_w = min(weight), 
            med_w = median(weight), 
            max_w = max(weight))


# consensus clustering ----------------------------------------------------

expo_frac_mat <- glass_cnsig_data %>% 
  pivot_wider(names_from = CNSig, values_from = weight) %>% 
  column_to_rownames('Samples') %>% 
  t()

cnsig_cluster_res <- ConsensusClusterPlus(expo_frac_mat, maxK = 12, seed = 123, distance = 'pearson', title = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/glass_cnsig_clustering_wes_tp', plot = 'pdf')

# save(cnsig_cluster_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/glass_cnsig_clustering_wes_tp/cnsig_cluster_res.Rdata')

# k = 4 is optimal number of cluster

sam_cluster <- cnsig_cluster_res[[4]]$consensusClass

# identical(names(sam_cluster), colnames(expo_frac_mat))

cnsig_group_cli_data <- tibble(sample = names(sam_cluster), cluster = as.character(sam_cluster), sample_barcode = str_sub(sample, 1, 15)) %>% 
  inner_join(glass_cli_sam, by = 'sample_barcode') %>% 
  left_join(cli_case, by = 'case_barcode') %>% 
  left_join(variants_titan_pp, by = 'sample_barcode')

# save(cnsig_group_cli_data, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/cnsig_group_cli_data.Rdata')



# Survival analysis -------------------------------------------------------


cnsig_group_sur_dat <- cnsig_group_cli_data %>% 
  rename(OS = case_vital_status, OS.time = case_overall_survival_mo) %>% 
  filter(!is.na(OS), !is.na(OS.time)) %>% 
  mutate(OS = ifelse(OS %in% 'dead', 1, 0), 
         cluster = factor(cluster, levels = c('1', '2', '3', '4')))


# summary(coxph(Surv(OS.time, OS) ~ cluster + case_age_diagnosis_years + case_sex + grade + idh_codel_subtype, data = cnsig_group_sur_dat))
# not significant

# median survival time
# surv_median(surv_fit(Surv(OS.time, OS) ~ cluster, data = cnsig_group_sur_dat))

# pairwise_survdiff(Surv(OS.time, OS) ~ cluster, data = cnsig_group_sur_dat, p.adjust.method = "none")

def_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ cluster, data = cnsig_group_sur_dat), 
                       pval = TRUE, 
                       legend.labs = c('1', '2', '3', '4'), 
#                       palette = c('#E64B35FF', '#4DBBD5FF', '#00A087FF', '#3C5488FF'), 
                       ggtheme = theme_bw())

final_plot <- def_plot$plot +
  labs(x = 'Survival time (Months)', y = 'OS probability', col = NULL) +
  theme(panel.grid = element_blank(), 
        legend.position = c(1, 1), 
        legend.justification = c(1, 1), 
        legend.background = element_rect(fill = NA), 
        legend.key = element_rect(fill = NA))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/cnsig_cluster_kmplot_wes_tp.pdf', width = 6, height = 6)

final_plot

dev.off()



# heatmap plotting --------------------------------------------------------

clust_col <- cnsig_cluster_res[[4]]$consensusTree

cli_mol_feat <- cnsig_group_cli_data %>% 
  select(sample_barcode, cluster, histology, grade, idh_codel_subtype, case_sex, case_age_diagnosis_years, purity, ploidy)

set.seed(122)

anno_col <- HeatmapAnnotation(sex = cli_mol_feat$case_sex, 
                              age = cli_mol_feat$case_age_diagnosis_years, 
                              purity = cli_mol_feat$purity, 
                              ploidy = cli_mol_feat$ploidy, 
                              Cluster = cli_mol_feat$cluster, 
                              annotation_legend_param = list(sex = list(nrow = 1), 
                                                             age = list(direction = "horizontal"), 
                                                             purity = list(direction = "horizontal"), 
                                                             ploidy = list(direction = "horizontal"), 
                                                             Cluster = list(nrow = 1)))

cnsig_heatmap <- Heatmap(expo_frac_mat, 
                         name = 'exposure', 
                         col = colorRamp2(c(0, 0.5, 1), c('#006CB1', 'white', '#C9372E')), 
                         cluster_rows = FALSE, 
                         cluster_columns = clust_col, 
                         show_column_names = FALSE, 
                         column_dend_height = unit(2, 'cm'), 
                         height = unit(6.5, 'cm'), 
                         top_annotation = anno_col, 
                         heatmap_legend_param = list(legend_direction = "horizontal"))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/cnsig_cluster_heatmap_wes_tp.pdf', width = 16, height = 10)

draw(cnsig_heatmap, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")

dev.off()




# cluster gt summary ------------------------------------------------------

require(gtsummary)
require(gt)

glass_cnsig_group_cli_data_gt_res <- cnsig_group_cli_data %>% 
  select(cluster, grade, mgmt_methylation, treatment_tmz, treatment_radiotherapy, 
         idh_codel_subtype, treatment_alkylating_agent, case_sex, case_age_diagnosis_years, purity:ploidy) %>% 
  tbl_summary(by = cluster, missing = "no", statistic = all_continuous() ~ "{median} ({p0}, {p100})") %>%
  add_n() %>%
  bold_labels() %>% 
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_q(method = 'BH')

# get absolute p values

glass_subtype_cli_mole_test_p <- tibble(variable = glass_cnsig_group_cli_data_gt_res$table_body$variable, 
                                  p_val = glass_cnsig_group_cli_data_gt_res$table_body$p.value, 
                                  q_val = glass_cnsig_group_cli_data_gt_res$table_body$q.value) %>% 
  distinct(variable, .keep_all = TRUE)

write_tsv(glass_subtype_cli_mole_test_p, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/glass_subtype_cli_mole_test_p.tsv')

# get summary result

glass_cnsig_subtype_cli_mole_summary <- glass_cnsig_group_cli_data_gt_res %>% as_flex_table()

save(glass_cnsig_subtype_cli_mole_summary, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/glass_cnsig_subtype_cli_mole_summary.Rdata')



### categorical variables ( barplot)

bar_plot <- function(var_name){
  
  # var_name <- 'histologic_type'
  
  sin_bar_plot <- cnsig_group_cli_data %>% 
    select(cluster, sin_var = all_of(var_name)) %>% 
    drop_na() %>% 
    ggplot(aes(cluster, fill = sin_var)) +
    geom_bar(position = 'fill') +
    labs(x = NULL, y = var_name, fill = NULL) +
    # ggsci::scale_fill_npg() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    # guides(fill = guide_legend(nrow = 2)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic() +
    theme(axis.ticks.x = element_blank())
  
}


### combine plots

comb_plots <- list(bar_plot('grade'), bar_plot('idh_codel_subtype'))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/GLASS_validation/cluster_grade_idh_plots.pdf', width = 9, height = 4)

wrap_plots(comb_plots)

dev.off()







