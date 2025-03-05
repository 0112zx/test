
require(tidyverse)
require(gtsummary)
require(gt)
require(patchwork)
require(ggalluvial)

# in windows --------------------------------------------------------------

# load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Data/Processed/glioma_cnsig_subtype_label.Rdata')
# tcga_immune_dat <- read_tsv('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Data/Original/TCGA_immune_related_score.tsv')


load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')
tcga_immune_dat <- read_tsv('/boot3/bio_liaojl/Common_data/TCGA_immune_related_score.tsv')

glioma_cnsig_subtype_immune_data <- glioma_cnsig_subtype_label %>% 
  select(patient, cluster) %>% 
  inner_join(tcga_immune_dat, by = c('patient' = 'TCGA Participant Barcode')) %>% 
  select(-`TCGA Study`, -`TCGA Subtype`, -`Intratumor Heterogeneity`, -(`Silent Mutation Rate`:`Homologous Recombination Defects`))

cnsig_subtype_immune_gt_res <- glioma_cnsig_subtype_immune_data %>% 
  select(-patient, -`TIL Regional Fraction`) %>% # TIL Regional Fraction, error: all observations are in the same group
  tbl_summary(by = cluster, missing = "no", statistic = all_continuous() ~ "{median} ({p0}, {p100})") %>%
  add_n() %>%
  bold_labels() %>% 
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_q(method = 'BH')


# get absolute p values

subtype_immune_test_p <- tibble(variable = cnsig_subtype_immune_gt_res$table_body$variable,
                                p_val = cnsig_subtype_immune_gt_res$table_body$p.value, 
                                q_val = p.adjust(p_val, method = 'BH')) %>% 
  distinct(variable, .keep_all = TRUE)

write_tsv(subtype_immune_test_p, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_test_p.tsv')

# get summary result

subtype_immune_summary <- cnsig_subtype_immune_gt_res %>% as_flex_table()

save(subtype_immune_summary, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_summary.Rdata')

# save to microsoft word

subtype_immune_summary %>% flextable::save_as_docx(path = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_summary.docx')

# load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_immune_difference/subtype_immune_summary.Rdata')


# representative immune feature plotting ----------------------------------

# 'Leukocyte Fraction', 'Stromal Fraction', 'IFN-gamma Response', 'TGF-beta Response', 'SNV Neoantigens', 'TCR Richness', 'CTA Score'

# immune_feature_plot <- function(immune_feature){
#   
#   immune_feature <- 'Leukocyte Fraction'
#   
#   sin_feature_plot <- glioma_cnsig_subtype_immune_data %>% 
#     select(cluster, sin_var = {{immune_feature}}) %>% 
#     ggplot(aes(cluster, sin_var, fill = cluster)) +
#     geom_boxplot(outlier.color = '#EAEAEA') +
#     ggpubr::stat_compare_means(label = "p.signif") +
#     labs(x = NULL, y = immune_feature) +
#     scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
#     theme_bw() +
#     theme(axis.text.x = element_blank(), 
#           axis.ticks.x = element_blank(), 
#           legend.position = 'none')
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/temp.pdf')
  # 
  # sin_feature_plot
  # 
  # dev.off()
  
  
# }


repre_immune_feature <- c('Leukocyte Fraction', 'Stromal Fraction', 'IFN-gamma Response', 'TGF-beta Response', 'TCR Richness', 'CTA Score')

immune_feature_plot_data <- glioma_cnsig_subtype_immune_data %>% 
  select(patient, cluster, {{repre_immune_feature}}) %>% 
  pivot_longer(-(patient:cluster), names_to = 'immune_feature', values_to = 'score') %>% 
  group_by(immune_feature) %>% 
  mutate(scaled_score = scale(score)[, 1]) %>% 
  ungroup()


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_feature_boxplot.pdf', width = 8, height = 4)

immune_feature_plot_data %>% 
  ggplot(aes(immune_feature, scaled_score, fill = cluster)) +
  geom_boxplot(outlier.color = '#EAEAEA') +
  ggpubr::stat_compare_means(label = "p.signif", label.y = 5) +
  labs(x = NULL, y = 'Scaled score', fill = NULL) +
  coord_cartesian(ylim = c(-3, 5)) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()



# alluvial plot (TCGA immune subtype) -----------------------------------------------------------


alluvial_plot_data <- glioma_cnsig_subtype_immune_data %>% 
  select(patient, cluster, `Immune Subtype`) %>% 
  drop_na() %>% 
  pivot_longer(-patient, names_to = 'classification', values_to = 'type') %>% 
  mutate(freq = 1)


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_alluvial_plot.pdf')

alluvial_plot_data %>% 
  ggplot(aes(classification, freq, stratum = type, alluvium = patient, fill = type, label = type)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH', 'C4', 'C5', 'C3', 'C1', 'C6'), values = c("#D7191C", "#FDAE61", "#ABDDA4", "#2B83BA", '#BC80BD', '#B15928', '#FCCDE5', '#FDBF6F', '#CCEBC5', '#FFED6F')) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(), 
        panel.background = element_blank(), 
        panel.grid = element_blank())


dev.off()







# immune filtration (CIBERSORT) -------------------------------------------------------

tcga_cibersort_data <- read_tsv('/boot3/bio_liaojl/Common_data/TCGA.Kallisto.fullIDs.cibersort.relative.tsv')

subtype_cibersort_test_data <- tcga_cibersort_data %>% 
  select(-CancerType, -(P.value:RMSE)) %>% 
  mutate(SampleID = str_replace_all(str_sub(SampleID, 1, 15), '\\.', '-')) %>% 
  pivot_longer(-SampleID, names_to = 'cell_type', values_to = 'fraction') %>% 
  mutate(cell_type = str_replace_all(cell_type, '\\.', ' ')) %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = c('SampleID' = 'sample'))
# 626 samples

cibersort_test_res <- subtype_cibersort_test_data %>% 
  group_by(cell_type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(fraction ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data) %>% 
  arrange(q_val)
# 12 significant

write_tsv(cibersort_test_res, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/cibersort_test_res.tsv')


# plotting

cibersort_signif_label <- cibersort_test_res %>% 
  mutate(sig_lab = case_when(
    q_val < 0.001 ~ '***', 
    q_val >= 0.001 & q_val < 0.01 ~ '**', 
    q_val >= 0.01 & q_val < 0.05 ~ '*', 
    q_val >= 0.05 ~ ''
  ))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_cibersort_boxplot_.pdf', width = 15, height = 3.5)

subtype_cibersort_test_data %>% 
  group_by(cell_type) %>% 
  mutate(fraction = scale(fraction)[, 1]) %>% # z-score
  ungroup() %>% 
  ggplot() +
  geom_boxplot(aes(cell_type, fraction, fill = cluster), outlier.color = '#EAEAEA') +
  geom_text(aes(cell_type, 1, label = sig_lab), data = cibersort_signif_label) +
  labs(x = NULL, y = 'Proportion', fill = NULL) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()




# immune filtration (TIMER) -------------------------------------------------------

subtype_timer_data <- read_csv('/boot3/bio_liaojl/Common_data/infiltration_estimation_for_tcga.csv') %>% 
  rename(sample = cell_type) %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample') %>% 
  select(sample, cluster, ends_with('_TIMER')) %>% 
  pivot_longer(-(sample:cluster), names_to = 'cell_type', values_to = 'score') %>% 
  mutate(cell_type = str_remove_all(cell_type, '_TIMER'))

subtype_timer_test_res <- subtype_timer_data %>% 
  group_by(cell_type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(score ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)
# all(6) significant

subtype_timer_stat <- subtype_timer_data %>% 
  group_by(cell_type, cluster) %>% 
  summarise(med_score = signif(median(score), digits = 3)) %>% 
  ungroup()
  
write_tsv(subtype_timer_stat, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_timer_stat.tsv')


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_timer_boxplot.pdf', width = 6, height = 3.5)

subtype_timer_data %>% 
  ggplot() +
  geom_boxplot(aes(cell_type, score, fill = cluster), outlier.color = '#EAEAEA') +
  labs(x = NULL, y = 'Timer Score', fill = NULL) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


# representative plots
# T cell CD4+, T cell CD8+, B cell

my_comparisons = list(c('Diploid', 'Diploid CIN'), c('Diploid', 'Tetraploid'), c('Diploid', 'Tetraploid CTH'), c('Diploid CIN', 'Tetraploid'), c('Diploid CIN', 'Tetraploid CTH'), c('Tetraploid', 'Tetraploid CTH'))

cell_infiltration_plot <- function(cell_type){
  
  # cell_type <- 'T cell CD8+'
  
  sin_ci_plot <- subtype_timer_data %>% 
    filter(cell_type %in% {{cell_type}}) %>% 
    ggplot(aes(cluster, score, col = cluster)) +
    geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
    geom_boxplot(outlier.color = NA, fill = NA) +
    # ggpubr::stat_compare_means(comparisons = my_comparisons, label = "p.signif", hide.ns = TRUE) +
    ggpubr::geom_pwc(label = "p.signif", hide.ns = TRUE) + 
    ggpubr::stat_compare_means() +
    labs(x = NULL, y = 'Infiltration level', title = cell_type) +
    scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    # scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank(), 
          plot.title = element_text(hjust = 0.5))
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/temp.pdf')
  # 
  # sin_ci_plot
  # 
  # dev.off()
  
}


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_timer_b_tcell_boxplot.pdf', width = 14, height = 6)

cell_infiltration_plot('T cell CD4+') + cell_infiltration_plot('T cell CD8+') + cell_infiltration_plot('B cell')

dev.off()


# 'Neutrophil', 'Macrophage', 'Myeloid dendritic cell'

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_timer_other_cell_boxplot.pdf', width = 14, height = 6)

cell_infiltration_plot('Neutrophil') + cell_infiltration_plot('Macrophage') + cell_infiltration_plot('Myeloid dendritic cell')

dev.off()












# xCell differences -----------------------------------------------------------

# load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_cli_data.Rdata')

require(ComplexHeatmap)
require(circlize)

xCell_TCGA_RSEM <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/xCell_TCGA_RSEM.txt')

xCell_glioma_RSEM <- xCell_TCGA_RSEM %>% 
  set_names(~str_replace_all(., '\\.', '-')) %>% 
  select(cell_type, any_of(glioma_cnsig_subtype_label$sample))
# 617 samples

# statistical test among cnsig subtypes

xCell_glioma_RSEM_test_res <- xCell_glioma_RSEM %>% 
  pivot_longer(-cell_type, names_to = 'sample', values_to = 'score') %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample') %>% 
  group_by(cell_type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(score ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)

# heatmap

xcell_mat <- xCell_glioma_RSEM %>% column_to_rownames('cell_type') %>% as.matrix()
glioma_cnsig_subtype_label_sort <- glioma_cnsig_subtype_label %>% filter(sample %in% colnames(xcell_mat)) %>% arrange(cluster)

xcell_mat_scaled <- xCell_glioma_RSEM %>% 
  pivot_longer(-cell_type, names_to = 'sample', values_to = 'score') %>% 
  group_by(cell_type) %>% 
  mutate(score = scale(score)[, 1]) %>% 
  ungroup() %>% 
  pivot_wider(names_from = sample, values_from = score) %>% 
  column_to_rownames('cell_type') %>% 
  as.matrix()


anno_col <- HeatmapAnnotation(Cluster = glioma_cnsig_subtype_label_sort$cluster)

set.seed(1234)

cluster_heatmap <- Heatmap(xcell_mat_scaled[, glioma_cnsig_subtype_label_sort$sample], 
                           name = 'score', 
                           col = colorRamp2(c(-1.5, 0, 1.5), c('#006CB1', 'white', '#C9372E')),
                           cluster_rows = TRUE, 
                           show_row_names = TRUE, 
                           cluster_columns = FALSE, 
                           show_column_names = FALSE, 
                           top_annotation = anno_col)



pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/xcell_cluster_heatmap_scaled.pdf', width = 13, height = 9)

cluster_heatmap

dev.off()





# neoantigens -------------------------------------------------------------

glioma_neoantigen_dat <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/TCGA_immune_related_score.tsv') %>% 
  filter(`TCGA Study` %in% c('LGG', 'GBM')) %>% 
  select(patient = `TCGA Participant Barcode`, snv_neoantigen = `SNV Neoantigens`, indel_neoantigen = `Indel Neoantigens`) %>% 
  pivot_longer(-patient, names_to = 'type', values_to = 'value') %>% 
  drop_na(value) %>% 
  inner_join(select(glioma_cnsig_subtype_label, patient, cluster), by = 'patient')

glioma_neoantigen_test <- glioma_neoantigen_dat %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(value ~ cluster, data = .)$p.value)) %>% 
  select(-data)
# type                p_val
# <chr>               <dbl>
# snv_neoantigen   4.74e-26
# indel_neoantigen 1.46e- 1

# boxplot

snv_neo_plot <- glioma_neoantigen_dat %>% 
  filter(type %in% 'snv_neoantigen', value < 5000) %>% # remove TCGA-DU-6392 of which snv neoantigen was 5122
  mutate(log_value = log10(value + 1)) %>% 
  ggplot(aes(cluster, log_value, col = cluster)) +
  geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  ggpubr::geom_pwc(label = "p.signif", hide.ns = TRUE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'log10(SNV Neoantigens)') +
  scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))

indel_neo_plot <- glioma_neoantigen_dat %>% 
  filter(type %in% 'indel_neoantigen') %>%
  mutate(log_value = log10(value + 1)) %>% 
  ggplot(aes(cluster, log_value, col = cluster)) +
  geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  ggpubr::geom_pwc(label = "p.signif", hide.ns = TRUE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'log10(Indel Neoantigens)') +
  scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(), 
        plot.title = element_text(hjust = 0.5))



pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_neoantigen_boxplot.pdf', width = 10, height = 6)

snv_neo_plot + indel_neo_plot

dev.off()



# selected immune gene signature ------------------------------------------

glioma_exp_data <- vroom::vroom('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_exp_data.tsv') %>%
  mutate_if(is.numeric, ~log2(.+1))

glioma_exp_mat <- glioma_exp_data %>% column_to_rownames('Hugo_Symbol') %>% as.matrix()


# TCGA immune-related score

glioma_immune_related_score <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/TCGA_immune_related_score.tsv') %>% 
  filter(`TCGA Study` %in% c('LGG', 'GBM')) %>% 
  select(patient = `TCGA Participant Barcode`, Proliferation:`TGF-beta Response`) %>% 
  rename(`Lymphocyte Infiltration` = `Lymphocyte Infiltration Signature Score`)

# immune function

require(GSVA)

# function: ssGSEA score computation

ssGSEA_compute <- function(expr_mat, gene_list){
  
  # expr_mat <- crc_exp_mat
  # gene_list <- c2_cp_gene
  
  ssgsea_res <- gsva(expr_mat, gene_list, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
  
  ssgsea_res_dat <- ssgsea_res %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('Sample') %>% 
    as_tibble() %>% 
    pivot_longer(-Sample, names_to = 'gene_set', values_to = 'ssGSEA_score')
  
}

immune_gene_signature <- qusage::read.gmt('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/IMM.gmt')

immune_function_signature <- immune_gene_signature[c('APC_co_inhibition', 'APC_co_stimulation', 'Cytolytic_activity', 'HLA', 'Inflammation-promoting', 'MHC_class_I', 'Parainflammation')]

immune_function_ssgsea_score <- ssGSEA_compute(glioma_exp_mat, immune_function_signature) %>% rename(sample = Sample, gene_signature = gene_set, score = ssGSEA_score)


# combine all immune gene signature

glioma_imm_genesig_comb <- glioma_immune_related_score %>% 
  pivot_longer(-patient, names_to = 'gene_signature', values_to = 'score') %>% 
  drop_na(score) %>% 
  inner_join(glioma_cnsig_subtype_label, by = 'patient') %>% 
  select(-patient, -cluster) %>% 
  bind_rows(immune_function_ssgsea_score) %>% 
  inner_join(glioma_cnsig_subtype_label, by = 'sample')
  
imm_genesig_test_res <- glioma_imm_genesig_comb %>% 
  group_by(gene_signature) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(score ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)

save(imm_genesig_test_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/imm_genesig_test_res.Rdata')


# dotplot

img_signature_plot_dat <- glioma_imm_genesig_comb %>% 
  group_by(gene_signature) %>% 
  mutate(score = scale(score)[, 1]) %>% # z-score
  ungroup() %>% 
  group_by(gene_signature, cluster) %>% 
  summarise(score = median(score)) %>% 
  ungroup() %>% 
  inner_join(imm_genesig_test_res, by = 'gene_signature') %>% 
  mutate(log_q = -log10(q_val))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_gene_signature_dotplot.pdf', width = 5, height = 4)

img_signature_plot_dat %>% 
  mutate(cluster = factor(cluster, levels = c('Diploid', 'Tetraploid', 'Diploid CIN', 'Tetraploid CTH')), 
         gene_signature = str_replace_all(gene_signature, '_', ' ')) %>% 
  ggplot(aes(cluster, gene_signature, fill = score), color = 'grey60') +
  geom_point(aes(size = log_q), shape = 21) +
  labs(x = NULL, y = NULL, fill = 'Median score (z-scored)', size = '-log10(q)') +
  scale_fill_gradientn(colors = c('#006CB1', 'white', '#C9372E')) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()




# immune checkpoint expression --------------------------------------------

immune_check_expr <- glioma_exp_data %>% 
  filter(Hugo_Symbol %in% c('BTLA', 'CD244', 'CD274', 'CD276', 'CD40', 'CD40LG', 'CD48', 'CD80', 'CD86', 'CEACAM1', 'CLEC4G', 'CTLA4', 'FGL1', 'HAVCR2', 'HLA-A', 'HLA-B', 'HLA-C', 'HMGB1', 'KIR3DL1', 'LAG3', 'LGALS3', 'LGALS9', 'NECTIN2', 'PDCD1', 'PDCD1LG2', 'PVR', 'SNCA', 'TIGIT', 'TNFRSF14', 'TNFRSF18', 'TNFRSF4', 'TNFRSF9', 'TNFSF18', 'TNFSF4', 'TNFSF9', 'VSIR')) %>% 
  pivot_longer(-Hugo_Symbol, names_to = 'sample', values_to = 'geneExp') %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample')

# 34 genes

immune_check_test_res <- immune_check_expr %>% 
  group_by(Hugo_Symbol) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(geneExp ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data) %>% 
  arrange(q_val)
# 32 significant

write_tsv(immune_check_test_res, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/immune_check_test_res.tsv')

# representative gene plotting
# CD274(PDL1), PDCD1(PD1), PDCD1LG2(PDL2), CD152(CTLA4)

rep_immune_check_data <- immune_check_expr %>% filter(Hugo_Symbol %in% c('CD274', 'PDCD1', 'PDCD1LG2', 'CTLA4'))


immune_check_syn <- c(
  'CD274' = 'CD274(PDL1)', 
  'PDCD1' = 'PDCD1(PD1)', 
  'PDCD1LG2' = 'PDCD1LG2(PDL2)', 
  'CTLA4' = 'CTLA4(CD152)'
)


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune_check_violin_plot_.pdf', width = 12, height = 4)

rep_immune_check_data %>% 
  ggplot(aes(cluster, geneExp, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
  ggpubr::geom_pwc(label = "p.signif", hide.ns = TRUE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'log2(geneExp+1)', fill = NULL) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  facet_grid(.~Hugo_Symbol, labeller = as_labeller(immune_check_syn)) +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid = element_blank(), 
        strip.text.x = element_text(margin = margin(0.2, 0, 0.2, 0, "cm"), size = 10))

dev.off()


# immune filtration (28 immune signature) -------------------------------------------------------


immune28_signature <- qusage::read.gmt('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/Maker28_PMID_28052254.gmt')

immune28_ssgsea_score <- ssGSEA_compute(glioma_exp_mat, immune28_signature) %>% rename(sample = Sample)

immune28_ssgsea_score_test <- immune28_ssgsea_score %>% inner_join(glioma_cnsig_subtype_label, by = 'sample')

immune28_ssgsea_score_test_res <- immune28_ssgsea_score_test %>% 
  group_by(gene_set) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(ssGSEA_score ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)

immune28_signif_label <- immune28_ssgsea_score_test_res %>% 
  mutate(sig_lab = case_when(
    q_val < 0.0001 ~ '****', 
    q_val >= 0.0001 & q_val < 0.001 ~ '***', 
    q_val >= 0.001 & q_val < 0.01 ~ '**', 
    q_val >= 0.01 & q_val < 0.05 ~ '*', 
    q_val >= 0.05 ~ ''
  ))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_immune_difference/subtype_immune28_boxplot.pdf', width = 15, height = 5)

immune28_ssgsea_score_test %>% 
  ggplot() +
  geom_boxplot(aes(gene_set, ssGSEA_score, fill = cluster), outlier.color = '#EAEAEA') +
  geom_text(aes(gene_set, 0.7, label = sig_lab), data = immune28_signif_label) +
  labs(x = NULL, y = 'Infiltration score', fill = NULL) +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'top')

dev.off()
















