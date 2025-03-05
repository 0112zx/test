
library(tidyverse)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(circlize)
library(ComplexHeatmap)
library(ggh4x)
library(patchwork)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_weight.Rdata')
glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')

# glioma_cnsig_weight_addsubtype <- glioma_cnsig_weight %>% inner_join(select(glioma_cli_data, bcr_patient_barcode, IDH_codel_subtype), by = c('sample' = 'bcr_patient_barcode'))

# consensus clustering ----------------------------------------------------

expo_frac_mat <- glioma_cnsig_weight %>% 
  pivot_wider(names_from = CNSig, values_from = weight) %>% 
  column_to_rownames('sample') %>% 
  t()

cnsig_cluster_res_euclidean <- ConsensusClusterPlus(expo_frac_mat, maxK = 12, seed = 123, distance = 'euclidean', title = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_clustering', plot = 'pdf')

# save(cnsig_cluster_res_euclidean, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_clustering/cnsig_cluster_res_euclidean.Rdata')

# k = 6 is optimal number of cluster

sam_cluster_euclidean <- cnsig_cluster_res_euclidean[[6]]$consensusClass

# identical(names(sam_cluster_euclidean), colnames(expo_frac_mat))

cnsig_group_cli_data <- tibble(patient = names(sam_cluster_euclidean), cluster = as.character(sam_cluster_euclidean)) %>% 
  inner_join(glioma_cli_data, by = c('patient' = 'bcr_patient_barcode'))

save(cnsig_group_cli_data, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/cnsig_group_cli_data.Rdata')


glioma_cnsig_subtype_cli_data <- cnsig_group_cli_data %>% 
  mutate(kps_score = case_when(
    kps_score < 80 ~ '<80', 
    # kps_score >= 50 & kps_score <= 70 ~ '50-70', 
    kps_score >= 80 ~ '80-100', 
    TRUE ~ NA_character_
  ), race = case_when(
    race %in% 'WHITE' ~ 'white', 
    race %in% c('AMERICAN INDIAN OR ALASKA NATIVE', 'ASIAN', 'BLACK OR AFRICAN AMERICAN') ~ 'others', 
    TRUE ~ NA_character_
  ), WGD = ifelse(WGD %in% '0', '0', '1/1+'), 
  transcriptome_subtype = case_when(
    transcriptome_subtype %in% 'PN' ~ 'proneural', 
    transcriptome_subtype %in% 'NE' ~ 'neural', 
    transcriptome_subtype %in% 'ME' ~ 'mesenchymal', 
    transcriptome_subtype %in% 'CL' ~ 'classical'
  )) %>% 
  filter(cluster %in% c('1', '2', '3', '4')) %>% 
  mutate(cluster = case_when(
    cluster %in% '1' ~ 'Diploid CIN', 
    cluster %in% '2' ~ 'Diploid', 
    cluster %in% '3' ~ 'Tetraploid CTH', 
    cluster %in% '4' ~ 'Tetraploid'
  )) %>% # factor order
  mutate(histologic_type = factor(histologic_type, levels = c('oligodendroglioma', 'oligoastrocytoma', 'astrocytoma', 'glioblastoma')), 
         grade = factor(grade, levels = c('G2', 'G3', 'G4')), 
         IDH_codel_subtype = factor(IDH_codel_subtype, levels = c('IDHmut-codel', 'IDHmut-non-codel', 'IDHwt')), 
         transcriptome_subtype = factor(transcriptome_subtype, levels = c('proneural', 'neural', 'mesenchymal', 'classical')))


save(glioma_cnsig_subtype_cli_data, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_cli_data.Rdata')

glioma_cnsig_subtype_label <- glioma_cnsig_subtype_cli_data %>% select(patient, sample, cluster)

save(glioma_cnsig_subtype_label, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')


# heatmap plotting --------------------------------------------------------

# annotation

# clust_col <- cnsig_cluster_res_euclidean[[6]]$consensusTree
# 
# anno_col <- cnsig_group_cli_data %>% 
#   select(patient, cluster, TP53mut, ATRXmut, CNH, aneuploidy_score, genome_doublings, `1p`:`22q`) %>% 
#   mutate_all(as.factor) %>% 
#   column_to_rownames('patient')
# 
# pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_clustering/CNSig_consensus_clustering_heatmap.pdf', width = 12, height = 5)
# 
# print(pheatmap(expo_frac_mat,
#                # color = colorRamp2(c(0, 1), c("white", "#ff6150")),
#               annotation_col = anno_col,
#                # annotation_colors = annotCol,
#                cluster_rows = FALSE,
#                cluster_cols = clust_col,
#                show_rownames = TRUE,
#                show_colnames = FALSE))
# 
# dev.off()


clust_col <- cnsig_cluster_res_euclidean[[6]]$consensusTree

# identical(cnsig_group_cli_data$patient, colnames(expo_frac_mat))

cli_mol_feat <- cnsig_group_cli_data %>% 
  select(patient, cluster, TP53mut, ATRXmut, CNH, aneuploidy_score, purity, ploidy, WGD) %>% 
  mutate(TP53mut = case_when(
    TP53mut == 1 ~ 'Mutant', 
    TP53mut == 0 ~ 'Wild type'
  ), ATRXmut = case_when(
    ATRXmut == 1 ~ 'Mutant', 
    ATRXmut == 0 ~ 'Wild type'
  ))

arm_cnv_mat <- cnsig_group_cli_data %>% 
  select(patient, `1p`:`22q`) %>% 
  mutate_if(is.numeric, ~as.character(.)) %>% 
  column_to_rownames('patient') %>% 
  t()


anno_col <- HeatmapAnnotation(`TP53 status` = cli_mol_feat$TP53mut, 
                              `ATRX status` = cli_mol_feat$ATRXmut, 
                              `Tumor purity` = cli_mol_feat$purity, 
                              `Tumor ploidy` = cli_mol_feat$ploidy, 
                              CNH = cli_mol_feat$CNH, 
                              WGD = cli_mol_feat$WGD, 
                              Aneuploidy = cli_mol_feat$aneuploidy_score, 
                              Cluster = cli_mol_feat$cluster)

set.seed(1234)

cnsig_heatmap <- Heatmap(expo_frac_mat, 
                         name = 'exposure', 
                         col = colorRamp2(c(0, 0.5, 1), c('#006CB1', 'white', '#C9372E')), 
                         cluster_rows = FALSE, 
                         cluster_columns = clust_col, 
                         show_column_names = FALSE, 
                         column_dend_height = unit(2, 'cm'), 
                         height = unit(6.5, 'cm'), 
                         top_annotation = anno_col)


arm_cnv_heatmap <- Heatmap(arm_cnv_mat, 
                           name = 'cnv', 
                           col = structure(c('#383584', 'white', '#B0191D'), names = c('-1', '0', '1')), 
                           cluster_rows = FALSE, 
                           cluster_columns = FALSE, 
                           show_column_names = FALSE, 
                           height = unit(12, 'cm'))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_clustering/CNSig_consensus_clustering_heatmap.pdf', width = 16, height = 10)

cnsig_heatmap %v% arm_cnv_heatmap

dev.off()


# CN Signature definition plotting ----------------------------------------

glioma_cnsig_def_tib <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_def_tib.tsv')

cnsig_plot_dat <- glioma_cnsig_def_tib %>% 
  pivot_longer(-CNclass, names_to = 'CNSig', values_to = 'Weight') %>% 
  separate(CNclass, into = c('Total_CN', 'Hetero_State', 'Size'), sep = ':') %>% 
  mutate(Hetero_State = case_when(
    Hetero_State %in% 'het' ~ 'Heterozygous', 
    Hetero_State %in% 'homdel' ~ 'HD', 
    TRUE ~ Hetero_State
  ), 
  Hetero_State = factor(Hetero_State, levels = c('HD', 'LOH', 'Heterozygous')), 
  Total_CN = factor(Total_CN, levels = c('0', '1', '2', '3-4', '5-8', '9+')))

tcn_color <- tibble(tcn = c('0', '1', '2', '3-4', '5-8', '9+'), 
                    color = c('#0000CD', '#545454', '#228B22', '#7D26CD', '#CD8500', '#8B0A50'))

# Hetero_State color: #A9A9A9

CNSig1_plot <- cnsig_plot_dat %>% 
  filter(CNSig %in% 'CNSig1') %>% 
  ggplot(aes(Size, Weight, fill = Total_CN)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, y = NULL, fill = NULL) +
  scale_fill_manual(breaks = tcn_color$tcn, values = tcn_color$color) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = 'none', 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        strip.background.x = element_rect(color = NA), 
        strip.background.y = element_blank(), 
        strip.text.y = element_text(angle = 0), 
        panel.spacing = unit(0.05, 'cm')) +
  facet_nested(CNSig ~ Hetero_State + Total_CN)

cnsig_plot_list <- list(CNSig1_plot, cnsig_plot('CNSig2'), cnsig_plot('CNSig3'), cnsig_plot('CNSig4'), cnsig_plot('CNSig5'), cnsig_plot('CNSig6'))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/glioma_cnsig_def_barplot.pdf')

wrap_plots(cnsig_plot_list, ncol = 1)

dev.off()


# plot function for cn signatures excluding CNSig1


cnsig_plot <- function(cnsig_name){
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/temp1.pdf')
  
  cnsig_plot_dat %>% 
    filter(CNSig %in% cnsig_name) %>% 
    ggplot(aes(Size, Weight, fill = Total_CN)) +
    geom_bar(stat = 'identity') +
    labs(x = NULL, y = NULL, fill = NULL) +
    scale_fill_manual(breaks = tcn_color$tcn, values = tcn_color$color) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = 'none', 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_blank(), 
          strip.text.y = element_text(angle = 0), 
          panel.spacing = unit(0.05, 'cm')) +
    facet_nested(CNSig ~ Hetero_State + Total_CN)
  
  # dev.off()
  
}









