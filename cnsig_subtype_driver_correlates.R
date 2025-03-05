
library(tidyverse)

source('/boot3/bio_liaojl/pub6/Temp/Liaojl/Code/common_R_functions/getGeneSamMat.R')

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')
glioma_mut_data <- vroom::vroom('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_mut_data.tsv')
glioma_cancer_gene <- read_tsv('/boot3/bio_liaojl/Common_data/tcga_cancer_gene.tsv') %>% filter(Cancer %in% c('LGG', 'GBM')) %>% distinct(Gene)
glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')
# 32 genes

# glioma_SMG_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/glioma_SMG_cell2016.tsv')
# get same results when using glioma_SMG_data


# cancer gene mutation frequency ------------------------------------------

# get cancer gene mutation status

cancer_gene_mut_data <- glioma_mut_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  semi_join(glioma_cancer_gene, by = c('Hugo_Symbol' = 'Gene'))
# 32 cancer genes

total_mut_sam <- glioma_mut_data %>% distinct(sample) %>% pull(sample)

driver_gs_mat <- getGeneSamMat(cancer_gene_mut_data, 'Hugo_Symbol', 'sample', total_mut_sam = total_mut_sam, status = 'ToF')
driver_gs_status <- driver_gs_mat %>% 
  t() %>% 
  as_tibble(rownames = 'sample')
# 900 samples, 32 genes

# test

glioma_cancer_gene_test_data <- glioma_cnsig_subtype_label %>% 
  select(-patient) %>% 
  inner_join(driver_gs_status, by = 'sample') %>% 
  pivot_longer(-(sample:cluster), names_to = 'Hugo_Symbol', values_to = 'status')

save(glioma_cancer_gene_test_data, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/glioma_cancer_gene_test_data.Rdata')

glioma_cancer_gene_stat <- glioma_cancer_gene_test_data %>% 
  group_by(Hugo_Symbol, cluster) %>% 
  summarise(mut_num = sum(status), 
            sam_num = n(), 
            mut_frac = mean(status)) %>% 
  ungroup()

glioma_cancer_gene_test_res <- glioma_cancer_gene_test_data %>% 
  group_by(Hugo_Symbol) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~fisher.test(.$cluster, .$status)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  filter(q_val < 0.05) %>% 
  select(-data)
# 10 significant genes

# write_tsv(glioma_cancer_gene_test_res, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/glioma_cancer_gene_test_res.tsv')

glioma_cancer_gene_stat_signif <- glioma_cancer_gene_stat %>% semi_join(glioma_cancer_gene_test_res, by = 'Hugo_Symbol')

# write_tsv(glioma_cancer_gene_stat_signif, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/glioma_cancer_gene_stat_signif.tsv')

# Oncogenic pathway mutated frequency -------------------------------------

# 10 Oncogenic pathways

onco_pathway_gene <- read_tsv('/boot3/bio_liaojl/Common_data/curated_pathways.txt') %>% select(Gene, Pathway)

onco_pathway_mut_data <- glioma_mut_data %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  inner_join(onco_pathway_gene, by = c('Hugo_Symbol' = 'Gene'))

onco_pathway_sam_mat <- getGeneSamMat(onco_pathway_mut_data, 'Pathway', 'sample', total_mut_sam = total_mut_sam, status = 'ToF')
onco_pathway_sam_status <- onco_pathway_sam_mat %>% t() %>% as_tibble(rownames = 'sample')


# test

glioma_oncopath_test_data <- glioma_cnsig_subtype_label %>% 
  select(-patient) %>% 
  inner_join(onco_pathway_sam_status, by = 'sample') %>% 
  pivot_longer(-(sample:cluster), names_to = 'pathway', values_to = 'status') %>% 
  group_by(pathway) %>% 
  filter(mean(status) > 0.01) %>% # NRF2 was excluded, only 2 mutated samples
  ungroup()
# 849 samples, 10 pathways

glioma_oncopath_test_res <- glioma_oncopath_test_data %>% 
  group_by(pathway) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~fisher.test(.$cluster, .$status)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  filter(q_val < 0.05) %>% 
  select(-data)
# 7 significant pathways

# write_tsv(glioma_oncopath_test_res, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/glioma_oncopath_test_res.tsv')

# plotting

oncopath_text_data <- glioma_oncopath_test_data %>% 
  group_by(cluster, pathway) %>% 
  summarise(mut_num = sum(status), 
            sam_num = n(), 
            mut_frac = mean(status)) %>% 
  ungroup() %>% 
  mutate(text_lab = str_c(signif(mut_frac, digits = 3), '(', mut_num, '/', sam_num, ')'), 
         x_pos = mut_frac + 0.2)


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_oncopath_frac_barplot.pdf', width = 10, height = 3)

glioma_oncopath_test_data %>% 
  ggplot(aes(y = pathway)) +
  geom_bar(aes(fill = status), position = 'fill') +
  geom_text(aes(x = x_pos, label = text_lab), data = oncopath_text_data, size = 3) +
  scale_fill_manual(breaks = c('TRUE', 'FALSE'), values = c('#C0392B', '#ECF0F1')) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  facet_grid(. ~ cluster) +
  theme_bw() +
  theme(legend.position = 'none', 
        axis.ticks.y = element_blank(), 
        panel.grid = element_blank(), 
        strip.text.x = element_text(margin = margin(0.25, 0, 0.25, 0, "cm"), size = 12), 
        panel.spacing.x = unit(0.4, 'cm'))

dev.off()




# Cancer gene interaction -------------------------------------------------

library(maftools)

# remove low frequency glioma cancer genes

cnsig_subtype_mut_sam <- glioma_mut_data %>% semi_join(glioma_cnsig_subtype_label, by = 'sample') %>% distinct(sample)
# 849

cnsig_subtype_cancer_gene_freq <- glioma_mut_data %>% 
  semi_join(glioma_cnsig_subtype_label, by = 'sample') %>% 
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Frame_Shift_Del', 'Frame_Shift_Ins')) %>% 
  semi_join(glioma_cancer_gene, by = c('Hugo_Symbol' = 'Gene')) %>% 
  distinct(sample, Hugo_Symbol) %>% 
  count(Hugo_Symbol) %>% 
  mutate(total_sam = nrow(cnsig_subtype_mut_sam), frac = n/total_sam) %>% 
  arrange(desc(frac)) %>% 
  filter(n >= 20 & frac >= 0.02)
# 19 genes

glioma_mut_data_maf <- glioma_mut_data %>% 
  semi_join(glioma_cnsig_subtype_label, by = 'sample') %>% 
  select(Tumor_Sample_Barcode = sample, Hugo_Symbol:Tumor_Seq_Allele2)

# write_tsv(glioma_mut_data_maf, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_mut.maf')

glioma_mut_data_maf_nest <- glioma_mut_data_maf %>% 
  inner_join(select(glioma_cnsig_subtype_label, -patient), by = c('Tumor_Sample_Barcode' = 'sample')) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  ungroup()


# write_tsv(glioma_mut_data_maf_nest$data[[1]], '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_D_mut.maf')
# write_tsv(glioma_mut_data_maf_nest$data[[2]], '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_T_mut.maf')
# write_tsv(glioma_mut_data_maf_nest$data[[3]], '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_D_CIN_mut.maf')
# write_tsv(glioma_mut_data_maf_nest$data[[4]], '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_T_CTH_mut.maf')

# cancer gene mutation interactions (Co-occurence, Mutually exclusive)

D_maf <- read.maf(maf = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_D_mut.maf')
T_maf <- read.maf(maf = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_T_mut.maf')
D_CIN_maf <- read.maf(maf = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_D_CIN_mut.maf')
T_CTH_maf <- read.maf(maf = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_T_CTH_mut.maf')

D_dri_interation_res <- somaticInteractions(maf = D_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol) %>% as_tibble() %>% mutate(q_val = p.adjust(pValue, method = 'BH'))
T_dri_interation_res <- somaticInteractions(maf = T_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol) %>% as_tibble() %>% mutate(q_val = p.adjust(pValue, method = 'BH'))
D_CIN_dri_interation_res <- somaticInteractions(maf = D_CIN_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol) %>% as_tibble() %>% mutate(q_val = p.adjust(pValue, method = 'BH'))
T_CTH_dri_interation_res <- somaticInteractions(maf = T_CTH_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol) %>% as_tibble() %>% mutate(q_val = p.adjust(pValue, method = 'BH'))

cnsig_subtype_cancer_gene_inter <- list(Diploid = D_dri_interation_res, Tetraploid = T_dri_interation_res, `Diploid CIN` = D_CIN_dri_interation_res, `Tetraploid CTH` = T_CTH_dri_interation_res) %>% 
  map_df(bind_rows, .id = 'Subtype')

write_tsv(cnsig_subtype_cancer_gene_inter, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/cnsig_subtype_cancer_gene_inter.tsv')

# plotting


gene_order <- c('IDH1', 'TP53', 'ATRX', 'CIC', 'FUBP1', 'PIK3CA', 'NOTCH1', 'PTEN', 'NF1', 'ARID1A', 'PIK3R1', 'EGFR', 'IDH2', 'NIPBL', 'ZBTB20', 'SMARCA4', 'RB1', 'SPTA1', 'KEL')

# Diploid
pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/Diploid_cancer_gene_inter.pdf')
somaticInteractions(maf = D_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol, geneOrder = rev(gene_order), showSum = FALSE)
dev.off()

# Tetraploid
pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/Tetraploid_cancer_gene_inter.pdf')
somaticInteractions(maf = T_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol, geneOrder = rev(gene_order), showSum = FALSE)
dev.off()

# Diploid CIN
pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/Diploid_CIN_cancer_gene_inter.pdf')
somaticInteractions(maf = D_CIN_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol, geneOrder = rev(gene_order), showSum = FALSE)
dev.off()

# Tetraploid CTH
pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/Tetraploid_CTH_cancer_gene_inter.pdf')
somaticInteractions(maf = T_CTH_maf, top = 50, genes = cnsig_subtype_cancer_gene_freq$Hugo_Symbol, geneOrder = rev(gene_order), showSum = FALSE)
dev.off()

# c('#8C510A', '#F0F3F3', '#01665E')


# SBS signature frequency -------------------------------------------------

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_sbs_sig_weight.Rdata')

glioma_subtype_sbs_sig_data <- glioma_cnsig_subtype_label %>% inner_join(glioma_sbs_sig_weight, by = c('sample' = 'Samples'))
# 780 samples

glioma_sbs_sig_test_res <- glioma_subtype_sbs_sig_data %>% 
  group_by(Signature) %>% 
  filter(mean(Weight > 0) >= 0.01) %>% # signature frequency >= 1%
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~kruskal.test(Weight ~ cluster, data = .)$p.value, NA)), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  filter(q_val < 0.05) %>% 
  select(-data)

# write_tsv(glioma_sbs_sig_test_res, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/glioma_sbs_sig_test_res.tsv')

# oncoplot ----------------------------------------------------------------

require(maftools)

glioma_anno_cli <- glioma_cnsig_subtype_label %>% 
  inner_join(select(glioma_cli_data, sample, TMB), by = 'sample') %>% 
  mutate(logTMB = log10(TMB+1)) %>% 
  select(Tumor_Sample_Barcode = sample, cluster, logTMB)

# write_tsv(glioma_anno_cli, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_anno_cli.tsv')

glioma_maf <- read.maf(maf = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_mut.maf', 
                       clinicalData = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_maftools/glioma_anno_cli.tsv', 
                       verbose = FALSE)

cluster_colors <- RColorBrewer::brewer.pal(n = 4, name = 'Spectral')
names(cluster_colors) <- c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH')

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/cnsig_subtype_oncoplot.pdf', width = 14, height = 8)

onco_sam_order <- oncoplot(maf = glioma_maf, 
                        genes = glioma_cancer_gene$Gene, 
                        clinicalFeatures = 'cluster', 
                        topBarData = 'logTMB', 
                        borderCol = NA, 
                        drawRowBar = FALSE, 
                        bgCol = NA, 
                        sortByAnnotation = TRUE, 
                        annotationColor = list(cluster = cluster_colors), 
                        writeMatrix = TRUE, 
                        drawBox = TRUE, 
                        removeNonMutated = FALSE, 
                        anno_height = 0.5, 
                        showTitle = FALSE, 
                        showPct = FALSE)

dev.off()


# cancer gene variant classification among subtypes

cancer_gene_varclass_subtype <- glioma_mut_data %>% 
  semi_join(glioma_cancer_gene, by = c('Hugo_Symbol' = 'Gene')) %>% 
  filter(Variant_Classification %in% c('Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Translation_Start_Site')) %>% 
  inner_join(glioma_cnsig_subtype_label, by = 'sample') %>% 
  select(sample, Hugo_Symbol, Variant_Classification, cluster) %>% 
  group_by(sample, Hugo_Symbol) %>% 
  mutate(type = ifelse(n() >= 2, 'Multi_Hit', Variant_Classification)) %>% 
  ungroup() %>% 
  distinct(sample, Hugo_Symbol, cluster, type)

cancer_gene_varclass_subtype_test_res <- cancer_gene_varclass_subtype %>% 
  group_by(Hugo_Symbol) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~chisq.test(.$cluster, .$type)$p.value, NA)), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  filter(p_val < 0.05)
# non significant after BH adjustification

# Hugo_Symbol data                p_val q_val
# <chr>       <list>              <dbl> <dbl>
# TP53        <tibble [336 × 3]> 0.0320 0.298
# PTPN11      <tibble [19 × 3]>  0.0297 0.298
# CIC         <tibble [107 × 3]> 0.0260 0.298


cancer_gene_varclass_subtype_stat <- cancer_gene_varclass_subtype %>% 
  count(Hugo_Symbol, cluster, type) %>% 
  group_by(Hugo_Symbol, cluster) %>% 
  mutate(frac = signif(n/sum(n), digits = 3) , 
         label = str_c(n, frac, sep = '/')) %>% 
  ungroup() %>% 
  select(-(n:frac)) %>% 
  pivot_wider(names_from = cluster, values_from = label, values_fill = '0/0') %>% 
  semi_join(cancer_gene_varclass_subtype_test_res, by = 'Hugo_Symbol')

write_tsv(cancer_gene_varclass_subtype_stat, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/cancer_gene_varclass_subtype_stat.tsv')

# SBS signature plotting --------------------------------------------------

sbs_sig_freq <- glioma_subtype_sbs_sig_data %>% 
  group_by(Signature) %>% 
  summarise(sbs_sig_count = sum(Weight > 0), sam_count = n(), sbs_sig_freq = mean(Weight > 0)) %>% 
  arrange(desc(sbs_sig_freq))
# signature freq < 0.03
# Other: SBS2, 7c, 7d, 10a, 13, 14

sbs_sig_plot_data <- tibble(sample = onco_sam_order) %>% 
  left_join(select(glioma_subtype_sbs_sig_data, sample, Signature, Weight), by = 'sample') %>% 
  mutate(sample = factor(sample, levels = onco_sam_order), 
         Signature = ifelse(Signature %in% str_c('SBS', c(2, '7c', '7d', '10a', 13, 14)), 'Other', Signature)) %>% 
  group_by(sample, Signature) %>% 
  summarise(Weight = sum(Weight)) %>% 
  ungroup()

sbs_sig_order <- glioma_subtype_sbs_sig_data %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE) %>% c(., 'Other')

sbs_signature_colors <- tribble(
  ~Signature, ~Color, 
  'SBS1', '#8F6504', 
  'SBS5', '#927C64', 
  'SBS37', '#C41707', 
  'SBS40', '#FFB50E', 
  'SBS6', '#FF6444', 
  'SBS7b', '#FFFA00', 
  'SBS10b', '#FF5FB3', 
  'SBS11', '#0373C6', 
  'SBS15', '#9C8CC4', 
  'SBS19', '#90AD1C', 
  'SBS30', '#B10DA1', 
  'SBS42', '#CC6677', 
  'Other', '#EAEAEA'
)


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/sbs_sig_weight_barplot.pdf', width = 12, height = 2.5)

sbs_sig_plot_data %>% 
  mutate(Signature = factor(Signature, levels = sbs_sig_order)) %>% 
  ggplot(aes(sample, Weight, fill = Signature)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, y = 'Signature weight', fill = NULL) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(breaks = sbs_signature_colors$Signature, values = sbs_signature_colors$Color) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())
  

dev.off()


# plotting prevalence and weight of SBS signatures among individual subtypes

SBSsig_order <- glioma_subtype_sbs_sig_data %>% distinct(Signature) %>% pull(Signature) %>% str_sort(numeric = TRUE)

sbs_sig_pw_prep <- glioma_subtype_sbs_sig_data %>% 
  mutate(Signature = factor(Signature, levels = SBSsig_order)) %>% 
  select(sample, CNSig = Signature, weight = Weight, cluster) %>% 
  group_by(cluster) %>% 
  nest() %>% 
  ungroup() %>% 
  # mutate(cluster = factor(cluster, levels = c('Diploid CIN', 'Diploid', 'Tetraploid CTH', 'Tetraploid'))) %>% 
  arrange(cluster)

sbs_signature_colors_full <- tribble(
  ~Signature, ~Color, 
  'SBS1', '#8F6504', 
  'SBS5', '#927C64', 
  'SBS37', '#C41707', 
  'SBS40', '#FFB50E', 
  'SBS6', '#FF6444', 
  'SBS7b', '#FFFA00', 
  'SBS10b', '#FF5FB3', 
  'SBS11', '#0373C6', 
  'SBS15', '#9C8CC4', 
  'SBS19', '#90AD1C', 
  'SBS30', '#B10DA1', 
  'SBS42', '#CC6677', 
  'SBS2', '#FFD599', 
  'SBS7c', '#DEA0FD', 
  'SBS7d', '#FFB6C1', 
  'SBS10a', '#88CCEE', 
  'SBS13', '#525975', 
  'SBS14', '#7ED7D1'
)




# function: plotting prevalence and weight of SBS signatures among individual subtypes


pw_plot <- function(sig_weight, plot_title = NULL){
  
  # sig_weight <- sbs_sig_pw_prep$data[[1]]
  
  sig_weight_filter <- sig_weight %>% group_by(CNSig) %>% filter(sum(weight > 0) >= 1) %>% ungroup()
  
  bar_p_data <- sig_weight_filter %>% 
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
    scale_fill_manual(breaks = sbs_signature_colors_full$Signature, values = sbs_signature_colors_full$Color) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank())
  
  box_w <- sig_weight_filter %>% 
    filter(weight > 0) %>% 
    ggplot(aes(CNSig, weight, col = CNSig)) +
    geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
    geom_boxplot(outlier.color = NA, fill = NA) +
    labs(x = NULL, y = 'Weight', title = plot_title) +
    scale_color_manual(breaks = sbs_signature_colors_full$Signature, values = sbs_signature_colors_full$Color) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank(), 
          plot.title = element_text(hjust = 0.5))

  comb_plot <- box_w / bar_p
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbs_signature_prevalence_comb.pdf', width = 6, height = 6)
  # 
  # comb_plot
  # 
  # dev.off()
}


sbs_sig_pw_plots <- sbs_sig_pw_prep %>% mutate(plots = map2(data, cluster, ~pw_plot(.x, .y)))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbs_signature_prevalence_comb.pdf', width = 10, height = 10)

wrap_plots(sbs_sig_pw_plots$plots)

dev.off()

# statistics

sbssig_stat_prep_data <- glioma_subtype_sbs_sig_data %>% 
  mutate(Signature = factor(Signature, levels = SBSsig_order)) %>% 
  select(sample, CNSig = Signature, weight = Weight, cluster)

# prevalence

subtype_sbssig_stat_prev <- sbssig_stat_prep_data %>% 
  group_by(CNSig, cluster) %>% 
  summarise(frac = signif(mean(weight>0), digits = 3)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = cluster, values_from = frac)

write_tsv(subtype_sbssig_stat_prev, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbssig_stat_prev.tsv')

# median in samples where signatures are active

subtype_sbssig_stat_pww_med <- sbssig_stat_prep_data %>% 
  filter(weight > 0) %>% 
  group_by(CNSig, cluster) %>% 
  summarise(med_w = signif(median(weight), digits = 3)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = cluster, values_from = med_w)

write_tsv(subtype_sbssig_stat_pww_med, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbssig_stat_pww_med.tsv')

# median in all samples

subtype_sbssig_stat_all_med <- sbssig_stat_prep_data %>% 
  group_by(CNSig, cluster) %>% 
  summarise(med_w = signif(median(weight), digits = 3)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = cluster, values_from = med_w)

write_tsv(subtype_sbssig_stat_all_med, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbssig_stat_all_med.tsv')

# mean in all samples

subtype_sbssig_stat_all_mean <- sbssig_stat_prep_data %>% 
  group_by(CNSig, cluster) %>% 
  summarise(med_w = signif(mean(weight), digits = 3)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = cluster, values_from = med_w)

write_tsv(subtype_sbssig_stat_all_mean, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/subtype_sbssig_stat_all_mean.tsv')














# Supplementary Figure 2 --------------------------------------------------

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_cli_data.Rdata')

signif_cancer_gene_dat <- glioma_cancer_gene_test_data %>% 
  semi_join(glioma_cancer_gene_test_res, by = 'Hugo_Symbol') %>% 
  select(-cluster) %>% 
  mutate(status = ifelse(status == TRUE, 'MUT', 'WT'), 
         status = factor(status, levels = c('WT', 'MUT'))) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = status)
  
signif_sbs_signature_dat <- glioma_subtype_sbs_sig_data %>% 
  semi_join(glioma_sbs_sig_test_res, by = 'Signature') %>% 
  select(sample, Signature, Weight) %>% 
  pivot_wider(names_from = Signature, values_from = Weight)

supp_plot_data <- glioma_cnsig_subtype_cli_data %>% 
  select(cluster, sample, TMB) %>% 
  drop_na() %>% 
  full_join(signif_cancer_gene_dat, by = 'sample') %>% 
  full_join(signif_sbs_signature_dat, by = 'sample')


### continuous variables (boxplot + point plot)

pb_plot <- function(var_name){
  
  # var_name <- 'age'
  
  if(var_name %in% 'TMB'){ # log10(TMB)
    
    sin_pb_plot <- supp_plot_data %>% 
      select(cluster, sin_var = all_of(var_name)) %>% 
      drop_na() %>% 
      ggplot(aes(cluster, sin_var, col = cluster)) +
      geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
      geom_boxplot(outlier.color = NA, fill = NA) +
      ggpubr::stat_compare_means() +
      labs(x = NULL, y = var_name) +
      scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      scale_y_log10() +
      theme_classic() +
      theme(legend.position = 'none', 
            axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.ticks.x = element_blank())
    
  }else{
    
    sin_pb_plot <- supp_plot_data %>% 
      select(cluster, sin_var = all_of(var_name)) %>% 
      drop_na() %>% 
      ggplot(aes(cluster, sin_var, col = cluster)) +
      geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
      geom_boxplot(outlier.color = NA, fill = NA) +
      ggpubr::stat_compare_means() +
      labs(x = NULL, y = str_c(var_name, ' weight')) +
      scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      theme_classic() +
      theme(legend.position = 'none', 
            axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.ticks.x = element_blank())
    
  }
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/temp.pdf')
  # 
  # sin_pb_plot
  # 
  # dev.off()
  
}


### categorical variables ( barplot)

bar_plot <- function(var_name){
  
  # var_name <- 'histologic_type'
  
  sin_bar_plot <- supp_plot_data %>% 
    select(cluster, sin_var = all_of(var_name)) %>% 
    drop_na() %>% 
    ggplot(aes(cluster, fill = sin_var)) +
    geom_bar(position = 'fill') +
    labs(x = NULL, y = 'Proportion', fill = var_name) +
    scale_fill_grey(start = 0.8, end = 0.2) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    theme_classic() +
    theme(legend.position = 'top', 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.ticks.x = element_blank())
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/temp.pdf')
  # 
  # sin_bar_plot
  # 
  # dev.off()
  
}


### combine plots

comb_plots <- list(pb_plot('TMB'), bar_plot('TP53'), bar_plot('EGFR'), bar_plot('PTEN'), bar_plot('IDH1'), 
                   bar_plot('ZBTB20'), bar_plot('NOTCH1'), bar_plot('ARID1A'), bar_plot('IDH2'), bar_plot('CIC'), bar_plot('FUBP1'), 
                   pb_plot('SBS1'), pb_plot('SBS5'), pb_plot('SBS7d'), pb_plot('SBS37'), pb_plot('SBS40'))


# age, histologic_type, grade, purity, ploidy, aneuploidy_score, Intratumor_Heterogeneity, CNH, IDH_codel_subtype, transcriptome_subtype

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/Figure_S2.pdf', width = 18, height = 10)

wrap_plots(comb_plots) + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A')

dev.off()


### for thesis

comb_plots_thesis <- list(pb_plot('TMB'), bar_plot('TP53'), bar_plot('EGFR'), bar_plot('PTEN'), bar_plot('IDH1'), 
                   bar_plot('ZBTB20'), bar_plot('NOTCH1'), bar_plot('ARID1A'), bar_plot('IDH2'), bar_plot('CIC'), bar_plot('FUBP1'))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_mut_difference/CNSig_subtype_cancer_gene_diff_thesis.pdf', width = 10, height = 11)

wrap_plots(comb_plots_thesis) + plot_layout(nrow = 3) + plot_annotation(tag_levels = 'A')

dev.off()



# gistic2 recurrent scnv analysis (discarded) ------------------------------

# gistic2 cytoband status

tcga_gistic2_data <- vroom::vroom('/boot3/bio_liaojl/Common_data/all_thresholded.by_genes_whitelisted.tsv') %>% rename_all(~str_sub(., 1, 15))

glioma_gistic2_data <- tcga_gistic2_data %>% 
  select(Hugo_Symbol = `Gene Symbol`, cytoband = Cytoband, all_of(intersect(glioma_cnsig_subtype_label$sample, colnames(tcga_gistic2_data))))

glioma_gistic_cytoband_status <- glioma_gistic2_data %>% 
  select(-Hugo_Symbol) %>% 
  distinct(cytoband, .keep_all = TRUE) %>% 
  column_to_rownames('cytoband') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  as_tibble() %>% 
  pivot_longer(-sample, names_to = 'cytoband', values_to = 'status') %>% 
  mutate(status = case_when(
    status >= 1 ~ 'Amp', 
    status <= -1 ~ 'Del', 
    TRUE ~ 'Neutral'
  ))

recurrent_amp_status <- glioma_gistic_cytoband_status %>% 
  mutate(status = ifelse(status %in% 'Amp', 'Amp', 'non-Amp')) %>% 
  group_by(cytoband) %>% 
  filter(mean(status %in% 'Amp') >= 0.25) %>% 
  ungroup()


amp_test_res <- recurrent_amp_status %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample') %>% 
  group_by(cytoband) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(test_res = map_dbl(data, ~fisher.test(.$status, .$cluster)$p.value))
  



















