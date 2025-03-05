
library(tidyverse)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(limma)

# load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')
# glioma_exp_data <- vroom::vroom('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_exp_data.tsv') %>% 
#   mutate_if(is.numeric, ~log2(.+1))

load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Data/Processed/glioma_cnsig_subtype_label.Rdata')
glioma_exp_data <- vroom::vroom('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Data/Processed/glioma_exp_data.tsv') %>% 
  mutate_if(is.numeric, ~log2(.+1))

common_label <- glioma_cnsig_subtype_label %>% 
  filter(sample %in% colnames(glioma_exp_data)) %>% 
  dplyr::select(-patient)
# 625 glioma samples

glioma_exp_mat <- glioma_exp_data %>% column_to_rownames('Hugo_Symbol') %>% as.matrix()


# function: limma DEG identification in pairwise subtype (e.g., Diploid CIN vs Diploid)

get_DE_res_re_pair <- function(sam_group_data){
  
  # sam_group_data <- D_sam_group_dat

  sam_group_dat <- sam_group_data %>% arrange(cluster)
  
  design <- model.matrix(~ sam_group_dat$cluster)
  
  fit <- lmFit(glioma_exp_mat[, sam_group_dat$sample], design)
  fit <- eBayes(fit) 
  allDiff <- topTable(fit, number = Inf)
  
  allDiff_label <- allDiff %>% 
    rownames_to_column('Hugo_Symbol') %>% 
    as_tibble() %>% 
    mutate(status = case_when(
      logFC >= 1 & adj.P.Val < 0.05 ~ 'up', 
      logFC <= -1 & adj.P.Val < 0.05 ~ 'down', 
      TRUE ~ 'non_DEG'
    ), label = ifelse(status %in% 'non_DEG', status, 'DEG'))
  
}

# 1. Diploid CIN vs Diploid

D_sam_group_dat <- common_label %>% filter(cluster %in% c('Diploid', 'Diploid CIN')) %>% mutate(cluster = factor(cluster, levels = c('Diploid', 'Diploid CIN')))
D_within_DE_res <- get_DE_res_re_pair(D_sam_group_dat) %>% mutate(cluster = 'Diploid CIN vs Diploid') %>% select(cluster, everything()) # 1833 significant

# 2. Tetraploid CTH vs Tetraploid

T_sam_group_dat <- common_label %>% filter(cluster %in% c('Tetraploid', 'Tetraploid CTH')) %>% mutate(cluster = factor(cluster, levels = c('Tetraploid', 'Tetraploid CTH')))
T_within_DE_res <- get_DE_res_re_pair(T_sam_group_dat) %>% mutate(cluster = 'Tetraploid CTH vs Tetraploid') %>% select(cluster, everything()) # 666 significant

# 3. (Tetraploid + Tetraploid CTH) vs (Diploid + Diploid CIN)

DT_sam_group_dat <- common_label %>% mutate(cluster = ifelse(cluster %in% c('Diploid', 'Diploid CIN'), 'Diploid_all', 'Tetraploid_all'), cluster = factor(cluster, levels = c('Diploid_all', 'Tetraploid_all')))
DT_within_DE_res <- get_DE_res_re_pair(DT_sam_group_dat) %>% mutate(cluster = 'Tetraploid all vs Diploid all') %>% select(cluster, everything()) # 0 significant

# 4. Tetraploid vs Diploid

T_vDsam_group_dat <- common_label %>% filter(cluster %in% c('Diploid', 'Tetraploid')) %>% mutate(cluster = factor(cluster, levels = c('Diploid', 'Tetraploid')))
T_vD_DE_res <- get_DE_res_re_pair(T_vDsam_group_dat) %>% mutate(cluster = 'Tetraploid vs Diploid') %>% select(cluster, everything()) # 0 significant

# 5. Tetraploid CTH vs Diploid CIN

TC_vDCsam_group_dat <- common_label %>% filter(cluster %in% c('Diploid CIN', 'Tetraploid CTH')) %>% mutate(cluster = factor(cluster, levels = c('Diploid CIN', 'Tetraploid CTH')))
TC_vDC_DE_res <- get_DE_res_re_pair(TC_vDCsam_group_dat) %>% mutate(cluster = 'Tetraploid CTH vs Diploid CIN') %>% select(cluster, everything()) # 1 significant

# 6. Tetraploid vs Diploid CIN

T_vDCsam_group_dat <- common_label %>% filter(cluster %in% c('Diploid CIN', 'Tetraploid')) %>% mutate(cluster = factor(cluster, levels = c('Diploid CIN', 'Tetraploid')))
T_vDC_DE_res <- get_DE_res_re_pair(T_vDCsam_group_dat) %>% mutate(cluster = 'Tetraploid vs Diploid CIN') %>% select(cluster, everything()) # 1279 significant

# 7. Tetraploid CTH vs Diploid

TC_vDsam_group_dat <- common_label %>% filter(cluster %in% c('Diploid', 'Tetraploid CTH')) %>% mutate(cluster = factor(cluster, levels = c('Diploid', 'Tetraploid CTH')))
TC_vD_DE_res <- get_DE_res_re_pair(TC_vDsam_group_dat) %>% mutate(cluster = 'Tetraploid CTH vs Diploid') %>% select(cluster, everything()) # 1224 significant

# combine all pairs excluding 3, 4 and 5

pair_DE_res_comb <- list(D_within_DE_res, T_within_DE_res, T_vDC_DE_res, TC_vD_DE_res) %>% reduce(bind_rows)

save(pair_DE_res_comb, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_analysis/pair_DE_res_comb.Rdata')

pair_DE_signif_comb <- pair_DE_res_comb %>% filter(label %in% 'DEG')

save(pair_DE_signif_comb, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_analysis/pair_DE_signif_comb.Rdata')
# load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DE_signif_comb.Rdata')


library(ggvenn)

venn_dat <- list(`Diploid CIN vs Diploid` = filter(pair_DE_signif_comb, cluster %in% 'Diploid CIN vs Diploid')$Hugo_Symbol, 
                 `Tetraploid CTH vs Tetraploid` = filter(pair_DE_signif_comb, cluster %in% 'Tetraploid CTH vs Tetraploid')$Hugo_Symbol, 
                 `Tetraploid vs Diploid CIN` = filter(pair_DE_signif_comb, cluster %in% 'Tetraploid vs Diploid CIN')$Hugo_Symbol, 
                 `Tetraploid CTH vs Diploid` = filter(pair_DE_signif_comb, cluster %in% 'Tetraploid CTH vs Diploid')$Hugo_Symbol)

ggvenn(venn_dat, c('Diploid CIN vs Diploid', 'Tetraploid CTH vs Tetraploid', 'Tetraploid vs Diploid CIN', 'Tetraploid CTH vs Diploid'))
ggsave('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DEG_vennplot.pdf', width = 8, height = 8)

com_gene <- Reduce(intersect, venn_dat) # 580 common DEGs


# heatmap showing expression of common DEGs --------------------------------------

# load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Data/Processed/glioma_cnsig_subtype_cli_data.Rdata')

set.seed(123456)

heatmap_sam_order <- common_label %>% 
  mutate(cluster = factor(cluster, levels = c('Diploid', 'Tetraploid', 'Diploid CIN', 'Tetraploid CTH'))) %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  mutate(sample_m = sample(sample)) %>% 
  ungroup()

heatmap_sam_order_m <- tibble(sample = heatmap_sam_order$sample_m) %>% left_join(select(heatmap_sam_order, -sample_m), by = 'sample')

common_gene_expr_mat <- glioma_exp_data %>% 
  filter(Hugo_Symbol %in% com_gene) %>% 
  pivot_longer(-Hugo_Symbol, names_to = 'sample', values_to = 'geneExp') %>% 
  group_by(Hugo_Symbol) %>%
  mutate(geneExp = scale(geneExp)[, 1]) %>% # z-score
  ungroup() %>% 
  pivot_wider(names_from = sample, values_from = geneExp) %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  as.matrix()
  
anno_col <- HeatmapAnnotation(Cluster = heatmap_sam_order_m$cluster, 
                              col = list(Cluster = c('Diploid' = '#D7191C', 'Diploid CIN' = '#FDAE61', 'Tetraploid' = '#ABDDA4', 'Tetraploid CTH' = '#2B83BA')))


pdf('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/common_DEG_expr_heatmap.pdf', width = 9, height = 6)

Heatmap(common_gene_expr_mat[, heatmap_sam_order_m$sample], 
        name = ' ', 
        col = colorRamp2(c(-1, 0, 1), c('#006CB1', 'white', '#C9372E')),
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        show_column_names = FALSE, 
        show_row_names = FALSE, 
        show_row_dend = FALSE, 
        top_annotation = anno_col)

dev.off()



# Reactome pathway enrichment analysis (in windows) ------------------------------------


# load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DE_signif_comb.Rdata')
load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DE_signif_comb.Rdata')

# 1. for DEGs between pairwise subtypes

DEG_signif_rectome_res <- pair_DE_signif_comb %>% 
  group_by(cluster) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(gene_id = map(data, ~clusterProfiler::bitr(.$Hugo_Symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID), 
         rectome_res = map(gene_id, ~ReactomePA::enrichPathway(gene = ., pvalueCutoff = 0.05, readable = T)))

# clusterProfiler::enrichGO(.$Hugo_Symbol, OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'SYMBOL')
# ReactomePA::dotplot(DEG_signif_rectome_res$rectome_res[[1]], showCategory = 15) # Diploid CIN vs Diploid

pair_DEG_reactome_res <- DEG_signif_rectome_res %>% dplyr::select(cluster, rectome_res)

save(pair_DEG_reactome_res, file = 'E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DEG_reactome_res.Rdata')

load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DEG_reactome_res.Rdata')


pair_DEG_reactome_signif_tib <- pair_DEG_reactome_res %>% 
  mutate(rectome_tib = map(rectome_res, ~as_tibble(.@result))) %>% 
  select(-rectome_res) %>% 
  unnest(rectome_tib) %>% 
  filter(p.adjust < 0.05)


pair_DEG_reactome_plot_dat <- pair_DEG_reactome_signif_tib %>% 
  mutate(number_temp = as.numeric(str_extract(BgRatio, '\\d+(?=/)')), 
         log_p_adj = -log10(p.adjust)) %>% 
  filter(number_temp >= 95) %>% 
  select(cluster, Description, log_p_adj)

pathway_order <- pair_DEG_reactome_plot_dat %>% distinct(Description) %>% pull(Description) %>% str_sort(decreasing = TRUE)


pdf('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DEG_reactome_heatmap.pdf', width = 6, height = 12)

pair_DEG_reactome_plot_dat %>% 
  mutate(Description = factor(Description, levels = pathway_order)) %>% 
  ggplot(aes(cluster, Description, fill = log_p_adj)) +
  geom_tile(col = 'grey60') +
  labs(x = NULL, y = NULL, fill = '-log10(q)') +
  scale_fill_gradientn(colors = c('#FFE500', '#16898C', '#471C50')) +
  # scale_x_discrete(expand = c(0, 0), position = 'top') +
  scale_x_discrete(expand = c(0, 0)) +
  # scale_x_discrete(expand = c(0, 0), labels = function(x) str_wrap(x, width = 60)) +
  scale_y_discrete(expand = c(0, 0), labels = function(x) str_sub(x, 1, 55)) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# 2. for common DEGs among subtypes

com_gene_id <- clusterProfiler::bitr(com_gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')$ENTREZID
common_DEGs_reactome_res <- ReactomePA::enrichPathway(gene = com_gene_id, pvalueCutoff = 0.05, readable = T)

save(common_DEGs_reactome_res, file = 'E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/common_DEGs_reactome_res.Rdata')

# barplot showing some enriched pathways

com_rectome_plot_data <- common_DEGs_reactome_res@result %>% 
  as_tibble() %>% 
  filter(p.adjust < 0.05) %>% 
  mutate(number_temp = as.numeric(str_extract(BgRatio, '\\d+(?=/)'))) %>% 
  filter(number_temp >= 90) %>% 
  separate(GeneRatio, into = c('num1', 'num2'), sep = '/') %>% 
  mutate(num1 = as.numeric(num1), 
         num2 = as.numeric(num2), 
         GeneRatio = num1/num2, 
         log_p_adj = -log10(p.adjust)) %>% 
  dplyr::select(Description, GeneRatio, log_p_adj) %>% 
  .[1:20, ] %>% 
  arrange(GeneRatio)

com_rectome_barplot <- com_rectome_plot_data %>% 
  mutate(Description = fct_inorder(Description)) %>% 
  ggplot(aes(GeneRatio, Description, fill = log_p_adj)) +
  geom_bar(stat = 'identity') +
  labs(y = NULL, fill = '-log10(q)') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradientn(colors = c('#FFE500', '#16898C', '#471C50')) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank())


pdf('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/com_rectome_barplot.pdf', width = 8, height = 6)

com_rectome_barplot

dev.off()



# volcano plot -----------------------------------------

load('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DE_res_comb.Rdata')

pair_volcano_plot <- pair_DE_res_comb %>% 
  group_by(cluster) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(vol_plot = map2(cluster, data, ~volcano_plot(.x, .y)))


pdf('E:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_analysis/pair_DEG_volcano_plot.pdf', width = 10, height = 10)

wrap_plots(pair_volcano_plot$vol_plot) + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A')

dev.off()


# function: volcano plotting

volcano_plot <- function(t_name, DE_res, max_overlap = 10){
    
  vol_plot_res <- DE_res %>% 
    mutate(log_q = -log10(adj.P.Val), 
           t_label = ifelse(status %in% c('down', 'up'), Hugo_Symbol, '')) %>% 
    ggplot(aes(logFC, log_q, col = status)) +
    labs(x = 'log2 fold change', y = '-log10 q-value', title = t_name) +
    scale_color_manual(breaks = c('down', 'non_DEG', 'up'), values = c('#006CB1', '#EAEAEA', '#C9372E')) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = '#EAEAEA') +
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = '#EAEAEA') +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = t_label), max.overlaps = max_overlap) +
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid = element_blank(), 
          plot.title = element_text(hjust = 0.5))
      
}



# Hallmark gene set ssGSEA analysis ---------------------------------------

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

hallmark_gene <- qusage::read.gmt('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/h.all.v2023.1.Hs.symbols.gmt')

hallmark_ssgsea_score <- ssGSEA_compute(glioma_exp_mat, hallmark_gene)

hallmark_test_dat <- hallmark_ssgsea_score %>% inner_join(select(glioma_cnsig_subtype_label, Sample = sample, cluster), by = 'Sample')

hallmark_test_res <- hallmark_test_dat %>% 
  group_by(gene_set) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, ~kruskal.test(ssGSEA_score ~ cluster, data = .)$p.value), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)

hallmark_test_signif <- hallmark_test_res %>% filter(q_val < 0.05) # 47 significant

# dotplot

hallmark_plot_dat <- hallmark_test_dat %>% 
  group_by(gene_set) %>% 
  mutate(ssGSEA_score = scale(ssGSEA_score)[, 1]) %>% # z-score
  ungroup() %>% 
  group_by(gene_set, cluster) %>% 
  summarise(med_ssgsea_score = median(ssGSEA_score)) %>% 
  ungroup() %>% 
  inner_join(hallmark_test_res, by = 'gene_set') %>% 
  mutate(log_q = -log10(q_val))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_analysis/hallmark_subtype_dotplot.pdf', width = 13, height = 3.5)

hallmark_plot_dat %>% 
  mutate(gene_set = str_remove(gene_set, 'HALLMARK_'), 
         gene_set = str_replace_all(gene_set, '_', ' '), 
         gene_set = str_to_title(gene_set), 
         cluster = factor(cluster, levels = c('Diploid', 'Tetraploid', 'Diploid CIN', 'Tetraploid CTH'))) %>% 
  ggplot(aes(gene_set, cluster, fill = med_ssgsea_score), color = 'grey60') +
  geom_point(aes(size = log_q), shape = 21) +
  labs(x = NULL, y = NULL, fill = 'Median ssGSEA score (z-scored)', size = '-log10(q)') +
  scale_fill_gradientn(colors = c('#006CB1', 'white', '#C9372E')) +
  ggthemes::theme_few() +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = 'top')

dev.off()


