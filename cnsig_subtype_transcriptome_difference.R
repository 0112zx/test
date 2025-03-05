
library(tidyverse)
library(edgeR)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
# library(clusterProfiler)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')
glioma_exp_data <- vroom::vroom('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_exp_data.tsv')

common_label <- glioma_cnsig_subtype_label %>% 
  filter(sample %in% colnames(glioma_exp_data)) %>% 
  select(-patient)
# 625 glioma samples

glioma_exp_filter <- glioma_exp_data %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  as_tibble() %>% 
  select_if(~mean(. > 0) > 0.8) # 17613 -> 14919, where TTR was excluded

glioma_exp_longer <- glioma_exp_filter %>% pivot_longer(-sample, names_to = 'Hugo_Symbol', values_to = 'geneExp')

# DE gene using KW test ------------------------------------

kw_test_data <- glioma_exp_longer %>% 
  inner_join(common_label, by = 'sample') %>% 
  group_by(Hugo_Symbol) %>% 
  nest() %>% 
  ungroup()

DE_kw_test_res <- kw_test_data %>% 
  mutate(res = map(data, ~{
    
    sin_p_val <- kruskal.test(geneExp ~ cluster, data = .)$p.value
    sin_fc_data <- group_by(., cluster) %>% summarise(mean_Exp = mean(geneExp)) %>% arrange(mean_Exp) %>% pull(mean_Exp)
    sin_logfc <- log2(sin_fc_data[3]/sin_fc_data[2])

    tibble(abs_logFC = sin_logfc, p_val = sin_p_val)
    
  })) %>% 
  select(-data) %>% 
  unnest(res) %>% 
  mutate(q_val = p.adjust(p_val, method = 'BH'), 
         status = ifelse(abs_logFC >= 1 & q_val < 0.05, 'signif', 'not.signif'))
# 438 significant genes

save(DE_kw_test_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/DE_kw_test_res.Rdata')
# load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/DE_kw_test_res.Rdata')

# in windows for gene functional enrichment analysis

load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/DE_kw_test_res.Rdata')

DE_kw_test_res_signif <- DE_kw_test_res %>% filter(status %in% 'signif')

DE_kw_go_res <- enrichGO(DE_kw_test_res_signif$Hugo_Symbol, OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'SYMBOL')

save(DE_kw_go_res, file = 'G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/DE_kw_go_res.Rdata')

# when using clusterProfiler::enrichKEGG
# require('R.utils')
# R.utils::setOption("clusterProfiler.download.method",'auto')



# DE gene enrichment plot -------------------------------------------------

DE_enrich_go_plot_data <- DE_kw_go_res@result %>% 
  as_tibble() %>% 
  filter(p.adjust < 0.05) %>% 
  mutate(GeneRatio = Count/369, 
         log_p_adj = -log10(p.adjust)) %>% 
  select(Description, GeneRatio, log_p_adj) %>% 
  arrange(desc(log_p_adj)) %>% 
  .[1:30, ] %>% 
  arrange(GeneRatio)

pdf('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/DE_enrich_go_barplot.pdf', width = 6.5, height = 6.5)

DE_enrich_go_plot_data %>% 
  mutate(Description = fct_inorder(Description)) %>% 
  ggplot(aes(GeneRatio, Description, fill = log_p_adj)) +
  geom_bar(stat = 'identity') +
  labs(y = NULL, fill = '-log10(p.adjust)') +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradientn(colors = c('#FFE500', '#16898C', '#471C50')) +
  theme_bw() +
  theme(axis.ticks.y = element_blank(), 
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank())

dev.off()

# for thesis

pdf('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/DE_enrich_go_barplot_for_thesis.pdf', width = 10, height = 6.5)

DE_enrich_go_plot_data %>% 
  arrange(desc(GeneRatio)) %>% 
  mutate(Description = fct_inorder(Description)) %>% 
  ggplot(aes(Description, GeneRatio, fill = log_p_adj)) +
  geom_bar(stat = 'identity') +
  labs(x = NULL, fill = '-log10(p.adjust)') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradientn(colors = c('#FFE500', '#16898C', '#471C50')) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

dev.off()


# DE gene subtype expression heatmap --------------------------------------

DE_kw_test_res_signif_top <- DE_kw_test_res %>% 
  filter(status %in% 'signif') %>% 
  # arrange(desc(abs_logFC)) %>% 
  arrange(q_val) %>%
  .[1:30, ]

DE_exp_plot_data <- kw_test_data %>% 
  semi_join(DE_kw_test_res_signif_top, by = 'Hugo_Symbol') %>% 
  unnest(data) %>% 
  # group_by(Hugo_Symbol) %>% 
  # mutate(geneExp = scale(geneExp)[, 1]) %>% 
  # ungroup()
  mutate(geneExp = log2(geneExp + 1))

heatmap_sam_order <- DE_exp_plot_data %>% 
  distinct(sample, cluster) %>% 
  arrange(cluster)

DE_exp_plot_mat <- DE_exp_plot_data %>% 
  select(-cluster) %>% 
  pivot_wider(names_from = sample, values_from = geneExp) %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  as.matrix() %>% 
  .[, heatmap_sam_order$sample]


anno_col <- HeatmapAnnotation(Cluster = heatmap_sam_order$cluster, 
                              col = list(Cluster = c('Diploid' = '#D7191C', 'Diploid CIN' = '#FDAE61', 'Tetraploid' = '#ABDDA4', 'Tetraploid CTH' = '#2B83BA')))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/DE_subtype_heatmap.pdf', width = 10, height = 6)
# pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/DE_subtype_heatmap_all.pdf', width = 10, height = 6)

Heatmap(DE_exp_plot_mat, 
        name = 'expression', 
        col = colorRamp2(c(0, 7.5, 15), c('#006CB1', 'white', '#C9372E')),
        # col = colorRamp2(c(-4, 0, 4), c("#3E4E9A", "black", "yellow")),
        cluster_rows = TRUE, 
        cluster_columns = FALSE, 
        show_column_names = FALSE, 
        # row_names_side = 'left', 
        show_row_names = FALSE, 
        show_row_dend = FALSE, 
        top_annotation = anno_col)

dev.off()


# DE gene using Wilcoxon rank-sum test (within) ------------------------------------

glioma_exp_filter <- glioma_exp_data %>% 
  column_to_rownames('Hugo_Symbol') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  as_tibble() %>% 
  select_if(~mean(. > 0) > 0.8) # 17613 -> 14919, where TTR was excluded

glioma_exp_longer <- glioma_exp_filter %>% pivot_longer(-sample, names_to = 'Hugo_Symbol', values_to = 'geneExp')


# function: DE gene using Wilcoxon rank-sum test (within)

wilcox_DE_res <- function(sam_group_data){
  
  # sam_group_data <- D_sam_group_dat
  
  test_data <- glioma_exp_longer %>% 
    inner_join(sam_group_data, by = 'sample') %>% 
    group_by(Hugo_Symbol) %>% 
    nest() %>% 
    ungroup()
  
  test_res <- test_data %>% 
    mutate(res = map(data, ~{
      
      sin_p_val <- wilcox.test(geneExp ~ cluster, data = .)$p.value
      sin_fc_data <- group_by(., cluster) %>% summarise(mean_Exp = mean(geneExp)) %>% arrange(cluster) %>% pull(mean_Exp)
      sin_logfc <- log2(sin_fc_data[2]/sin_fc_data[1])
      
      tibble(logFC = sin_logfc, p_val = sin_p_val)
      
    })) %>% 
    select(-data) %>% 
    unnest(res) %>% 
    mutate(q_val = p.adjust(p_val, method = 'BH'), 
           status = case_when(
             logFC >= 1 & q_val < 0.05 ~ 'up', 
             logFC <= -1 & q_val < 0.05 ~ 'down', 
             TRUE ~ 'not.signif'
           ))
  
}


# 1. Diploid CIN vs Diploid

D_sam_group_dat <- common_label %>% filter(cluster %in% c('Diploid', 'Diploid CIN'))
D_within_DE_wilcox_res <- wilcox_DE_res(D_sam_group_dat) # 1807 significant

save(D_within_DE_wilcox_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/D_within_DE_wilcox_res.Rdata')

# in windows enrichment analysis

load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/D_within_DE_wilcox_res.Rdata')

D_within_DE_wilcox_res_signif <- D_within_DE_wilcox_res %>% filter(status != 'not.signif')
D_within_DE_wilcox_go_res <- enrichGO(D_within_DE_wilcox_res_signif$Hugo_Symbol, OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'SYMBOL')

save(D_within_DE_wilcox_go_res, file = 'G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/D_within_DE_wilcox_go_res.Rdata')

# dotplot(D_within_DE_wilcox_go_res, showCategory = 25)



# 2. Tetraploid CTH vs Tetraploid

T_sam_group_dat <- common_label %>% filter(cluster %in% c('Tetraploid', 'Tetraploid CTH'))
T_within_DE_wilcox_res <- wilcox_DE_res(T_sam_group_dat) # 534 significant

save(T_within_DE_wilcox_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/T_within_DE_wilcox_res.Rdata')

# in windows enrichment analysis

load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/T_within_DE_wilcox_res.Rdata')

T_within_DE_wilcox_res_signif <- T_within_DE_wilcox_res %>% filter(status != 'not.signif')
T_within_DE_wilcox_go_res <- enrichGO(T_within_DE_wilcox_res_signif$Hugo_Symbol, OrgDb = org.Hs.eg.db, ont = 'BP', pAdjustMethod = 'BH', pvalueCutoff = 0.01, qvalueCutoff = 0.05, keyType = 'SYMBOL')

save(T_within_DE_wilcox_go_res, file = 'G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/T_within_DE_wilcox_go_res.Rdata')

# dotplot(T_within_DE_wilcox_go_res, showCategory = 25)



# 3. (Tetraploid + Tetraploid CTH) vs (Diploid + Diploid CIN)

DT_sam_group_dat <- common_label %>% mutate(cluster = ifelse(cluster %in% c('Diploid', 'Diploid CIN'), 'Diploid_all', 'Tetraploid_all'))
DT_within_DE_wilcox_res <- wilcox_DE_res(DT_sam_group_dat) # 1807 significant

save(DT_within_DE_wilcox_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/DT_within_DE_wilcox_res.Rdata')

# !!!! no significant genes


# volcano and enrichment plotting -----------------------------------------


load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/D_within_DE_wilcox_go_res.Rdata')
load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/T_within_DE_wilcox_go_res.Rdata')



D_T_plots_res <- list(D_vol = volcano_plot(D_within_DE_wilcox_res, max_overlap = 20), D_enrich = enrich_barplot(D_within_DE_wilcox_go_res), 
                      T_vol = volcano_plot(T_within_DE_wilcox_res, max_overlap = 30), T_enrich = enrich_barplot(T_within_DE_wilcox_go_res))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/D_T_DE_volcano_enrich_comb.pdf', width = 10, height = 10)

wrap_plots(D_T_plots_res) + plot_layout(nrow = 2, widths = c(2.5, 1)) + plot_annotation(tag_levels = 'A')

dev.off()


# function: volcano plotting

volcano_plot <- function(DE_res, max_overlap = 10){
  
  # DE_res <- D_within_DE_wilcox_res
  
  vol_plot_res <- DE_res %>% 
    mutate(log_q = -log10(q_val), 
           t_label = ifelse(status %in% c('down', 'up'), Hugo_Symbol, '')) %>% 
    ggplot(aes(logFC, log_q, col = status)) +
    labs(x = 'log2 fold change', y = '-log10 q-value') +
    scale_color_manual(breaks = c('down', 'not.signif', 'up'), values = c('#006CB1', '#EAEAEA', '#C9372E')) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = '#EAEAEA') +
    geom_vline(xintercept = c(-1,1), linetype = "dashed", color = '#EAEAEA') +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = t_label), max.overlaps = max_overlap) +
    theme_bw() + 
    theme(legend.position = 'none', 
          panel.grid = element_blank())
      
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/temp.pdf')
  # 
  # vol_plot_res
  # 
  # dev.off()
  
}


# function: gene enrichment barplot plotting

enrich_barplot <- function(DE_go_res){
  
  # DE_go_res <- D_within_DE_wilcox_go_res
  
  DE_go_plot_data <- DE_go_res@result %>% 
    as_tibble() %>% 
    filter(p.adjust < 0.05) %>% 
    separate(GeneRatio, into = c('num1', 'num2'), sep = '/') %>% 
    mutate(num1 = as.numeric(num1), 
           num2 = as.numeric(num2), 
           GeneRatio = num1/num2, 
           log_p_adj = -log10(p.adjust)) %>% 
    select(Description, GeneRatio, log_p_adj) %>% 
    arrange(desc(log_p_adj)) %>% 
    .[1:30, ] %>% 
    arrange(GeneRatio)
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/temp.pdf', width = 6.5, height = 6.5)
  
  plot_res <- DE_go_plot_data %>% 
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
  
  # dev.off()
  
}

# DE gene using edgeR -----------------------------------------------

glioma_exp_mat <- glioma_exp_data %>% column_to_rownames('Hugo_Symbol') %>% as.matrix()

# function: using edgeR including filtering and normalization

get_DE_res <- function(subtype_name){
  
  # subtype_name <- 'Diploid'
  
  sam_group_dat <- common_label %>% 
    mutate(label = ifelse(cluster %in% subtype_name, subtype_name, 'other'), 
           label = factor(label, levels = c('other', subtype_name))) %>% 
    arrange(label)
  
  # DGEList object
  
  dgelist <- DGEList(counts = glioma_exp_mat[, sam_group_dat$sample], group = sam_group_dat$label)
  
  # filter low count data
  
  keep <- rowSums(cpm(dgelist) > 1 ) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  
  # normalization
  
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  
  # differential expression analysis
  
  design <- model.matrix(~sam_group_dat$label)
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE) # estimate expression value diviation
  fit <- glmQLFit(dge, design, robust = TRUE) # model fit
  lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
  
  gdiff_raw <- lrt %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Hugo_Symbol') %>% 
    as_tibble() %>% 
    mutate(status = case_when(
      logFC >= 1 & FDR < 0.01 ~ 'up', 
      logFC <= -1 & FDR < 0.01 ~ 'down', 
      TRUE ~ 'not.signif'
    ))
  
}

D_DE_gene_res <- get_DE_res('Diploid') # 1492 significant
D_CIN_DE_gene_res <- get_DE_res('Diploid CIN') # 1230 significant
T_DE_gene_res <- get_DE_res('Tetraploid') # 357 significant
T_CTH_DE_gene_res <- get_DE_res('Tetraploid CTH') # 474 significant

# subtype_DE_gene_res_list <- list(D_DE_gene_res, D_CIN_DE_gene_res, T_DE_gene_res, T_CTH_DE_gene_res)
# save(subtype_DE_gene_res_list, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/subtype_DE_gene_res_list.Rdata')


# function: using edgeR (filtering and normalization by self)

glioma_exp_filter_mat <- glioma_exp_filter %>% column_to_rownames('sample') %>% t() %>% as.matrix()

get_DE_res_re <- function(subtype_name){
  
  # subtype_name <- 'Tetraploid'
  
  sam_group_dat <- common_label %>% 
    mutate(label = ifelse(cluster %in% subtype_name, subtype_name, 'other'), 
           label = factor(label, levels = c('other', subtype_name))) %>% 
    arrange(label)
  
  # DGEList object
  
  dgelist <- DGEList(counts = glioma_exp_filter_mat[, sam_group_dat$sample], group = sam_group_dat$label)
  
  # differential expression analysis
  
  design <- model.matrix(~sam_group_dat$label)
  dge <- estimateDisp(dgelist, design, robust = TRUE) # estimate expression value diviation
  fit <- glmQLFit(dge, design, robust = TRUE) # model fit
  lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
  
  gdiff_raw <- lrt %>% 
    as.data.frame() %>% 
    rownames_to_column(var = 'Hugo_Symbol') %>% 
    as_tibble() %>% 
    mutate(status = case_when(
      logFC >= 1 & FDR < 0.01 ~ 'up', 
      logFC <= -1 & FDR < 0.01 ~ 'down', 
      TRUE ~ 'not.signif'
    ))
  
}


D_DE_gene_res_re <- get_DE_res_re('Diploid') # 1168 significant
D_CIN_DE_gene_res_re <- get_DE_res_re('Diploid CIN') # 929 significant
T_DE_gene_res_re <- get_DE_res_re('Tetraploid') # 146 significant
T_CTH_DE_gene_res_re <- get_DE_res_re('Tetraploid CTH') # 219 significant

# subtype_DE_gene_res_re_list <- list(D_DE_gene_res_re, D_CIN_DE_gene_res_re, T_DE_gene_res_re, T_CTH_DE_gene_res_re)
# save(subtype_DE_gene_res_re_list, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_transcriptome_difference/subtype_DE_gene_res_re_list.Rdata')
# load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_transcriptome_difference/subtype_DE_gene_res_re_list.Rdata')


# D_DE_gene_res_re_signif <- subtype_DE_gene_res_re_list[[1]] %>% filter(status != 'not.signif')
# D_CIN_DE_gene_res_re_signif <- subtype_DE_gene_res_re_list[[2]] %>% filter(status != 'not.signif')
# T_DE_gene_res_re_signif <- subtype_DE_gene_res_re_list[[3]] %>% filter(status != 'not.signif')
# T_CTH_DE_gene_res_re_signif <- subtype_DE_gene_res_re_list[[4]] %>% filter(status != 'not.signif')
# 
# length(intersect(D_DE_gene_res_re_signif$Hugo_Symbol, D_CIN_DE_gene_res_re_signif$Hugo_Symbol)) # 760
# length(intersect(T_DE_gene_res_re_signif$Hugo_Symbol, T_CTH_DE_gene_res_re_signif$Hugo_Symbol)) # 38
# length(intersect(D_DE_gene_res_re_signif$Hugo_Symbol, T_DE_gene_res_re_signif$Hugo_Symbol)) # 87
# length(intersect(D_DE_gene_res_re_signif$Hugo_Symbol, T_CTH_DE_gene_res_re_signif$Hugo_Symbol)) # 171



# DE gene using Wilcoxon rank-sum test (e.g., Diploid vs other) ------------------------------------


get_DE_res_wilcox <- function(subtype_name){
  
  # subtype_name <- 'Diploid'
  
  sam_group_dat <- common_label %>% 
    mutate(label = ifelse(cluster %in% subtype_name, subtype_name, 'other'), 
           label = factor(label, levels = c('other', subtype_name))) %>% 
    arrange(label)
  
  test_data <- glioma_exp_longer %>% 
    inner_join(select(sam_group_dat, sample, label), by = 'sample') %>% 
    group_by(Hugo_Symbol) %>% 
    nest() %>% 
    ungroup()
  
  test_res <- test_data %>% 
    mutate(res = map(data, ~{
      
      sin_p_val <- wilcox.test(geneExp ~ label, data = .)$p.value
      sin_fc_data <- group_by(., label) %>% summarise(mean_Exp = mean(geneExp)) %>% pull(mean_Exp)
      sin_logfc <- log2(sin_fc_data[2]/sin_fc_data[1])
      
      tibble(logFC = sin_logfc, p_val = sin_p_val)
      
    })) %>% 
    select(-data) %>% 
    unnest(res) %>% 
    mutate(q_val = p.adjust(p_val, method = 'BH'), 
           status = case_when(
             logFC >= 1 & q_val < 0.01 ~ 'up', 
             logFC <= -1 & q_val < 0.01 ~ 'down', 
             TRUE ~ 'not.signif'
           ))
  
}

D_DE_gene_res_wilcox <- get_DE_res_wilcox('Diploid') # 1087 significant
D_CIN_DE_gene_res_wilcox <- get_DE_res_wilcox('Diploid CIN') # 868
T_DE_gene_res_wilcox <- get_DE_res_wilcox('Tetraploid') # 1
T_CTH_DE_gene_res_wilcox <- get_DE_res_wilcox('Tetraploid CTH') # 125









