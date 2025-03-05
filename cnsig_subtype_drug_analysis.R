
# common drugs in glioma therapy (chatGPT)
# Temozolomide, Cisplatin, Carboplatin, Cyclophosphamide, Methotrexate, Irinotecan, Etoposide, Camptothecin, Bevacizumab, Carmustine

library(tidyverse)
library(patchwork)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')

# drug response (limited samples) -----------------------------------------------------------

drug_gbm <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/nationwidechildrens.org_clinical_drug_gbm.txt')
drug_lgg <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/nationwidechildrens.org_clinical_drug_lgg.txt')

drug_glioma <- drug_gbm %>% 
  bind_rows(drug_lgg) %>% 
  mutate(response = ifelse(measure_of_response %in% c('[Discrepancy]', '[Not Applicable]', '[Not Available]', '[Unknown]'), NA, measure_of_response)) %>% 
  select(Sample = bcr_patient_barcode, drug = drug_name, type = therapy_type, response) %>% 
  drop_na() %>% 
  mutate(response = case_when(
    response %in% 'Clinical Progressive Disease' ~ 'PD', 
    response %in% 'Complete Response' ~ 'CR', 
    response %in% 'Partial Response' ~ 'PR', 
    response %in% 'Stable Disease' ~ 'SD'
  ), response_binary = ifelse(response %in% c('PD', 'SD'), 'SD/PD', 'CR/PR'))

# common drugs in glioma therapy(from drug data)

glioma_drug_name <- union(drug_glioma$drug, c('Temozolomide', 'Cisplatin', 'Carboplatin', 'Cyclophosphamide', 'Methotrexate', 'Irinotecan', 'Etoposide', 'Camptothecin', 'Bevacizumab'))

save(glioma_drug_name, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_drug_name.Rdata')


# top drugs with sufficient samples

drug_sam_stat <- drug_glioma %>% distinct(Sample, drug) %>% count(drug) %>% top_n(5, n)
# drug             n
# <chr>        <int>
# Avastin         20
# CCNU            10
# Procarbazine    11
# Temodar         67
# Temozolomide    62

# Temozolomide or Temodar(Temozolomide的商品名)

tmz_response_dat <- drug_glioma %>% 
  filter(drug %in% c('Temodar', 'Temozolomide')) %>% 
  mutate(label = ifelse(response_binary %in% 'SD/PD', 0, 1)) %>% 
  select(Sample, label) %>% 
  group_by(Sample) %>% 
  summarise(new_label = sum(label)) %>% 
  mutate(response_binary = ifelse(new_label > 0, 'CR/PR', 'SD/PD')) %>% 
  select(-new_label)
# response_binary     n
# <chr>           <int>
# CR/PR              22
# SD/PD             107

# overall

drug_response_dat <- drug_glioma %>% 
  mutate(label = ifelse(response_binary %in% 'SD/PD', 0, 1)) %>% 
  select(Sample, label) %>% 
  group_by(Sample) %>% 
  summarise(new_label = sum(label)) %>% 
  mutate(response_binary = ifelse(new_label > 0, 'CR/PR', 'SD/PD')) %>% 
  select(-new_label)
# response_binary     n
# <chr>           <int>
# CR/PR              29
# SD/PD             127

save(tmz_response_dat, drug_response_dat, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_drug_response_dat.Rdata')



# comparision among cnsig subtypes

cnsig_subtype_tmz_response <- tmz_response_dat %>% inner_join(select(glioma_cnsig_subtype_label, Sample = patient, cluster), by = 'Sample')
cnsig_subtype_tmz_response_p <- fisher.test(cnsig_subtype_tmz_response$response_binary, cnsig_subtype_tmz_response$cluster)$p.value
# p = 0.79


cnsig_subtype_drug_response <- drug_response_dat %>% inner_join(select(glioma_cnsig_subtype_label, Sample = patient, cluster), by = 'Sample')
cnsig_subtype_drug_response_p <- fisher.test(cnsig_subtype_drug_response$response_binary, cnsig_subtype_drug_response$cluster)$p.value
# p = 0.44


# radiation response (limited samples) ------------------------------------------------------

radiation_gbm <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/nationwidechildrens.org_clinical_radiation_gbm.txt')
radiation_lgg <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/nationwidechildrens.org_clinical_radiation_lgg.txt')

radiation_glioma <- radiation_gbm %>% bind_rows(radiation_lgg)

radiation_response_dat <- radiation_glioma %>% 
  mutate(response = ifelse(measure_of_response %in% c('[Not Applicable]', '[Not Available]', '[Unknown]'), NA, measure_of_response)) %>% 
  drop_na(response) %>% 
  mutate(response_binary = ifelse(response %in% c('Complete Response', 'Partial Response'), 'CR/PR', 'SD/PD')) %>% 
  select(Sample = bcr_patient_barcode, response_binary) %>% 
  mutate(label = ifelse(response_binary %in% 'SD/PD', 0, 1)) %>% 
  select(Sample, label) %>% 
  group_by(Sample) %>% 
  summarise(new_label = sum(label)) %>% 
  mutate(response_binary = ifelse(new_label > 0, 'CR/PR', 'SD/PD')) %>% 
  select(-new_label)
# response_binary     n
# <chr>           <int>
# CR/PR              47
# SD/PD             133

save(radiation_response_dat, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/radiation_response_dat.Rdata')

cnsig_subtype_radiation_response <- radiation_response_dat %>% inner_join(select(glioma_cnsig_subtype_label, Sample = patient, cluster), by = 'Sample')
cnsig_subtype_radiation_response_p <- fisher.test(cnsig_subtype_radiation_response$response_binary, cnsig_subtype_radiation_response$cluster)$p.value
# p = 0.65


# drug sensitivity prediction ---------------------------------------------

library(oncoPredict)
library(sva)
library(preprocessCore)
library(ridge)
library(glmnet)
library(impute)

# 基于GDSC预测TCGA胶质瘤样本的药物敏感性

glioma_exp_data <- vroom::vroom('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_exp_data.tsv') %>% 
  mutate_if(is.numeric, ~log2(.+1))

glioma_exp_mat <- glioma_exp_data %>% column_to_rownames('Hugo_Symbol')

# GDSC表达谱和药物敏感数据(训练)

GDSC_cell_line_cancer_type <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/GDSC_cell_line_cancer_type.tsv')
GDSC_cell_line_expr_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/GDSC_Cell_line_RMA_proc_basalExp.txt')
GDSC_cell_line_logIC50 <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/GDSC_cell_lines_drug_logIC50.tsv') # log(IC50), log(exp(1))

glioma_identifier <- GDSC_cell_line_cancer_type %>% filter(`GDSC_Tissue_descriptor 2` %in% 'glioma') %>% pull(COSMIC_identifier)
# 53 glioma cell lines

GDSC_cell_line_expr_mat <- GDSC_cell_line_expr_data %>% 
  select(Hugo_Symbol = GENE_SYMBOLS, starts_with('DATA')) %>% 
  set_names(~str_remove(., 'DATA\\.')) %>% 
  drop_na(Hugo_Symbol) %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE) %>% 
  column_to_rownames('Hugo_Symbol')
# 17419 genes, 1018 cell lines

GDSC_cell_line_IC50_mat <- GDSC_cell_line_logIC50 %>% 
  filter(`Cell_line cosmic_identifiers` %in% glioma_identifier) %>% 
  select(-Sample_Names) %>% 
  column_to_rownames('Cell_line cosmic_identifiers') %>% 
  mutate_all(~exp(.))
# 行为细胞系，列为药物

GDSC_cell_line_IC50_mat <- GDSC_cell_line_IC50_mat[, apply(GDSC_cell_line_IC50_mat, 2, function(x) {mean(is.na(x)) < 0.2})] # 去除在超过20%样本中都缺失的药物
GDSC_cell_line_IC50_mat <- GDSC_cell_line_IC50_mat[apply(GDSC_cell_line_IC50_mat, 1, function(x) {mean(is.na(x)) < 0.2}), ] # 去除在超过20%药物中都缺失的样本
GDSC_cell_line_IC50_mat <- t(impute.knn(t(GDSC_cell_line_IC50_mat))$data) # KNN填补缺失值
# 46 cell lines, 216 drugs

setwd('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_drug_analysis')

calcPhenotype(trainingExprData = as.matrix(GDSC_cell_line_expr_mat), # 函数会自己匹配细胞系表达和药敏的共有样本
              trainingPtype = GDSC_cell_line_IC50_mat, # 但如果要选用特殊的细胞系类型，比如“生殖系统”、“消化系统”等，需要提前对细胞系做预筛选，函数不提供该操作
              testExprData = as.matrix(glioma_exp_mat),
              batchCorrect = 'standardize',
              powerTransformPhenotype = FALSE, # 是否进行对数转化，确保原始数据已经经过log转化后就为FALSE，否则默认为TRUE
              removeLowVaryingGenes = 0.2,
              removeLowVaringGenesFrom = 'homogenizeData', # 在过滤基因时候采用的数据类型选项：'homogenizeData' （批次效应消除后的数据）或 'rawData'（原始数据）
              minNumSamples = 10, # 确定最小训练样本数目（一般不存在这个问题因为GDSC，CTRP，PRISM都有很多细胞系）
              selection = 1, # 确定如何处理重复基因选项： -1是询问用户, 1 取均值, 以及 2 移除重复
              printOutput = TRUE,
              pcr = FALSE, # 确定是否采用PCA对基因表达进行降维选项：默认不进行
              report_pc = FALSE,
              percent = 80, # 确定在pcr=TRUE时，主成份分析需要达到的解释度，默认80，即所得到的主成份的解释度相加要大于80%
              cc = FALSE, # 确定是否需要列出相关性结果来确定潜在的药物biomarker
              rsq = FALSE)

# batchCorrect选项: "eb"对应ComBat, "qn"对应quantiles normalization, "standardize", 以及 "none"
# "eb" 在使用微阵列数据做训练集（细胞系），并且应用到微阵列测试集（临床样本）时比较有效
# "standardize" 在使用微阵列数据做训练集（细胞系），但应用在RNA-seq测试集（临床样本）时比较有效


# comparison of predicted drug sensitivity 

drug_pred <- read_csv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_drug_analysis/calcPhenotype_Output/DrugPredictions.csv')

drug_pred_test_data <- drug_pred %>% 
  pivot_longer(-sample, names_to = 'drug', values_to = 'IC50') %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample') %>% 
  mutate(logIC50 = log(IC50))

drug_pred_test_res <- drug_pred_test_data %>% 
  group_by(drug) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(p_val = map_dbl(data, possibly(~kruskal.test(IC50 ~ cluster, data = .)$p.value, NA)), 
         q_val = p.adjust(p_val, method = 'BH')) %>% 
  select(-data)

save(drug_pred_test_res, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_drug_analysis/drug_pred_test_res.Rdata')

# load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_drug_name.Rdata') # glioma_drug_name
# glioma_common_drug_res <- drug_pred_test_res %>% filter(drug %in% glioma_drug_name)


# function: violinplot

drug_violin <- function(drug_name){
  
  # drug_name <- '5-Fluorouracil'
  
  sin_violin <- drug_pred_test_data %>% 
    filter(drug %in% drug_name, !is.na(logIC50)) %>% 
    ggplot(aes(cluster, logIC50, fill = cluster)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
    ggpubr::geom_pwc(hide.ns = TRUE) + 
    ggpubr::stat_compare_means() +
    labs(x = NULL, y = 'predicted logIC50', title = drug_name) +
    # ggsci::scale_fill_jco() +
    scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
    theme_classic() +
    theme(legend.position = 'none', 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1), 
          plot.title = element_text(hjust = 0.5))
  
}

# final drugs (FDA approved for glioma treatment)(https://www.cancer.gov/about-cancer/treatment/drugs/brain)

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_drug_analysis/selected_drug_logIC50_violinplot.pdf', width = 12, height = 4)

drug_violin('Temozolomide') + drug_violin('Dabrafenib') + drug_violin('Trametinib') + plot_annotation(tag_levels = 'A')

dev.off()



# TIDE score --------------------------------------------------------------


gbm_tide_score <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/TCGA.GBM.RNASeq.norm_subtract.OS_base') %>% 
  mutate(TIDE_score = ifelse(Exclusion >= Dysfunction, Exclusion, Dysfunction)) %>% 
  select(Sample, TIDE_score, Exclusion, Dysfunction)

lgg_tide_score <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/TCGA.LGG.RNASeq.OS_base') %>% 
  mutate(TIDE_score = ifelse(Exclusion >= Dysfunction, Exclusion, Dysfunction)) %>% 
  select(Sample, TIDE_score, Exclusion, Dysfunction)

glioma_tide_score <- gbm_tide_score %>% bind_rows(lgg_tide_score)

save(glioma_tide_score, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_tide_score.Rdata')

cnsig_subtype_tide_score <- glioma_tide_score %>% inner_join(select(glioma_cnsig_subtype_label, Sample = patient, cluster), by = 'Sample')

# TIDE score

tide_score_plot <- cnsig_subtype_tide_score %>% 
  ggplot(aes(cluster, TIDE_score, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
  ggpubr::geom_pwc(hide.ns = FALSE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'TIDE score') +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

# Exclusion score

exclusion_score_plot <- cnsig_subtype_tide_score %>% 
  ggplot(aes(cluster, Exclusion, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
  ggpubr::geom_pwc(hide.ns = FALSE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'Exclusion score') +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

# Dysfunction score

dysfunction_score_plot <- cnsig_subtype_tide_score %>% 
  ggplot(aes(cluster, Dysfunction, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
  ggpubr::geom_pwc(hide.ns = FALSE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'Dysfunction score') +
  scale_fill_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_drug_analysis/cnsig_subtype_tide_score_boxplot.pdf', width = 12, height = 4)

tide_score_plot + exclusion_score_plot + dysfunction_score_plot

dev.off()





