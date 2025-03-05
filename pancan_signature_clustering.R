
library(tidyverse)
library(survival)
library(survminer)
library(ConsensusClusterPlus)
library(pheatmap)
library(circlize)
library(patchwork)

tcga_ttype_patient <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv') %>% 
  mutate(patient_code = str_sub(bcr_patient_barcode, 9, 12)) %>% 
  distinct(type, bcr_patient_barcode, patient_code)
# 11,160 patients

load('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/MutSig_comprehensive_analysis/Processed/comb_cli_data.Rdata')

tcga_cli_data <- comb_cli_data %>% 
  inner_join(tcga_ttype_patient, by = c('Ori_Cancer_Type' = 'type', 'Sample' = 'patient_code')) %>% 
  rename(type = Ori_Cancer_Type) %>% 
  select(type, bcr_patient_barcode, age:grade, OS, OS.time, purity:TMB)
# 11,160 patients

tcga_pancan_CNsig21_attributions <- read_tsv('/boot3/bio_liaojl/TCGA_CN_signature_transcriptome_work/Data/Original/tcga_pancan_CNsig21_attributions.tsv') %>% 
  mutate(Sample = str_sub(Sample, 1, 12)) %>% # important, some samples named as 'TCGA-02-0001-a'
  distinct(Sample, .keep_all = TRUE)
# 9,699 samples

tcga_cx_sig <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Original/tcga_cx_sig.tsv') %>% 
  mutate(Patient = str_sub(Patient, 1, 12)) %>% 
  distinct(Patient, .keep_all = TRUE)
# 6,335 samples

tcga_sbs_sig <- read_csv('/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Original/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv') %>% 
  mutate(patient = str_sub(`Sample Names`, 1, 12)) %>% 
  select(patient, starts_with('SBS')) %>% 
  distinct(patient, .keep_all = TRUE) %>% 
  pivot_longer(-patient, names_to = 'signature', values_to = 'value') %>% 
  group_by(patient) %>% 
  filter(sum(value) > 0) %>% 
  mutate(value = value/sum(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = signature, values_from = value)

# SBS and CN signature

comb_sbs_cn_data <- tcga_pancan_CNsig21_attributions %>% 
  inner_join(tcga_sbs_sig, by = c('Sample' = 'patient')) %>% 
  inner_join(select(tcga_cli_data, type, Sample = bcr_patient_barcode), by = 'Sample') %>% 
  rename(patient = Sample) %>% 
  select(type, patient, everything())
# 8,208 patients

# SBS and CX signature

comb_sbs_cx_data <- tcga_cx_sig %>% 
  rename(patient = Patient) %>% 
  inner_join(tcga_sbs_sig, by = 'patient') %>% 
  inner_join(select(tcga_cli_data, type, patient = bcr_patient_barcode), by = 'patient') %>% 
  select(type, patient, everything())
# 5,236 patients

# only CN signature

cn_alone_data <- tcga_pancan_CNsig21_attributions %>% 
  inner_join(select(tcga_cli_data, type, Sample = bcr_patient_barcode), by = 'Sample') %>% 
  rename(patient = Sample) %>% 
  select(type, patient, everything())
# 9,634 patients

# only CX signature

cx_alone_data <- tcga_cx_sig %>% 
  rename(patient = Patient) %>% 
  inner_join(select(tcga_cli_data, type, patient = bcr_patient_barcode), by = 'patient') %>% 
  select(type, patient, everything())
# 6,297 patients


# clustering

comb_sbs_cn_clus_res <- ttype_sig_clustering(comb_sbs_cn_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/comb_sbs_cn')
comb_sbs_cx_clus_res <- ttype_sig_clustering(comb_sbs_cx_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/comb_sbs_cx')
cn_alone_clus_res <- ttype_sig_clustering(cn_alone_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cn_alone')
cx_alone_clus_res <- ttype_sig_clustering(cx_alone_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cx_alone')


# functions ---------------------------------------------------------------

# consensus clustering based on signature data

ttype_sig_clustering <- function(pancan_sig_data, dir_path){
  
  # pancan_sig_data <- comb_sbs_cn_data
  # dir_path <- '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/comb_sbs_cn'
  
  # matrix
  
  expo_frac_mat <- pancan_sig_data %>% 
    group_by(type) %>% 
    nest() %>% 
    ungroup() %>% 
    mutate(mat_data = map(data, ~{
      select_if(., ~mean(. > 0) > 0.05) %>% column_to_rownames('patient') %>% t()
    })) %>% 
    select(-data)

  # create directory
  
  for(i in expo_frac_mat$type){
    dir.create(str_c(dir_path, '/', i), showWarnings = FALSE)
  } # dir.create() does not crash if the directory already exists, it just prints out a warning
  
  # consensus clustering
  
  clustering_res <- expo_frac_mat %>% 
    mutate(path = str_c(dir_path, '/', type), 
           clus_res = map2(mat_data, path, possibly(~ConsensusClusterPlus(.x, maxK = 12, distance = 'euclidean', seed = 123, title = .y, plot = 'pdf'), NULL)))
  
}

# extract significant cluster cox cofficients from summary data

get_signif_coff <- function(cox_summary_res){
  
  # cox_summary_res <- comb_sbs_cn_cox
  
  # extract cox efficient information
  
  cox_efficients <- cox_summary_res %>% 
    select(-log_rank_p) %>% 
    pivot_longer(-type, names_to = 'covar_comb', values_to = 'summary_res') %>% 
    filter(!map_dbl(summary_res, ~is.null(.))) %>%  # remove null rows
    mutate(efficient_res = map(summary_res, ~{
      tibble(variable = names(.$coefficients[, 1]),
             coef = .$coefficients[, 1],
             HR = .$coefficients[, 2],
             HR_confint_lower = .$conf.int[, 3],
             HR_confint_upper = .$conf.int[, 4],
             p_val = .$coefficients[, 5])
    })) %>% 
    select(-summary_res) %>% 
    unnest(cols = efficient_res) %>% 
    mutate(covar_comb = str_replace_all(str_remove(covar_comb, '_res'), '_', '+'))
  
  # significant cox result of cluster
  
  cox_signif <- cox_efficients %>% 
    mutate(cluster = str_detect(variable, 'cluster'), 
           signif = ifelse(p_val < 0.05, TRUE, FALSE)) %>% 
    group_by(type, covar_comb) %>% 
    filter(any(cluster == TRUE & signif == TRUE)) %>% 
    ungroup() %>% 
    select(type:p_val)
  
}

# km-plots plotting

kmplot_fun <- function(ttype, sur_dat){
  
  # ttype <- 'ACC'
  # sur_dat <- tcga_cnsig_sur_res$data[[1]]
  
  ggsurvplot(surv_fit(Surv(OS.time, OS) ~ cluster, data = sur_dat),
             pval = TRUE,
             xlab = 'Survival time (Months)',
             # legend.labs = c('Low', 'High'),
             # palette = c('#5082AF', '#C9372E'),
             risk.table = TRUE,
             title = ttype, 
             ggtheme = theme_bw(),
             data = sur_dat)$plot
  
}



# manually select optimal number of clusters per cancer type --------------

optimal_cluster_num <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/optimal_cluster_num.tsv')

# merge cluster and clinical data -----------------------------------------

# SBS and CN signature

comb_sbs_cn_clus_cli <- comb_sbs_cn_clus_res %>% 
  inner_join(optimal_cluster_num, by = 'type') %>% 
  select(type, clus_res, num = comb_sbs_cn_num) %>% 
  drop_na() %>% 
  mutate(optimal_res = map2(clus_res, num, ~{
    tibble(patient = names(.x[[.y]]$consensusClass), cluster = as.character(.x[[.y]]$consensusClass))
  })) %>% 
  select(type, optimal_res) %>% 
  unnest(cols = c(optimal_res)) %>% 
  inner_join(tcga_cli_data, by = c('patient' = 'bcr_patient_barcode', 'type'))

save(comb_sbs_cn_clus_cli, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/comb_sbs_cn_clus_cli.Rdata')

# SBS and CX signature

comb_sbs_cx_clus_cli <- comb_sbs_cx_clus_res %>% 
  inner_join(optimal_cluster_num, by = 'type') %>% 
  select(type, clus_res, num = comb_sbs_cx_num) %>% 
  drop_na() %>% 
  mutate(optimal_res = map2(clus_res, num, ~{
    tibble(patient = names(.x[[.y]]$consensusClass), cluster = as.character(.x[[.y]]$consensusClass))
  })) %>% 
  select(type, optimal_res) %>% 
  unnest(cols = c(optimal_res)) %>% 
  inner_join(tcga_cli_data, by = c('patient' = 'bcr_patient_barcode', 'type'))

save(comb_sbs_cx_clus_cli, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/comb_sbs_cx_clus_cli.Rdata')

# CN signature alone

cn_alone_clus_cli <- cn_alone_clus_res %>% 
  inner_join(optimal_cluster_num, by = 'type') %>% 
  select(type, clus_res, num = cn_alone_num) %>% 
  drop_na() %>% 
  mutate(optimal_res = map2(clus_res, num, ~{
    tibble(patient = names(.x[[.y]]$consensusClass), cluster = as.character(.x[[.y]]$consensusClass))
  })) %>% 
  select(type, optimal_res) %>% 
  unnest(cols = c(optimal_res)) %>% 
  inner_join(tcga_cli_data, by = c('patient' = 'bcr_patient_barcode', 'type'))

save(cn_alone_clus_cli, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/cn_alone_clus_cli.Rdata')

# CX signature alone

cx_alone_clus_cli <- cx_alone_clus_res %>% 
  inner_join(optimal_cluster_num, by = 'type') %>% 
  select(type, clus_res, num = cx_alone_num) %>% 
  drop_na() %>% 
  mutate(optimal_res = map2(clus_res, num, ~{
    tibble(patient = names(.x[[.y]]$consensusClass), cluster = as.character(.x[[.y]]$consensusClass))
  })) %>% 
  select(type, optimal_res) %>% 
  unnest(cols = c(optimal_res)) %>% 
  inner_join(tcga_cli_data, by = c('patient' = 'bcr_patient_barcode', 'type'))

save(cx_alone_clus_cli, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/cx_alone_clus_cli.Rdata')

# survival analysis ------------------------------------------------------------

### SBS and CN signature

comb_sbs_cn_cox <- comb_sbs_cn_clus_cli %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(log_rank_p = map_dbl(data, ~survdiff(Surv(OS.time, OS) ~ cluster, data = .)$pvalue)) %>% 
  filter(log_rank_p < 0.05) %>% # log-rank test must significant, 9 cancer types
  mutate(age_sex_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + grade, data = .)), NULL)), 
         age_sex_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + stage, data = .)), NULL)), 
         age_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + grade, data = .)), NULL)), 
         age_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + stage, data = .)), NULL)), 
         age_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age, data = .)), NULL))) %>% 
  select(-data)

# save(comb_sbs_cn_cox, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/comb_sbs_cn_cox.Rdata')

comb_sbs_cn_cox_signif <- get_signif_coff(comb_sbs_cn_cox)
# 8 cancer types

comb_sbs_cn_cox_signif_clus_num <- comb_sbs_cn_clus_cli %>% 
  semi_join(comb_sbs_cn_cox_signif, by = 'type') %>% 
  count(type, cluster)

# write_tsv(comb_sbs_cn_cox_signif, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cn_cox_signif.tsv')

# km-plots

comb_sbs_cn_kmplot <- comb_sbs_cn_clus_cli %>% 
  semi_join(comb_sbs_cn_cox_signif, by = 'type') %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cn_signif_kmplots.pdf', width = 16, height = 12)

wrap_plots(comb_sbs_cn_kmplot$km_plots)

dev.off()


# km-plots (remove cluster containing less than 10 patients)

comb_sbs_cn_kmplot_p10 <- comb_sbs_cn_clus_cli %>% 
  semi_join(comb_sbs_cn_cox_signif, by = 'type') %>% 
  group_by(type, cluster) %>% 
  filter(n() >= 10) %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cn_signif_kmplots_p10.pdf', width = 16, height = 12)

wrap_plots(comb_sbs_cn_kmplot_p10$km_plots)

dev.off()


### SBS and CX signature

comb_sbs_cx_cox <- comb_sbs_cx_clus_cli %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(log_rank_p = map_dbl(data, ~survdiff(Surv(OS.time, OS) ~ cluster, data = .)$pvalue)) %>% 
  filter(log_rank_p < 0.05) %>% # 12 cancer types
  mutate(age_sex_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + grade, data = .)), NULL)), 
         age_sex_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + stage, data = .)), NULL)), 
         age_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + grade, data = .)), NULL)), 
         age_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + stage, data = .)), NULL)), 
         age_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age, data = .)), NULL))) %>% 
  select(-data)

# save(comb_sbs_cx_cox, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/comb_sbs_cx_cox.Rdata')

comb_sbs_cx_cox_signif <- get_signif_coff(comb_sbs_cx_cox)

write_tsv(comb_sbs_cx_cox_signif, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cx_cox_signif.tsv')

# km-plots

comb_sbs_cx_kmplot <- comb_sbs_cx_clus_cli %>% 
  semi_join(comb_sbs_cx_cox_signif, by = 'type') %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cx_signif_kmplots.pdf', width = 16, height = 12)

wrap_plots(comb_sbs_cx_kmplot$km_plots)

dev.off()

# km-plots (remove cluster containing less than 10 patients)

comb_sbs_cx_kmplot_p10 <- comb_sbs_cx_clus_cli %>% 
  semi_join(comb_sbs_cx_cox_signif, by = 'type') %>% 
  group_by(type, cluster) %>% 
  filter(n() >= 10) %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/comb_sbs_cx_signif_kmplots_p10.pdf', width = 16, height = 12)

wrap_plots(comb_sbs_cx_kmplot_p10$km_plots)

dev.off()


# CN signature alone

cn_alone_cox <- cn_alone_clus_cli %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(log_rank_p = map_dbl(data, ~survdiff(Surv(OS.time, OS) ~ cluster, data = .)$pvalue)) %>% 
  filter(log_rank_p < 0.05) %>% # 10 cancer types
  mutate(age_sex_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + grade, data = .)), NULL)), 
         age_sex_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + stage, data = .)), NULL)), 
         age_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + grade, data = .)), NULL)), 
         age_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + stage, data = .)), NULL)), 
         age_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age, data = .)), NULL))) %>% 
  select(-data)

save(cn_alone_cox, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/cn_alone_cox.Rdata')

cn_alone_cox_signif <- get_signif_coff(cn_alone_cox)

write_tsv(cn_alone_cox_signif, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cn_alone_cox_signif.tsv')

# km-plots

cn_alone_kmplot <- cn_alone_clus_cli %>% 
  semi_join(cn_alone_cox_signif, by = 'type') %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cn_alone_signif_kmplots.pdf', width = 16, height = 12)

wrap_plots(cn_alone_kmplot$km_plots)

dev.off()

# km-plots (remove cluster containing less than 10 patients)

cn_alone_kmplot_p10 <- cn_alone_clus_cli %>% 
  semi_join(cn_alone_cox_signif, by = 'type') %>% 
  group_by(type, cluster) %>% 
  filter(n() >= 10) %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cn_alone_signif_kmplots_p10.pdf', width = 16, height = 12)

wrap_plots(cn_alone_kmplot_p10$km_plots)

dev.off()


# CX signature alone

cx_alone_cox <- cx_alone_clus_cli %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(log_rank_p = map_dbl(data, ~survdiff(Surv(OS.time, OS) ~ cluster, data = .)$pvalue)) %>% 
  filter(log_rank_p < 0.05) %>% # 13 cancer types
  mutate(age_sex_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + grade, data = .)), NULL)), 
         age_sex_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + stage, data = .)), NULL)), 
         age_grade_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + grade, data = .)), NULL)), 
         age_stage_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age + stage, data = .)), NULL)), 
         age_res = map(data, possibly(~summary(coxph(Surv(OS.time, OS) ~ cluster + age, data = .)), NULL))) %>% 
  select(-data)

save(cx_alone_cox, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Data/Processed/signature_clustering/cx_alone_cox.Rdata')

cx_alone_cox_signif <- get_signif_coff(cx_alone_cox)

write_tsv(cx_alone_cox_signif, '/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cx_alone_cox_signif.tsv')

# km-plots

cx_alone_kmplot <- cx_alone_clus_cli %>% 
  semi_join(cx_alone_cox_signif, by = 'type') %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cx_alone_signif_kmplots.pdf', width = 16, height = 12)

wrap_plots(cx_alone_kmplot$km_plots)

dev.off()


# km-plots (remove cluster containing less than 10 patients)

cx_alone_kmplot_p10 <- cx_alone_clus_cli %>% 
  semi_join(cx_alone_cox_signif, by = 'type') %>% 
  group_by(type, cluster) %>% 
  filter(n() >= 10) %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup() %>% 
  mutate(km_plots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_sig_clustering/cluster_cox_result/cx_alone_signif_kmplots_p10.pdf', width = 16, height = 12)

wrap_plots(cx_alone_kmplot_p10$km_plots)

dev.off()






