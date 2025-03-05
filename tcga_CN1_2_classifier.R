
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

tcga_cli <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv') %>% 
  select(bcr_patient_barcode, type, age = age_at_initial_pathologic_diagnosis, gender, race, OS:DSS.time, PFI, PFI.time)

tcga_pancan_CNsig21_attributions <- read_tsv('/boot3/bio_liaojl/TCGA_CN_signature_transcriptome_work/Data/Original/tcga_pancan_CNsig21_attributions.tsv')


cnsig_sur_data <- tcga_pancan_CNsig21_attributions %>% 
  mutate(Sample = str_sub(Sample, 1, 12)) %>% 
  select(Sample, CN1, CN2) %>% 
  inner_join(select(tcga_cli, bcr_patient_barcode, type, OS, OS.time), by = c('Sample' = 'bcr_patient_barcode')) %>% 
  group_by(type) %>% 
  mutate(CN1_ntile = ntile(CN1, 3), 
         CN2_ntile = ntile(CN2, 3), 
         CN1_class = case_when(
           CN1_ntile == 1 ~ 'CN1L', 
           CN1_ntile == 2 ~ 'CN1M', 
           CN1_ntile == 3 ~ 'CN1H'
         ), 
         CN2_class = case_when(
           CN2_ntile == 1 ~ 'CN2L', 
           CN2_ntile == 2 ~ 'CN2M', 
           CN2_ntile == 3 ~ 'CN2H'
         ), 
         class = str_c(CN1_class, CN2_class, sep = '-')) %>% 
  ungroup()


cnsig_sur_data_med <- tcga_pancan_CNsig21_attributions %>% 
  mutate(Sample = str_sub(Sample, 1, 12)) %>% 
  select(Sample, CN1, CN2) %>% 
  inner_join(select(tcga_cli, bcr_patient_barcode, type, OS, OS.time), by = c('Sample' = 'bcr_patient_barcode')) %>% 
  group_by(type) %>% 
  mutate(CN1_class = ifelse(CN1 > median(CN1), 'CN1H', 'CN1L'), 
         CN2_class = ifelse(CN2 > median(CN2), 'CN2H', 'CN2L'), 
         class = str_c(CN1_class, CN2_class, sep = '-')) %>% 
  ungroup()




# survival analysis -------------------------------------------------------

# combine CN1 and CN2 tiles

tcga_cnsig_sur_res <- cnsig_sur_data %>% 
  select(type, class, OS, OS.time) %>% 
  group_by(type) %>% 
  nest() %>% 
  mutate(kmplots = map2(type, data, ~kmplot_fun(.x, .y)))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_cn1_2_kmplots.pdf', width = 25, height = 18)

wrap_plots(tcga_cnsig_sur_res$kmplots)

dev.off()


# combine CN1 and CN2 median

tcga_cnsig_sur_res_med <- cnsig_sur_data_med %>% 
  select(type, class, OS, OS.time) %>% 
  group_by(type) %>% 
  nest() %>% 
  ungroup()
  mutate(kmplots = map2(type, data, ~kmplot_fun(.x, .y)))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_cn1_2_kmplots_med.pdf', width = 25, height = 18)

wrap_plots(tcga_cnsig_sur_res_med$kmplots)

dev.off()


# single CN1 3tiles

tcga_cnsig_sur_CN1_res <- cnsig_sur_data %>% 
  select(type, class = CN1_class, OS, OS.time) %>% 
  group_by(type) %>% 
  nest() %>% 
  mutate(kmplots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_cn1_kmplots.pdf', width = 25, height = 18)

wrap_plots(tcga_cnsig_sur_CN1_res$kmplots)

dev.off()

# single CN2 3tiles

tcga_cnsig_sur_CN2_res <- cnsig_sur_data %>% 
  select(type, class = CN2_class, OS, OS.time) %>% 
  group_by(type) %>% 
  nest() %>% 
  mutate(kmplots = map2(type, data, ~kmplot_fun(.x, .y)))

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/tcga_cn2_kmplots.pdf', width = 25, height = 18)

wrap_plots(tcga_cnsig_sur_CN2_res$kmplots)

dev.off()





kmplot_fun <- function(ttype, sur_dat){
  
  # ttype <- 'ACC'
  # sur_dat <- tcga_cnsig_sur_res$data[[1]]
  
  # pdf('/boot3/bio_liaojl/temp.pdf')
  
  ggsurvplot(surv_fit(Surv(OS.time, OS) ~ class, data = sur_dat),
             pval = TRUE,
             xlab = 'Survival time (Months)',
             # legend.labs = c('Low', 'High'),
             # palette = c('#5082AF', '#C9372E'),
             risk.table = TRUE,
             title = ttype, 
             ggtheme = theme_bw(),
             data = sur_dat)$plot

  # dev.off()
  
}












