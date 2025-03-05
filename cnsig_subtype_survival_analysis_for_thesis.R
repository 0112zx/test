
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_weight.Rdata')
load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')
glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')


# survival difference among subtypes --------------------------------------


subtype_sur_data <- glioma_cli_data %>% 
  select(sample, age, gender, grade, IDH_codel_subtype, OS:PFI.time) %>% 
  inner_join(select(glioma_cnsig_subtype_label, sample, cluster), by = 'sample') %>% 
  rename(group = cluster) %>% 
  mutate(group = factor(group, levels = c('Diploid', 'Tetraploid', 'Tetraploid CTH', 'Diploid CIN')))

# pairwise comparison of subgroups survival

# pairwise_survdiff(Surv(OS.time, OS) ~ group, data = subtype_sur_data, p.adjust.method = "none") # OS
# pairwise_survdiff(Surv(DSS.time, DSS) ~ group, data = subtype_sur_data, p.adjust.method = "none") # DSS
# pairwise_survdiff(Surv(PFI.time, PFI) ~ group, data = subtype_sur_data, p.adjust.method = "none") # PFI

# median survival time

# surv_median(surv_fit(Surv(OS.time, OS) ~ group, data = subtype_sur_data)) # OS
# surv_median(surv_fit(Surv(DSS.time, DSS) ~ group, data = subtype_sur_data)) # DSS
# surv_median(surv_fit(Surv(PFI.time, PFI) ~ group, data = subtype_sur_data)) # PFI

# kmplots

subtype_os_kmplot <- kmplot_fun('OS', subtype_sur_data, c('Diploid', 'Tetraploid', 'Tetraploid CTH', 'Diploid CIN'), c("#D7191C", "#ABDDA4", "#2B83BA", "#FDAE61"))
subtype_dss_kmplot <- kmplot_fun('DSS', subtype_sur_data, c('Diploid', 'Tetraploid', 'Tetraploid CTH', 'Diploid CIN'), c("#D7191C", "#ABDDA4", "#2B83BA", "#FDAE61"))
subtype_pfi_kmplot <- kmplot_fun('PFI', subtype_sur_data, c('Diploid', 'Tetraploid', 'Tetraploid CTH', 'Diploid CIN'), c("#D7191C", "#ABDDA4", "#2B83BA", "#FDAE61"))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_prognosis_for_thesis/subtype_kmplots.pdf', width = 12, height = 4)

subtype_os_kmplot + subtype_dss_kmplot + subtype_pfi_kmplot

dev.off()


# functions

kmplot_fun <- function(type, sur_dat, curve_labs, curve_colors, title = NULL){
  
  # type <- 'OS'
  # sur_dat <- subtype_sur_data
  # curve_labs <- c('Diploid', 'Tetraploid', 'Tetraploid CTH', 'Diploid CIN')
  # curve_colors <- c("#D7191C", "#ABDDA4", "#2B83BA", "#FDAE61")
  
  fit_formu <- as.formula(str_c('Surv(', type, '.time, ', type, ') ~ group'))
  
  fit_res <- surv_fit(fit_formu, data = sur_dat)
  
  def_plot <- ggsurvplot(fit_res, pval = TRUE, legend.labs = curve_labs, palette = curve_colors, ggtheme = theme_bw())
  
  final_plot <- def_plot$plot +
    labs(x = 'Survival time (Days)', y = str_c(type, ' probability'), col = NULL, title = title) +
    theme(panel.grid = element_blank(), 
          legend.position = c(1, 1), 
          legend.justification = c(1, 1), 
          legend.background = element_rect(fill = NA), 
          legend.key = element_rect(fill = NA), 
          plot.title = element_text(hjust = 0.5))
  
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_prognosis_for_thesis/subtype_os_kmplot.pdf')
  # 
  # final_plot
  # 
  # dev.off()
  
}


# prognosis potential of individual copy number signatures --------------------------------------


comb_cnsig_sur_data <- glioma_cnsig_weight %>% 
  inner_join(select(glioma_cli_data, bcr_patient_barcode, age, gender, grade, IDH_codel_subtype, OS, OS.time), by = c('sample' = 'bcr_patient_barcode')) %>% 
  group_by(CNSig) %>% 
  mutate(group = ifelse(weight > median(weight), 'High', 'Low'), 
         group = factor(group, levels = c('Low', 'High'))) %>% 
  nest() %>% 
  ungroup()


comb_cnsig_kmplots <- comb_cnsig_sur_data %>% 
  mutate(kmplots = map2(data, CNSig, ~kmplot_fun(type = 'OS', sur_dat = .x, title = .y, curve_labs = c('Low', 'High'), curve_colors = c('#006CB1', '#C9372E'))))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_prognosis_for_thesis/cnsig_kmplots.pdf', width = 11, height = 7)

wrap_plots(comb_cnsig_kmplots$kmplots)

dev.off()


# median survival time

comb_cnsig_surv_median <- comb_cnsig_sur_data %>% 
  mutate(surv_med_res = map(data, ~surv_median(surv_fit(Surv(OS.time, OS) ~ group, data = .))))


# COX analysis ------------------------------------------------------------


# among subtypes

subtype_cox_res <- summary(coxph(Surv(OS.time, OS) ~ group + age + gender + grade + IDH_codel_subtype, data = subtype_sur_data))
# not significant

# individual copy number signatures

comb_cnsig_cox_res <- comb_cnsig_sur_data %>% 
  mutate(cox_res = map(data, ~summary(coxph(Surv(OS.time, OS) ~ group + age + gender + grade + IDH_codel_subtype, data = .))))
# not significant


















































