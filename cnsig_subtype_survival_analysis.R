
library(tidyverse)
library(survival)
library(survminer)
library(patchwork)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_weight.Rdata')
glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')

comb_cnsig_sur_data <- glioma_cnsig_weight %>% 
  inner_join(select(glioma_cli_data, bcr_patient_barcode, IDH_codel_subtype, OS, OS.time), by = c('sample' = 'bcr_patient_barcode')) %>% 
  mutate(CNSig_lab = str_remove_all(CNSig, '-')) %>% 
  filter(!is.na(IDH_codel_subtype)) %>% 
  group_by(IDH_codel_subtype, CNSig) %>% 
  mutate(group = ifelse(weight > median(weight), 'High', 'Low'), 
         group = str_c(CNSig_lab, group, sep = '-')) %>% 
  nest() %>% 
  ungroup()


# kmplots


comb_cnsig_kmplots <- comb_cnsig_sur_data %>% mutate(kmplots = map2(IDH_codel_subtype, data, ~kmplot_fun(.x, .y)))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_prognosis/glima_subtype_cnsig_prognosis.pdf', width = 25, height = 18)

wrap_plots(comb_cnsig_kmplots$kmplots)

dev.off()


# functions

kmplot_fun <- function(ttype, sur_dat){
  
  # ttype <- 'ACC'
  # sur_dat <- tcga_cnsig_sur_res$data[[1]]
  
  # pdf('/boot3/bio_liaojl/temp.pdf')
  
  ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group, data = sur_dat),
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

# IDHmut-non-codel: CNSig1-High vs CNSig1-Low classified by median of CNSig1 contribution, p = 0.00058


# CNSig1 in IDHmut-non-codel by different thresholds ----------------------------------------------


IDHmut_noncodel_sur_data <- comb_cnsig_sur_data$data[[7]] %>% 
  mutate(group_0 = ifelse(weight > 0, 'CNSig1-present', 'CNSig1-absent'), 
         group_1 = ifelse(weight == 1, 'CNSig1-whole', 'CNSig1-frac'), 
         group_3tile = as.character(ntile(weight, 3)), 
         group_4tile = as.character(ntile(weight, 4)), 
         group_10tile = as.character(ntile(weight, 10)))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_prognosis/IDHmut_noncodel_CNSig1_kmplots.pdf', width = 15, height = 10)

group_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group, data = IDHmut_noncodel_sur_data),
                         pval = TRUE,
                         xlab = 'Survival time (Months)',
                         risk.table = TRUE,
                         ggtheme = theme_bw(),
                         data = IDHmut_noncodel_sur_data)$plot

group_0_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group_0, data = IDHmut_noncodel_sur_data),
                           pval = TRUE,
                           xlab = 'Survival time (Months)',
                           risk.table = TRUE,
                           ggtheme = theme_bw(),
                           data = IDHmut_noncodel_sur_data)$plot

group_1_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group_1, data = IDHmut_noncodel_sur_data),
                           pval = TRUE,
                           xlab = 'Survival time (Months)',
                           risk.table = TRUE,
                           ggtheme = theme_bw(),
                           data = IDHmut_noncodel_sur_data)$plot

group_3tile_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group_3tile, data = IDHmut_noncodel_sur_data),
                               pval = TRUE,
                               xlab = 'Survival time (Months)',
                               risk.table = TRUE,
                               ggtheme = theme_bw(),
                               data = IDHmut_noncodel_sur_data)$plot

group_4tile_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group_4tile, data = IDHmut_noncodel_sur_data),
                               pval = TRUE,
                               xlab = 'Survival time (Months)',
                               risk.table = TRUE,
                               ggtheme = theme_bw(),
                               data = IDHmut_noncodel_sur_data)$plot

group_10tile_plot <- ggsurvplot(surv_fit(Surv(OS.time, OS) ~ group_10tile, data = IDHmut_noncodel_sur_data),
                                pval = TRUE,
                                xlab = 'Survival time (Months)',
                                risk.table = TRUE,
                                ggtheme = theme_bw(),
                                data = IDHmut_noncodel_sur_data)$plot


group_plot + group_0_plot + group_1_plot + group_3tile_plot + group_4tile_plot + group_10tile_plot


dev.off()


CNSig1_group_idh_noncodel_cli <- comb_cnsig_sur_data$data[[7]] %>% 
  select(sample, CNSig1_weight = weight, CNSig1_group = group) %>% 
  inner_join(glioma_cli_data, by = c('sample' = 'bcr_patient_barcode')) %>% 
  mutate(CNSig1_group = factor(CNSig1_group, levels = c('CNSig1-Low', 'CNSig1-High')))

save(CNSig1_group_idh_noncodel_cli, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/CNSig1_group_idh_noncodel_cli.Rdata')


# CNSig1 in IDHmut-non-codel, Cox analysis ----------------------------------------------

# univariate cox analysis

uni_cox_res <- summary(coxph(Surv(OS.time, OS) ~ CNSig1_group, data = CNSig1_group_idh_noncodel_cli))
# p = 0.000837

uni_cox_res_continuous <- summary(coxph(Surv(OS.time, OS) ~ CNSig1_weight, data = CNSig1_group_idh_noncodel_cli))
# p = 0.0078

# multivariate cox analysis

multi_cox_res <- summary(coxph(Surv(OS.time, OS) ~ CNSig1_group + age + gender + race + grade + histologic_type + TMB, data = CNSig1_group_idh_noncodel_cli))


multi_cox_res <- summary(coxph(Surv(OS.time, OS) ~ CNSig1_group + age + gender + race, data = CNSig1_group_idh_noncodel_cli))


multi_cox_res <- summary(coxph(Surv(OS.time, OS) ~ CNSig1_group + age + grade, data = CNSig1_group_idh_noncodel_cli))



multicox_res_ori <- coxph(sur_formu, data = sin_cancer_data)
multicox_res <- summary(multicox_res_ori)







# survival among subgroups

sur_group_dat <- cnsig_group_cli_data %>% 
  group_by(cluster) %>% 
  filter(n() > 10) %>% 
  ungroup()

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_clustering/subgroup_kmplot.pdf')

ggsurvplot(surv_fit(Surv(OS.time, OS) ~ cluster, data = sur_group_dat),
           pval = TRUE,
           xlab = 'Survival time (Months)',
           # legend.labs = c('Low', 'High'),
           # palette = c('#5082AF', '#C9372E'),
           risk.table = TRUE,
           ggtheme = theme_bw(),
           data = sur_group_dat)

dev.off()

pairwise_res <- pairwise_survdiff(Surv(OS.time, OS) ~ cluster, data = sur_group_dat)




# cox analysis



summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + IDH_codel_subtype, data = cnsig_group_cli_data))

summary(coxph(Surv(OS.time, OS) ~ cluster + age + gender + grade, data = cnsig_group_cli_data))

















































