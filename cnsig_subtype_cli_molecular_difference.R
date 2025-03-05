
require(tidyverse)
require(gtsummary)
require(gt)
require(patchwork)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_cli_data.Rdata')

glioma_cnsig_subtype_cli_mole_gt_res <- glioma_cnsig_subtype_cli_data %>% 
  mutate(race = factor(race, levels = c('white', 'others')), 
         histologic_type = factor(histologic_type, levels = c('oligodendroglioma', 'oligoastrocytoma', 'astrocytoma', 'glioblastoma')), 
         transcriptome_subtype = factor(transcriptome_subtype, levels = c('proneural', 'neural', 'mesenchymal', 'classical')), 
         kps_score = factor(kps_score, levels = c('<80', '80-100'))) %>% 
  select(cluster, age:race, KPS = kps_score, histologic_type, grade, TP53mut, ATRXmut, chr19_20_co_gain, chr7_gain_or_chr10_loss, 
         MGMT_promoter_status, TERT_promoter_status, IDH_codel_subtype, transcriptome_subtype, 
         TelomeraseSignatureScore, CNH:TMB, Intratumor_Heterogeneity:Fraction_Altered) %>% 
  tbl_summary(by = cluster, missing = "no", statistic = all_continuous() ~ "{median} ({p0}, {p100})") %>%
  add_n() %>%
  bold_labels() %>% 
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_q(method = 'BH')


# get absolute p values

subtype_cli_mole_test_p <- tibble(variable = glioma_cnsig_subtype_cli_mole_gt_res$table_body$variable, 
                                  p_val = glioma_cnsig_subtype_cli_mole_gt_res$table_body$p.value, 
                                  q_val = glioma_cnsig_subtype_cli_mole_gt_res$table_body$q.value) %>% 
  distinct(variable, .keep_all = TRUE)

write_tsv(subtype_cli_mole_test_p, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/subtype_cli_mole_test_p.tsv')

# get summary result

glioma_cnsig_subtype_cli_mole_summary <- glioma_cnsig_subtype_cli_mole_gt_res %>% as_flex_table()

save(glioma_cnsig_subtype_cli_mole_summary, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/glioma_cnsig_subtype_cli_mole_summary.Rdata')
# load('G:/拷贝数突变signature课题/拷贝数signature应用课题/拷贝数signature分型/Result/CNsig_subtype_cli_molecular_difference/glioma_cnsig_subtype_cli_mole_summary.Rdata')

# save to microsoft word

glioma_cnsig_subtype_cli_mole_summary %>% flextable::save_as_docx(path = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/glioma_cnsig_subtype_cli_molecular_summary.docx')


# representative plots ----------------------------------------------------------------

# age, histologic_type, grade, purity, ploidy, aneuploidy_score, Intratumor_Heterogeneity, CNH, IDH_codel_subtype, transcriptome_subtype

### continuous variables (boxplot + point plot)

pb_plot <- function(var_name){
  
  # var_name <- 'age'
  
    sin_pb_plot <- glioma_cnsig_subtype_cli_data %>% 
      select(cluster, sin_var = all_of(var_name)) %>% 
      ggplot(aes(cluster, sin_var, col = cluster)) +
      geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
      geom_boxplot(outlier.color = NA, fill = NA) +
      ggpubr::stat_compare_means() +
      labs(x = NULL, y = var_name) +
      scale_color_manual(breaks = c('Diploid', 'Diploid CIN', 'Tetraploid', 'Tetraploid CTH'), values = RColorBrewer::brewer.pal(n = 4, name = 'Spectral')) +
      # ggsci::scale_color_npg() +
      scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
      theme_classic() +
      theme(legend.position = 'none', 
            axis.text.x = element_text(angle = 45, hjust = 1), 
            axis.ticks.x = element_blank())
    
  # pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/temp.pdf')
  # 
  # sin_pb_plot
  # 
  # dev.off()
  
}


### categorical variables ( barplot)

bar_plot <- function(var_name){
  
  # var_name <- 'histologic_type'
  
  sin_bar_plot <- glioma_cnsig_subtype_cli_data %>% 
    select(cluster, sin_var = all_of(var_name)) %>% 
    drop_na() %>% 
    ggplot(aes(cluster, fill = sin_var)) +
    geom_bar(position = 'fill') +
    labs(x = NULL, y = var_name, fill = NULL) +
    # ggsci::scale_fill_npg() +
    scale_fill_grey(start = 0.8, end = 0.2) +
    guides(fill = guide_legend(nrow = 2)) +
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

comb_plots <- list(pb_plot('age'), bar_plot('histologic_type'), bar_plot('grade'), pb_plot('purity'), pb_plot('ploidy'), 
                   pb_plot('aneuploidy_score'), pb_plot('Intratumor_Heterogeneity'), pb_plot('CNH'), bar_plot('IDH_codel_subtype'), bar_plot('transcriptome_subtype'))


# age, histologic_type, grade, purity, ploidy, aneuploidy_score, Intratumor_Heterogeneity, CNH, IDH_codel_subtype, transcriptome_subtype

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_cli_molecular_difference/subtype_cli_mole_plots.pdf', width = 18, height = 10)

wrap_plots(comb_plots) + plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A')

dev.off()









