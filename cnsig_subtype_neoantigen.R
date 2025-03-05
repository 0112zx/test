
library(tidyverse)

load('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_subtype_label.Rdata')

neoantigen_dat <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Original/TCIA-NeoantigensData.tsv') %>% 
  filter(disease %in% 'GBM')
# learn from PMID: Figure 3C: Neoantigen load across cancers
# neoantigen load was calculated as number of all peptides of each gene

glioma_mut_dat <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_mut_data.tsv')


sam_neo_count <- neoantigen_dat %>% 
  count(patientBarcode, name = 'n_neo') %>% 
  mutate(log_neo = log10(n_neo)) %>% 
  inner_join(select(glioma_cnsig_subtype_label, patient, cluster), by = c('patientBarcode' = 'patient'))


neo_count_plot <- sam_neo_count %>% 
  ggplot(aes(cluster, log_neo, fill = cluster)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, fill = 'white', outlier.color = NA) +
  ggpubr::geom_pwc(hide.ns = FALSE) + 
  ggpubr::stat_compare_means() +
  labs(x = NULL, y = 'log10(neoantigens)') +
  ggsci::scale_fill_jco() +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        plot.title = element_text(hjust = 0.5))


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CNsig_subtype_neoantigen/neo_count_plot.pdf')

neo_count_plot

dev.off()







