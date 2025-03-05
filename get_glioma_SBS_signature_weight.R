
library(tidyverse)

glioma_sig_qc <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/SBS_signature_refitting/Assignment_Solution/Solution_Stats/Assignment_Solution_Samples_Stats.txt')
# 900 samples

# remove low quality samples

glioma_sig_qc_filter <- glioma_sig_qc %>% 
  filter(`Total Mutations` >= 20, `Cosine Similarity` >= 0.6)
# 828 samples

glioma_sbs_sig_weight <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/SBS_signature_refitting/Assignment_Solution/Activities/Assignment_Solution_Activities.txt') %>% 
  semi_join(glioma_sig_qc_filter, by = c('Samples' = 'Sample Names')) %>% 
  pivot_longer(-Samples, names_to = 'Signature', values_to = 'Exposure') %>% 
  group_by(Samples) %>% 
  mutate(Weight = Exposure/sum(Exposure)) %>% 
  ungroup()

save(glioma_sbs_sig_weight, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_sbs_sig_weight.Rdata')







