
library(tidyverse)

cosmic_sig <- read_tsv('/boot3/bio_liaojl/Common_data/COSMIC_v3.3_SBS_GRCh37_exome.txt')

glioma_cosmic_sig_ref <- cosmic_sig %>% 
  select(Type, all_of(str_c('SBS', c(1, 2, 5, 6, '7b', '7c', '7d', '10a', '10b', 11, 13, 14, 15, 19, 30, 37, 40, 42))))

write_tsv(glioma_cosmic_sig_ref, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cosmic_sig_ref.txt')










