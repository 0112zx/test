
library(tidyverse)

tcga_mut_data <- vroom::vroom('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_mut.tsv')

glioma_mut_data <- tcga_mut_data %>% 
  filter(ttype %in% c('LGG', 'GBM')) %>% 
  mutate(sample = str_sub(Tumor_Sample_Barcode, 1, 15)) %>% 
  select(sample, Hugo_Symbol:Tumor_Seq_Allele2, HGVSp_Short:IMPACT)

write_tsv(glioma_mut_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_mut_data.tsv')

















