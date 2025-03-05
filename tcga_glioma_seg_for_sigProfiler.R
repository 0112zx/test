
library(tidyverse)

glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')
lgg_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/lgg_cli_data.tsv')

# ABSOLUTE ----------------------------------------------------------------

tcga_cn_dat <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Original/TCGA_mastercalls.abs_segtabs.fixed.txt')

glioma_cn_data <- tcga_cn_dat %>% semi_join(glioma_cli_data, by = c('Sample' = 'sample'))
# 1074 samples

glioma_seg_ascat_form <- glioma_cn_data %>% 
  group_by(Sample) %>% 
  filter(n() > 1) %>% # remove 7 samples with only 1 segment, others all more than 20 segments
  ungroup() %>% 
  select(sample = Sample, chr = Chromosome, startpos = Start, endpos = End, nMajor = Modal_HSCN_2, nMinor = Modal_HSCN_1)
# 1067 samples

write_tsv(glioma_seg_ascat_form, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/glioma_seg_ascat_form.tsv')


# ASCAT -------------------------------------------------------------------

tcga_cn_ascat <- tibble(f_name = list.files('/boot3/bio_liaojl/Common_data/segments/')) %>% # 10,674 samples
  mutate(dat = map(f_name, ~read_tsv(str_c('/boot3/bio_liaojl/Common_data/segments/', .)))) %>% 
  unnest(dat) %>% 
  select(-f_name)
# bbb <- tcga_cn_ascat %>% distinct(sample) %>% mutate(patient = str_sub(sample, 1, 12)) # no duplicated patients

# vroom::vroom_write(tcga_cn_ascat, '/boot3/bio_liaojl/Common_data/tcga_cn_ascat.tsv')

glioma_cn_ascat <- tcga_cn_ascat %>% semi_join(glioma_cli_data, by = c('sample' = 'bcr_patient_barcode'))
# 1004 samples

write_tsv(glioma_cn_ascat, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/glioma_cn_ascat.tsv')


lgg_cn_ascat <- tcga_cn_ascat %>% semi_join(lgg_cli_data, by = c('sample' = 'bcr_patient_barcode'))
# 498 samples

write_tsv(lgg_cn_ascat, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler/lgg_cn_ascat.tsv')




