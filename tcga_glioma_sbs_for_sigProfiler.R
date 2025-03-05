
library(tidyverse)

glioma_mut_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_mut_data.tsv')

glioma_snv_data <- glioma_mut_data %>% 
  filter(Variant_Type %in% 'SNP') %>% 
  rename(Sample = sample, mut_type = Variant_Type, chrom = Chromosome, pos_start = Start_Position, pos_end = End_Position, ref = Reference_Allele, alt = Tumor_Seq_Allele2) %>% 
  mutate(Project = 'Glioma', ID = '.', Genome = 'GRCh37', Type = 'SOMATIC') %>% 
  select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type) %>% 
  distinct_all()

write_tsv(glioma_snv_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/for_SigProfiler_SBS/glioma_snv_data.tsv')







