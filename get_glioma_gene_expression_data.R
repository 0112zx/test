
library(tidyverse)

glioma_cli_data <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')

genecode_v19_protein_coding_gene <- read_tsv('/boot3/bio_liaojl/Common_data/processed_data/genecode_v19_protein_coding_gene.tsv')

tcga_exp_data <- vroom::vroom('/boot3/bio_liaojl/TCGA_CN_signature_transcriptome_work/Data/Original/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv')

glioma_exp_data_pre <- tcga_exp_data %>% 
  column_to_rownames('gene_id') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('sample') %>% 
  as_tibble() %>% 
  mutate(sample = str_sub(sample, 1, 15)) %>% 
  distinct(sample, .keep_all = TRUE) %>% 
  semi_join(glioma_cli_data, by = 'sample') %>% 
  column_to_rownames('sample') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('gene_id') %>% 
  as_tibble()
# 668 glioma samples


glioma_exp_data <- glioma_exp_data_pre %>% 
  separate(gene_id, into = c('Hugo_Symbol', NA), sep = '\\|') %>% 
  semi_join(genecode_v19_protein_coding_gene, by = 'Hugo_Symbol') %>% 
  distinct(Hugo_Symbol, .keep_all = TRUE)
# 668 glioma samples, 17613 genes


vroom::vroom_write(glioma_exp_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_exp_data.tsv')












