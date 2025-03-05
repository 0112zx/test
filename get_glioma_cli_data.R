
library(tidyverse)

# clinical features

tcga_glioma_cli <- read_tsv('/boot3/bio_liaojl/pub6/Temp/Liaojl/Data/Mutation_burden_research/Processed/tcga_pan_cli.tsv') %>% 
  filter(type %in% c('LGG', 'GBM')) %>% 
  select(bcr_patient_barcode, type, age = age_at_initial_pathologic_diagnosis, gender, race, OS:DSS.time, PFI, PFI.time)
# 1111 patients

# molecular features

mol_fea_dat <- read_tsv('/boot3/bio_liaojl/Common_data/lgggbm_tcga_pub_clinical_data.tsv') %>% 
  select(bcr_patient_barcode = `Patient ID`, 
         sample = `Sample ID`, 
#         BRAF_KIAA1549_fusion = `BRAF-KIAA1549 fusion`, # only 1 fusion
#         BRAF_V600E_status = `BRAF V600E status`, # only 5 mutant
         chr19_20_co_gain = `Chr 19/20 co-gain`, 
         chr7_gain_or_chr10_loss = `Chr 7 gain/Chr 10 loss`, 
         grade = `Neoplasm Histologic Grade`, 
         histologic_type = `Neoplasm Histologic Type Name`, 
         IDH_codel_subtype = `IDH/codel subtype`, 
         IDH_status = `IDH status`, 
         kps_score = `Karnofsky Performance Score`, 
         MGMT_promoter_status = `MGMT promoter status`, 
         TERT_promoter_status = `TERT promoter status`, 
         transcriptome_subtype = `Transcriptome Subtype`)
# 1122 patients/samples


# Telomere data

telo_dat <- read_tsv('/boot3/bio_liaojl/Common_data/Barthel2017_SuppTab1_TelomereLength_TCGA.csv') %>% 
  filter(Disease %in% c('LGG', 'GBM')) %>% 
  select(sample = SampleID, TP53mut, ATRXmut, TelomeraseSignatureScore, ATRXstatus, DAXXstatus)
# 822 samples

# copy number heterogeneity

CNH_dat <- read_tsv('/boot3/bio_liaojl/Common_data/Erik_van_Dijk_2021nc_CNH.tsv') %>% 
  filter(Type %in% c('LGG', 'GBM')) %>% 
  mutate(sample = str_sub(Samplename, 1, 15)) %>% 
  select(sample, CNH)
# 1085 samples

# purity, ploidy, WGD

pp_dat <- read_tsv('/boot3/bio_liaojl/Common_data/summary.ascatv3TCGA.penalty70.hg19.tsv') %>% 
  filter(cancer_type %in% c('LGG', 'GBM')) %>% 
  mutate(sample = str_sub(barcodeTumour, 1, 15)) %>% 
  select(sample, purity, ploidy, WGD) %>% 
  distinct(sample, .keep_all = TRUE)
# 1010 samples

# aneuploidy, TMB, (arm cnv)

aneu_dat <- read_tsv('/boot3/bio_liaojl/Common_data/Taylor2018_S2.csv') %>% 
  filter(Type %in% c('LGG', 'GBM')) %>% 
  select(sample = Sample, aneuploidy_score = `AneuploidyScore(AS)`, TMB = `Non-silentMutationsperMb`, `1p`:`22q`)
# 1086 samples

# Intratumor Heterogeneity, Fraction Altered

hetero_frac_alt_dat <- read_tsv('/boot3/bio_liaojl/Common_data/TCGA_immune_related_score.tsv') %>% 
  filter(`TCGA Study` %in% c('LGG', 'GBM')) %>% 
  select(bcr_patient_barcode = `TCGA Participant Barcode`, Intratumor_Heterogeneity = `Intratumor Heterogeneity`, Fraction_Altered = `Fraction Altered`)

# merge data

all_mol_dat <- list(mol_fea_dat, telo_dat, CNH_dat, pp_dat, aneu_dat) %>% reduce(left_join, by = 'sample')

glioma_cli_data <- tcga_glioma_cli %>% 
  left_join(all_mol_dat, by = 'bcr_patient_barcode') %>% 
  select(sample, bcr_patient_barcode, everything())
# 1111 patients/samples

write_tsv(glioma_cli_data, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cli_data.tsv')

# glioma_cnsig_subtype_cli_data originated from glioma_cnsig_clustering.R

