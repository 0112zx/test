
library(tidyverse)
library(patchwork)

glioma_cnsig_def <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/ASCAT/sigprofiler_output/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Signatures/CNV48_De-Novo_Signatures.txt') %>% 
  column_to_rownames('MutationType') %>% 
  t()

# Similarity to 21 CN signatures -----------------------------------------

CNsig21_definitions <- read_tsv('/boot3/bio_liaojl/Common_data/CNsig21_definitions.tsv') %>% 
  column_to_rownames('CNclass') %>% 
  t()

# pairwise cosine similarity

cn_sig_sim <- c()

for(i in rownames(glioma_cnsig_def)){
  
  i_cos_vec <- c()
  
  for(j in rownames(CNsig21_definitions)){
    
    i_cos_vec <- c(i_cos_vec, lsa::cosine(glioma_cnsig_def[i, ], CNsig21_definitions[j, ]))
    
  }
  
  cn_sig_sim <- rbind(cn_sig_sim, i_cos_vec)
  
}

# rownames(cn_sig_sim) <- str_c('CN-Sig-', 1:6)
rownames(cn_sig_sim) <- rownames(glioma_cnsig_def)
colnames(cn_sig_sim) <- rownames(CNsig21_definitions)

cn_sig_sim <- cn_sig_sim %>% as.data.frame() %>% rownames_to_column('Signature')

write_tsv(cn_sig_sim, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/cn_sig_similarity.tsv')

cn_sig_sim <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/cn_sig_similarity.tsv')

# plot

CN_sig_order <- rownames(CNsig21_definitions) %>% str_sort(numeric = TRUE)

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/cn_sig_similarity.pdf', width = 12, height = 6)

cn_sig_sim %>% 
  pivot_longer(-Signature, names_to = 'CN_signature', values_to = 'cosine_similarity') %>% 
  mutate(cosine_similarity = round(cosine_similarity, 2), CN_signature = factor(CN_signature, levels = CN_sig_order)) %>% 
  ggplot(aes(CN_signature, Signature, fill = cosine_similarity)) +
  geom_tile() +
  geom_text(aes(label = cosine_similarity)) +
  labs(x = NULL, y = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()


rownames(glioma_cnsig_def) <- c(str_c('CNSig', 1:4), 'CNSig6', 'CNSig5')
# save(glioma_cnsig_def, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_def.Rdata')

glioma_cnsig_def_tib <- glioma_cnsig_def[str_c('CNSig', 1:6), ] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('CNclass') %>% 
  as_tibble()

write_tsv(glioma_cnsig_def_tib, '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_def_tib.tsv')

# replot 

cn_sig_sim_re <- cn_sig_sim %>% mutate(Signature = c(str_c('CNSig', 1:4), 'CNSig6', 'CNSig5')) %>% arrange(Signature)

pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/cn_sig_similarity_re.pdf', width = 12, height = 4)

cn_sig_sim_re %>% 
  pivot_longer(-Signature, names_to = 'CN_signature', values_to = 'cosine_similarity') %>% 
  mutate(cosine_similarity = round(cosine_similarity, 2), CN_signature = factor(CN_signature, levels = CN_sig_order)) %>% 
  ggplot(aes(CN_signature, Signature, fill = cosine_similarity)) +
  geom_tile() +
  geom_text(aes(label = cosine_similarity)) +
  labs(x = NULL, y = NULL) +
  scale_fill_gradientn(colors = c('white', '#C9372E')) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()



# CN signataure exposure data ---------------------------------------------

glioma_cnsig_activity <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/ASCAT/sigprofiler_output/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Activities/CNV48_De-Novo_Activities_refit.txt')

colnames(glioma_cnsig_activity) <- c('sample', str_c('CNSig', 1:4), 'CNSig6', 'CNSig5')

glioma_cnsig_weight <- glioma_cnsig_activity %>% 
  pivot_longer(-sample, names_to = 'CNSig', values_to = 'count') %>% 
  group_by(sample) %>% 
  mutate(weight = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  arrange(sample, CNSig)

save(glioma_cnsig_weight, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/glioma_cnsig_weight.Rdata')


# plot barplot(prevalence) and boxplot(weight)

bar_p_data <- glioma_cnsig_weight %>% 
  group_by(CNSig) %>% 
  summarise(frac = mean(weight > 0), 
            exp1 = sum(weight > 0), 
            exp2 = n())

bar_p <- bar_p_data %>% 
  mutate(expression = map2(exp1, exp2, ~bquote(over(.(.x), .(.y))))) %>% 
  ggplot(aes(CNSig, frac, fill = CNSig)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = expression), parse = TRUE, vjust = -0.15, size = 3) +
  labs(x = NULL, y = 'Prevalence') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  ggsci::scale_fill_jco() +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank())

box_w <- glioma_cnsig_weight %>% 
  filter(weight > 0) %>% 
  ggplot(aes(CNSig, weight, col = CNSig)) +
  geom_jitter(width = 0.25, size = 1, color = '#EAEAEA') +
  geom_boxplot(outlier.color = NA, fill = NA) +
  labs(x = NULL, y = 'Weight') +
  ggsci::scale_color_jco() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())


pdf('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result/CN_signature_identification/CNSig_prevalence_comb.pdf', width = 6, height = 6)

box_w / bar_p

dev.off()


cnsig_stat <- glioma_cnsig_weight %>% 
  filter(weight > 0) %>% 
  group_by(CNSig) %>% 
  summarise(min_w = min(weight), 
            med_w = median(weight), 
            max_w = max(weight))



# CN signataure exposure data (LGG) ---------------------------------------------

lgg_cnsig_activity <- read_tsv('/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Result_LGG/CN_signature_identification/sigprofiler_output/CNV48/Suggested_Solution/CNV48_De-Novo_Solution/Activities/CNV48_De-Novo_Activities_refit.txt')

# A: CN1, B: CN20, C: CN2, D: CN9, E: CN24

colnames(lgg_cnsig_activity) <- c('sample', 'CNSig1', 'CNSig4', 'CNSig2', 'CNSig3', 'CNSig5')

lgg_cnsig_weight <- lgg_cnsig_activity %>% 
  pivot_longer(-sample, names_to = 'CNSig', values_to = 'count') %>% 
  group_by(sample) %>% 
  mutate(weight = count/sum(count)) %>% 
  ungroup() %>% 
  select(-count) %>% 
  arrange(sample, CNSig)
# 498 samples

save(lgg_cnsig_weight, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/lgg_cnsig_weight.Rdata')

lgg_cnsig_activity_scale <- lgg_cnsig_activity %>% 
  pivot_longer(-sample, names_to = 'CNSig', values_to = 'count') %>% 
  mutate(count_ln = log(count + 1)) %>% 
  group_by(CNSig) %>% 
  mutate(count_ln_scale = scale(count_ln)[, 1]) %>% 
  ungroup() %>% 
  select(sample, CNSig, count_ln_scale) %>% 
  arrange(sample, CNSig)

save(lgg_cnsig_activity_scale, file = '/boot3/bio_liaojl/Copy_number_signature_application_work/Glioma_classification/Data/Processed/lgg_cnsig_activity_scale.Rdata')










