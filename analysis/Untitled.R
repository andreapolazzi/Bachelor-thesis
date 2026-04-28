library(tidyverse)
library(readxl)
library(here)
source('scripts/clean_columns.R')

data <- read_xlsx(here('data', 'raw', 'drivers.xlsx'))

data2 <- readRDS(here('data/processed/drivers_telins.rds'))

data3 <- data %>% 
  left_join(data2 %>% select(patient_id, coverage, tumor_purity, tel_ins, tel_ins_norm),
            by = c('Patient_ID_IF_mod' = 'patient_id'))
data3 <- data3 %>% 
  relocate(coverage, tumor_purity, tel_ins, tel_ins_norm, .after = Mixed_def...4)

data3 <- clean_drivers_data(data3)

data3 <- data3 %>% 
  select(patient_id, cancer_type, cancer_type_code, mixed_def, coverage, tumor_purity, tel_ins, 
         tel_ins_norm, tumor_tf_rate, blood_tf_rate, gene, category, driver, primaryTumorLocation, 
         `Intratel_reads/total_reads`, `tel_reads/total reads`, TCAGGG_singletons_norm_by_all_reads, 
         TGAGGG_singletons_norm_by_all_reads, TTGGGG_singletons_norm_by_all_reads, TTCGGG_singletons_norm_by_all_reads,
         TTTGGG_singletons_norm_by_all_reads, ATAGGG_singletons_norm_by_all_reads, CATGGG_singletons_norm_by_all_reads,
         CTAGGG_singletons_norm_by_all_reads, GTAGGG_singletons_norm_by_all_reads, TAAGGG_singletons_norm_by_all_reads,
         tvr_all, tvr_intratel)

saveRDS(data3, file = 'drivers_full.RDS')

pat_summary <- data3 %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  skim()


library(skimr)
library(openxlsx)
report <- skim(data3)
write.xlsx(report, 'dataset_summary2.xlsx')
write.xlsx(pat_summary, "Summary_Pazienti_Unici.xlsx")
