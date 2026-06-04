library(tidyverse)
library(here)
source('scripts/clean_columns.R')

# Drivers telins ####
# Load datasets ####
drivers <- readxl::read_xlsx('data/raw/drivers.xlsx')
clean_drivers <- clean_drivers_data(drivers)
clean_drivers <- clean_drivers %>% 
  distinct(patient_id, .keep_all = TRUE) %>% 
  rename(intratel_reads_total = `Intratel_reads/total_reads`,
         intratel_reads_total_order = `Intratel_reads/total_reads_ORDER`,
         tel_reads_total = `tel_reads/total reads`,
         tel_reads_total_order = `tel_reads/total reads_ORDER`)

code_conversion <- readxl::read_xlsx('data/raw/code_conversion_patient_number.xlsx')

coverage <- readxl::read_xlsx('data/raw/Samples_hartwig_coverage.xlsx')

TTAGGG_counts <- read.csv('../telomerehunter_results/spectrum.080426.maps.TTAGGG.csv', sep = ';')

# Check structures and matches
glimpse(code_conversion)
summary(code_conversion)

code_conversion %>% 
  filter(!complete.cases(.)) %>% 
  print(n = Inf)

code_conv_clean <- code_conversion %>% 
  drop_na(Patient_ID_IF, Patient_ID_IF_mod) %>% 
  distinct(Patient_ID_IF_mod, .keep_all = TRUE)

ids_code <- code_conv_clean$Patient_ID_IF_mod
ids_drivers <- clean_drivers$patient_id

length(ids_code)
length(ids_drivers)
# They differ

common_ids <- intersect(ids_code, ids_drivers)

only_in_code  <- setdiff(ids_code, ids_drivers)
only_in_driv  <- setdiff(ids_drivers, ids_code)

# Let's see which ones are only in code_conversion
code_conv_clean %>% 
  filter(Patient_ID_IF_mod %in% only_in_code) %>% 
  print(n = Inf)

# Join the datasets, adding patient_code column
drivers_joined <- clean_drivers %>% 
  left_join(code_conv_clean %>% 
              select(Patient_ID_IF_mod, ...1),
            by = join_by(patient_id == Patient_ID_IF_mod)
            ) %>% 
  select(patient_id, ...1, everything()) %>% 
  rename(patient_code = ...1)

# Now let's do the same for the coverage
coverage <- coverage %>% 
  rename(patient_code = sample)

drivers_full <- drivers_joined %>% 
  left_join(
    coverage %>% select(patient_code, coverage),
    by = 'patient_code'
  ) %>% 
  relocate(coverage, .after = patient_code)

sum(!is.na(drivers_full$coverage))
sum(is.na(drivers_full$coverage)) 

drivers_full %>% 
  filter(!complete.cases(coverage)) %>% 
  print(n = Inf)

# Tumor purity ####
metadata<- readxl::read_xlsx('data/raw/Metadatos_2_to_ALV.xlsx')

drivers_final <- drivers_full %>% 
  left_join(
    metadata %>% select(Patient_ID_IF_mod, tumorPurity),
    by = join_by(patient_id == Patient_ID_IF_mod)
  ) %>% 
  relocate(tumorPurity, .after = coverage) %>% 
  rename(tumor_purity = tumorPurity)

summary(drivers_final)


# Clean TTAGGG dataset ####
TTAGGG_counts_tum_no_junction <- TTAGGG_counts %>% 
  filter(str_starts(Sample, "tumor")) %>% 
  select(-ncol(TTAGGG_counts), -contains("junction"))


# Average TTAGGG counts for each band - visualization ####
band_means <- TTAGGG_counts_tum_no_junction %>% 
  select(-Sample) %>% 
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "band",
    values_to = "mean_counts"
  )

# tabella delle bande più alte
band_means %>% 
  mutate(is_outlier = mean_counts > quantile(mean_counts, 0.95)) %>% 
  arrange(desc(mean_counts)) %>% 
  slice(1:50) %>% 
  View()

# plot delle 50 bande con media più alta
band_means %>%
  arrange(desc(mean_counts)) %>%
  slice(1:50) %>%
  ggplot(aes(x = reorder(band, log(mean_counts)), y = log(mean_counts))) +
  geom_col() +
  coord_flip()


# Adding telomere insertion information ####

# 1. colonne genomiche
genomic_cols <- colnames(TTAGGG_counts_tum_no_junction) %>%
  setdiff("Sample")

# 2. parse nome colonna -> chr + arm + banda
col_info <- tibble(col = genomic_cols) %>%
  separate(col, into = c("chr", "arm_band"), sep = "_", remove = FALSE) %>%
  mutate(
    arm = substr(arm_band, 1, 1),
    band = as.numeric(str_extract(arm_band, "\\d+\\.?\\d*"))
  )

# 3. identifica bande terminali
terminal_cols <- col_info %>%
  group_by(chr, arm) %>%
  filter(band == max(band)) %>%
  pull(col)

# 4. identifica bande outlier forti in base alla media
outlier_cols <- band_means %>%
  filter(mean_counts > quantile(mean_counts, 0.95)) %>%
  pull(band)

# highlight = colonne interessanti: no terminali MA outlier
band_means_flag <- band_means %>%
  mutate(
    is_outlier = mean_counts > quantile(mean_counts, 0.95),
    is_terminal = band %in% terminal_cols,
    highlight = is_outlier & !is_terminal
  )

# distribution of log(counts) highlight + outlier
band_means_flag %>%
  arrange(desc(mean_counts)) %>%
  slice(1:50) %>%
  ggplot(aes(x = reorder(band, mean_counts), y = log(mean_counts), fill = highlight, colour = is_outlier)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red"))+
  scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = 'blue'))

# distribution of log(counts) highlight + terminal
band_means_flag %>%
  arrange(desc(mean_counts)) %>%
  slice(1:60) %>%
  ggplot(aes(x = reorder(band, mean_counts), y = log(mean_counts), fill = highlight, colour = is_terminal)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "red"))+
  scale_color_manual(values = c('FALSE' = 'grey70', 'TRUE' = 'green'))

# 5. scelta colonne
# bande outlier non terminali
highlighted_cols <- band_means_flag %>%
  filter(highlight) %>%
  pull(band)

# highlighted ma da tenere (selezione manuale)
keep_cols <- c(
  "X9_q13", "X2_q21.2", "X9_p24.2", "X2_q13", "X2_p22.3",
  "Y_q11.223", "X9_q34.11", "X_p22.11", "X5_p14.3",
  "X8_q12.3", "X22_q11.21"
)

# colonne finali da escludere: tutte le terminali + tutti gli highlighted NON presenti in keep_cols
exclude_cols <- c(
  terminal_cols,
  setdiff(highlighted_cols, keep_cols)
)

# colonne da tenere per la somma finale
internal_cols <- setdiff(genomic_cols, exclude_cols)

# 6. somma finale
tel_ins_dataset <- TTAGGG_counts_tum_no_junction %>% 
  mutate(tel_ins = rowSums(across(all_of(internal_cols)), na.rm = TRUE)) %>% 
  select(Sample, tel_ins, everything())



# New tel_ins column ####
tel_ins_dataset_join <- tel_ins_dataset %>% 
  mutate(patient_code = str_remove(Sample, "^tumor_"),
         patient_code = str_remove(patient_code, "-.*$")
         ) %>% 
  select(patient_code, tel_ins)

drivers_new <- drivers_final %>% 
  left_join(tel_ins_dataset_join, by = 'patient_code') %>% 
  mutate(
    tel_ins_norm = case_when(
      is.na(tel_ins) ~ NA_real_,
      is.na(coverage) | is.na(tumor_purity) ~ NA_real_,
      coverage == 0 | tumor_purity == 0 ~ NA_real_,
      TRUE ~ tel_ins / (coverage * tumor_purity)
    )
  ) %>% 
  relocate(tel_ins, tel_ins_norm, .after = mixed_def)


saveRDS(drivers_new, 'data/processed/drivers_telins.rds')

plot(drivers_new$tumor_tf_rate, drivers_new$tel_ins_norm)

drivers_new %>% 
  filter(tumor_tf_rate<50, tel_ins_norm<7500) %>% 
  ggplot(aes(tumor_tf_rate, tel_ins_norm)) +
  geom_point()


# Blood TF rate primary ####
## Load datasets####
library(here)
primary_data <- readxl::read_xlsx(here('data', 'raw', 'PCAWG_variables_PCA_woK.xlsx'))
clean_pr_data <- clean_pcawg_data(primary_data)
blood_data <- readxl::read_xlsx(here('data', 'raw', 'PCAWG_TF_rate_tumor_blood.xlsx'))

primary_blood_data <- clean_pr_data %>% 
  left_join(blood_data,
            by = c('donor_id' = 'icgc_donor_id')) %>% 
  relocate(`TF rate_blood`, .after = tf_rate) %>%
  rename(tf_blood_rate = `TF rate_blood`, tf_primary_rate = tf_rate) %>% 
  select(-c(ncol(.)-2, ncol(.)-1))

write_xlsx(primary_blood_data, here('data', 'processed', 'PCAWG_primary.xlsx'))  


# Metastatic dataset + metadata ####
library(writexl)
library(here)
library(tidyverse)
df1 <- readxl::read_xlsx('data/raw/drivers.xlsx')
source(here('scripts', 'clean_columns.R'))
source(here('scripts', 'useful_functions.R'))
df1 <- clean_drivers_data(df1)
df2 <- readxl::read_xlsx(here('data', 'raw', 'Tumor_blood_to_ALV_2.xlsx'), sheet = 'Sheet1')
df3 <- readxl::read_xlsx(here('data', 'raw', 'Metadatos_to_ALV.xlsx'))
df4 <- readxl::read_xlsx(here('data', 'raw', 'Metadatos_2_to_ALV.xlsx'))

df3 <- df3[-1,]

metastatic_full <- df1 %>% 
  left_join(df2, join_by('patient_id' == 'Patient_ID_IF_mod'), suffix = c('', '_df2')) %>% 
  left_join(df3, join_by('patient_id' == 'Patient_ID_IF_mod'), suffix = c('', '_df3')) %>% 
  left_join(df4, join_by('patient_id' == 'Patient_ID_IF_mod'), suffix = c('', '_df4'))

metastatic_full <- janitor::clean_names(metastatic_full)

metastatic_full <- metastatic_full %>% 
  rename_with(~ str_replace(.x, "^([a-z]{6})(?=_singleton)", toupper))

# Group similar columns together
metastatic_full <- metastatic_full %>% 
  relocate(cancer_type2, cancer_type_16, cancer_type_23, .after = cancer_type) %>% 
  relocate(cancer_type_code2, cancer_type_code_17, cancer_type_code_24, cancer_type_code_df3, .after = cancer_type_code) %>% 
  relocate(mixed, mixed_def2, mixed_def_2, .after = mixed_def) %>% 
  relocate(tumor_telomere_fusion_rate_gt_end_coverage_purity, .after = tumor_tf_rate) %>% 
  relocate(primary_tumor_location_df2, primary_tumor_location_df4, .after = primary_tumor_location) %>% 
  relocate(tumor_tf_blood_tf, .after = blood_tf_rate_outwards) %>% 
  relocate(primary_tumor_type, primary_tumor_sub_type, primary_tumor_sub_location, primary_tumor_extra_details, .after = cancer_type_23) %>% 
  relocate(cancer_subtype, .after = primary_tumor_sub_type) %>% 
  relocate(tumor_purity_2, .after = tumor_purity) %>% 
  relocate(biopsy_site_2, biopsy_site_df4, simplified_biopsiy_site, .after = biopsy_site) %>% 
  relocate(biopsy_location_df4, .after = biopsy_location) %>% 
  relocate(primary_tumor_type_df4, .after = primary_tumor_type) %>% 
  relocate(primary_tumor_sub_type_df4, .after = primary_tumor_sub_type) %>% 
  relocate(primary_tumor_extra_details_df4, .after = primary_tumor_extra_details) %>% 
  relocate(has_systemic_pre_treatment_df4, .after = has_systemic_pre_treatment) %>% 
  relocate(has_radiotherapy_pre_treatment_df4, .after = has_radiotherapy_pre_treatment) %>% 
  relocate(treatment_given_df4, .after = treatment_given) %>% 
  relocate(treatment_start_date_df4, .after = treatment_start_date) %>% 
  relocate(treatment_end_date_df4, .after = treatment_end_date) %>% 
  relocate(consolidated_treatment_type_df4, .after = consolidated_treatment_type) %>% 
  relocate(response_date_df4, .after = response_date) %>% 
  relocate(response_measured_df4, .after = response_measured) %>% 
  relocate(first_response_df4, .after = first_response) %>% 
  relocate(treatment_df4, .after = treatment) %>% 
  relocate(primary_tumor_sub_location_df4, .after = primary_tumor_sub_location) %>% 
  relocate(tumor_purity_2, .after = mixed_def2) %>% 
  relocate(equal, .after = mixed_def) %>% 
  relocate(tissue_type, .before = tissue_type_14) %>% 
  relocate(tissue_group, .before = tissue_group_15)

  
# Drop real duplicates or order columns
metastatic_full <- metastatic_full %>% 
                      select(-cancer_type_code_df3, -cancer_type_code_17) %>%
                      select(-cancer_type_16, -cancer_type_df3) %>% 
                      select(-mixed_def_2, -cancer_subtype_df3) %>% 
                      select(-tumor_telomere_fusion_rate_gt_end_coverage_purity) %>% 
                      select(-blood_telomere_fusion_rate_gt_end_coverage_purity) %>% 
                      select(-tumor_tf_rate_inwards_df2) %>% 
                      select(-tumor_tf_rate_outwards_df2) %>% 
                      select(-blood_tf_rate_inwards_df2) %>% 
                      select(-blood_tf_rate_outwards_df2) %>% 
                      select(-primary_tumor_location_df2) %>% 
                      select(-ends_with('_order')) %>%
                      select(-ends_with('_rank')) %>% 
                      select(-tissue_type_21, -tissue_group_22) %>%
                      select(-tumor_purity) %>% 
                      select(-gender_df3, -gender_df4) %>% 
                      select(-birth_year_df4, -death_date_df4, -icgc_sample_id, -icgc_specimen_id) %>% 
                      select(-biopsy_date_df4)

write_xlsx(metastatic_full, here('data','processed','metastatic_full.xlsx'))

# Another filtering
data <- readxl::read_xlsx(here('data', 'processed', 'metastatic_full.xlsx'), na = c('', 'NA', 'null', 'NULL', 'unknown', 'Unknown'))
data_no_red <- data %>% 
  select(-c(cancer_type2, cancer_type_23, primary_tumor_sub_location, primary_tumor_type_df4, primary_tumor_sub_location_df4, primary_tumor_extra_details, 
            cancer_type_code, cancer_type_code2, cancer_type_code_24, equal, mixed, mixed_def2, chromosome, chromosome_band,
            transcript_id, canonical_transcript, missense, nonsense, splice, frameshift, inframe, biallelic,
            primary_tumor_location_df4, TCAGGG_singletons_norm_by_intratel_reads, TGAGGG_singletons_norm_by_intratel,
            TTGGGG_singletons_norm_by_intratel, TTCGGG_singletons_norm_by_intratel_reads, TTTGGG_singletons_norm_by_intratel_reads,
            ATAGGG_singletons_norm_by_intratel_reads, CATGGG_singletons_norm_by_intratel_reads, CTAGGG_singletons_norm_by_intratel_reads,
            GTAGGG_singletons_norm_by_intratel_reads, TAAGGG_singletons_norm_by_intratel_reads, tvr_intratel, tissue_type, tissue_type_14,
            biopsy_location, biopsy_location_df4, has_systemic_pre_treatment, has_radiotherapy_pre_treatment, treatment_given, treatment_start_date,
            treatment_end_date, response_date, response_measured, first_response, treatment_df4, progression_status_code,
            is_blacklisted_sample, x12, x16, is_blacklisted, is_blacklisted_cohort, blacklist_comment, n_biopsies_in_patient, n_cancer_types_in_patient,
            is_selected_biopsy, doids, consolidated_treatment_type))

# Handle 0's that should be NA's
no_zero_cols <- c('primary_tumor_type', 'primary_tumor_sub_type')
data_no_red_na <- data_no_red %>% 
  mutate(across(all_of(no_zero_cols), ~na_if(.x, '0')))

write_xlsx(data_no_red_na, here('data','processed','metastatic_reduced.xlsx'))
# Data has been reordered manually afterwards



# Metastatic intrachromosomal reads TH2 ####
met_red_edit <- readxl::read_xlsx(here('data', 'processed', 'metastatic_reduced_edited.xlsx'))
th2_results <- read.csv(here('data','raw','telomerehunter2_r3.results','telomerehunter2_r3.summaries.tsv'), sep = '\t')
th_singleton_results <- readr::read_tsv(
  here('..', 'telomerehunter_results', 'telomere-hunter-all-summary-080426_nolog.tsv.gz'),
  show_col_types = FALSE
)
code_conv_df <- readRDS(here('data', 'processed', 'drivers_telins.rds'))
total_reads_used_df <- readxl::read_xlsx(here('data', 'raw', 'Samples_hartwig_coverage_A.xlsx'))
total_reads_used_missing_df <- readr::read_csv(
  here('data', 'raw', 'total_reads_used_missing.csv'),
  show_col_types = FALSE
) %>%
  transmute(sample = patient_code, Total_reads_used = TOTAL_READS_USED)

total_reads_used_df <- total_reads_used_df %>%
  rows_patch(total_reads_used_missing_df, by = "sample", unmatched = "ignore") %>%
  bind_rows(
    total_reads_used_missing_df %>%
      anti_join(total_reads_used_df, by = "sample")
  )

singleton_patterns <- c(
  "TCAGGG", "TGAGGG", "TTGGGG", "TTCGGG", "TTTGGG",
  "ATAGGG", "CATGGG", "CTAGGG", "GTAGGG", "TAAGGG"
)
singleton_norm_cols <- paste0(singleton_patterns, "_singletons_norm_by_all_reads")
singleton_dist_cols <- paste0(singleton_patterns, "_singleton_dist")
singleton_rel_tol <- 5e-3

# add patient code
met_red_edit_pat_code <- met_red_edit %>% 
  left_join(code_conv_df %>% select(patient_id, patient_code), by = 'patient_id') %>% 
  relocate(patient_code, .after = patient_id)

# TelomereHunter2's singleton distance is:
#   log2(singleton_tumor_norm / singleton_control_norm) -
#   log2(tel_content_tumor / tel_content_control)
# where singleton_norm = singleton_count / total_reads and tel_content are BOTH on the
# SAME TelomereHunter-native total_reads of the input BAM.
#
# IMPORTANT: TH was run on TelFusDetector-filtered BAMs, so total_reads is the filtered
# read count (~1e5), not WGS depth (~1e9). That is fine here because the same total_reads
# normalizes both the singleton term and the tel_content term, so the filtered scale
# cancels and the metric is preserved.

# Step 1: tumor & control singleton norms + tel_content per matched pair (keyed by tumor_code).
singleton_dist_long <- th_singleton_results %>%
  mutate(tumor_code = str_remove(PID, "-.*$")) %>%
  select(PID, tumor_code, sample, total_reads, tel_content, all_of(singleton_norm_cols)) %>%
  pivot_longer(
    cols = all_of(singleton_norm_cols),
    names_to = "singleton_col",
    values_to = "singleton_norm"
  ) %>%
  mutate(pattern = str_remove(singleton_col, "_singletons_norm_by_all_reads$")) %>%
  group_by(tumor_code, pattern) %>%
  summarise(
    norm_tumor          = first(na.omit(singleton_norm[sample == "tumor"])),
    norm_control        = first(na.omit(singleton_norm[sample == "control"])),
    total_reads_tumor   = first(na.omit(total_reads[sample == "tumor"])),
    tel_content_tumor   = first(na.omit(tel_content[sample == "tumor"])),
    tel_content_control = first(na.omit(tel_content[sample == "control"])),
    .groups = "drop"
  )

# Step 2: sanity check the join at the RAW-COUNT level. The metastatic table stores tumor
# singletons normalized by WGS reads, while the TH summary normalizes by filtered total_reads;
# recovering the raw count from each (norm * its own denominator) must agree if the right
# samples are joined.
tumor_singleton_check <- met_red_edit_pat_code %>%
  select(patient_id, patient_code, all_of(singleton_norm_cols)) %>%
  pivot_longer(
    cols = all_of(singleton_norm_cols),
    names_to = "singleton_col",
    values_to = "singleton_norm_existing"
  ) %>%
  mutate(pattern = str_remove(singleton_col, "_singletons_norm_by_all_reads$")) %>%
  left_join(total_reads_used_df %>% select(patient_code = sample, wgs = Total_reads_used),
            by = "patient_code") %>%
  left_join(
    singleton_dist_long %>% select(tumor_code, pattern, norm_tumor, total_reads_tumor),
    by = join_by(patient_code == tumor_code, pattern)
  ) %>%
  filter(!is.na(singleton_norm_existing), !is.na(norm_tumor), !is.na(wgs)) %>%
  mutate(
    raw_met = singleton_norm_existing * wgs,
    raw_th  = norm_tumor * total_reads_tumor,
    rel_diff = abs(raw_met - raw_th) /
      pmax(abs(raw_met), abs(raw_th), .Machine$double.eps),
    within_tolerance = (raw_met == 0 & raw_th == 0) | rel_diff <= singleton_rel_tol
  )

failed_singleton_patients <- tumor_singleton_check %>%
  filter(!within_tolerance) %>%
  summarise(
    n_failed_singleton_checks = n_distinct(pattern),
    .by = c(patient_id, patient_code)
  )

if (nrow(failed_singleton_patients) > 0) {
  warning("TH-native tumor singleton norms disagree with the metastatic table for ",
          nrow(failed_singleton_patients), " patient(s); inspect failed_singleton_patients.")
  print(failed_singleton_patients)
}

# Step 3: singleton_dist on the TH-native scale (filtered total_reads cancels).
# Drop the few patients whose tumor counts are inconsistent between the metastatic table and
# the TH summary (flagged above).
singleton_dist_df <- singleton_dist_long %>%
  filter(!tumor_code %in% failed_singleton_patients$patient_code) %>%
  mutate(
    singleton_dist = if_else(
      norm_tumor == 0 | norm_control == 0 |
        tel_content_tumor == 0 | tel_content_control == 0,
      NA_real_,
      log2(norm_tumor / norm_control) - log2(tel_content_tumor / tel_content_control)
    ),
    singleton_dist_col = paste0(pattern, "_singleton_dist")
  ) %>%
  select(tumor_code, singleton_dist_col, singleton_dist) %>%
  pivot_wider(names_from = singleton_dist_col, values_from = singleton_dist)

# Add Total_reads_used to the main table for intrachromosomal normalization
total_reads_used <- met_red_edit_pat_code %>%
  left_join(total_reads_used_df %>% select(sample, Total_reads_used),
            join_by(patient_code == sample))

# Extract tumor intrachromosomal reads from TH2, keyed by patient code
intrachrom_df <- th2_results %>%
  filter(sample == 'tumor') %>%
  mutate(patient_code = str_remove(PID, "-.*$"))

# Assemble final table: intrachromosomal reads + singleton distances + normalization
met_intrachrom <- total_reads_used %>%
  left_join(intrachrom_df %>% select(patient_code, intrachromosomal_reads), by = 'patient_code'
  ) %>% 
  left_join(singleton_dist_df, by = join_by(patient_code == tumor_code)) %>%
  mutate(intrachrom_reads_total_reads = intrachromosomal_reads/Total_reads_used) %>% 
  relocate(intrachrom_reads_total_reads, .before = intratel_reads_total_reads) %>% 
  relocate(any_of(singleton_dist_cols), .after = TAAGGG_singletons_norm_by_all_reads) %>%
  select(where(~ !all(is.na(.))), -patient_code, -Total_reads_used, -intrachromosomal_reads, -CATGGG_singletons_norm_by_all_reads, -CATGGG_singleton_dist)
library(writexl)
write_xlsx(met_intrachrom, here('data', 'processed', 'metastatic_red_edit_singleton_dist.xlsx'))



# Cancer grouping ####
## Primary ####
data <- readxl::read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))

mes_ori <- c('Bone-Epith', 'Bone-Osteosarc', 'CNS-LGG', 'CNS-PiloAstro', 'Panc-Endocrine', 'SoftTissue-Leiomyo', 'SoftTissue-Liposarc')
nonmes_ori <- setdiff(unique(data$cancer_type), mes_ori)

mesenchymal_grouping_prim <- data %>% 
  mutate(cancer_group =
    case_when(
      cancer_type %in% mes_ori ~ 'Mesenchymal_origin',
      cancer_type %in% nonmes_ori ~ 'Non_Mesenchymal_origin',
      TRUE ~ NA_character_
    )
  )

write_xlsx(mesenchymal_grouping, here('data', 'processed', 'PCAWG_primary.xlsx'))


## Metastatic ####
data <- readxl::read_xlsx(here('data', 'processed', 'metastatic_red_edit_singleton_dist.xlsx'))

ori_data <- readxl::read_xlsx(here('data','raw', 'Hartwig_mesenchimal.xlsx'))

mesenchymal_grouping_met <- data %>% 
  left_join(ori_data, join_by(mixed_def == Cancer_type)) %>% 
  rename(tumor_origin = 'Mesenchimal_origin?') %>% 
  relocate(tumor_origin, .after = cancer_type) %>%
  mutate(tumor_origin = case_when(
    tumor_origin == 'Mesenchimal' ~ 'Mesenchymal',
    tumor_origin == 'Non_mesenchimal' ~ 'Non_mesenchymal',
    TRUE ~ tumor_origin
  ))

writexl::write_xlsx(mesenchymal_grouping_met, here('data', 'processed', 'metastatic_red_edit_mesenchymal.xlsx'))



# Gene mutations primary ####
genes_raw <- readxl::read_xlsx(
  here('data', 'raw', 'TableS3_panorama_driver_mutations_pcawg_v2_18042018_IF_mod.xlsx')
)
genes_raw <- genes_raw[-1, ]  # remove placeholder header row

# POT1 and DLG2 are absent from this dataset and are dropped.
genes_of_interest <- c('ATRX', 'DAXX', 'TERT', 'TSC2', 'MEN1', 'MET', 'KRAS', 'VHL')
gene_cols         <- c('ATRX_DAXX_trunc', 'TERT_mod', 'TSC2', 'MEN1', 'MET', 'KRAS', 'VHL')

# ATRX/DAXX and TERT: flag only rows where TMM_associated_mut_summary matches the
# expected class label. All other genes: flag by presence in the gene column.
genes_long <- genes_raw %>%
  filter(gene %in% genes_of_interest) %>%
  select(sample = `sample...1`, gene, tmm = TMM_associated_mut_summary) %>%
  mutate(
    mut_col = case_when(
      gene %in% c('ATRX', 'DAXX') ~ 'ATRX_DAXX_trunc',
      gene == 'TERT'               ~ 'TERT_mod',
      TRUE                         ~ gene
    ),
    present = case_when(
      gene %in% c('ATRX', 'DAXX') ~ tmm == 'ATRX_DAXX_trunc',
      gene == 'TERT'               ~ tmm == 'TERT_mod',
      TRUE                         ~ TRUE
    )
  ) %>%
  filter(present) %>%
  distinct(sample, mut_col)

# All samples present in the genes dataset: used to distinguish "no mutation found"
# (FALSE) from "sample not in genes dataset" (NA) after the join.
genes_samples <- genes_raw %>% distinct(sample = `sample...1`)

genes_wide <- genes_long %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    id_cols     = sample,
    names_from  = mut_col,
    values_from = present,
    values_fill = FALSE
  ) %>%
  right_join(genes_samples, by = 'sample') %>%
  mutate(across(all_of(gene_cols), ~ replace_na(.x, FALSE))) %>%
  select(sample, all_of(gene_cols))

# Join gene mutation columns to primary dataset.
# Samples absent from the genes dataset retain NA (unknown), not FALSE (wild-type).
data_primary <- readxl::read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))

data_primary_genes <- data_primary %>%
  left_join(genes_wide, by = join_by(donor_id == sample))

# Samples in primary but not in genes dataset
missing_from_genes <- data_primary_genes %>%
  filter(is.na(ATRX_DAXX_trunc)) %>%
  pull(donor_id)
cat('Primary samples absent from genes dataset (n =', length(missing_from_genes), '):\n')
print(missing_from_genes)

library(writexl)
write_xlsx(data_primary_genes, here('data', 'processed', 'PCAWG_primary+genes.xlsx'))


