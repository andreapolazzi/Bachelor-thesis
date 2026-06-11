library(tidyverse)
library(readxl)
library(here)

met <- read_xlsx(here('data', 'processed', 'metastatic_red_edit_singleton_dist.xlsx')) %>%
  distinct(patient_id, .keep_all = TRUE) %>%
  select(patient_id, cancer_type, primary_tumor_type)

# Group derived from cancer_type (free-text, regex-based)
derive_from_cancer_type <- function(x) {
  case_when(
    str_detect(x, regex('melanoma', ignore_case = TRUE)) ~ 'Melanoma',
    str_detect(x, regex('sarcoma|gastrointestinal stromal', ignore_case = TRUE)) ~ 'Sarcoma',
    str_detect(x, regex('carcinoma', ignore_case = TRUE)) ~ 'Carcinoma',
    is.na(x) ~ NA_character_,
    TRUE ~ 'Other'
  )
}

# Group derived from primary_tumor_type (exact logic used in the plot chunk)
derive_from_ptt <- function(x) {
  case_when(
    x == 'Carcinoma' ~ 'Carcinoma',
    x == 'Melanoma'  ~ 'Melanoma',
    str_detect(x, regex('sarcoma|gastrointestinal stromal', ignore_case = TRUE)) ~ 'Sarcoma',
    is.na(x) ~ NA_character_,
    TRUE ~ 'Other'
  )
}

discordance <- met %>%
  mutate(
    group_from_cancer_type = derive_from_cancer_type(cancer_type),
    group_from_ptt         = derive_from_ptt(primary_tumor_type)
  ) %>%
  filter(
    is.na(group_from_cancer_type) != is.na(group_from_ptt) |
    (!is.na(group_from_cancer_type) & !is.na(group_from_ptt) &
       group_from_cancer_type != group_from_ptt)
  ) %>%
  select(patient_id, cancer_type, group_from_cancer_type, primary_tumor_type, group_from_ptt)

View(discordance)
cat(nrow(discordance), "discordant patients\n")
