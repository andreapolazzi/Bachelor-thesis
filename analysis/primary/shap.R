library(shapr)
library(shapviz)
library(treeshap)
library(randomForest)

data <- readxl::read_xlsx(here('data', 'processed', 'PCAWG_primary.xlsx'))
clean_data <- data %>% 
  drop_na() %>% 
  filter(tf_primary_rate<10)
clean_data <- clean_data %>% 
  mutate(
    tf_primary_rate = log1p(tf_primary_rate),
    tf_blood_rate = log1p(tf_blood_rate),
    telomere_insertion_rate = log1p(telomere_insertion_rate),
    TERT_FPKM = log1p(TERT_FPKM),
    alt_status = factor(
      alt_status,
      levels = c("ALT-low", "ALT-high")
    )
  ) %>%  
  select(-donor_id, -cancer_type, -Specimen.Type.Summary)

