# Standardize the two source schemas before downstream analyses refer to columns
# by name. These functions intentionally preserve values and column order except
# for the explicit removal of a confirmed duplicate in `clean_drivers_data()`.

clean_drivers_data <- function(data) {
  
  colnames(data)[colnames(data) == 'Patient_ID_IF_mod'] <- 'patient_id'
  
  colnames(data)[colnames(data) == 'cancer_type...2'] <- 'cancer_type'
  
  colnames(data)[colnames(data) ==
                   'Tumor_Telomere fusion rate (Gt-end/coverage*purity)...5'] <- 'tumor_tf_rate'
  
  colnames(data)[colnames(data) ==
                   'Total TVR_normalized_by_all-reads'] <- 'tvr_all'
  
  colnames(data)[colnames(data) ==
                   'Total TVR_normalized_by_intratel-reads'] <- 'tvr_intratel'
  
  colnames(data)[colnames(data) == 'cancer_type_code...3'] <- 'cancer_type_code'
  
  colnames(data)[colnames(data) == 'Mixed_def...4'] <- 'mixed_def'
  
  colnames(data)[colnames(data) == 'Blood_Telomere fusion rate (Gt-end/coverage*purity)'] <- 'blood_tf_rate'
  
  colnames(data)[colnames(data) == 'cancer_type...28'] <- 'cancer_type2'
  
  colnames(data)[colnames(data) == 'cancer_type_code...29'] <- 'cancer_type_code2'
  
  colnames(data)[colnames(data) == 'Mixed_def...30'] <- 'mixed_def2'
  
  data$`Tumor_Telomere fusion rate (Gt-end/coverage*purity)...32` <- NULL
  
  colnames(data)[colnames(data) == 'TTGGGG_singletons_norm_by_all_reads_ORDER...42'] <- 'TTGGGG_singletons_norm_by_all_reads_ORDER'
  
  colnames(data)[colnames(data) == 'TTGGGG_singletons_norm_by_all_reads_ORDER...62'] <- 'TTGGGG_singletons_norm_by_intratel_ORDER'
  
  colnames(data)[colnames(data) == 'CATGGG_singletons_norm_by_all_reads_ORDER...50'] <- 'CATGGG_singletons_norm_by_all_reads_ORDER'
  
  colnames(data)[colnames(data) == 'CATGGG_singletons_norm_by_all_reads_ORDER...70'] <- 'CATGGG_singletons_norm_by_intratel_ORDER'
  
  return(data)
}

clean_pcawg_data <- function(data) {
  names(data) <- c(
    'donor_id',
    'cancer_type',
    'alt_status',
    'tf_rate',
    'telomere_insertion_rate',
    'telomere_content_log2',
    'TERT_FPKM',
    'ATAGGG_singleton_dist',
    'CTAGGG_singleton_dist',
    'GTAGGG_singleton_dist',
    'TAAGGG_singleton_dist',
    'TCAGGG_singleton_dist',
    'TGAGGG_singleton_dist',
    'TTCGGG_singleton_dist',
    'TTGGGG_singleton_dist',
    'TTTGGG_singleton_dist'
  )

  # Coerce measurement columns explicitly because Excel imports can mix numeric
  # values with literal missing-value strings.
  numeric_cols <- c(
    'tf_rate',
    'telomere_insertion_rate',
    'telomere_content_log2',
    'TERT_FPKM',
    'ATAGGG_singleton_dist',
    'CTAGGG_singleton_dist',
    'GTAGGG_singleton_dist',
    'TAAGGG_singleton_dist',
    'TCAGGG_singleton_dist',
    'TGAGGG_singleton_dist',
    'TTCGGG_singleton_dist',
    'TTGGGG_singleton_dist',
    'TTTGGG_singleton_dist'
  )

  # Treat source-level missing-value tokens as missing before numeric conversion.
  data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
    x <- trimws(as.character(x))
    x[x %in% c('NA', '', 'NaN')] <- NA
    as.numeric(x)
  })

  # Fix the ALT reference level so model coefficients consistently describe
  # ALT-high relative to ALT-low.
  data$alt_status  <- factor(data$alt_status, levels = c('ALT-low', 'ALT-high'))
  data$cancer_type <- factor(data$cancer_type)

  data
}
