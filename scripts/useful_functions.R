# Compare two potentially duplicated columns while handling numeric coercion,
# missing values, and floating-point tolerance. Discordant rows are returned for
# inspection rather than resolved automatically.
compare_columns <- function(data, col1, col2, id_col = 'patient_id', tol = 1e-8) {
  
  if (!all(c(col1, col2, id_col) %in% names(data))) {
    stop('One or more of the specified columns do not exist in the dataset.')
  }
  
  x <- data[[col1]]
  y <- data[[col2]]
  
  if (is.numeric(x) && !is.numeric(y)) {
    y_num <- suppressWarnings(as.numeric(y))
    if (!all(is.na(y) == is.na(y_num))) {
      warning(sprintf("'%s' is not numeric and cannot be cleanly converted: comparing as text.", col2))
      x_comp <- as.character(x)
      y_comp <- as.character(y)
      mode_type <- 'character'
    } else {
      x_comp <- x
      y_comp <- y_num
      mode_type <- 'numeric'
    }
    
  } else if (!is.numeric(x) && is.numeric(y)) {
    x_num <- suppressWarnings(as.numeric(x))
    if (!all(is.na(x) == is.na(x_num))) {
      warning(sprintf("'%s' is not numeric and cannot be cleanly converted: comparing as text.", col1))
      x_comp <- as.character(x)
      y_comp <- as.character(y)
      mode_type <- 'character'
    } else {
      x_comp <- x_num
      y_comp <- y
      mode_type <- 'numeric'
    }
    
  } else if (is.numeric(x) && is.numeric(y)) {
    x_comp <- x
    y_comp <- y
    mode_type <- 'numeric'
    
  } else {
    x_comp <- as.character(x)
    y_comp <- as.character(y)
    mode_type <- 'character'
  }
  
  if (mode_type == 'numeric') {
    equal_vec <- ifelse(
      is.na(x_comp) & is.na(y_comp),
      TRUE,
      ifelse(
        is.na(x_comp) | is.na(y_comp),
        FALSE,
        abs(x_comp - y_comp) <= tol
      )
    )
  } else {
    equal_vec <- ifelse(
      is.na(x_comp) & is.na(y_comp),
      TRUE,
      ifelse(
        is.na(x_comp) | is.na(y_comp),
        FALSE,
        trimws(x_comp) == trimws(y_comp)
      )
    )
  }
  
  diff_rows <- data[!equal_vec, c(id_col, col1, col2), drop = FALSE]
  
  if (all(equal_vec)) {
    print(TRUE)
    return(invisible(TRUE))
  } else {
    print(FALSE)
    View(diff_rows)
    return(invisible(diff_rows))
  }
}
