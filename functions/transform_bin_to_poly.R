transform_bin_to_poly <- function(test, poly_col_indices, testname) {
  res <- matrix(NA, nrow = nrow(test), ncol = 0)
  # Add poly items to res
  for (poly in poly_col_indices) {
    res <- cbind(res, rowSums(test[, poly]))
  }
  # Add remaining bin items to res
  res <- cbind(res, test[, -unlist(poly_col_indices)])
  colnames(res) <- paste0(testname, 1:ncol(res))
  return(res)
}
