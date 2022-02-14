from_b_to_d <- function(a, b) {
  # d needs an extra column sof 0's
  d <- matrix(0, nrow = nrow(b), ncol = ncol(b) + 1)
  for (i in 2:ncol(d)) {
    d[, i] <- d[, i - 1] - b[, i - 1] * a
  }
  return(d)
}
