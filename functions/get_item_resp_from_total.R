get_item_resp_from_total <- function(totals, w_vec) {
  res <- matrix(0, nrow = length(totals), ncol = length(w_vec))
  for (i in seq_along(totals)) {
    correct_indices <- sample(seq_along(w_vec), totals[i], prob = w_vec)
    answer_set <- rep(0, length(w_vec))
    answer_set[correct_indices] <- 1
    res[i, ] <- answer_set
  }
  return(res)
}
