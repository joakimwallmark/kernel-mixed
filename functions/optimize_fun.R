source("functions/get_item_resp_from_total.R")
optimize_fun <- function(weights, true_probs, totals, threshold = 0.01) {
  curr_w <- weights
  step <- 0.3
  
  while (T) {
    resp_mat <- colMeans(get_item_resp_from_total(totals, curr_w))
    for (w in 1:length(curr_w)) {
      if (resp_mat[w] < true_probs[w]) {
        curr_w[w] <- curr_w[w] + step
      }
      else {
        curr_w[w] <- max(curr_w[w] - step, 0.01)
      }
    }
    step <- step * 0.75
    if (step < threshold) break
  }
  
  return(curr_w)
}