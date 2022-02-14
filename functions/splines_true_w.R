library(gss)

# computes "true" equating based on KE and EE from IRT model
splines_true_w <- function(scores, x_den, y_den, xa_den = NULL, ya_den = NULL) {
  p_x <- pssden(x_den, scores)
  true <- qssden(y_den, p_x)
  # weights from X for weighted global indicies
  W <- pssden(x_den, scores + 0.5) - pssden(x_den, scores - 0.5)
  return(list(true = true, W = W))
}
