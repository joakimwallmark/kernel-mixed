library(gss)
source("functions/get_item_resp_from_total.R")
source("functions/transform_bin_to_poly.R")

generate_data_splines <- function(i, arg_list) {
  n <- arg_list[[1]]
  cor_xa <- arg_list[[2]][1] # correlation between X and A
  cor_ya <- arg_list[[2]][2]
  x_den <- arg_list[[3]]$X
  y_den <- arg_list[[3]]$Y
  xa_den <- arg_list[[3]]$XA
  ya_den <- arg_list[[3]]$YA
  poly_indices <- arg_list[[4]]
  w_x <- arg_list$w_x
  w_ax <- arg_list$w_ax
  w_y <- arg_list$w_y
  w_ay <- arg_list$w_ay
  filename <- arg_list$filename

  # Cov matrix
  # X, Y, AP, AQ is order in matrix
  cov_mat_xa <- matrix(c(
    1, cor_xa,
    cor_xa, 1
  ), ncol = 2, byrow = TRUE)
  cov_mat_ya <- matrix(c(
    1, cor_ya,
    cor_ya, 1
  ), ncol = 2, byrow = TRUE)
  # random from correlated std Normals
  norm_xa <- MASS::mvrnorm(n, c(0, 0), cov_mat_xa)
  norm_ya <- MASS::mvrnorm(n, c(0, 0), cov_mat_ya)
  p_norm_xa <- pnorm(norm_xa)
  p_norm_ya <- pnorm(norm_ya)
  # get corresponding values in splines densities
  x <- round(qssden(x_den, p_norm_xa[, 1]))
  xa <- round(qssden(xa_den, p_norm_xa[, 2]))
  ya <- round(qssden(ya_den, p_norm_ya[, 2]))
  if (identical(xa_den, ya_den)) {
    y <- round(qssden(y_den, p_norm_ya[, 1]))
  } else { # see paper for explanation. P!=Q
    y <- round(qssden(y_den, pssden(xa_den, qssden(ya_den, p_norm_ya[, 1]))))
  }

  x <- get_item_resp_from_total(x, w_x)
  xa <- get_item_resp_from_total(xa, w_ax)
  y <- get_item_resp_from_total(y, w_y)
  ya <- get_item_resp_from_total(ya, w_ay)

  x <- transform_bin_to_poly(x, poly_indices$XY, "X")
  xa <- transform_bin_to_poly(xa, poly_indices$A, "A")
  y <- transform_bin_to_poly(y, poly_indices$XY, "Y")
  ya <- transform_bin_to_poly(ya, poly_indices$A, "A")

  colnames(x) <- paste(rep("X", ncol(x)), seq_along(x), sep = "")
  colnames(xa) <- paste(rep("A", ncol(xa)), seq_along(xa), sep = "")
  colnames(y) <- paste(rep("Y", ncol(y)), seq_along(y), sep = "")
  colnames(ya) <- paste(rep("A", ncol(ya)), seq_along(ya), sep = "")

  result <- list(X = x, XA = xa, Y = y, YA = ya)
  save(result, file = paste("data/spline datasets/splines i", i, " ", filename, sep = ""))
  return(result)
}
