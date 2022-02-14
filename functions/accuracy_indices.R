# compute bias, SEE, RMSE, and global measures of the same
# sim_res is list of matrices from sim. iterations
accuracyIndices <- function(sim_res, true_eqs, W, true_eq_names, method_names) {
  no_methods <- length(sim_res)
  no_scores <- length(W)
  scores <- 0:(length(W) - 1)

  true_eq_list <- vector("list", length(true_eqs))
  for (true_eq in 1:length(true_eqs)) {
    local_list <- vector("list", no_methods)
    bias_global <- vector("numeric", no_methods)
    SEE_global <- vector("numeric", no_methods)
    RMSE_global <- vector("numeric", no_methods)
    bias_wglobal <- vector("numeric", no_methods)
    SEE_wglobal <- vector("numeric", no_methods)
    RMSE_wglobal <- vector("numeric", no_methods)
    names(bias_global) <- names(SEE_global) <- names(RMSE_global) <-
      names(bias_wglobal) <- names(SEE_wglobal) <- names(RMSE_wglobal) <- method_names
    for (i in 1:no_methods) {
      mateq <- sim_res[[i]][, scores + 1]
      local_bias <- colMeans(mateq) - true_eqs[[true_eq]]
      local_SEE <- apply(mateq, 2, sd)
      local_RMSE <- apply(
        t(t(mateq) - true_eqs[[true_eq]]), 2,
        function(diff) {
          sqrt(sum((diff)^2) / nrow(mateq))
        }
      )
      names(local_bias) <- names(local_SEE) <- names(local_RMSE) <- scores
      local_list[[i]] <- list("bias" = local_bias, "SEE" = local_SEE, "RMSE" = local_RMSE)

      bias_global[i] <- mean(abs(local_bias))
      SEE_global[i] <- mean(local_SEE)
      RMSE_global[i] <- mean(local_RMSE)
      bias_wglobal[i] <- sum(W * abs(local_bias))
      SEE_wglobal[i] <- sum(W * local_SEE)
      RMSE_wglobal[i] <- sum(W * local_RMSE)
    }
    global_list <- list(
      "bias_global" = bias_global, "SEE_global" = SEE_global,
      "RMSE_global" = RMSE_global, "bias_wglobal" = bias_wglobal,
      "SEE_wglobal" = SEE_wglobal, "RMSE_wglobal" = RMSE_wglobal
    )
    names(local_list) <- method_names

    true_eq_list[[true_eq]] <- list("local" = local_list, "global" = global_list)
  }
  names(true_eq_list) <- true_eq_names
  return(true_eq_list)
}
