# compute bias, SEE, RMSE, and global measures of the same
# sim_res is list of matrices from sim. iterations
accuracy_indices <- function(sim_res, true_eqs, w, true_eq_names, method_names) {
  no_methods <- length(sim_res)
  scores <- 0:(length(w) - 1)

  true_eq_list <- vector("list", length(true_eqs))
  for (true_eq in seq_along(true_eqs)) {
    local_list <- vector("list", no_methods)
    bias_global <- vector("numeric", no_methods)
    see_global <- vector("numeric", no_methods)
    rmse_global <- vector("numeric", no_methods)
    bias_wglobal <- vector("numeric", no_methods)
    see_wglobal <- vector("numeric", no_methods)
    rmse_wglobal <- vector("numeric", no_methods)
    names(bias_global) <- names(see_global) <- names(rmse_global) <-
      names(bias_wglobal) <- names(see_wglobal) <- names(rmse_wglobal) <- method_names
    for (i in 1:no_methods) {
      mateq <- sim_res[[i]][, scores + 1]
      local_bias <- colMeans(mateq) - true_eqs[[true_eq]]
      local_see <- apply(mateq, 2, sd)
      local_rmse <- apply(
        t(t(mateq) - true_eqs[[true_eq]]), 2,
        function(diff) {
          sqrt(sum((diff)^2) / nrow(mateq))
        }
      )
      names(local_bias) <- names(local_see) <- names(local_rmse) <- scores
      local_list[[i]] <- list("bias" = local_bias, "SEE" = local_see, "RMSE" = local_rmse)

      bias_global[i] <- mean(abs(local_bias))
      see_global[i] <- mean(local_see)
      rmse_global[i] <- mean(local_rmse)
      bias_wglobal[i] <- sum(w * abs(local_bias))
      see_wglobal[i] <- sum(w * local_see)
      rmse_wglobal[i] <- sum(w * local_rmse)
    }
    global_list <- list(
      "bias_global" = bias_global, "SEE_global" = see_global,
      "RMSE_global" = rmse_global, "bias_wglobal" = bias_wglobal,
      "SEE_wglobal" = see_wglobal, "RMSE_wglobal" = rmse_wglobal
    )
    names(local_list) <- method_names

    true_eq_list[[true_eq]] <- list("local" = local_list, "global" = global_list)
  }
  names(true_eq_list) <- true_eq_names
  return(true_eq_list)
}
