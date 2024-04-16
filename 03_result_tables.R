# Splines -----------------------------------------------------------------
load(file = "fitted_splines.RData")
n <- 1500
R <- 1000
bin_items <- c(50)
poly_items <- c(10)
bin_items_a <- c(25)
poly_items_a <- c(5)

harder_y <- c(FALSE, TRUE) # for filenames
diff_pop <- c(FALSE, TRUE) # for filenames

method_names <- c("IRTOSE", "IRTOSENEAT", "IRTKE", "KE", "IRTKECE", "KECE", "IRTKEPSE", "KEPSE")
mat_rows <- length(harder_y) * length(diff_pop) * length(method_names)
res_mat <- data.frame(matrix("", nrow = mat_rows, ncol = 9))
res_mat[, 4:9] <- sapply(res_mat[, 4:9], as.numeric)
rowi <- 1:8
# loop through scenarios
for (item_scen_A in seq_along(bin_items_a)) {
  for (item_scen in seq_along(bin_items)) {
    for (pop_scen in seq_along(diff_pop)) {
      for (test_scen in harder_y) {
        no_bin <- bin_items[item_scen]
        no_poly <- poly_items[item_scen]
        no_bin_a <- bin_items_a[item_scen_A]
        no_poly_a <- poly_items_a[item_scen_A]
        filename_err <- paste("data/spline/spline Error R", R, " n", n, " b", no_bin, "p", no_poly,
          " bA", no_bin_a, "pA", no_poly_a,
          " diffy ", test_scen, " diffpop ", diff_pop[pop_scen], ".RData",
          sep = ""
        )
        load(filename_err)
        res_mat[rowi, 1] <- method_names
        if (test_scen) {
          res_mat[rowi, 2] <- "x"
        }
        if (diff_pop[pop_scen]) {
          res_mat[rowi, 3] <- "x"
        }
        res_mat[rowi, 4] <- abs(res$spline$global$bias_global)
        res_mat[rowi, 5] <- abs(res$spline$global$bias_wglobal)
        res_mat[rowi, 6] <- res$spline$global$SEE_global
        res_mat[rowi, 7] <- res$spline$global$SEE_wglobal
        res_mat[rowi, 8] <- res$spline$global$RMSE_global
        res_mat[rowi, 9] <- res$spline$global$RMSE_wglobal

        rowi <- rowi + 8
      }
    }
  }
}

library(xtable)
print(xtable(res_mat), include.rownames = FALSE, include.colnames = FALSE)


# IRT ---------------------------------------------------------------------
n <- 1500
R <- 1000

no_bin_poly <- c(50, 10, 25, 5)
easier_y <- c(FALSE, TRUE)
diff_pop <- c(FALSE, TRUE)

method_names <- c("IRTOSE", "IRTOSENEAT", "IRTKE", "KE", "IRTKECE", "KECE", "IRTKEPSE", "KEPSE")
mat_rows <- length(easier_y) * length(diff_pop) * length(method_names)
res_mat <- data.frame(matrix("", nrow = mat_rows, ncol = 9))
res_mat[, 4:9] <- sapply(res_mat[, 4:9], as.numeric)
res_mat_ee <- res_mat
rowi <- 1:8

for (pop_scen in 1:length(diff_pop)) {
  for (test_scen in easier_y) {
    no_bin <- no_bin_poly[1]
    no_poly <- no_bin_poly[2]
    no_bin_a <- no_bin_poly[3]
    no_poly_a <- no_bin_poly[4]
    filename_err <- paste("data/irt/irt Error R", R, " n", n, " b", no_bin, "p", no_poly,
      " bA", no_bin_a, "pA", no_poly_a,
      " diffy ", test_scen, " diffpop ", diff_pop[pop_scen], ".RData",
      sep = ""
    )
    load(filename_err)
    res_mat[rowi, 1] <- method_names
    if (test_scen) {
      res_mat[rowi, 2] <- "x"
    }
    if (diff_pop[pop_scen]) {
      res_mat[rowi, 3] <- "x"
    }
    res_mat[rowi, 4] <- abs(res$KEtrue$global$bias_global)
    res_mat[rowi, 5] <- abs(res$KEtrue$global$bias_wglobal)
    res_mat[rowi, 6] <- res$KEtrue$global$SEE_global
    res_mat[rowi, 7] <- res$KEtrue$global$SEE_wglobal
    res_mat[rowi, 8] <- res$KEtrue$global$RMSE_global
    res_mat[rowi, 9] <- res$KEtrue$global$RMSE_wglobal

    res_mat_ee[rowi, 1] <- method_names
    if (test_scen) {
      res_mat_ee[rowi, 2] <- "x"
    }
    if (diff_pop[pop_scen]) {
      res_mat_ee[rowi, 3] <- "x"
    }
    res_mat_ee[rowi, 4] <- abs(res$EEtrue$global$bias_global)
    res_mat_ee[rowi, 5] <- abs(res$EEtrue$global$bias_wglobal)
    res_mat_ee[rowi, 6] <- res$EEtrue$global$SEE_global
    res_mat_ee[rowi, 7] <- res$EEtrue$global$SEE_wglobal
    res_mat_ee[rowi, 8] <- res$EEtrue$global$RMSE_global
    res_mat_ee[rowi, 9] <- res$EEtrue$global$RMSE_wglobal

    rowi <- rowi + 8
  }
}

library(xtable)
print(xtable(res_mat), include.rownames = FALSE, include.colnames = FALSE)
print(xtable(res_mat_ee), include.rownames = FALSE, include.colnames = FALSE)
