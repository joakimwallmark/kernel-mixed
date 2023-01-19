library(gss)
source("functions/sim_parallel.R")
source("functions/get_item_resp_from_total.R")
source("functions/irtke_eg.R")
source("functions/irtke_neat.R")
source("functions/ke_eg.R")
source("functions/ke_neat.R")
source("functions/accuracy_indices.R")
source("functions/generate_data_splines.R")
source("functions/transform_bin_to_poly.R")
source("functions/splines_true_w.R")

options(scipen = 999)

load(file = "fitted_splines.RData")
load(file = "weights.RData") # weights for generating item scores from test scores
n <- 1500
r <- 1000
bin_items <- c(50, 35, 20)
poly_items <- c(10, 15, 20)
bin_items_a <- c(25, 10)
poly_items_a <- c(5, 10)

different_y <- c(F, T) # for filenames
diff_pop <- c(F, T) # for filenames
w_x <- list(opt_w_x)
w_y <- list(opt_w_y, opt_w_ey)
w_a <- list(opt_w_a, opt_w_ea)
method_names <- c("IRTKE", "KE", "IRTKECE", "KECE", "IRTKEPSE", "KEPSE")
cor_list <- list(c(0.85, 0.85)) # correlations between XA and YA

# loop through scenarios
set.seed(123)
start.time <- Sys.time()
for (item_scen_a in 1:length(bin_items_a)) {
  max_score_a <- a_list[[1]]$domain[2, ] # get max from spline densities
  no_bin_a <- bin_items_a[item_scen_a]
  no_poly_a <- poly_items_a[item_scen_a]
  # randomly assign which binary item indices should be used for each poly item
  # (used in generateDataSplinesNEAT for anchor items)
  pol_ind_a <- base::split(
    sample(1:max_score_a, no_poly_a * 3),
    1:no_poly_a
  )
  for (item_scen in 1:length(bin_items)) {
    max_score_xy <- x_list[[1]]$domain[2, ]
    no_bin <- bin_items[item_scen]
    no_poly <- poly_items[item_scen]
    # randomly assign which binary item indices should be used for each poly item
    # (used in generateDataSplinesNEAT for X and Y)
    pol_indices <- base::split(
      sample(1:max_score_xy, no_poly * 3),
      1:no_poly
    )
    for (pop_scen in 1:length(a_list)) { # a_list elements are various densities for anchors
      for (test_scen in 1:length(y_list)) { # a_list elements are various densities for anchors
        filename <- paste("R", r, " n", n, " b", no_bin, "p", no_poly, " bA", no_bin_a, "pA", no_poly_a,
          " diffy ", different_y[test_scen], " diffpop ", diff_pop[pop_scen], ".RData",
          sep = ""
        )
        filename_iter <- paste("data/spline/spline ", filename, sep = "")
        filename_err <- paste("data/spline/spline Error ", filename, sep = "")
        filename_true <- paste("data/spline/spline True ", filename, sep = "")
        print(filename)
        true <- splines_true_w(
          scores = 0:max_score_xy, x_list[[1]], y_list[[test_scen]],
          a_list[[1]], a_list[[pop_scen]]
        )
        save(true, file = filename_true)
        cat_xy <- c(rep(4, no_poly), rep(2, no_bin))
        cat_a <- c(rep(4, no_poly_a), rep(2, no_bin_a))
        simres <- sim_parallel(
          iter = r, data_gen_FUN = generate_data_splines,
          data_gen_FUN_args = list(n,
            max = list(
              "X" = max_score_xy,
              "Y" = max_score_xy,
              "A" = max_score_a
            ),
            correl = cor_list[[1]],
            densities = list(
              X = x_list[[1]],
              Y = y_list[[test_scen]], # easier Y is in 2
              XA = a_list[[1]],
              YA = a_list[[pop_scen]]
            ), # easier A is in 2
            polind = list(XY = pol_indices, A = pol_ind_a),
            w_x = w_x[[1]], w_y = w_y[[test_scen]], w_ax = w_a[[1]], w_ay = w_a[[pop_scen]],
            filename = filename
          ),
          method_FUNs = list(irtke_eg, ke_eg, irtke_neat, ke_neat, irtke_neat, ke_neat),
          method_FUNs_args = list(
            list(max_score_xy, cat_x = cat_xy, cat_y = cat_xy),
            max_score_xy,
            list("CE", max_score_xy, max_score_a, cat_x = cat_xy, cat_y = cat_xy, cat_a = cat_a),
            list("CE", max_score_xy, max_score_a),
            list("PSE", max_score_xy, max_score_a, cat_x = cat_xy, cat_y = cat_xy, cat_a = cat_a),
            list("PSE", max_score_xy, max_score_a)
          ),

          packages_vec = c("kequate", "mirt", "gss"),
          export_vec = c("get_item_resp_from_total", "transform_bin_to_poly"),
          res_filename = filename_iter
        )
        res <- accuracyIndices(simres, list(true$true), true$W, "spline", method_names)
        save(res, file = filename_err)
      }
    }
  }
}
time <- Sys.time() - start.time
save(time, file = paste("data/spline/time splines R", r, ".RData", sep = ""))
print(time)
