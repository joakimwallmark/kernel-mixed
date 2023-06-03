source("functions/sim_parallel.R")
source("functions/generate_data_irt.R")
source("functions/accuracy_indices.R")
source("functions/irt_true_w.R")
source("functions/irtke_neat.R")
source("functions/irtke_eg.R")
source("functions/ke_neat.R")
source("functions/ke_eg.R")
source("functions/irtose_eg.R")
source("functions/irtose_neat.R")
source("functions/irt_probs.R")
options(scipen = 999)

load("sim_real_item_parameters.RData")

# functions for generating theta scores
theta_funs <- list(
  list(
    P = function(n) {
      rnorm(n)
    },
    Q = function(n) {
      rnorm(n)
    }
  ),
  list(
    P = function(n) {
      rnorm(n)
    },
    Q = function(n) {
      rnorm(n, mean = 0.5, sd = 1.2)
    }
  )
)
# scenarios with population thetas with corresponding densities for true IRT computation
thetas <- list(
  list(
    P = seq(-6, 6, length = 1e4),
    Q = seq(-6, 6, length = 1e4)
  ),
  list(
    P = seq(-6, 6, length = 1e4),
    Q = seq(0.5 - 1.2 * 6, 0.5 + 1.2 * 6, length = 1e4)
  )
)
theta_den <- list(
  list(
    P = dnorm(thetas[[1]]$P),
    Q = dnorm(thetas[[1]]$Q)
  ),
  list(
    P = dnorm(thetas[[2]]$P),
    Q = dnorm(thetas[[2]]$Q, mean = 0.5, sd = 1.2)
  )
)

# vector indicating whether pop differ in given scenario
diff_pop <- c(FALSE, TRUE)

true_eq_names <- c("KEtrue", "EEtrue")
method_names <- c("IRTOSE", "IRTOSENEAT", "IRTKE", "KE", "IRTKECE", "KECE", "IRTKEPSE", "KEPSE")

n <- 1500
r <- 1000

# loop through scenarios
set.seed(12)
start_time <- Sys.time()
for (pop_scen in seq_along(thetas)) {
  for (item_scen in seq_along(item_pars)) {
    no_bin <- no_bin_poly[[item_scen]][1]
    no_poly <- no_bin_poly[[item_scen]][2]
    no_bin_a <- no_bin_poly[[item_scen]][3]
    no_poly_a <- no_bin_poly[[item_scen]][4]
    filename <- paste("R", r, " n", n, " b", no_bin, "p", no_poly, " bA", no_bin_a, "pA", no_poly_a,
      " diffy ", easier_y[item_scen], " diffpop ", diff_pop[pop_scen], ".RData",
      sep = ""
    )
    filename_iter <- paste("data/irt/irt ", filename, sep = "")
    filename_err <- paste("data/irt/irt Error ", filename, sep = "")
    filename_true <- paste("data/irt/irt True ", filename, sep = "")
    print(filename)

    x_pars <- item_pars[[item_scen]]$x_pars
    y_pars <- item_pars[[item_scen]]$y_pars
    a_pars <- item_pars[[item_scen]]$a_pars
    true <- irt_true_w(x_pars, y_pars, thetas[[pop_scen]], theta_den[[pop_scen]])
    save(true, file = filename_true)

    names <- paste("I", 1:(no_bin + no_poly), sep = "")
    names_a <- paste("A", 1:(no_bin_a + no_poly_a), sep = "")
    cat_x <- sapply(apply(x_pars[, 2:ncol(x_pars)], MARGIN = 1, na.omit), length) + 1
    cat_y <- sapply(apply(y_pars[, 2:ncol(y_pars)], MARGIN = 1, na.omit), length) + 1
    cat_a <- sapply(apply(a_pars[, 2:ncol(a_pars)], MARGIN = 1, na.omit), length) + 1
    max_score_xy <- sum(cat_x - 1)
    max_score_a <- sum(cat_a - 1)
    simres <- sim_parallel(
      iter = r, 
      data_gen_fun = generate_data_irt, 
      data_gen_fun_args = list(
         n, 
         n,
         theta_funs = theta_funs[[pop_scen]],
         x_pars = x_pars, 
         y_pars = y_pars, 
         a_pars = a_pars,
         names = names, names_a = names_a,
         filename = filename
       ), 
      method_funs = list(irtose_eg, irtose_neat, irtke_eg, ke_eg, irtke_neat, ke_neat, irtke_neat, ke_neat),
      methods_funs_args = list(
       list(max_score_xy, cat_x = cat_x, cat_y = cat_y),
       list(max_score_xy, max_score_a, cat_x = cat_x, cat_y = cat_y, cat_a = cat_a),
       list(max_score_xy, cat_x = cat_x, cat_y = cat_y),
       max_score_xy,
       list("CE", max_score_xy, max_score_a, cat_x = cat_x, cat_y = cat_y, cat_a = cat_a),
       list("CE", max_score_xy, max_score_a),
       list("PSE", max_score_xy, max_score_a, cat_x = cat_x, cat_y = cat_y, cat_a = cat_a),
       list("PSE", max_score_xy, max_score_a)
      ),
      packages_vec = c("kequate", "mirt"),
      # export_vec = c("from_b_to_d"),
      export_vec = c(
        "from_b_to_d",
        "cmnom_kequate",
        "probpl_kequate",
        "polyprob_kequate",
        "eqcpoly_kequate",
        "adjgpcmmirt_kequate",
        "irtinput_kequate"
      ), 
      res_filename = filename_iter
    )

    res <- accuracy_indices(simres, true[-3], true$W, true_eq_names, method_names)
    save(res, file = filename_err)
  }
}
