library(mirt)

# actual simulations ------------------------------------------------------
load("data/mock_data.RData")
# get IRT coefficients from real data
real_x_m <- mirt(data = X_14, model = 1, itemtype = "gpcm", SE = T)
real_y_m <- mirt(data = Y_13, model = 1, itemtype = "gpcm", SE = T)

x_irt_coefs <- coef(real_x_m, IRTpars = T, simplify = T)$items
y_irt_coefs <- coef(real_y_m, IRTpars = T, simplify = T)$items

# One scenario with original 2014 data for X, Y and A
# Last 30 items are anchor items
x_coef <- y_coef <- x_irt_coefs[1:60, ]
a_coef <- x_irt_coefs[61:90, ]

# Add polytomous from Y to make more polytomous ---------------------------
# 30 points, 15 on anchor
sum(sapply(apply(y_poly_coefs[1:10, 2:ncol(y_poly_coefs[1:10, ])], MARGIN = 1, na.omit), length))
sum(sapply(apply(y_poly_coefs[11:15, 2:ncol(y_poly_coefs[11:15, ])], MARGIN = 1, na.omit), length))
# remove 30/15 random items (keep 20/10) from the dichotomous on X/A to keep test length
set.seed(1243)
bin_xy <- x_dich_coefs[sample(1:50, 20), ]
bin_a <- x_dich_coefs[sample(51:75, 10), ]
mean(bin_a[, 2])
mean(x_dich_coefs[51:75, 2])
mean(bin_xy[, 2])
mean(x_dich_coefs[1:50, 2])
x_poly_coef <- y_poly_coef <- rbind(bin_xy, x_poly_coefs[1:10, ], y_poly_coefs[1:10, ])
a_poly_coef <- rbind(bin_a, x_poly_coefs[11:15, ], y_poly_coefs[11:15, ])


# Create easier Y ---------------------------------------------------------
# get easier items all forms sorted by diffculty
all_binary <- rbind(x_dich_coefs, y_dich_coefs)
all_polytomous <- rbind(x_poly_coefs, y_poly_coefs)
# get 50% easiest
bin_easiest <- sort(all_binary[, 2])[1:100]
# get mean b for each typ of polytomous items (difficulty)
# sort and find from 67% easiest items
pol_3cat <- sort(rowMeans(all_polytomous[c(1:3, 6:8, 11:13, 16:18, 21:23, 26:28), 2:3]))[1:12]
pol_5cat <- sort(rowMeans(all_polytomous[c(4, 9, 14, 19, 24, 29), 2:5]))[1:4]
pol_6cat <- sort(rowMeans(all_polytomous[c(5, 10, 15, 20, 25, 30), 2:6]))[1:4]
# sample from the easiest to create easier test form
set.seed(1242)
easy_y_bin <- sample(bin_easiest, 50)
easy_y_pol_bin <- sample(bin_easiest, 20)
easy_y_pol3 <- sample(pol_3cat, 6)
easy_y_pol_pol3 <- pol_3cat
easy_y_pol5 <- sample(pol_5cat, 2)
easy_y_pol_pol5 <- pol_5cat
easy_y_pol6 <- sample(pol_6cat, 2)
easy_y_pol_pol6 <- pol_6cat
easy_y_coef <- rbind(
  all_binary[names(easy_y_bin), ],
  all_polytomous[c(names(easy_y_pol3), names(easy_y_pol5), names(easy_y_pol6)), ]
)
easy_y_poly_coef <- rbind(
  all_binary[names(easy_y_pol_bin), ],
  all_polytomous[c(names(easy_y_pol_pol3), names(easy_y_pol_pol5), names(easy_y_pol_pol6)), ]
)

# Save parameters in scenario list for simulations ------------------------
# one list item per scenario with scenario parameters
item_pars <- list(
  list(x_pars = x_coef, y_pars = y_coef, a_pars = a_coef),
  list(x_pars = x_coef, y_pars = y_coef, a_pars = a_poly_coef),
  list(x_pars = x_coef, y_pars = easy_y_coef, a_pars = a_coef),
  list(x_pars = x_coef, y_pars = easy_y_coef, a_pars = a_poly_coef),
  list(x_pars = x_poly_coef, y_pars = y_poly_coef, a_pars = a_coef),
  list(x_pars = x_poly_coef, y_pars = y_poly_coef, a_pars = a_poly_coef),
  list(x_pars = x_poly_coef, y_pars = easy_y_poly_coef, a_pars = a_coef),
  list(x_pars = x_poly_coef, y_pars = easy_y_poly_coef, a_pars = a_poly_coef)
)
no_bin_poly <- list(
  c(50, 10, 25, 5),
  c(50, 10, 10, 10),
  c(50, 10, 25, 5),
  c(50, 10, 10, 10),
  c(20, 20, 25, 5),
  c(20, 20, 10, 10),
  c(20, 20, 25, 5),
  c(20, 20, 10, 10)
)
# vector indicating whether harder Y test is used in given scenario
easier_y <- c(F, F, T, T, F, F, T, T)

save(item_pars, no_bin_poly, easier_y, file = "sim_real_item_parameters2.RData")
