# select log-linear model based on AIC or BIC
select_ll_model <- function(neatk_x, neatk_y, max_power = 6, max_interaction = 3, use_bic = FALSE) {
  max_x <- max_power
  max_a <- max_power
  max_xa <- max_interaction
  model <- c("frequency~I(X)+I(X^2)+I(A)+I(A^2)")
  lowestx <- glm(model, family = "poisson", data = neatk_x, x = TRUE)
  lowesty <- glm(model, family = "poisson", data = neatk_y, x = TRUE)
  for (exp_A in 3:max_a) {
    uniques_a <- paste0("I(A^", 1:exp_A, ")")
    for (exp_X in 3:max_x) {
      uniques_x <- paste0("I(X^", 1:exp_X, ")")
      for (exp_XA_X in 1:max_xa) {
        for (exp_XA_A in 1:max_xa) {
          combos_xa <- expand.grid(uniques_x[1:exp_XA_X], uniques_a[1:exp_XA_A])
          combos_xa <- paste(combos_xa$Var1, combos_xa$Var2, sep = ":")
          model <- reformulate(c(uniques_x, uniques_a, combos_xa), response = "frequency")

          xmod <- glm(model, family = "poisson", data = neatk_x, x = TRUE)
          ymod <- glm(model, family = "poisson", data = neatk_y, x = TRUE)
          if (use_bic) {
            if (BIC(lowestx) > BIC(xmod)) {
              lowestx <- xmod
            }
            if (BIC(lowesty) > BIC(ymod)) {
              lowesty <- ymod
            }
          } else {
            if (AIC(lowestx) > AIC(xmod)) {
              lowestx <- xmod
            }
            if (AIC(lowesty) > AIC(ymod)) {
              lowesty <- ymod
            }
          }
        }
      }
    }
  }
  return(list(lowestx, lowesty))
}
