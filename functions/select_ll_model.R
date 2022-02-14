# select log-linear model based on AIC or BIC
select_ll_model <- function(neatk_X, neatk_Y, max_power = 6, max_interaction = 3, use_BIC = F) {
  max_x <- max_power
  max_a <- max_power
  max_xa <- max_interaction
  model <- c("frequency~I(X)+I(X^2)+I(A)+I(A^2)")
  lowestx <- glm(model, family = "poisson", data = neatk_X, x = TRUE)
  lowesty <- glm(model, family = "poisson", data = neatk_Y, x = TRUE)
  for (exp_A in 3:max_a) {
    uniquesA <- paste0("I(A^", 1:exp_A, ")")
    for (exp_X in 3:max_x) {
      uniquesX <- paste0("I(X^", 1:exp_X, ")")
      for (exp_XA_X in 1:max_xa) {
        for (exp_XA_A in 1:max_xa) {
          combosXA <- expand.grid(uniquesX[1:exp_XA_X], uniquesA[1:exp_XA_A])
          combosXA <- paste(combosXA$Var1, combosXA$Var2, sep = ":")
          model <- reformulate(c(uniquesX, uniquesA, combosXA), response = "frequency")

          xmod <- glm(model, family = "poisson", data = neatk_X, x = TRUE)
          ymod <- glm(model, family = "poisson", data = neatk_Y, x = TRUE)
          if (use_BIC) {
            if (BIC(lowestx) > BIC(xmod)) {
              lowestx <- xmod
            }
            if (BIC(lowesty) > BIC(ymod)) {
              lowesty <- ymod
            }
          }
          else {
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
