# IRT kernel equating from X to Y
# i is iteration index from simulations for errors
# 2 items in argList: max_score, categories vector
ke_eg <- function(data, i, max_score) {
  x <- data$X
  y <- data$Y

  ke <- tryCatch(
    {
      ke_eg_x <- kefreq(in1 = rowSums(x), xscores = 0:max_score)
      ke_eg_y <- kefreq(in1 = rowSums(y), xscores = 0:max_score)
      # select log-linear model based on AIC
      model <- c("frequency~I(X)+I(X^2)")
      lowestx <- glm(model, family = "poisson", data = ke_eg_x, x = TRUE)
      lowesty <- glm(model, family = "poisson", data = ke_eg_y, x = TRUE)
      for (expon in 3:5) {
        model <- paste(model, "+I(X^", expon, ")", sep = "")
        xmod <- glm(model, family = "poisson", data = ke_eg_x, x = TRUE)
        ymod <- glm(model, family = "poisson", data = ke_eg_y, x = TRUE)
        if (lowestx$aic > xmod$aic) {
          lowestx <- xmod
        }
        if (lowesty$aic > ymod$aic) {
          lowesty <- ymod
        }
      }

      ke <- kequate(design = "EG", x = 0:max_score, y = 0:max_score, r = lowestx, s = lowesty, kernel = "gaussian")
    },
    error = function(e) {
      print(e)
      # save dataset in case of error
      save(e, x, y, max_score,
        file = paste("data/KE error ", i, ".RData", sep = "")
      )
    }
  )
  result <- c(ke@equating$eqYx, ke@equating$SEEYx, ke@PRE$PREYx)
  names(result) <- c(
    paste("eqYx", 0:max_score), paste("SEEYx", 0:max_score),
    paste("PREYx", 1:10)
  )
  return(result)
}
