# IRT kernel equating from X to Y
# i is iteration index from simulations for errors
ke_neat <- function(data, i, arg_list) {
  source("functions/select_ll_model.R") # foreach doesn't find this function in parallel so this needs to be here??...
  x <- cbind(data[[1]], data[[2]])
  y <- cbind(data[[3]], data[[4]])
  method <- arg_list[[1]]
  max_score <- arg_list[[2]]
  max_score_a <- arg_list[[3]]

  if (method == "CE") {
    res_names <- c(
      paste("eqYx", 0:max_score), paste("SEEYx", 0:max_score),
      paste("PREAx", 1:10), paste("PREYa", 1:10)
    )
  } else {
    res_names <- c(
      paste("eqYx", 0:max_score), paste("SEEYx", 0:max_score),
      paste("PREYx", 1:10)
    )
  }
  a_items_x <- grep("A", colnames(x))
  a_items_y <- grep("A", colnames(y))
  ke_eg_x <- kefreq(rowSums(x[, -a_items_x]), 0:max_score, rowSums(x[, a_items_x]), 0:max_score_a)
  ke_eg_y <- kefreq(rowSums(y[, -a_items_x]), 0:max_score, rowSums(y[, a_items_x]), 0:max_score_a)

  lowbic <- select_ll_model(ke_eg_x, ke_eg_y, use_BIC = T)
  glm_x <- lowbic[[1]]
  glm_y <- lowbic[[2]]

  ke <- tryCatch(
    {
      if (method == "CE") {
        ke <- kequate(
          design = paste("NEAT_", method, sep = ""),
          x = 0:max_score, y = 0:max_score, a = 0:max_score_a, P = glm_x, Q = glm_y, kernel = "gaussian"
        )
      }
      else {
        ke <- kequate(
          design = paste("NEAT_", method, sep = ""),
          x = 0:max_score, y = 0:max_score, P = glm_x, Q = glm_y, kernel = "gaussian"
        )
      }
    },
    error = function(e) {
      print(e)
      # save dataset in case of error, return NAs
      save(e, method, x, y, max_score, max_score_a, glm_x, glm_y,
        file = paste("Data/KE error ", i, ".RData", sep = "")
      )
      result <- rep(0, length(res_names))
      names(result) <- res_names
      return(result)
    }
  )
  if (method == "CE") {
    result <- c(ke@equating$eqYx, ke@equating$SEEYx, ke@PRE$PREAx, ke@PRE$PREYa)
    names(result) <- res_names
  }
  else {
    result <- c(ke@equating$eqYx, ke@equating$SEEYx, ke@PRE$PREYx)
    names(result) <- res_names
  }
  return(result)
}
