# IRT kernel equating from X to Y
# i is iteration index from simulations for errors
# 2 items in argList: max_score, categories vector
irtke_neat <- function(data, i, arg_list) {
  x <- cbind(data[[1]], data[[2]])
  y <- cbind(data[[3]], data[[4]])
  method <- arg_list[[1]]
  max_score <- arg_list[[2]]
  max_score_a <- arg_list[[3]]
  cat_x <- arg_list$cat_x
  cat_y <- arg_list$cat_y
  cat_a <- arg_list$cat_a

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

  irtke <- tryCatch(
    {
      gpc_x <- mirt(data = x, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)
      gpc_y <- mirt(data = y, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)

      irtke <- irtose(
        design = method, gpc_x, gpc_y,
        0:max_score, 0:max_score, 0:max_score_a,
        catsX = cat_x, catsY = cat_y,
        catsA = cat_a, model = "GPCM"
      )
    },
    error = function(e) {
      print(e)
      # save dataset in case of error
      save(e, method, x, y, max_score, max_score_a, cat_x, cat_y, cat_a,
        file = paste("Data/IRTKE error ", i, ".RData", sep = "")
      )
      result <- rep(0, length(res_names))
      names(result) <- res_names
      return(result)
    }
  )
  if (method == "CE") {
    result <- c(irtke@equating$eqYx, irtke@equating$SEEYx, irtke@PRE$PREAx, irtke@PRE$PREYa)
    names(result) <- res_names
  } else {
    result <- c(irtke@equating$eqYx, irtke@equating$SEEYx, irtke@PRE$PREYx)
    names(result) <- res_names
  }
  return(result)
}
