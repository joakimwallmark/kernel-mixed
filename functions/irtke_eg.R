source("functions/IRTOSE.R")
# IRT kernel equating from X to Y
# i is iteration index from simulations for errors
# 2 items in argList: max_score, categories vector
irtke_eg <- function(data, i, arg_list) {
  x <- data$X
  y <- data$Y
  max_score <- arg_list[[1]]
  cat_x <- arg_list$cat_x
  cat_y <- arg_list$cat_y

  irtke <- tryCatch(
    {
      gpc_x <- mirt(data = x, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)
      gpc_y <- mirt(data = y, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)

      irtke <- irtose(
        design = "EG", gpc_x, gpc_y, 0:max_score, 0:max_score,
        catsX = cat_x, catsY = cat_y, model = "GPCM"
      )
    },
    error = function(e) {
      print(e)
      # save dataset in case of error
      save(e, x, y, max_score, cat_x, cat_y,
        file = paste("Data/IRTKE error ", i, ".RData", sep = "")
      )
    }
  )
  result <- c(irtke@equating$eqYx, irtke@equating$SEEYx, irtke@PRE$PREYx)
  names(result) <- c(
    paste("eqYx", 0:max_score), paste("SEEYx", 0:max_score),
    paste("PREYx", 1:10)
  )
  return(result)
}
