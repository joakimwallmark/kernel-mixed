# IRTOSE from X to Y
# i is iteration index from simulations for errors
# 2 items in argList: max_score, categories vector
irtose_neat <- function(data, i, arg_list) {
  source("functions/IRTOSE.R")
  x <- cbind(data[[1]], data[[2]])
  y <- cbind(data[[3]], data[[4]])
  max_score <- arg_list[[1]]
  max_score_a <- arg_list[[2]]
  cat_x <- arg_list$cat_x
  cat_y <- arg_list$cat_y
  cat_a <- arg_list$cat_a
  

  res_names <- c(paste("eqYx", 0:max_score), paste("SEEYx", 0:max_score), paste("PREYx", 1:10))
  
  irtose <- tryCatch(
    {
      gpc_x <- mirt(data = x, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)
      gpc_y <- mirt(data = y, model = 1, itemtype = "gpcm", quadpts = 61, SE = T)
      
      irtose <- IRT_OSE_NEAT(
        gpc_x, gpc_y,
        0:max_score, 0:max_score, 0:max_score_a,
        catsX = cat_x, catsY = cat_y,
        catsA = cat_a, model = "GPCM"
      )
    },
    error = function(e) {
      print(e)
      # save dataset in case of error
      save(e, method, x, y, max_score, max_score_a, cat_x, cat_y, cat_a,
           file = paste("Data/IRTOSE error ", i, ".RData", sep = "")
      )
      result <- rep(0, length(res_names))
      names(result) <- res_names
      return(result)
    }
  )
  result <- c(irtose, rep(0, length(irtose)), rep(0, 10))
  names(result) <- res_names
  return(result)
}
