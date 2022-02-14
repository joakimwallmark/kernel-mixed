source("functions/from_b_to_d.R")
generate_data_irt <- function(i, arg_list) {
  nP <- arg_list[[1]]
  nQ <- arg_list[[2]]
  theta_p <- arg_list$theta_funs$P(nP)
  theta_q <- arg_list$theta_funs$Q(nQ)
  # generate response for given scenario
  x_pars <- arg_list$x_pars
  y_pars <- arg_list$y_pars
  a_pars <- arg_list$a_pars
  filename <- arg_list$filename

  # transform to difficulty to d format to use simdata()
  dx <- from_b_to_d(x_pars[, 1], x_pars[, 2:ncol(x_pars), drop = F])
  dy <- from_b_to_d(y_pars[, 1], y_pars[, 2:ncol(x_pars), drop = F])
  da <- from_b_to_d(a_pars[, 1], a_pars[, 2:ncol(x_pars), drop = F])

  x_gpc <- simdata(x_pars[, 1], dx,
    itemtype = "gpcm", Theta = matrix(theta_p, ncol = 1)
  )
  y_gpc <- simdata(y_pars[, 1], dy,
    itemtype = "gpcm", Theta = matrix(theta_q, ncol = 1)
  )
  ap_gpc <- simdata(a_pars[, 1], da,
    itemtype = "gpcm", Theta = matrix(theta_p, ncol = 1)
  )
  aq_gpc <- simdata(a_pars[, 1], da,
    itemtype = "gpcm", Theta = matrix(theta_q, ncol = 1)
  )

  colnames(x_gpc) <- colnames(y_gpc) <- arg_list$names
  colnames(ap_gpc) <- colnames(aq_gpc) <- arg_list$names_a

  result <- list(X = x_gpc, AX = ap_gpc, Y = y_gpc, AY = aq_gpc)
  save(result, file = paste("data/irt datasets/irt i", i, " ", filename, sep = ""))

  return(result)
}
