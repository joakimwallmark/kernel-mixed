library(kequate)
library(equate)


# computes "true" equating based on KE and EE from IRT model
# thetas is list of theta values for each pop. thetasdens is list of corresponding dens
irt_true_w <- function(x_pars, y_pars, thetas, thetasdens) {
  x_pars_b <- apply(x_pars[, 2:ncol(x_pars), drop = FALSE], MARGIN = 1, na.omit)
  y_pars_b <- apply(y_pars[, 2:ncol(y_pars), drop = FALSE], MARGIN = 1, na.omit)
  if (class(x_pars_b) != "list") {
    if (class(x_pars_b) != "vector") {
      x_pars_b <- as.list(x_pars_b)
      y_pars_b <- as.list(y_pars_b)
    }
    else {
      x_pars_b <- split(t(x_pars_b), 1:ncol(x_pars_b))
      y_pars_b <- split(t(y_pars_b), 1:ncol(y_pars_b))
    }
  }
  cat_x <- sapply(x_pars_b, length) + 1
  cat_y <- sapply(y_pars_b, length) + 1

  r <- list()
  s <- list()
  for (pop_i in 1:2) { # get P(X=x) and P(Y=y) for both populations
    p_poly_x <- polyprob(x_pars[, 1], x_pars_b,
      cat_x,
      model = "GPCM", thetas[[pop_i]]
    )
    p_poly_y <- polyprob(y_pars[, 1], y_pars_b,
      cat_y,
      model = "GPCM", thetas[[pop_i]]
    )

    # Average over pop.
    r[[pop_i]] <- rowSums(cmnom2(cat_x, p_poly_x, thetas[[pop_i]], thetasdens[[pop_i]]))
    s[[pop_i]] <- rowSums(cmnom2(cat_y, p_poly_y, thetas[[pop_i]], thetasdens[[pop_i]]))
  }
  # average to get P(X=x) for T
  r_t <- 0.5 * r[[1]] + 0.5 * r[[2]]
  s_t <- 0.5 * s[[1]] + 0.5 * s[[2]]

  # Kernel equating r och s över T
  keq <- kequate(
    design = "EG", x = 0:(length(r_t) - 1), y = 0:(length(s_t) - 1),
    r = r_t, s = s_t, DMP = diag(length(r_t)), DMQ = diag(length(s_t)),
    N = 1000, M = 1000, kernel = "gaussian"
  )

  # EE using true probs (multipled to avoid errors because of small probs)
  eg_x <- freqtab(0, 0:(length(r_t) - 1))
  eg_y <- freqtab(0, 0:(length(s_t) - 1))
  eg_x[] <- r_t * 1e10
  eg_y[] <- s_t * 1e10
  eg_eq <- equate(eg_x, eg_y, type = "equipercentile", method = "none", smoothmethod = "none")

  return(list(IRTKE = keq@equating$eqYx, EE = eg_eq$con[[2]], W = r_t))
}
