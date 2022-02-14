# get list of P(x|theta) for all thetas in qpoints
# list has length of all possible x
binprob <- function(a, b, c, qpoints) {
  items <- length(b)
  probs <- vector("list", items)
  for (i in 1:items) {
    out <- matrix(0, nrow = 2, ncol = length(qpoints))
    out[2, ] <- c[i] + (1 - c[i]) / (1 + exp(-a[i] * (qpoints - b[i])))
    out[1, ] <- 1 - out[2, ]
    probs[[i]] <- out
  }
  return(probs)
}

# get list of P(X_i|theta) for all thetas in qpoints
# each list element corresponds to one item and is a matrix with col=theta and row = item category
polyprob <- function(a, b, cats, model, qpoints) {
  JX <- length(cats)
  kX <- sum(cats) - JX
  probs <- vector("list", JX)
  if (model == "GPCM") {
    for (i in 1:JX) {
      out <- matrix(0, nrow = cats[i], ncol = length(qpoints))
      denom <- 0
      for (j in 1:(cats[i] - 1)) {
        temp <- 0
        for (l in 1:j) {
          temp <- a[i] * (qpoints - b[[i]][l]) + temp
        }
        denom <- exp(temp) + denom
      }
      out[1, ] <- 1 / (1 + denom)
      for (j in 1:(cats[i] - 1)) {
        numer <- exp(j * a[i] * qpoints - a[i] * sum(b[[i]][1:j]))
        out[j + 1, ] <- numer / (1 + denom)
      }
      probs[[i]] <- out
    }
  }
  if (model == "GRM") {
    for (i in 1:JX) {
      out <- matrix(0, nrow = cats[i], ncol = length(qpoints))
      out[1, ] <- 1 - 1 / (1 + exp(-a[i] * (qpoints - b[[i]][1])))
      out[cats[i], ] <- 1 / (1 + exp(-a[i] * (qpoints -
        b[[i]][cats[i] - 1])))
      for (j in 2:(cats[i] - 1)) {
        out[j, ] <- 1 / (1 + exp(-a[i] *
          (qpoints - b[[i]][j - 1]))) - 1 / (1 + exp(-a[i] *
          (qpoints - b[[i]][j])))
      }
      probs[[i]] <- out
    }
  }
  return(probs)
}

# compound multinomial, row sums of resulting matrix are P(x=row_index)
cmnom <- function(cats, probs, qpoints, tmu = 0, tsig = 1) {
  nsprobs <- probs_given_each_qpoint(cats, probs, qpoints)
  t(t(nsprobs) * (dnorm(qpoints, tmu, tsig) / sum(dnorm(qpoints, tmu, tsig))))
}

cmnom2 <- function(cats, probs, qpoints, dens_val) {
  nsprobs <- probs_given_each_qpoint(cats, probs, qpoints)
  t(t(nsprobs) * (dens_val / sum(dens_val)))
}

# Resulting matrix: P(x=row_index|theta=qpoints[col_index])
probs_given_each_qpoint <- function(cats, probs, qpoints) {
  JX <- length(cats)
  xK <- sum(cats) - JX
  nsprobs <- matrix(0, nrow = xK + 1, ncol = length(qpoints))
  nsprobs[1:(cats[1]), ] <- probs[[1]]
  if (JX > 1) {
    for (i in 2:JX) {
      maxsc <- sum(cats[1:i]) - length(cats[1:i])
      sprobs <- nsprobs
      nsprobs <- matrix(0, nrow = xK + 1, ncol = length(qpoints))
      for (j in 1:(maxsc - cats[i] + 2)) {
        nsprobs[j:(cats[i] + j - 1), ] <- t(sprobs[j, ] * t(probs[[i]])) +
          nsprobs[j:(cats[i] + j - 1), ]
      }
    }
  }
  nsprobs
}
