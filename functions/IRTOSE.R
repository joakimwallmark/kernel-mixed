library(equate)

# Parameter linking part is taken from kequate package PSE method
IRT_OSE_EG <- function(P, Q, x, y, a = 0, qpoints = seq(-6, 6, by = 0.1),
                       model = "2pl", catsX = 0, catsY = 0, catsA = 0,
                       wS = 0.5, eqcoef = "mean-mean",
                       robust = FALSE, distribution = list("normal", par = data.frame(mu = 0, sigma = 1))) {
  if (model == "2pl") {
    nX <- length(x) - 1
    nY <- length(y) - 1
    if (inherits(P, "ltm")) {
      bx <- as.numeric(coef.ltm(P)[, 1])[1:(length(x) - 
                                              1)]
      ax <- as.numeric(coef.ltm(P)[, 2])[1:(length(x) - 
                                              1)]
      if (P$IRT.param) {
        bxltm <- -ax * bx
      }
      else {
        bxltm <- bx
        bx <- -bxltm/ax
      }
      N <- dim(P$X)[1]
      ltmP <- P
      P <- ltmP$X
      covalphaP <- vcov.ltm(ltmP, robust = robust)
    }
    else {
      if (inherits(P, c("SingleGroupClass", "ConfirmatoryClass"))) {
        myspP <- extract.mirt(P, "parvec")
        bx <- myspP[seq(2, 2 * nX, by = 2)]
        ax <- myspP[seq(1, 2 * nX, by = 2)]
        mycovP <- extract.mirt(P, "vcov")
        upmat <- cbind(mycovP[seq(2, 2 * nX, by = 2), 
                              seq(2, 2 * nX, by = 2)], mycovP[seq(2, 2 * 
                                                                    nX, by = 2), seq(1, 2 * nX, by = 2)])
        lowmat <- cbind(mycovP[seq(1, 2 * nX, by = 2), 
                               seq(2, 2 * nX, by = 2)], mycovP[seq(1, 2 * 
                                                                     nX, by = 2), seq(1, 2 * nX, by = 2)])
        covP <- rbind(upmat, lowmat)
        bxltm <- bx
        bx <- -bxltm/ax
        dataP <- extract.mirt(P, "data")
        N <- nrow(dataP)
        ltmP <- P
        P <- dataP
        covalphaP <- covP
      }
      else {
        if (is.matrix(P)) {
          if ((length(x) - 1) != ncol(P)) 
            return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
          ltmP <- ltm(P ~ z1, IRT.param = FALSE)
          bx <- as.numeric(coef.ltm(ltmP)[, 1])[1:(length(x) - 
                                                     1)]
          ax <- as.numeric(coef.ltm(ltmP)[, 2])[1:(length(x) - 
                                                     1)]
          bxltm <- bx
          bx <- -bxltm/ax
          N <- dim(P)[1]
          covalphaP <- vcov.ltm(ltmP, robust = robust)
        }
        else {
          return("Unsupported input. P must be either an object created by the package ltm or a matrix of responses.")
        }
      }
    }
    if (inherits(Q, "ltm")) {
      by <- as.numeric(coef.ltm(Q)[, 1])[1:(length(y) - 
                                              1)]
      ay <- as.numeric(coef.ltm(Q)[, 2])[1:(length(y) - 
                                              1)]
      if (Q$IRT.param) {
        byltm <- -ay * by
      }
      else {
        byltm <- by
        by <- -byltm/ay
      }
      M <- dim(Q$X)[1]
      ltmQ <- Q
      Q <- ltmQ$X
      covalphaQ <- vcov.ltm(ltmQ, robust = robust)
    }
    else {
      if (inherits(Q, c("SingleGroupClass", "ConfirmatoryClass"))) {
        myspQ <- extract.mirt(Q, "parvec")
        by <- myspQ[seq(2, 2 * nY, by = 2)]
        ay <- myspQ[seq(1, 2 * nY, by = 2)]
        mycovQ <- extract.mirt(Q, "vcov")
        upmat <- cbind(mycovQ[seq(2, 2 * nY, by = 2), 
                              seq(2, 2 * nY, by = 2)], mycovQ[seq(2, 2 * 
                                                                    nY, by = 2), seq(1, 2 * nY, by = 2)])
        lowmat <- cbind(mycovQ[seq(1, 2 * nY, by = 2), 
                               seq(2, 2 * nY, by = 2)], mycovQ[seq(1, 2 * 
                                                                     nY, by = 2), seq(1, 2 * nY, by = 2)])
        covQ <- rbind(upmat, lowmat)
        byltm <- by
        by <- -byltm/ay
        dataQ <- extract.mirt(Q, "data")
        M <- nrow(dataQ)
        ltmQ <- Q
        Q <- dataQ
        covalphaQ <- covQ
      }
      else {
        return("Unsupported input. P must be either an object created by the package ltm or a matrix of responses.")
      }
    }
    irtx <- probpl(qpoints, bx, model, a = ax)
    irty <- probpl(qpoints, by, model, a = ay)
  }
  if (model == "3pl") {
    nX <- length(x) - 1
    nY <- length(y) - 1
    if (inherits(P, "tpm")) {
      cx <- as.numeric(coef.tpm(P)[, 1])[1:(length(x) - 
                                              1)]
      bx <- as.numeric(coef.tpm(P)[, 2])[1:(length(x) - 
                                              1)]
      ax <- as.numeric(coef.tpm(P)[, 3])[1:(length(x) - 
                                              1)]
      if (P$IRT.param) {
        bxltm <- -ax * bx
      }
      else {
        bxltm <- bx
        bx <- -bxltm/ax
      }
      N <- dim(P$X)[1]
      ltmP <- P
      P <- ltmP$X
      covalphaP <- vcov.ltm(ltmP, robust = robust)
    }
    else {
      if (inherits(P, c("SingleGroupClass", "ConfirmatoryClass"))) {
        myspP <- extract.mirt(P, "parvec")
        ax <- myspP[seq(1, 3 * nX, by = 3)]
        bx <- myspP[seq(2, 3 * nX, by = 3)]
        cx <- myspP[seq(3, 3 * nX, by = 3)]
        mycovP <- extract.mirt(P, "vcov")
        upmat <- cbind(mycovP[seq(3, 3 * nX, by = 3), 
                              seq(3, 3 * nX, by = 3)], mycovP[seq(3, 3 * 
                                                                    nX, by = 3), seq(2, 3 * nX, by = 3)], mycovP[seq(3, 
                                                                                                                     3 * nX, by = 3), seq(1, 3 * nX, by = 3)])
        midmat <- cbind(mycovP[seq(2, 3 * nX, by = 3), 
                               seq(3, 3 * nX, by = 3)], mycovP[seq(2, 3 * 
                                                                     nX, by = 3), seq(2, 3 * nX, by = 3)], mycovP[seq(2, 
                                                                                                                      3 * nX, by = 3), seq(1, 3 * nX, by = 3)])
        lowmat <- cbind(mycovP[seq(1, 3 * nX, by = 3), 
                               seq(3, 3 * nX, by = 3)], mycovP[seq(1, 3 * 
                                                                     nX, by = 3), seq(2, 3 * nX, by = 3)], mycovP[seq(1, 
                                                                                                                      3 * nX, by = 3), seq(1, 3 * nX, by = 3)])
        covP <- rbind(upmat, midmat, lowmat)
        bxltm <- bx
        bx <- -bxltm/ax
        cxltm <- cx
        cx <- exp(cx)/(1 + exp(cx))
        dataP <- extract.mirt(P, "data")
        N <- nrow(dataP)
        ltmP <- P
        P <- dataP
        covalphaP <- covP
      }
      else {
        if (is.matrix(P)) {
          if ((length(x) - 1) != ncol(P)) 
            return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
          ltmP <- ltm(P ~ z1, IRT.param = FALSE)
          cx <- as.numeric(coef.tpm(ltmP)[, 1])[1:(length(x) - 
                                                     1)]
          bx <- as.numeric(coef.tpm(ltmP)[, 2])[1:(length(x) - 
                                                     1)]
          ax <- as.numeric(coef.tpm(ltmP)[, 3])[1:(length(x) - 
                                                     1)]
          bxltm <- bx
          bx <- -bxltm/ax
          N <- dim(P)[1]
          covalphaP <- vcov.ltm(ltmP, robust = robust)
        }
        else {
          return("Unsupported input. P must be either an object created by the package ltm or a matrix of responses.")
        }
      }
    }
    if (inherits(Q, "tpm")) {
      cy <- as.numeric(coef.tpm(Q)[, 1])[1:(length(y) - 
                                              1)]
      by <- as.numeric(coef.tpm(Q)[, 2])[1:(length(y) - 
                                              1)]
      ay <- as.numeric(coef.tpm(Q)[, 3])[1:(length(y) - 
                                              1)]
      if (Q$IRT.param) {
        byltm <- -ay * by
      }
      else {
        byltm <- by
        by <- -byltm/ay
      }
      M <- dim(Q$X)[1]
      ltmP <- Q
      Q <- ltmP$X
      covalphaQ <- vcov.ltm(ltmQ, robust = robust)
    }
    else {
      if (inherits(Q, c("SingleGroupClass", "ConfirmatoryClass"))) {
        myspQ <- extract.mirt(Q, "parvec")
        ay <- myspQ[seq(1, 3 * nY, by = 3)]
        by <- myspQ[seq(2, 3 * nY, by = 3)]
        cy <- myspQ[seq(3, 3 * nY, by = 3)]
        mycovQ <- extract.mirt(Q, "vcov")
        upmat <- cbind(mycovQ[seq(3, 3 * nY, by = 3), 
                              seq(3, 3 * nY, by = 3)], mycovQ[seq(3, 3 * 
                                                                    nY, by = 3), seq(2, 3 * nY, by = 3)], mycovQ[seq(3, 
                                                                                                                     3 * nY, by = 3), seq(1, 3 * nY, by = 3)])
        midmat <- cbind(mycovQ[seq(2, 3 * nY, by = 3), 
                               seq(3, 3 * nY, by = 3)], mycovQ[seq(2, 3 * 
                                                                     nY, by = 3), seq(2, 3 * nY, by = 3)], mycovQ[seq(2, 
                                                                                                                      3 * nY, by = 3), seq(1, 3 * nY, by = 3)])
        lowmat <- cbind(mycovQ[seq(1, 3 * nY, by = 3), 
                               seq(3, 3 * nY, by = 3)], mycovQ[seq(1, 3 * 
                                                                     nY, by = 3), seq(2, 3 * nY, by = 3)], mycovQ[seq(1, 
                                                                                                                      3 * nY, by = 3), seq(1, 3 * nY, by = 3)])
        covQ <- rbind(upmat, midmat, lowmat)
        byltm <- by
        by <- -byltm/ay
        cyltm <- cy
        cy <- exp(cy)/(1 + exp(cy))
        dataQ <- extract.mirt(Q, "data")
        M <- nrow(dataQ)
        ltmQ <- Q
        Q <- dataQ
        covalphaQ <- covQ
      }
      else {
        if (is.matrix(Q)) {
          if ((length(y) - 1) != ncol(Q)) 
            return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
          ltmQ <- ltm(Q ~ z1, IRT.param = FALSE)
          cy <- as.numeric(coef.tpm(ltmQ)[, 1])[1:(length(y) - 1)]
          by <- as.numeric(coef.tpm(ltmQ)[, 2])[1:(length(y) - 1)]
          ay <- as.numeric(coef.tpm(ltmQ)[, 3])[1:(length(y) - 1)]
          byltm <- by
          by <- -byltm/ay
          M <- dim(Q)[1]
          covalphaQ <- vcov.ltm(ltmQ, robust = robust)
        }
        else {
          return("Unsupported input. P must be either an object created by the package ltm or a matrix of responses.")
        }
      }
    }
    irtx <- probpl(qpoints, bx, model, a = ax)
    irty <- probpl(qpoints, by, model, a = ay)
  }
  if (model == "GPCM") {
    JX <- length(catsX)
    JY <- length(catsY)
    kX <- sum(catsX) - JX
    kY <- sum(catsY) - JY
    gpcmP <- P
    gpcmQ <- Q
    ax <- numeric(JX)
    ay <- numeric(JY)
    bx <- vector("list", JX)
    by <- vector("list", JY)
    if (inherits(gpcmP, "gpcm") && inherits(gpcmQ, "gpcm")) {
      if (!P$IRT.param || !Q$IRT.param) 
        return("Please fit the IRT models using IRT.param = TRUE.")
      if (length(unique(catsX)) == 1) {
        ax <- coef(gpcmP)[, catsX[1]]
        for (i in 1:JX) bx[[i]] <- coef(gpcmP)[i, 
                                               1:(catsX[1] - 1)]
      }
      else {
        for (i in 1:JX) ax[i] <- coef(gpcmP)[[i]][catsX[i]]
        for (i in 1:JX) bx[[i]] <- coef(gpcmP)[[i]][1:(catsX[i] - 
                                                         1)]
      }
      if (length(unique(catsY)) == 1) {
        ay <- coef(gpcmQ)[, catsY[1]]
        for (i in 1:JY) by[[i]] <- coef(gpcmQ)[i, 
                                               1:(catsY[1] - 1)]
      }
      else {
        for (i in 1:JY) ay[i] <- coef(gpcmQ)[[i]][catsY[i]]
        for (i in 1:JY) by[[i]] <- coef(gpcmQ)[[i]][1:(catsY[i] - 
                                                         1)]
      }
      N <- nrow(gpcmP$X)
      M <- nrow(gpcmQ$X)
      P <- gpcmP$X
      Q <- gpcmQ$X
    }
    if (inherits(gpcmP, "SingleGroupClass") && inherits(gpcmQ, 
                                                        "SingleGroupClass")) {
      for (i in 1:JX) {
        ttpar <- extract.item(gpcmP, i)
        ax[i] <- ttpar@par[1]
        bx[[i]] <- -(tail(ttpar@par, catsX[i] - 1) - 
                       c(0, tail(ttpar@par, catsX[i] - 1)[-(catsX[i] - 
                                                              1)]))/ax[i]
      }
      for (i in 1:JY) {
        ttpar <- extract.item(gpcmQ, i)
        ay[i] <- ttpar@par[1]
        by[[i]] <- -(tail(ttpar@par, catsY[i] - 1) - 
                       c(0, tail(ttpar@par, catsY[i] - 1)[-(catsY[i] - 
                                                              1)]))/ay[i]
      }
      dataP <- extract.mirt(gpcmP, "data")
      dataQ <- extract.mirt(gpcmQ, "data")
      N <- nrow(dataP)
      M <- nrow(dataQ)
      P <- dataP
      Q <- dataQ
    }
    ltmP <- gpcmP
    ltmQ <- gpcmQ
    irtx <- polyprob_kequate(ax, bx, catsX, model, qpoints)
    irty <- polyprob_kequate(ay, by, catsY, model, qpoints)
  }
  if (model == "GRM") {
    JX <- length(catsX)
    JY <- length(catsY)
    kX <- sum(catsX) - JX
    kY <- sum(catsY) - JY
    grmP <- P
    grmQ <- Q
    ax <- numeric(JX)
    ay <- numeric(JY)
    bx <- vector("list", JX)
    by <- vector("list", JY)
    for (i in 1:JX) {
      ttpar <- extract.item(grmP, i)
      ax[i] <- ttpar@par[1]
      bx[[i]] <- -ttpar@par[-1]/ax[i]
    }
    for (i in 1:JY) {
      ttpar <- extract.item(grmQ, i)
      ay[i] <- ttpar@par[1]
      by[[i]] <- -ttpar@par[-1]/ay[i]
    }
    dataP <- extract.mirt(grmP, "data")
    dataQ <- extract.mirt(grmQ, "data")
    N <- nrow(dataP)
    M <- nrow(dataQ)
    P <- dataP
    Q <- dataQ
    ltmP <- grmP
    ltmQ <- grmQ
    irtx <- polyprob_kequate(ax, bx, catsX, model, qpoints)
    irty <- polyprob_kequate(ay, by, catsY, model, qpoints)
  }
  
  if (model %in% c("1pl", "2pl", "3pl")) {
    r <- LordWW(irtx, qpoints)
    s <- LordWW(irty, qpoints)
  }
  if (model %in% c("GPCM", "GRM")) {
    r <- rowSums(cmnom_kequate(catsX, irtx, qpoints))
    s <- rowSums(cmnom_kequate(catsY, irty, qpoints))
  }
  
  # EE using estimated probabilties
  eg_x <- freqtab(0, 0:(length(r) - 1))
  eg_y <- freqtab(0, 0:(length(s) - 1))
  eg_x[] <- r * 1e10  # (multipled to avoid errors because of small probs)
  eg_y[] <- s * 1e10
  eg_eq <- equate(eg_x, eg_y, type = "equipercentile", method = "none", smoothmethod = "none")
  
  return(eg_eq$con[[2]])
}


# Parameter linking part is taken from kequate package PSE method
IRT_OSE_NEAT <- function(P, Q, x, y, a = 0, qpoints = seq(-6, 6, by = 0.1),
                         model = "2pl", catsX = 0, catsY = 0, catsA = 0,
                         wS = 0.5, eqcoef = "mean-mean",
                         robust = FALSE, distribution = list("normal", par = data.frame(mu = 0, sigma = 1))) {
  # get parameters from model
  if (model == "2pl") {
    if (inherits(P, c("ConfirmatoryClass", "SingleGroupClass", "ltm"))) {
      Plist <- irtinput_kequate(P, x, a, robust, model,
                                catsX = catsX,
                                catsA = catsA
      )
    } else if (is.matrix(P)) {
      if ((length(a) + length(x) - 2) != ncol(P)) {
        return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
      }
      ltmP <- ltm(P ~ z1, IRT.param = FALSE)
      Plist <- irtinput_kequate(ltmP, x, a, robust, model)
    }
    if (inherits(Q, c("ConfirmatoryClass", "SingleGroupClass", "ltm"))) {
      Qlist <- irtinput_kequate(Q, y, a, robust, model,
                                catsX = catsY,
                                catsA = catsA
      )
    } else if (is.matrix(Q)) {
      if ((length(a) + length(x) - 2) != ncol(Q)) {
        return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
      }
      ltmQ <- ltm(Q ~ z1, IRT.param = FALSE)
      Qlist <- irtinput_kequate(ltmQ, x, a, robust, model)
    }
    ax <- Plist$ax
    aaP <- Plist$aaP
    bxltm <- Plist$bxltm
    baPltm <- Plist$baPltm
    bx <- Plist$bx
    baP <- Plist$baP
    ltmP <- Plist$ltmP
    P <- Plist$P
    N <- Plist$N
    covalphaP <- Plist$covP
    ay <- Qlist$ax
    aaQ <- Qlist$aaP
    byltm <- Qlist$bxltm
    baQltm <- Qlist$baPltm
    by <- Qlist$bx
    baQ <- Qlist$baP
    ltmQ <- Qlist$ltmP
    Q <- Qlist$P
    M <- Qlist$N
    covalphaQ <- Qlist$covP
  }
  if (model == "3pl") {
    if (inherits(P, c("ConfirmatoryClass", "SingleGroupClass", "tpm"))) {
      Plist <- irtinput_kequate(P, x, a, robust, model, catsX = catsX, catsA = catsA)
    } else if (is.matrix(P)) {
      if ((length(a) + length(x) - 2) != ncol(P)) {
        return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
      }
      ltmP <- tpm(P, IRT.param = FALSE)
      Plist <- irtinput_kequate(ltmP, x, a, robust, model)
    }
    if (inherits(Q, c("ConfirmatoryClass", "SingleGroupClass", "tpm"))) {
      Qlist <- irtinput_kequate(Q, y, a, robust, model, catsX = catsY, catsA = catsA)
    } else if (is.matrix(Q)) {
      if ((length(a) + length(x) - 2) != ncol(Q)) {
        return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
      }
      ltmQ <- tpm(Q, IRT.param = FALSE)
      Qlist <- irtinput_kequate(ltmQ, x, a, robust, model)
    }
    ax <- Plist$ax
    aaP <- Plist$aaP
    bxltm <- Plist$bxltm
    baPltm <- Plist$baPltm
    bx <- Plist$bx
    baP <- Plist$baP
    cx <- Plist$cx
    caP <- Plist$caP
    ltmP <- Plist$ltmP
    if (inherits(ltmP, c("ConfirmatoryClass", "SingleGroupClass"))) {
      cxmirt <- cx
      caPmirt <- caP
      cx <- exp(cx) / (1 + exp(cx))
      caP <- exp(caP) / (1 + exp(caP))
    }
    P <- Plist$P
    N <- Plist$N
    covalphaP <- Plist$covP
    ay <- Qlist$ax
    aaQ <- Qlist$aaP
    byltm <- Qlist$bxltm
    baQltm <- Qlist$baPltm
    by <- Qlist$bx
    baQ <- Qlist$baP
    cy <- Qlist$cx
    caQ <- Qlist$caP
    ltmQ <- Qlist$ltmP
    if (inherits(ltmP, c("ConfirmatoryClass", "SingleGroupClass"))) {
      cymirt <- cy
      caQmirt <- caQ
      cy <- exp(cy) / (1 + exp(cy))
      caQ <- exp(caQ) / (1 + exp(caQ))
    }
    Q <- Qlist$P
    M <- Qlist$N
    covalphaQ <- Qlist$covP
  }
  if (model == "GPCM") {
    JX <- length(catsX)
    JY <- length(catsY)
    JA <- length(catsA)
    kX <- sum(catsX) - JX
    kY <- sum(catsY) - JY
    kA <- sum(catsA) - JA
    gpcmP <- P
    gpcmQ <- Q
    parsP <- irtinput_kequate(gpcmP, x, a, robust, model, catsX = catsX, catsA = catsA)
    parsQ <- irtinput_kequate(gpcmQ, y, a, robust, model, catsX = catsY, catsA = catsA)
    ax <- parsP$ax
    bx <- parsP$bx
    aaP <- parsP$aaP
    baP <- parsP$baP
    ay <- parsQ$ax
    by <- parsQ$bx
    aaQ <- parsQ$aaP
    baQ <- parsQ$baP
    N <- parsP$N
    M <- parsQ$N
    P <- parsP$P
    Q <- parsQ$P
    ltmP <- gpcmP
    ltmQ <- gpcmQ
  }
  if (model == "GRM") {
    JX <- length(catsX)
    JY <- length(catsY)
    JA <- length(catsA)
    kX <- sum(catsX) - JX
    kY <- sum(catsY) - JY
    kA <- sum(catsA) - JA
    grmP <- P
    grmQ <- Q
    parsP <- irtinput_kequate(grmP, x, a, robust, model, catsX = catsX, catsA = catsA)
    parsQ <- irtinput_kequate(grmQ, y, a, robust, model, catsX = catsY, catsA = catsA)
    ax <- parsP$ax
    bx <- parsP$bx
    aaP <- parsP$aaP
    baP <- parsP$baP
    ay <- parsQ$ax
    by <- parsQ$bx
    aaQ <- parsQ$aaP
    baQ <- parsQ$baP
    N <- parsP$N
    M <- parsQ$N
    P <- parsP$P
    Q <- parsQ$P
    ltmP <- grmP
    ltmQ <- grmQ
  }
  if (model == "1pl") {
    irtpars <- list(bx = bx, by = by, baP = baP, baQ = baQ)
    irtrP <- probpl_kequate(qpoints, bx, model)
    irtrQ <- probpl_kequate(betaA * qpoints + betaB, bx, model)
    irtsQ <- probpl_kequate(qpoints, by, model)
    irtsP <- probpl_kequate(betaA * qpoints + betaB, by, model)
  }
  if (model == "2pl") {
    ltmPcf <- matrix(c(bxltm, baPltm, ax, aaP), nrow = length(c(bxltm, baPltm)))
    ltmQcf <- matrix(c(byltm, baQltm, ay, aaQ), nrow = length(c(byltm, baQltm)))
    dimnames(ltmPcf)[[1]] <- c(sprintf("X%03d", 1:(length(x) - 1)), sprintf("A%03d", 1:(length(a) - 1)))
    dimnames(ltmQcf)[[1]] <- c(sprintf("Y%03d", 1:(length(y) - 1)), sprintf("A%03d", 1:(length(a) - 1)))
    avoidprint <- capture.output({
      modPQ <- modIRT(
        coef = list(test1 = ltmPcf, test2 = ltmQcf),
        var = list(test1 = covalphaP, test2 = covalphaQ), names = paste("test", 1:2, sep = ""), ltparam = TRUE
      )
    })
    eqc <- direc(mods = modPQ, which = c(2, 1), method = eqcoef)
    betaA <- eqc$A
    betaB <- eqc$B
    irtpars <- list(
      ax = ax, ay = ay, aaP = aaP, aaQ = aaQ,
      bx = bx, by = by, baP = baP, baQ = baQ
    )
    irtrP <- probpl_kequate(qpoints, bx, model, a = ax)
    irtrQ <- probpl_kequate(betaA * qpoints + betaB, bx, model,
                            a = ax
    )
    irtsQ <- probpl_kequate(qpoints, by, model, a = ay)
    irtsP <- probpl_kequate((qpoints - betaB) / betaA, by, model, a = ay)
    adjcovalphaP <- adjltm(covalphaP, pars = list(ax = ax, aa = aaP, bxltm = bxltm, baltm = baPltm), "PSE", model = "2pl")
    adjcovalphaQ <- adjltm(covalphaQ, pars = list(ax = ay, aa = aaQ, bxltm = byltm, baltm = baQltm), "PSE", model = "2pl")
  }
  if (model == "3pl") {
    if (inherits(ltmP, "tpm")) {
      ltmPcf <- matrix(c(
        cx, caP, bxltm, baPltm, ax,
        aaP
      ), nrow = length(c(ax, aaP)))
      ltmQcf <- matrix(c(
        cy, caQ, byltm, baQltm, ay,
        aaQ
      ), nrow = length(c(ay, aaQ)))
      dimnames(ltmPcf)[[1]] <- c(sprintf("X%03d", 1:(length(x) - 1)), sprintf("A%03d", 1:(length(a) - 1)))
      dimnames(ltmQcf)[[1]] <- c(sprintf("Y%03d", 1:(length(y) - 1)), sprintf("A%03d", 1:(length(a) - 1)))
      avoidprint <- capture.output({
        modPQ <- modIRT(coef = list(test1 = ltmPcf, test2 = ltmQcf), var = list(test1 = covalphaP, test2 = covalphaQ), names = paste("test", 1:2, sep = ""), ltparam = TRUE, lparam = FALSE)
      })
    }
    if (inherits(ltmP, c("ConfirmatoryClass", "SingleGroupClass"))) {
      mirtPcf <- matrix(c(
        cxmirt, caPmirt, bxltm,
        baPltm, ax, aaP
      ), nrow = length(c(ax, aaP)))
      mirtQcf <- matrix(c(
        cymirt, caQmirt, byltm,
        baQltm, ay, aaQ
      ), nrow = length(c(ay, aaQ)))
      dimnames(mirtPcf)[[1]] <- c(sprintf(
        "X%03d",
        1:(length(x) - 1)
      ), sprintf("A%03d", 1:(length(a) -
                               1)))
      dimnames(mirtQcf)[[1]] <- c(sprintf(
        "Y%03d",
        1:(length(y) - 1)
      ), sprintf("A%03d", 1:(length(a) -
                               1)))
      avoidprint <- capture.output({
        modPQ <- modIRT(coef = list(
          test1 = mirtPcf,
          test2 = mirtQcf
        ), var = list(
          test1 = covalphaP,
          test2 = covalphaQ
        ), names = paste("test",
                         1:2,
                         sep = ""
        ), ltparam = TRUE, lparam = TRUE)
      })
    }
    eqc <- direc(mods = modPQ, which = c(2, 1), method = eqcoef)
    betaA <- eqc$A
    betaB <- eqc$B
    adjcovalphaP <- modPQ$test1$var
    adjcovalphaQ <- modPQ$test2$var
    irtpars <- list(
      ax = ax, ay = ay, aaP = aaP, aaQ = aaQ,
      bx = bx, by = by, baP = baP, baQ = baQ, cx = cx,
      cy = cy, caP = caP, caQ = caQ
    )
    irtrP <- probpl_kequate(qpoints, bx, a = ax, c = cx)
    irtrQ <- probpl_kequate(betaA * qpoints + betaB, bx,
                            a = ax,
                            c = cx
    )
    irtsQ <- probpl_kequate(qpoints, by, a = ay, c = cy)
    irtsP <- probpl_kequate((qpoints - betaB) / betaA, by,
                            a = ay,
                            c = cy
    )
  }
  if (model %in% c("1pl", "2pl", "3pl")) {
    rP <- LordWW_kequate(irtrP, qpoints)
    rQ <- LordWW_kequate(irtrQ, qpoints)
    sQ <- LordWW_kequate(irtsQ, qpoints)
    sP <- LordWW_kequate(irtsP, qpoints)
  }
  if (model %in% c("GPCM", "GRM")) {
    estAB <- eqcpoly_kequate(
      aaP, baP, aaQ, baQ, catsA, model,
      eqcoef, qpoints, distribution
    )
    betaA <- estAB[1]
    betaB <- estAB[2]
    irtrP <- polyprob_kequate(ax, bx, catsX, model, qpoints)
    irtrQ <- polyprob_kequate(ax, bx, catsX, model, betaA *
                                qpoints + betaB)
    irtsP <- polyprob_kequate(ay, by, catsY, model, (qpoints -
                                                       betaB) / betaA)
    irtsQ <- polyprob_kequate(ay, by, catsY, model, qpoints)
    covalphaP <- extract.mirt(ltmP, "vcov")
    covalphaQ <- extract.mirt(ltmQ, "vcov")
    if (model == "GPCM") {
      adjderP <- adjgpcmmirt_kequate(ltmP)
      adjderQ <- adjgpcmmirt_kequate(ltmQ)
    }
    if (model == "GRM") {
      adjderP <- adjgrmmirt(ltmP)
      adjderQ <- adjgrmmirt(ltmQ)
    }
    adjcovalphaP <- t(adjderP) %*% covalphaP %*% adjderP
    adjcovalphaQ <- t(adjderQ) %*% covalphaQ %*% adjderQ
  }
  if (model %in% c("GPCM", "GRM")) {
    rP <- rowSums(cmnom_kequate(catsX, irtrP, qpoints))
    rQ <- rowSums(cmnom_kequate(catsX, irtrQ, qpoints))
    sQ <- rowSums(cmnom_kequate(catsY, irtsQ, qpoints))
    sP <- rowSums(cmnom_kequate(catsY, irtsP, qpoints))
  }
  r <- wS * rP + (1 - wS) * rQ
  s <- wS * sP + (1 - wS) * sQ
  
  # EE using estimated probabilties
  eg_x <- freqtab(0, 0:(length(r) - 1))
  eg_y <- freqtab(0, 0:(length(s) - 1))
  eg_x[] <- r * 1e10  # (multipled to avoid errors because of small probs)
  eg_y[] <- s * 1e10
  eg_eq <- equate(eg_x, eg_y, type = "equipercentile", method = "none", smoothmethod = "none")
  
  return(eg_eq$con[[2]])
}

cmnom_kequate <- function(cats, probs, qpoints) 
{
  JX <- length(cats)
  xK <- sum(cats) - JX
  nsprobs <- matrix(0, nrow = xK + 1, ncol = length(qpoints))
  nsprobs[1:(cats[1]), ] <- probs[[1]]
  if (JX > 1) {
    for (i in 2:JX) {
      maxsc <- sum(cats[1:i]) - length(cats[1:i])
      sprobs <- nsprobs
      nsprobs <- matrix(0, nrow = xK + 1, ncol = length(qpoints))
      for (j in 1:(maxsc - cats[i] + 2)) nsprobs[j:(cats[i] + 
                                                      j - 1), ] <- t(sprobs[j, ] * t(probs[[i]])) + 
        nsprobs[j:(cats[i] + j - 1), ]
    }
  }
  t(t(nsprobs) * (dnorm(qpoints)/sum(dnorm(qpoints))))
}


probpl_kequate <- function(qpoints, b, model = "3pl", a = 0, c = 0) {
  out <- matrix(0, ncol = length(qpoints), nrow = length(b))
  if (model == "1pl") {
    tmat <- matrix(rep(b, length(qpoints)),
                   nrow = length(b),
                   ncol = length(qpoints)
    )
    tmat <- t(t(tmat) - qpoints)
    out <- 1 / (1 + exp(tmat))
  }
  if (model == "2pl") {
    tmat <- matrix(rep(a, length(qpoints)),
                   nrow = length(a),
                   ncol = length(qpoints)
    )
    tmat <- -t(t(tmat) * qpoints) + tmat * b
    out <- 1 / (1 + exp(tmat))
  }
  if (model == "3pl") {
    tmat <- matrix(rep(a, length(qpoints)),
                   nrow = length(a),
                   ncol = length(qpoints)
    )
    tmat <- -t(t(tmat) * qpoints) + tmat * b
    out <- c + (1 - c) / (1 + exp(tmat))
  }
  return(out)
}

polyprob_kequate <- function(a, b, cats, model, qpoints) {
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
      out[cats[i], ] <- 1 / (1 + exp(-a[i] * (qpoints - b[[i]][cats[i] - 1])))
      for (j in 2:(cats[i] - 1)) {
        out[j, ] <- 1 / (1 + exp(-a[i] * (qpoints - b[[i]][j - 1]))) - 1 / (1 + exp(-a[i] * (qpoints - b[[i]][j])))
      }
      probs[[i]] <- out
    }
  }
  return(probs)
}

LordWW_kequate <- function(p, qpoints) {
  if (!is.matrix(p)) {
    p <- as.matrix(p)
  }
  n <- dim(p)[1]
  N <- dim(p)[2]
  test <- array(0, c(n + 1, N, n))
  q <- 1 - p
  l <- 1
  test[1, , l] <- q[1, ]
  test[2, , l] <- p[1, ]
  testR <- test[1, , l]
  testR2 <- test[2, , l]
  for (i in 2:n) {
    for (r in 1:(i + 1)) {
      if (r == 1) {
        testR <- testR * q[i, ]
      }
      if (r > 1 && r < (i + 1)) {
        test[r, , l + 1] <- test[r, , l] * q[i, ] +
          test[r - 1, , l] * p[i, ]
      }
      if (r == (i + 1)) {
        testR2 <- testR2 * p[i, ]
        test[r, , l + 1] <- testR2
      }
    }
    test[1, , l + 1] <- testR
    l <- l + 1
  }
  test[, , n] <- t(t(test[, , n]) * (dnorm(qpoints) / sum(dnorm(qpoints))))
  if (N == 1) {
    return(test[, , n])
  } else {
    return(rowSums(test[, , n]))
  }
}

eqcpoly_kequate <- function(aP, bP, aQ, bQ, cats, model, eqcoef, qpoints, distribution) {
  J <- length(cats)
  sA <- sum(aQ) / sum(aP)
  sB <- do.call(sum, bP) / (sum(cats) - J) - (sA * do.call(
    sum,
    bQ
  )) / (sum(cats) - J)
  if (eqcoef == "mean-mean") {
    outAB <- c(sA, sB)
  }
  return(outAB)
}

adjgpcmmirt_kequate <- function(input) {
  ncats <- extract.mirt(input, "K")
  nitem <- length(ncats)
  pdermat <- matrix(0, ncol = sum(ncats), nrow = sum(ncats))
  k <- 1
  for (i in 1:nitem) {
    ttpar <- extract.item(input, i)
    tempmat <- matrix(0, nrow = ncats[i], ncol = ncats[i])
    tempmat[1, 1] <- ttpar@par[2 + ncats[i] + 1] / ttpar@par[1]^2
    for (j in 2:(ncats[i] - 1)) tempmat[1, j] <- (ttpar@par[j + ncats[i] + 2] - ttpar@par[j + ncats[i] + 1]) / ttpar@par[1]^2
    tempmat[1, ncats[i]] <- 1
    for (j in 2:(ncats[i] - 1)) {
      tempmat[j, j - 1] <- -1 / ttpar@par[1]
      tempmat[j, j] <- 1 / ttpar@par[1]
    }
    tempmat[ncats[i], ncats[i] - 1] <- -1 / ttpar@par[1]
    pdermat[k:(k + ncats[i] - 1), k:(k + ncats[i] - 1)] <- tempmat
    k <- k + ncats[i]
  }
  pdermat
}

irtinput_kequate <- function(P, x, a, robust, model, SE = T, catsX = 0, catsA = 0) {
  covP <- NULL
  if (model == "2pl") {
    nX <- length(x) - 1
    nA <- length(a) - 1
    if (inherits(P, "ltm")) {
      bx <- as.numeric(coef.ltm(P)[, 1])[1:nX]
      baP <- as.numeric(coef.ltm(P)[, 1])[(nX + 1):(nX +
                                                      nA)]
      ax <- as.numeric(coef.ltm(P)[, 2])[1:nX]
      aaP <- as.numeric(coef.ltm(P)[, 2])[(nX + 1):(nX +
                                                      nA)]
      if (P$IRT.param) {
        bxltm <- -ax * bx
        baPltm <- -aaP * baP
      } else {
        bxltm <- bx
        baPltm <- baP
        bx <- -bxltm / ax
        baP <- -baPltm / aaP
      }
      N <- dim(P$X)[1]
      ltmP <- P
      P <- ltmP$X
      if (SE) {
        covP <- vcov.ltm(ltmP, robust = robust)
      }
    } else if (inherits(P, c("SingleGroupClass", "ConfirmatoryClass"))) {
      myspP <- extract.mirt(P, "parvec")
      b <- myspP[seq(2, 2 * (nX + nA), by = 2)]
      a <- myspP[seq(1, 2 * (nX + nA), by = 2)]
      ax <- a[1:nX]
      bx <- b[1:nX]
      baP <- b[(nX + 1):(nX + nA)]
      aaP <- a[(nX + 1):(nX + nA)]
      if (SE) {
        mycovP <- extract.mirt(P, "vcov")
        upmat <- cbind(mycovP[seq(2, 2 * (nX + nA), by = 2), seq(2, 2 * (nX + nA), by = 2)], mycovP[seq(2, 2 * (nX + nA), by = 2), seq(1, 2 * (nX + nA), by = 2)])
        lowmat <- cbind(mycovP[seq(1, 2 * (nX + nA), by = 2), seq(2, 2 * (nX + nA), by = 2)], mycovP[seq(1, 2 * (nX + nA), by = 2), seq(1, 2 * (nX + nA), by = 2)])
        covP <- rbind(upmat, lowmat)
      }
      bxltm <- bx
      baPltm <- baP
      bx <- -bxltm / ax
      baP <- -baPltm / aaP
      dataP <- extract.mirt(P, "data")
      N <- nrow(dataP)
      ltmP <- P
      P <- dataP
    } else if (is.matrix(P)) {
      if ((length(a) + length(x) - 2) != ncol(P)) {
        return("Unsupported input. Input matrices must have rows denoting individuals and columns denoting items.")
      }
      ltmP <- ltm(P ~ z1, IRT.param = FALSE)
      bx <- as.numeric(coef.ltm(ltmP)[, 1])[1:nX]
      baP <- as.numeric(coef.ltm(ltmP)[, 1])[(nX + 1):(nX +
                                                         nA)]
      ax <- as.numeric(coef.ltm(ltmP)[, 2])[1:nX]
      aaP <- as.numeric(coef.ltm(ltmP)[, 2])[(nX + 1):(nX +
                                                         nA)]
      bxltm <- bx
      baPltm <- baP
      bx <- -bxltm / ax
      baP <- -baPltm / aaP
      N <- dim(P)[1]
      if (SE) {
        covP <- vcov.ltm(ltmP, robust = robust)
      }
    } else if (is.list(P)) {
      if (P$IRT.param == T) {
        ax <- P$pars[1:nX, 1]
        bx <- P$pars[1:nX, 2]
        aaP <- P$pars[(nX + 1):(nX + nA), 1]
        baP <- P$pars[(nX + 1):(nX + nA), 2]
        baPltm <- -ax * bx
        bxltm <- -aaP * baP
        N <- P$N
        if (SE) {
          covP <- P$cov
        }
        ltmP <- NULL
        P <- matrix(0)
      } else {
        bxltm <- P$pars[1:nX, 2]
        baPltm <- P$pars[(nX + 1):(nX + nA), 2]
        ax <- P$pars[1:nX, 1]
        aaP <- P$pars[(nX + 1):(nX + nA), 1]
        bx <- -bxltm / ax
        baP <- -baPltm / aaP
        N <- P$N
        if (SE) {
          covP <- P$cov
        }
        ltmP <- NULL
        P <- matrix(0)
      }
    } else {
      return("Unsupported input. P must be either an object created by the packages ltm or mirt, a matrix of responses or a list with parameters, covariance matrix and sample size.")
    }
    return(list(
      bxltm = bxltm, baPltm = baPltm, ax = ax,
      aaP = aaP, bx = bx, baP = baP, N = N, covP = covP,
      ltmP = ltmP, P = P, IRT.param = T
    ))
  }
  if (model == "3pl") {
    nX <- length(x) - 1
    nA <- length(a) - 1
    if (inherits(P, "tpm")) {
      cx <- as.numeric(coef.tpm(P)[, 1])[1:(length(x) -
                                              1)]
      caP <- as.numeric(coef.tpm(P)[, 1])[(length(x)):((length(x) +
                                                          length(a) - 2))]
      bx <- as.numeric(coef.tpm(P)[, 2])[1:(length(x) -
                                              1)]
      baP <- as.numeric(coef.tpm(P)[, 2])[(length(x)):((length(x) +
                                                          length(a) - 2))]
      ax <- as.numeric(coef.tpm(P)[, 3])[1:(length(x) -
                                              1)]
      aaP <- as.numeric(coef.tpm(P)[, 3])[(length(x)):((length(x) +
                                                          length(a) - 2))]
      if (P$IRT.param) {
        bxltm <- -ax * bx
        baPltm <- -aaP * baP
      } else {
        bxltm <- bx
        baPltm <- baP
        bx <- -bxltm / ax
        baP <- -baPltm / aaP
      }
      N <- dim(P$X)[1]
      ltmP <- P
      P <- ltmP$X
      if (SE) {
        covP <- vcov.tpm(ltmP)
      }
    } else if (inherits(P, c("SingleGroupClass", "ConfirmatoryClass"))) {
      myspP <- extract.mirt(P, "parvec")
      a <- myspP[seq(1, 3 * (nX + nA), by = 3)]
      b <- myspP[seq(2, 3 * (nX + nA), by = 3)]
      cpar <- myspP[seq(3, 3 * (nX + nA), by = 3)]
      ax <- a[1:nX]
      bx <- b[1:nX]
      cx <- cpar[1:nX]
      aaP <- a[(nX + 1):(nX + nA)]
      baP <- b[(nX + 1):(nX + nA)]
      caP <- cpar[(nX + 1):(nX + nA)]
      if (SE) {
        mycovP <- extract.mirt(P, "vcov")
        upmat <- cbind(mycovP[seq(3, 3 * (nX + nA),
                                  by = 3
        ), seq(3, 3 * (nX + nA), by = 3)], mycovP[seq(3,
                                                      3 * (nX + nA),
                                                      by = 3
        ), seq(2, 3 * (nX + nA),
               by = 3
        )], mycovP[
          seq(3, 3 * (nX + nA), by = 3),
          seq(1, 3 * (nX + nA), by = 3)
        ])
        midmat <- cbind(mycovP[seq(2, 3 * (nX + nA),
                                   by = 3
        ), seq(3, 3 * (nX + nA), by = 3)], mycovP[seq(2,
                                                      3 * (nX + nA),
                                                      by = 3
        ), seq(2, 3 * (nX + nA),
               by = 3
        )], mycovP[
          seq(2, 3 * (nX + nA), by = 3),
          seq(1, 3 * (nX + nA), by = 3)
        ])
        lowmat <- cbind(mycovP[seq(1, 3 * (nX + nA),
                                   by = 3
        ), seq(3, 3 * (nX + nA), by = 3)], mycovP[seq(1,
                                                      3 * (nX + nA),
                                                      by = 3
        ), seq(2, 3 * (nX + nA),
               by = 3
        )], mycovP[
          seq(1, 3 * (nX + nA), by = 3),
          seq(1, 3 * (nX + nA), by = 3)
        ])
        covP <- rbind(upmat, midmat, lowmat)
      }
      bxltm <- bx
      baPltm <- baP
      bx <- -bxltm / ax
      baP <- -baPltm / aaP
      dataP <- extract.mirt(P, "data")
      N <- nrow(dataP)
      ltmP <- P
      P <- dataP
    } else {
      return("Unsupported input. P must be either an object created by the packages ltm or mirt, a matrix of responses or a list with parameters, covariance matrix and sample size.")
    }
    return(list(
      bxltm = bxltm, baPltm = baPltm, ax = ax,
      aaP = aaP, bx = bx, baP = baP, cx = cx, caP = caP,
      N = N, covP = covP, ltmP = ltmP, P = P, IRT.param = T
    ))
  }
  if (model %in% c("GPCM")) {
    JX <- length(catsX)
    JA <- length(catsA)
    ax <- numeric(JX)
    bx <- vector("list", JX)
    aaP <- numeric(JA)
    baP <- vector("list", JA)
    if (inherits(P, "SingleGroupClass")) {
      for (i in 1:JX) {
        ttpar <- extract.item(P, i)
        ax[i] <- ttpar@par[1]
        bx[[i]] <- -(tail(ttpar@par, catsX[i] - 1) -
                       c(0, tail(ttpar@par, catsX[i] - 1)[-(catsX[i] -
                                                              1)])) / ax[i]
      }
      for (i in 1:JA) {
        ttpar <- extract.item(P, i + JX)
        aaP[i] <- ttpar@par[1]
        baP[[i]] <- -(tail(ttpar@par, catsA[i] - 1) -
                        c(0, tail(ttpar@par, catsA[i] - 1)[-(catsA[i] -
                                                               1)])) / aaP[i]
      }
      dataP <- extract.mirt(P, "data")
      N <- nrow(dataP)
      if (SE) {
        covP <- extract.mirt(P, "vcov")
      }
      ltmP <- P
      P <- dataP
    }
    return(list(
      ax = ax, aaP = aaP, bx = bx, baP = baP,
      N = N, covP = covP, ltmP = ltmP, P = P, IRT.param = F
    ))
  }
  if (model %in% c("GRM")) {
    JX <- length(catsX)
    JA <- length(catsA)
    ax <- numeric(JX)
    bx <- vector("list", JX)
    aaP <- numeric(JA)
    baP <- vector("list", JA)
    if (inherits(P, "SingleGroupClass")) {
      for (i in 1:JX) {
        ttpar <- extract.item(P, i)
        ax[i] <- ttpar@par[1]
        bx[[i]] <- -ttpar@par[-1] / ax[i]
      }
      for (i in 1:JA) {
        ttpar <- extract.item(P, i + JX)
        aaP[i] <- ttpar@par[1]
        baP[[i]] <- -ttpar@par[-1] / aaP[i]
      }
      dataP <- extract.mirt(P, "tabdata")
      N <- nrow(dataP)
      if (SE) {
        covP <- extract.mirt(P, "vcov")
      }
      ltmP <- P
      P <- dataP
    }
    return(list(
      ax = ax, aaP = aaP, bx = bx, baP = baP,
      N = N, covP = covP, ltmP = ltmP, P = P, IRT.param = F
    ))
  }
}
