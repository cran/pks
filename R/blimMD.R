## Fitting BLIM by minimum discrepancy
blimMD <- function(K, N.R, R = as.binmat(N.R),
                   betafix = rep(NA, nitems), etafix = rep(NA, nitems),
                   incrule = c("minimum", "hypblc1", "hypblc2"), m = 1) {

  K       <- as.matrix(K)
  N.R     <- setNames(as.integer(N.R), names(N.R))  # convert to named int
  N       <- sum(N.R)
  nitems  <- ncol(K)
  npat    <- nrow(R)
  nstates <- nrow(K)

  ## Assigning state K given response R
  d.RK <- if(length(which(c(betafix, etafix) == 0))) {
    apply(K, 1, function(k) {
      RwoK <- t(R) & !k
      idx <- which(RwoK, arr.ind=TRUE)
      RwoK[idx[idx[, "row"] %in% which(etafix == 0), ]] <- NA
      
      KwoR <- k & !t(R)
      idx <- which(KwoR, arr.ind=TRUE)
      KwoR[idx[idx[, "row"] %in% which(betafix == 0), ]] <- NA
      colSums(RwoK) + colSums(KwoR)
    })
  } else {
    apply(K, 1, function(k) colSums(xor(t(R), k)))
  }
  d.min <- apply(d.RK, 1, min, na.rm=TRUE)            # minimum discrepancy

  i.RK  <- switch(match.arg(incrule),                 # inclusion rule
             minimum = (d.RK == d.min) &              !is.na(d.RK),
             hypblc1 = replace(1/(1 + d.RK - d.min)^m, is.na(d.RK), 0),
             hypblc2 = replace(1/(1 + d.RK)^m,         is.na(d.RK), 0))
  m.RK  <- i.RK/rowSums(i.RK) * N.R                   # P(K|R) * N(R)

  ## Minimum discrepancy distribution 
  disc.tab <- xtabs(N.R ~ d.min)
  disc     <- as.numeric(names(disc.tab)) %*% disc.tab / N

  ## Distribution of knowledge states
  P.K <- colSums(m.RK)/N
  names(P.K) <- if(is.null(rownames(K))) as.pattern(K) else rownames(K)

  ## Careless error and guessing parameters
  beta <- eta <- numeric(nitems)
  names(beta) <- names(eta) <-
    if(is.null(colnames(K))){
      make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
        sep="")
    }else colnames(K)
  for(j in seq_len(nitems)){
    beta[j] <- sum(m.RK[which(R[,j] == 0), which(K[,j] == 1)]) /
               sum(m.RK[,which(K[,j] == 1)])

    eta[j]  <- sum(m.RK[which(R[,j] == 1), which(K[,j] == 0)]) /
               sum(m.RK[,which(K[,j] == 0)])
  }
  beta[is.na(beta)] <- 0
   eta[is.na( eta)] <- 0
  beta[!is.na(betafix)] <- betafix[!is.na(betafix)]  # reset fixed parameters
   eta[!is.na( etafix)] <-  etafix[!is.na( etafix)]

  z <- list(discrepancy=c(disc), P.K=P.K, beta=beta, eta=eta,
    disc.tab=disc.tab, nstates=nstates, npatterns=npat, ntotal=N,
    nerror=NA, method="MD", iter=NA, loglik=NA)
  class(z) <- "blim"
  z
}

