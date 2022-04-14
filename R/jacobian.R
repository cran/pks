jacobian <- function(object, P.K = rep(1/nstates, nstates),
                     beta = rep(0.1, nitems), eta = rep(0.1, nitems),
                     betafix = rep(NA, nitems), etafix = rep(NA, nitems))
{
  K <- as.matrix(object$K)
# N <- object$ntotal
  nitems <- ncol(K)
  nstates <- nrow(K)
  R <- as.matrix(expand.grid(rep(list(0:1), nitems)))
  rownames(R) <- as.pattern(R)
  npatterns <- nrow(R)

  PRKfun <- if(length(which(c(betafix, etafix) == 0))) {
    getPRK[["apply"]] 
  } else {
    getPRK[["matmult"]] 
  }
  betanew <- beta
   etanew <-  eta
  betanew[!is.na(betafix)] <- betafix[!is.na(betafix)]
   etanew[!is.na( etafix)] <-  etafix[!is.na( etafix)]
  P.R.K <- do.call(PRKfun, list(betanew, etanew, K, R))

# P.R <- as.numeric(P.R.K %*% P.K)
  K.star <- K[-nstates, ]
  R.star <- R[-npatterns, ]
  P.K.star <- P.K[-nstates]
  P.R.K.star <- P.R.K[-npatterns, -nstates]
  P <- P.R.K.star - outer(P.R.K[-npatterns, nstates], rep(1, nstates - 1))
  colnames(P) <- paste("pi", as.pattern(K.star), sep="_")
  B <- sapply(1:nitems, function(q) (as.matrix(P.R.K.star[, K.star[, q] == 1]) %*% P.K.star[K.star[, q] == 1] +
      P.R.K[-npatterns, nstates] * P.K[nstates]) * (-1/(1 - beta[q]))^R.star[, q] * (1/beta[q])^(1-R.star[, q]))
  colnames(B) <- paste("beta", colnames(K), sep="_")
  E <- sapply(1:nitems, function(q) (as.matrix(P.R.K.star[, K.star[, q] == 0]) %*% P.K.star[K.star[, q] == 0]) *
      (1/eta[q])^R.star[, q] * (-1/(1-eta[q]))^(1-R.star[, q]))
  colnames(E) <- paste("eta", colnames(K), sep="_")
  J <- cbind(P, B, E)
  rownames(J) <- rownames(R.star)
  J
}

