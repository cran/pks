## Fitting the basic local independence model (BLIM) by MDML
blim <- function(K, N.R, method = c("MD", "ML", "MDML"), R = as.binmat(N.R),
  P.K = rep(1/nstates, nstates), beta = rep(0.1, nitems),
  eta = rep(0.1, nitems), errtype = c("both", "error", "guessing"),
  errequal = FALSE, incradius = 0, tol=0.0000001, maxiter = 10000) {

  K       <- as.matrix(K)
  N.R     <- setNames(as.integer(N.R), names(N.R))  # convert to named int
  N       <- sum(N.R)
  nitems  <- ncol(K)
  npat    <- nrow(R)
  nstates <- nrow(K)

  names(P.K) <- if(is.null(rownames(K))) as.pattern(K) else rownames(K)

  ## Assigning state K given response R
  em    <- switch(method <- match.arg(method), MD = 0, ML = 1, MDML = 1)
  md    <- switch(method, MD = 1, ML = 0, MDML = 1)
  d.RK  <- switch(errtype <- match.arg(errtype),
             both = t(apply(R, 1, function(r) apply(K, 1, function(q)
                      sum(xor(q, r))))),
            error = t(apply(R, 1, function(r) apply(K, 1, function(q)
                      if(any(q - r < 0)) NA else sum(q - r)))),
         guessing = t(apply(R, 1, function(r) apply(K, 1, function(q)
                      if(any(r - q < 0)) NA else sum(r - q)))))
  d.min <- apply(d.RK, 1, min, na.rm=TRUE)            # minimum discrepancy

  i.RK  <- (d.RK <= (d.min + incradius)) & !is.na(d.RK)

  ## Minimum discrepancy distribution 
  disc.tab <- xtabs(N.R ~ d.min)
  disc     <- as.numeric(names(disc.tab)) %*% disc.tab / N
  
  iter     <- 0
  maxdiff  <- 2 * tol

  while ((maxdiff > tol) && (iter < maxiter) &&
         ((md*(1 - em) != 1) || (iter == 0))) {
    pi.old   <- P.K
    beta.old <- beta
    eta.old  <- eta
    
    P.R.K  <- switch(errtype,
           both = t(apply(R, 1, function(r) apply(K, 1, function(q)
              prod(beta^((1-r)*q) * (1-beta)^(r*q) * eta^(r*(1-q)) * (1-eta)^((1-r)*(1-q)))))),
          error = t(apply(R, 1, function(r) apply(K, 1, function(q)
              prod(beta^((1-r)*q) * (1-beta)^(r*q) * 0^(r*(1-q)) * 1^((1-r)*(1-q)))))),
       guessing = t(apply(R, 1, function(r) apply(K, 1, function(q)
              prod(0^((1-r)*q) * 1^(r*q) * eta^(r*(1-q)) * (1-eta)^((1-r)*(1-q)))))))
    P.R    <- as.numeric(P.R.K %*% P.K)
    P.K.R  <- P.R.K * outer(1/P.R, P.K)         # prediction of P(K|R)
    mat.RK <- i.RK^md * P.K.R^em
    m.RK   <- (mat.RK / rowSums(mat.RK)) * N.R  # m.RK = E(M.RK) = P(K|R) * N(R)
    loglik <- sum(log(P.R) * N.R)

    ## Distribution of knowledge states
    P.K <- colSums(m.RK) / N

    ## Careless error and guessing parameters
    for(j in seq_len(nitems)) {
      beta[j] <- sum(m.RK[which(R[,j] == 0), which(K[,j] == 1)]) /
                 sum(m.RK[,which(K[,j] == 1)])

      eta[j]  <- sum(m.RK[which(R[,j] == 1), which(K[,j] == 0)]) /
                 sum(m.RK[,which(K[,j] == 0)])
    }

    maxdiff <- max(abs(c(P.K, beta, eta) - c(pi.old, beta.old, eta.old)))
    iter <- iter + 1
  }
  if(iter >= maxiter) warning("iteration maximum has been exceeded")

  ## Mean number of errors
  P.Kq <- numeric(nitems)
  for(j in seq_len(nitems))
    P.Kq[j] <- sum(P.K[which(K[,j] == 1)])
  nerror <- c("careless error" = sum(beta * P.Kq),
                 "lucky guess" = sum( eta * (1 - P.Kq)))

  ## Global beta and eta parameters
  if (errequal) {
    beta <- rep(nerror[1] / sum(P.Kq),            nitems)
     eta <- rep(nerror[2] / (nitems - sum(P.Kq)), nitems)
  }

  names(beta) <- names(eta) <-
    if(is.null(colnames(K))) {
      make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
        sep="")
    } else
      colnames(K)

  ## Number of parameters
  npar <- length(P.K) - 1 +
    (if(errtype == "both") 2 else 1) * (if(errequal) 1 else nitems)

  ## Goodness of fit
  fitted <- setNames(N*P.R, names(N.R))
  G2     <- 2*sum(N.R*log(N.R/fitted), na.rm=TRUE)
  df     <- (2^nitems - 1) - npar
  gof    <- c(G2=G2, df=df, pval = 1 - pchisq(G2, df))

  z <- list(discrepancy=c(disc), P.K=P.K, beta=beta, eta=eta,
    disc.tab=disc.tab, K=K, N.R=N.R, nitems=nitems, nstates=nstates,
    npatterns=npat, ntotal=N, nerror=nerror, npar=npar, errtype=errtype,
    method=method, iter=iter, loglik=loglik, fitted.values=fitted,
    goodness.of.fit=gof)
  class(z) <- "blim"
  z
}


print.blim <- function(x, P.Kshow = FALSE, errshow = TRUE,
  digits=max(3, getOption("digits") - 2), ...){
  cat("\nBasic local independence models (BLIMs)\n")
  cat("\nNumber of knowledge states:", x$nstates)
  cat("\nNumber of response patterns:", x$npatterns)
  cat("\nNumber of respondents:", x$ntotal)

  method <- switch(x$method,
            MD = "Minimum discrepancy",
            ML = "Maximum likelihood",
          MDML = "Minimum discrepancy maximum likelihood")
  cat("\n\nMethod:", method)
  cat("\nNumber of iterations:", x$iter)
  G2   <- x$goodness.of.fit[1]
  df   <- x$goodness.of.fit[2]
  pval <- x$goodness.of.fit[3]
  cat("\nGoodness of fit (2 log likelihood ratio):\n")
  cat("\tG2(", df, ") = ", format(G2, digits=digits), ", p = ",
      format(pval, digits=digits), "\n", sep="")

  cat("\nMinimum discrepancy distribution (mean = ",
    round(x$discrepancy, digits=digits), ")\n", sep="")
  disc.tab <- x$disc.tab
  names(dimnames(disc.tab)) <- NULL
  print(disc.tab)
  cat("\nMean number or errors (total = ",
    round(sum(x$nerror), digits=digits), ")\n", sep="")
  print(x$nerror)
  if(P.Kshow){
    cat("\nDistribution of knowledge states\n")
    printCoefmat(cbind("P(K)"=x$P.K), digits=digits, cs.ind=1, tst.ind=NULL)
  }
  if(errshow){
    cat("\nError and guessing parameters\n")
    printCoefmat(cbind(beta=x$beta, eta=x$eta), digits=digits, cs.ind=1:2,
      tst.ind=NULL)
  }
  cat("\n")
  invisible(x)
}


## Log-likelihood for blim objects
logLik.blim <- function(object, ...){
  if(length(list(...)))
    warning("extra arguments discarded")
  p <- object$npar
  val <- object$loglik
  attr(val, "df") <- p
  attr(val, "nobs") <- object$npatterns
  class(val) <- "logLik"
  val
}


## Number of obsevations
nobs.blim <- function(object, ...) object$npatterns


## Residuals for BLIMs
residuals.blim <- function(object, type=c("deviance", "pearson"), ...){

  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))

  type <- match.arg(type)
  wts <- object$ntotal
  y <- object$N.R / wts
  mu <- object$fitted.values/wts
  res <- switch(type,
    deviance = if(object$goodness['df'] > 0){
        d.res <- sqrt(pmax(dev.resids(y, mu, wts), 0))
        ifelse(y > mu, d.res, -d.res)  # sign
      }
      else rep.int(0, length(mu)),
    pearson = (y - mu) * sqrt(wts)/sqrt(mu)
  )
  if(!is.null(object$na.action)) res <- naresid(object$na.action, res)
  res
}


## Diagnostic plot for BLIMs
plot.blim <- function(x,
  xlab="Predicted response probabilities", ylab="Deviance residuals", ...){

  xres <- resid(x)
  mu   <- x$fitted.values/x$ntotal
  plot(mu, xres, xlab = xlab, ylab = ylab, type="n", ...)
  abline(h = 0, lty = 2)
  panel.smooth(mu, xres)
}


## Simulate responses from BLIM
simulate.blim <- function(object, nsim = 1, seed = NULL, ...){
     P.K <- object$P.K
    beta <- object$beta
     eta <- object$eta
       K <- as.matrix(object$K)
  nitems <- object$nitems
       N <- object$ntotal

  state.id <- sample(seq_along(P.K), N, replace=TRUE, prob=P.K)  # draw states

  P.1.K <- K*(1 - beta) + (1 - K)*eta                # P(resp = 1 | K)
  R     <- matrix(N * nitems, N, nitems)             # response matrix
  for(i in seq_len(N))
    R[i,] <- rbinom(nitems, 1, P.1.K[state.id[i],])  # draw a response

  as.pattern(R, freq = TRUE)
}


## Binary matrix to pattern
as.pattern <- function(R, freq = FALSE){
  if(freq){
    N.R <- table(apply(R, 1, paste, collapse=""))
    setNames(as.integer(N.R), names(N.R))          # convert to named int
  }else
    unname(apply(R, 1, paste, collapse=""))
}


## Pattern to named binary matrix
as.binmat <- function(N.R, uniq = TRUE, col.names = NULL){
  pat <- if(is.null(names(N.R))) N.R else names(N.R)
  R   <- if(uniq) strsplit(pat, "") else strsplit(rep(pat, N.R), "")
  R   <- do.call(rbind, R)
  storage.mode(R) <- "numeric"

  colnames(R) <- 
    if(is.null(col.names)){
      nitems <- ncol(R)
      make.unique(c("a", letters[(seq_len(nitems) %% 26) + 1])[-(nitems + 1)],
        sep="")
    }else
      col.names

  R
}

