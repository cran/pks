## Adapted from Stefanutti et al. (2012), DOI 10.3758/s13428-012-0187-z
blimit <- function(K, beta=NULL, eta=NULL, pi=NULL,
                   file_name=NULL) {
  nitems  <- ncol(K)
  nstates <- nrow(K)
  npar    <- nitems + nitems + nstates - 1
  
  if(is.null(beta)) beta <- .49*runif(nitems)
  if(is.null(eta))   eta <- .49*runif(nitems)
  if(is.null(pi)) {
    pi <- runif(nstates)  # FIXME: dirichlet sampling
    pi <- pi/sum(pi)
  }

  ## Initialization beta0new and eta0new
  beta0_new <- seq_len(nitems)
   eta0_new <- seq_len(nitems)
  
  ta   <- rep(1, nitems)
  tb   <- rep(1, nitems)
  tab  <- matrix(rep(1, 2*nitems), nrow=2)
  tap  <- rep(1, nitems)
  tbp  <- rep(1, nitems)
  tabp <- matrix(rep(1, 2*nitems), nrow=2)

  jacbas <- pksjacbas(K, beta, eta, pi)
  jac <- jacbas$jac
  
  ja <- jac[, seq_len(nitems)]
  jb <- jac[, (nitems + 1):(2*nitems)]
  jp <- jac[, 2*nitems + seq_len(nstates - 1)]
  
  r <- qr(jac)[["rank"]]
  
  # Print in a file called "file_name"
  if(!is.null(file_name)) sink(file_name)
  
  cat("\n\nB L I M I T\n")
  cat("BASIC LOCAL INDEPENDENCE MODEL IDENTIFICATION ANALYSIS\n\n")
  cat("Number of items:                 ", nitems, "\n")
  cat("Number of knowledge states:      ", nstates, "\n\n")
  cat("Total number of parameters:      ", npar, "\n")
  cat("Jacobian matrix rank:            ", r, "\n")
  cat("Null space dimension (NSD):      ", npar - r, "\n")
  
  info <- list(
     NItems = nitems,
    NStates = nstates,
       NPar = npar,
       Rank = r,
        NSD = npar - r
  )
  
  if(2^nitems - 1 <= r) {
    cat("\nWARNING: The model is not quantitatively testable.\n")
    cat("Jacobian matrix rank (", r, ") is not less than the number\n")
    cat("of independent observables (", 2^nitems-1, ").\n\n")
  }
  
  if(r < npar) {
    cat("Identification problems detected:\n")
    cat("Jacobian matrix is not full rank.\n\n")
    
    # ranks of the three submatrices
    ra <- qr(ja)[["rank"]]
    rb <- qr(jb)[["rank"]]
    rp <- qr(jp)[["rank"]]
    
    # Here all pairs of submatrices are assembled and their respective
    # ranks are computed
    rab <- qr(cbind(ja, jb))[["rank"]]
    rap <- qr(cbind(ja, jp))[["rank"]]
    rbp <- qr(cbind(jb, jp))[["rank"]]
    
    # Null space dimensions
    nsda <- length(beta) - ra
    nsdb <- length(eta) - rb
    nsdp <- nstates - rp - 1
    nsdab <- length(beta) + length(eta) - rab
    nsdap <- length(beta) + nstates - rap - 1
    nsdbp <- length(eta) + nstates - rbp - 1
    
    # Tradeoff dimensions
    tdab <- ra + rb - rab
    tdap <- ra + rp - rap
    tdbp <- rb + rp - rbp
    tdabp <- rab + rap + rbp - ra - rb - rp - r
    
    SRATable <- data.frame(
      SUBMATRIX = c("[BETA]", "[ETA]", "[PI]", "[BETA ETA]", "[BETA PI]",
                    "[ETA PI]", "[BETA ETA PI]"),
      NPAR = c(length(beta), length(eta), nstates - 1, 2*nitems,
               nitems + nstates - 1, nitems + nstates - 1, npar),
      RANK = c(ra, rb, rp, rab, rap, rbp, r),
      NSD = c(nsda, nsdb, nsdp, nsdab, nsdap, nsdbp, npar - r),
      TRADEOFF_DIM = c(nsda, nsdb, nsdp, tdab, tdap, tdbp, tdabp)
    )
    
    cat("Submatrix rank analysis table\n");
    cat("[BETA] = submatrix of the careless error parameters\n");
    cat("[ETA]  = submatrix of the lucky guess parameters\n");
    cat("[PI]   = submatrix of the state probabilities\n");
    print(SRATable)
    
    ITEMS_BETA <- paste0("BETA", seq_len(nitems))
    ITEMS_ETA  <- paste0("ETA",  seq_len(nitems))
   
    BETA <- beta
     ETA <- eta
    
    if(nsda > 0) {
      # what the following FOR loop does is: remove from [BETA] one
      # column at the time in a stepwise fashion. Each time compute the
      # rank of the matrix and see if it remains unchanged. If so, a
      # dependent column has been found and, thus, a dependent parameter
      # has been detected
      for(i in seq_len(nitems)) {
        index <- rep(1, nitems)
        index[i] <- 0
        if(qr(ja[, index])[["rank"]] == ra) {
          ta[i] <- 0
        }
      }
      ndep <- nitems - sum(ta)
      
      cat("\n\nItem diagnostics for [BETA] submatrix\n")
      cat("First-order tradeoff dimensions:   ", nsda, "\n")
      cat("0 = parameter is not independent\n")
      cat("1 = parameter is independent\n")
      cat("-------------------------------\n")
      Beta_diagnostics <- data.frame(
        ITEMS = seq_len(nitems),
        BETA = beta,
        INDEPENDENT = ta
      )
      print(Beta_diagnostics)
      cat("-------------------------------\n")
    }
    
    if(nsdb > 0) {
      # what the following FOR loop does is: remove from [ETA] one
      # column at the time in a stepwise fashion. Each time compute the
      # rank of the matrix and see if it remains unchanged. If so, a
      # dependent column has been found and, thus, a dependent parameter
      # has been detected
      for(i in seq_len(nitems)) {
        index <- rep(1, nitems)
        index[i] <- 0
        if(qr(jb[, index])[["rank"]] == rb) {
          tb[i] <- 0
        }
      }
      ndep <- nitems - sum(tb)
      
      cat("\n\nItem diagnostics for [ETA] submatrix\n")
      cat("First-order tradeoff dimensions:   ", nsdb, "\n")
      cat("0 = parameter is not independent\n")
      cat("1 = parameter is independent\n")
      cat("-------------------------------\n")
      Eta_diagnostics <- data.frame(
        ITEMS = seq_len(nitems),
        ETA = eta,
        INDEPENDENT = tb
      )
      print(Eta_diagnostics)
      cat("-------------------------------\n");
    }    

    if(tdab > 0) {
      # what the following FOR loop does is: remove from [BETA ETA] one
      # column at the time in a stepwise fashion. Each time compute the
      # rank of the matrix and see if it remains unchanged. If so, a
      # dependent column has been found and, thus, a dependent parameter
      # has been detected
      for(i in seq_len(nitems)) {
        index <- rep(1, nitems)
        index[i] <- 0
        if(qr(cbind(ja[, index], jb))[["rank"]] == rab) {
          tab[1, i] <- 0
        }
      }
      ndep <- nitems - sum(tab[1, ])
      for(i in seq_len(nitems)) {
        index <- rep(1, nitems)
        index[i] <- 0
        if(qr(cbind(ja, jb[, index]))[["rank"]] == rab) {
          tab[2, i] <- 0
        }
      }
      ndep <- ndep + nitems - sum(tab[2, ])  # FIXME: what is ndep for?
      
      cat("\n\nItem diagnostics for [BETA ETA] submatrix\n")
      cat("First-order tradeoff dimensions:   ", nsdab, "\n")
      cat("0 = parameter is not independent\n")
      cat("1 = parameter is independent\n")
      cat("-------------------------------\n")
      BetaEta_diagnostics <- data.frame(
        ITEMS = seq_len(nitems),
        BETA = beta,
        BETA_INDEPENDENT = tab[1, ],
        ETA = eta,
        ETA_INDEPENDENT = tab[2, ]
      )
      print(BetaEta_diagnostics)
      cat("-------------------------------\n")
    }

    if(tdap > 0) {
      nsb <- nullR(cbind(jp, ja))
      cat("\n\nItem diagnostics for [BETA PI] submatrix\n")
      cat("Second-order tradeoff dimensions:   ", tdap, "\n")
      cat("-------------------------------\n")
      
      PiBeta_diagnostics <- data.frame(ITEMS_BETA, round(BETA, 3),
                                       round(nsb[nstates:nrow(nsb), ], 2))
      if(ncol(PiBeta_diagnostics) - 2 != tdap) {
        stop("\nIt isn't possible obtain a rational base of the null space",
             "\n  for this set of parameters.",
             "\nRerun the function with a different set of parameters.")
      }
      
      colnames(PiBeta_diagnostics) <-
        c("PARAMS", "BETA", paste0("DIM", 3:ncol(PiBeta_diagnostics) - 2))
      print(PiBeta_diagnostics)
      cat("-------------------------------\n")
      tap <- nsb
      beta0_new <-
        if(is.matrix(nsb[nstates:nrow(nsb), ])) {
          which(rowSums(nsb[nstates:nrow(nsb), ] == 1) == 0)
        } else {
          which((nsb[nstates:nrow(nsb), ] == 1) == 0)
        }
    }
     
    if(tdbp > 0) {
      nsb <- nullR(cbind(jp, jb))
      cat("\n\nItem diagnostics for [ETA PI] submatrix\n")
      cat("Second-order tradeoff dimensions:   ", tdbp, "\n")
      cat("-------------------------------\n")
      
      PiEta_diagnostics <- data.frame(ITEMS_ETA, round(ETA, 3),
                                      round(nsb[nstates:nrow(nsb),], 2))
      if(ncol(PiEta_diagnostics) - 2 != tdbp) {
        stop("\nIt isn't possible obtain a rational base of the null space",
             "\n  for this set of parameters.",
             "\nRerun the function with a different set of parameters.")
      }
      
      colnames(PiEta_diagnostics) <-
        c("PARAMS", "ETA", paste0("DIM", 3:ncol(PiEta_diagnostics) - 2))
      print(PiEta_diagnostics)
      cat("-------------------------------\n")
      tbp <- nsb
      eta0_new <-
        if(is.matrix(nsb[nstates:nrow(nsb), ])) {
          which(rowSums(nsb[nstates:nrow(nsb), ] == 1) == 0)
        } else {
          which((nsb[nstates:nrow(nsb), ] == 1) == 0);
        }
    }

    # So far only tradeoffs between pairs of parameters have been checked.
    # Here we check whether there are tradeoffs involving all three types
    # of parameters. This happens when all three submatrices [BETA ETA],
    # [BETA PI] and [ETA PI] are full rank while [BETA ETA PI] is not.
    if(tdabp > 0) {
      for(i in seq_len(nitems)) {
        index <- rep(1, nitems)
        index[i] <- 0
        if(qr(cbind(ja[, index], jb, jp))[["rank"]] == r) {
          tabp[1, i] <- 0
        }
        if(qr(cbind(ja, jb[, index], jp))[["rank"]] == r){
          tabp[2, i] <- 0
        } 
      }         
      ja0_new <- ja[, beta0_new]
      jb0_new <- jb[, eta0_new]
      nsb <- nullR(cbind(jp, ja0_new, jb0_new))
      
      cat("\n\nItem diagnostics for [BETA ETA PI] submatrix\n")
      cat("Third-order tradeoff dimensions:   ", tdbp, "\n")
      cat("-------------------------------\n")
      
      PiBetaEta_diagnostics <-
        data.frame(c(ITEMS_BETA[beta0_new], ITEMS_ETA[eta0_new]),
                   round(c(BETA[beta0_new], ETA[eta0_new]), 3),
                   round(nsb[nstates:nrow(nsb), ], 2)
      )
      
      colnames(PiBetaEta_diagnostics) <-
        c("PARAMS", "VALUES", paste0("DIM", 3:ncol(PiBetaEta_diagnostics) - 2))
      print(PiBetaEta_diagnostics)
      cat("-------------------------------\n")
      tabp <- nsb
    }
    info$RankBeta <- ra
    info$RankEta <- rb
    info$RankPi <- rp
    info$RankBetaEta <- rab
    info$RankBetaPi <- rap
    info$RankEtaPi <- rbp
         
  } else {
    cat("No identification problems detected.\n")
  }
  info$DiagBetaEta <- tab
  info$DiagBetaPi <- tap
  info$DiagEtaPi <- tbp
  info$DiagBetaEtaPi <- tabp
  info$Jacobian <- jac
  info$beta <- beta
  info$eta <- eta
  info$pi <- pi
  
  if(!is.null(file_name)) sink()
  return(info)
}


pksjac <- function(pat, K, beta, eta, pi) {
  prk <- rho(beta, eta, K, pat)
  xa <- prk %*% diag(pi, length(pi), length(pi)) %*% K
  xb <- prk %*% diag(pi, length(pi), length(pi)) %*% (1 - K)
  da <- (xa*(1 - pat)) %*% diag(1/beta) -       (xa*pat) %*% diag(1/(1 - beta))
  db <-       (xb*pat) %*% diag(1/eta)  - (xb*(1 - pat)) %*% diag(1/(1 - eta)) 
  dp <- apply(prk, 2, function(x, y) x - y, y = prk[, ncol(prk)])
  dp <- t(as.matrix(dp[1:(ncol(prk) - 1)]))
  cbind(da, db, dp)
}


pksjacbas <- function(K, beta, eta, pi) {
  nitems <- ncol(K)
  nstates <- nrow(K)
  npar <- 2*nitems + nstates - 1
    
  r <- 0
  basis <- NULL
  jac <- NULL
  
  for(d in 0:(2^nitems - 1)) {
    x <- d
    k <- nitems
    p <- matrix(0, 1, nitems)
    while(k > 0 && x > 0) {
      p[k] <- x %% 2
      x <- x %/% 2
      k <- k - 1
    }
    
    basnew <- rbind(basis, p)
    jacrow <- pksjac(p, K, beta, eta, pi)
    jacnew <- rbind(jac, jacrow)
    rnew <- qr(jacnew)[["rank"]]
    if(rnew > r) {
      basis <- basnew
      jac <- jacnew
      r <- rnew
#     cat(paste0("\n", r, "/", npar, " ", p, collapse = ""))  # FIXME:
      cat(paste0(r, "/", npar, " ", paste(p, collapse = "")), "\n")
    }
    if(nrow(basis) == npar) break
  }
  list(jac=jac, basis=basis)
}


nullR <- function(A){
  # The results could be inaccurate for some value of the matrix A, because
  # the rational solution is less accurate to generate the base of the null
  # space.
  m <- nrow(A)
  n <- ncol(A)
  R <- rref(A)
  pivcol <- which(colSums(abs(R)) == 1)
  r <- length(pivcol)
  nopiv <- seq_len(n)
  nopiv <- nopiv[-pivcol]
  Z <- matrix(0, ncol = n - r, nrow = n)
  if(n > r){
    Z[nopiv, ] <- diag(n - r)
    if(r > 0) {
      Z[pivcol, ] = -R[seq_len(r), nopiv]
    }
  }
  Z
}


rho <- function(a, b, K, pat, miss = array(0, dim = dim(pat))) {
  m <- dim(pat)[1]
  eps <- 1e-9
  a <- a + (a < eps)*eps - (1 - a < eps)*eps  # prevents log of zero
  b <- b + (b < eps)*eps - (1 - b < eps)*eps
  u <- rep(1, m)
  xa <- u %*% t(log(a))
  ya <- u %*% t(log(1 - a))
  xb <- u %*% t(log(b))
  yb <- u %*% t(log(1 - b))
  p <-       (pat*(1 - miss)*ya) %*% t(K)
  q <-       (pat*(1 - miss)*xb) %*% t(1 - K)
  r <- ((1 - pat)*(1 - miss)*xa) %*% t(K)
  s <- ((1 - pat)*(1 - miss)*yb) %*% t(1 - K)
  c <- exp(p + q + r + s)
  return(c)
}


# Adapted from pracma package
# See also:
# https://stackoverflow.com/questions/3126759/reduced-row-echelon-form
rref <- function(A) {
  stopifnot(is.numeric(A))
  if(!is.matrix(A)) 
    stop("Input parameter 'A' must be a matrix.")
  nr <- nrow(A)
  nc <- ncol(A)
  tol <- .Machine$double.eps * max(nr, nc) * max(abs(A))
  r <- 1
  for(i in 1:nc) {
    pivot <- which.max(abs(A[r:nr, i]))
    pivot <- r + pivot - 1
    m <- abs(A[pivot, i])
    if(m <= tol) {
      A[r:nr, i] <- 0
    } else {
      A[c(pivot, r), i:nc] <- A[c(r, pivot), i:nc]
      A[r, i:nc] <- A[r, i:nc]/A[r, i]
      if(r == 1) {
        ridx <- (r + 1):nr
      } else if (r == nr) {
        ridx <- 1:(r - 1)
      } else {
        ridx <- c(1:(r - 1), (r + 1):nr)
      }
      A[ridx, i:nc] <-
        A[ridx, i:nc] - A[ridx, i, drop = FALSE] %*% A[r, i:nc, drop = FALSE]
      if (r == nr) 
        break
      r <- r + 1
    }
  }
  A[abs(A) < tol] <- 0
  A
}

