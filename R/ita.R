## Item tree analysis
ita <- function(R, L = NULL, makeK = FALSE) {

  ## Helpers
  isTransitive <- function(x) {
    P <- bb <= x
    all(P[P %*% P > 0])
  }

  getKfromP <- function(x) rbind(0L, sets::binary_closure(t(x)))

  getdmin <- function(x) mean(apply(x, 1, min, na.rm = TRUE))

  getdisc <- function(R, K, total = TRUE) {
    d.RK <- apply(K, 1, function(k) colSums(xor(t(R), k)))
    disc <- sapply(list(fit = d.RK, complexity = t(d.RK)), getdmin)
    if(total) sum(disc) else disc
  }

  getoptimalL <- function(bb, R, transitiveL) {
    LL <- rev(transitiveL)
    if(length(LL) < 2)
      return(LL)
    P <- bb <= LL[1]          # start from largest threshold
    dsum1 <- getdisc(R, getKfromP(P))
    for(i in 2:length(LL)) {  # reduce threshold value
      P <- bb <= LL[i]
      dsum <- getdisc(R, getKfromP(P))
      if(dsum1 < dsum) break  # stop when dsum starts getting worse
      dsum1 <- dsum
    }
    LL[i - 1]
  }

  stopifnot(is.matrix(R))
  nitems <- ncol(R)
  cn <- colnames(R)
  cn <- if(!is.null(cn)) cn else letters[seq_len(nitems)]

  ## Count all pairwise counterexamples to p < q
  i <- expand.grid(x = seq_len(nitems), y = seq_len(nitems))
  bb <- matrix(colSums((1 - R[, i$x]) * R[, i$y]),
               ncol = nitems, dimnames = list("<" = cn, ">" = cn))
  tL <- Filter(isTransitive, sort(unique(as.vector(bb))))

  ## Select optimal threshold value
  searchL <- is.null(L)
  if(searchL) L <- getoptimalL(bb, R, transitiveL = tL)

  ## Define relation according to threshold
  if(!isTransitive(L)) stop("relation not transitive for threshold L")
  P <- bb <= L

  ## Make K by closing its base (requires transitivity)
  K <- disc <- NULL
  if(makeK) {
    K <- getKfromP(P)
    disc <- getdisc(R, K, total = FALSE)
    disc <- c(disc, total = sum(disc))
  }
  retval <- list(K = K, discrepancy = disc, transitiveL = tL,
                 searchL = searchL, L = L, P = bb)
  class(retval) <- "ita"
  retval
}


print.ita <- function(x, digits = max(3, getOption("digits") - 3),
                      ...){
  cat("\nItem tree analysis (ITA)\n")
  cat("\nViolations of precedence relation:\n")
  print(x$P)
  cat(if(x$searchL) "\nOptimal threshold" else "\nThreshold", "(L):", x$L)
  cat("\nTransitive thresholds:\n")
  print(x$transitiveL)
  if(!is.null(x$K)) {
    cat("\nNumber of knowledge states:", nrow(x$K))
    cat("\nDiscrepancy:\n")
    print(x$discrepancy, digits = digits)
  }
  cat("\n")
  invisible(x)
}

