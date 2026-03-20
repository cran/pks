# - S3 methods
#   - predict.blim()
#   - slm(...) |> predict() by inheritance
#
# - newdata argument
#   - How to deal with non-observed response patterns? Just predict them.
#   - Do we need the full (normalized) P.R? No, renormalize for observed R.
#
# - User interface
#   predict(m, newdata = c("1000", "0110"))  # predict state for provided R
#   predict(m)        # predict for each unique pattern m$N.R
#   predict(m, type = "state")  # predict single state
#   predict(m, type = "probs")  # predict state distribution for each R
#
# - How to deal with a posteriori equally likely states?
#   - min, max, random
#
# - How about predictions by minimum discrepancy?
#   - Refactor blim(): getiRK(), getdRK(), or rather getMD()
#
# - Maybe for inspiration
#   https://github.com/cran/nnet/blob/master/R/multinom.R
#
# data(endm)
# m <- blim(endm$K, endm$N.R, method = "ML", etafix = c(0, NA, NA, NA))
# predict(m, type = "probs")
# predict(m)
# predict(m, newdata = "1100")
# predict(m, newdata = "1100", type = "probs")
# iRK <- pks:::getMD(endm$K, as.binmat(endm$N.R))$i.RK
# blim(endm$K, endm$N.R, method = "ML") |> predict(i.RK = iRK)


## Recompute P(K|R) for each R
getPKR <- function(beta, eta, K, R, P.K, PRKfun, i.RK = NULL, method = "ML",
                   betafix = NULL, etafix = NULL,
                   incradius = 0, normalize = FALSE) {
  if(method %in% c("ML", "MDML")) {
    P.K.R <- matrix(0,
      nrow = nrow(R), ncol = nrow(K),            # else apply() will drop
      dimnames = list(rownames(R), rownames(K))  #   dim if R is a singleton
    )
    P.K.R[] <- t(t(PRKfun(beta, eta, K, R)) * P.K)   # P(R|K) * P(K)
  }
  if(method %in% c("MD", "MDML") && is.null(i.RK)) {
    i.RK <- getMD(K = K, R = R, betafix = betafix,
                  etafix = etafix, incradius = incradius)[["i.RK"]]
  }
  if(method == "MDML") {
    P.K.R <- P.K.R/rowSums(P.K.R)
    P.K.R <- i.RK * P.K.R
  }
  if(method == "MD") P.K.R <- i.RK
  if(normalize) P.K.R <- P.K.R/rowSums(P.K.R)
  P.K.R
}


# Predict single K for each pattern (R is implied by P.K.R)
# - Note: even with min and max, there could be ties
predict1K <- function(P.K.R, K, method) {
  get1Kidx <- list(
    random = function(x) {
      sample(x, 1)
    },
    min = function(x) {
      x[order(rowSums(K[x, ]))][1]  # order by |K|
    },
    max = function(x) {
      x[order(rowSums(K[x, ]), decreasing = TRUE)][1]
    }
  )
  Klist <- apply(P.K.R, 1, function(x) which(x == max(x)), simplify = FALSE)
  Kidx <- integer(nrow(P.K.R))
  n_pred_states <- lengths(Klist)
  Kidx[n_pred_states == 1L] <-    # if there is a single max-probable state
    Klist[n_pred_states == 1L] |>
    unlist()
  if(any(n_pred_states > 1L)) {   # if there are several, select one
    Kidx[n_pred_states > 1L] <-
      sapply(Klist[n_pred_states > 1L], get1Kidx[[method]])
  }
  K[Kidx, , drop = FALSE]
}


predict.blim <- function(object, newdata = NULL, type = c("state", "probs"),
                         method = c("ML", "MD", "MDML"), quiet = FALSE,
                         ties.method = c("min", "max", "random"), i.RK = NULL,
                         incradius = object$incradius, as.pattern = TRUE,
                         ...) {
  stopifnot(inherits(object, "blim"))
  type <- match.arg(type)
  method <- match.arg(method)
  if(method != object$method & !quiet)
    message(
      sprintf("estimation method is %s, prediction method is %s",
              object$method, method)
    )
  ties.method <- match.arg(ties.method)
  R <- if(is.null(newdata)) {
    as.binmat(fitted(object))  # internal resp. patterns, possibly zeropadded
  } else {
    as.binmat(newdata)         # explicitly provided response patterns
  }
  PRKfun <-
    if(length(which(c(object$betafix, object$etafix) == 0)))
      getPRK[["apply"]]
    else
      getPRK[["matmult"]]
  if(type == "probs") {
    P.K.R <- getPKR(beta = object$beta, eta = object$eta, K = object$K, R = R,
                    P.K = object$P.K, PRKfun = PRKfun, i.RK = i.RK,
                    method = method, incradius = incradius,
                    betafix = object$betafix, etafix = object$etafix,
                    normalize = TRUE)
    rownames(P.K.R) <- seq_len(nrow(R))
    P.K.R
  } else {
    P.K.R <- getPKR(beta = object$beta, eta = object$eta, K = object$K, R = R,
                    P.K = object$P.K, PRKfun = PRKfun, i.RK = i.RK,
                    method = method, incradius = incradius,
                    betafix = object$betafix, etafix = object$etafix)
    Khat <- predict1K(P.K.R, object$K, ties.method)
    if(as.pattern) as.pattern(Khat, ...) else Khat
  }
}

