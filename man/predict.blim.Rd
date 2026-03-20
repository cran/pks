\name{predict.blim}
\alias{predict.blim}
\alias{predict.slm}
\title{Predict Knowledge States from Basic Local Independence Models (BLIMs)}
\description{
  Predict knowledge state or state probabilities based on a fitted \code{blim}
  object given a response pattern.
}
\usage{
\method{predict}{blim}(object, newdata = NULL, type = c("state", "probs"),
        method = c("ML", "MD", "MDML"), quiet = FALSE,
        ties.method = c("min", "max", "random"), i.RK = NULL,
        incradius = object$incradius, as.pattern = TRUE, \dots)
}
\arguments{
  \item{object}{an object of class \code{blim}, typically the result of a
    call to \code{\link{blim}}.}
  \item{newdata}{a character vector of response patterns with which to
    predict.}
  \item{type}{for each pattern, predict a knowledge state or the state
    probabilities.}
  \item{method}{how to compute the posterior state probabilities.
    (See Details.)}
  \item{quiet}{silence message when estimation and prediction methods differ.}
  \item{ties.method}{how to deal with a posteriori equally probable states.}
  \item{i.RK}{optional indicator matrix of states at minimum distance from
    response patterns.}
  \item{incradius}{see \code{\link{blim}}.}
  \item{as.pattern}{return a character vector via \code{\link{as.pattern}} or
    an indicator matrix.}
  \item{\dots}{further arguments passed to \code{\link{as.pattern}}.}
}
\details{
  Predicted is the modal posterior state (\code{type = "state"}) or the
  posterior distribution (\code{type = "probs"}) of knowledge states given a
  response pattern. Depending on the \code{method} argument, the posterior
  distribution is defined as:
  \itemize{
    \item Maximum likelihood (\code{method = "ML"})
      \deqn{P(K | R)_{ML} =
        \frac{P(R | K) P(K)}{\sum_K P(R | K) P(K)} =
        \frac{P(R | K) P(K)}{P(R)}}{%
      P(K|R)_ML = P(R|K) P(K) / sum_K(P(R|K) P(K))
                = P(R|K) P(K) / P(R)}
    \item Minimum discrepancy (\code{method = "MD"})
      \deqn{P(K | R)_{MD} =
        \frac{i_{RK}}{\sum_K i_{RK}}}{%
      P(K|R)_MD = i_RK / sum_K i_RK}
    \item Minimum discrepancy ML (\code{method = "MDML"})
      \deqn{P(K | R)_{MDML} =
        \frac{i_{RK} \cdot P(K | R)_{ML}}
             {\sum_K i_{RK} \cdot P(K | R)_{ML}}}{%
      P(K|R)_MDML = i_RK * P(K|R)_ML / sum_K(i_RK * P(K|R)_ML)}
  }
  where \eqn{i_{RK}}{i_RK} is a pattern-by-state indicator matrix that is one
  for each state \eqn{K} that is at minimum distance from pattern \eqn{R}.
}
\value{
  If \code{type = "state"}, a character vector of knowledge states (if
  \code{as.pattern = TRUE}) or a state-by-problem indicator matrix (if
  \code{as.pattern = FALSE}).

  If \code{type = "probs"}, a matrix of posterior state probabilities.
}
\seealso{
  \code{\link{blim}}, \code{\link{slm}}.
}
\examples{
data(DoignonFalmagne7)
m <- blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R)
predict(m)
predict(m, newdata = c("00100", "10100", "01111"), as.pattern = FALSE)

m <- slm(DoignonFalmagne7$K, DoignonFalmagne7$N.R)
predict(m, newdata = "00100")

data(endm)
m <- blim(endm$K, endm$N.R)
predict(m, type = "probs", method = "MD")
cbind(
  observed = names(m$N.R),
  min = predict(m, method = "MD", ties.method = "min"),
  max = predict(m, method = "MD", ties.method = "max"),
  rnd = predict(m, method = "MD", ties.method = "random")
) |> print(quote = FALSE)
}
\keyword{models}
