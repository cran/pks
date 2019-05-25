\name{jacobian}
\alias{jacobian}

\title{Jacobian Matrix for Basic Local Independence Model}
\description{Computes the Jacobian matrix for a basic local independence model
  (BLIM).
}
\usage{
jacobian(object, P.K = rep(1/nstates, nstates),
         beta = rep(0.1, nitems), eta = rep(0.1, nitems),
         errtype = c("both", "error", "guessing"))
}

\arguments{
  \item{object}{an object of class \code{blim}, typically the result of a call
    to \code{\link{blim}}.}
  \item{P.K}{the vector of parameter values for probabilities of knowledge
    states.}
  \item{beta}{the vector of parameter values for probabilities of a careless
    error.}
  \item{eta}{the vector of parameter values for probabilities of a lucky
    guess.}
  \item{errtype}{type of response errors that can occur: \code{error} for
    careless errors only, \code{guessing} for lucky guesses only, and
    \code{both} for both error types.}
}

\details{
  This is a draft version.  It may change in future releases.
}

\value{
  The Jacobian matrix.  The number of rows equals 2^(number of items) - 1,
  the number of columns equals the number of independent parameters in the
  model.
}

\references{
  Heller, J. (2017).
  Identifiability in probabilistic knowledge structures.
  \emph{Journal of Mathematical Psychology}, \bold{77}, 46--57.
  \doi{10.1016/j.jmp.2016.07.008}

  Stefanutti, L., Heller, J., Anselmi, P., & Robusto, E. (2012).
  Assessing the local identifiability of probabilistic knowledge structures.
  \emph{Behavior Research Methods}, \bold{44}, 1197--1211.
  \doi{10.3758/s13428-012-0187-z}
}

\seealso{
  \code{\link{blim}}, \code{\link{simulate.blim}},
  \code{\link{gradedness}}.
}

\examples{
data(endm)
blim1 <- blim(endm$K2, endm$N.R)

## Test of identifiability
J <- jacobian(blim1)
dim(J)
qr(J)$rank
}
\keyword{models}
