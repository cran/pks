\name{jacobian}
\alias{jacobian}

\title{Jacobian Matrix for Basic Local Independence Model}
\description{Computes the Jacobian matrix for a basic local independence model
  (BLIM).
}
\usage{
jacobian(object, P.K = rep(1/nstates, nstates),
         beta = rep(0.1, nitems), eta = rep(0.1, nitems),
         betafix = rep(NA, nitems), etafix = rep(NA, nitems))
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
  \item{betafix, etafix}{vectors of fixed error and guessing parameter values;
    \code{NA} indicates a free parameter.}
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
  \emph{Behavior Research Methods}, \bold{44}(4), 1197--1211.
  \doi{10.3758/s13428-012-0187-z}
}

\seealso{
  \code{\link{blim}}, \code{\link{simulate.blim}},
  \code{\link{gradedness}}.
}

\examples{
data(endm)
m <- blim(endm$K2, endm$N.R)

## Test of identifiability
J <- jacobian(m)
dim(J)
qr(J)$rank
}
\keyword{models}
