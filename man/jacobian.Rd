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
  To do ...
}

\value{
  The Jacobian matrix.  The number of rows equals 2^(number of items) - 1,
  the number of columns equals the number of independent parameters in the
  model.
}

\seealso{
  \code{\link{blim}}, \code{\link{simulate.blim}}.
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
