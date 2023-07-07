\name{blimit}
\alias{blimit}

\title{Basic Local Independence Model Identification Analysis}
\description{Tests the local identifiability of a basic local independence
  model (BLIM).
}
\usage{
blimit(K, beta = NULL, eta = NULL, pi = NULL, file_name = NULL)
}

\arguments{
  \item{K}{a state-by-problem indicator matrix representing the knowledge
    structure.  An element is one if the problem is contained in the state,
    and else zero.}
  \item{beta, eta, pi}{vectors of parameter values for probabilities of
    careless errors, lucky guesses, and knowledge states, respectively.}
  \item{file_name}{name of an output file.}
}

\details{
  See Stefanutti et al. (2012) for details.

  The \code{blimit} function has been adapted from code provided by Andrea
  Brancaccio, Debora de Chiusole, and Luca Stefanutti. It contains a function
  to compute the reduced row echelon form based on an implementation in the
  pracma package.
}

\value{
  A list having the following components:
  \item{NItems}{the number of items.}
  \item{NStates}{the number of knowledge states.}
  \item{NPar}{the number of parameters.}
  \item{Rank}{the rank of the Jacobian matrix.}
  \item{NSD}{the null space dimension.}
  \item{RankBeta, RankEta, RankPi, RankBetaEta, RankBetaPi, RankEtaPi}{the
    rank of submatrices of the Jacobian.}
  \item{DiagBetaEta, DiagBetaPi, DiagEtaPi, DiagBetaEtaPi}{diagnostic
    information about specific parameter trade-offs.}
  \item{Jacobian}{the Jacobian matrix.}
  \item{beta, eta, pi}{the parameter values used in the analysis.}
}

\references{
  Stefanutti, L., Heller, J., Anselmi, P., & Robusto, E. (2012).
  Assessing the local identifiability of probabilistic knowledge structures.
  \emph{Behavior Research Methods}, \bold{44}(4), 1197--1211.
  \doi{10.3758/s13428-012-0187-z}
}

\seealso{
  \code{\link{blim}}, \code{\link{jacobian}}.
}

\examples{
K <- as.binmat(c("0000", "1000", "0100", "1110", "1101", "1111"))

set.seed(1234)
info <- blimit(K)
}
\keyword{models}
