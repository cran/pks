\name{blim}
\alias{blim}
\alias{blimMD}
\alias{logLik.blim}
\alias{nobs.blim}

\title{Basic Local Independence Models (BLIMs)}
\description{Fits a basic local independence model (BLIM) for probabilistic
  knowledge structures by minimum discrepancy maximum likelihood estimation.
}
\usage{
blim(K, N.R, method = c("MD", "ML", "MDML"), R = as.binmat(N.R),
  P.K = rep(1/nstates, nstates),
  beta = if(errequal) 0.1 else rep(0.1, nitems),
   eta = if(errequal) 0.1 else rep(0.1, nitems),
  errtype = c("both", "error", "guessing"),
  errequal = FALSE, incradius = 0, tol = 1e-7, maxiter = 10000)

blimMD(K, N.R, R = as.binmat(N.R), errtype = c("both", "error",
  "guessing"), incrule = c("minimum", "hypblc1", "hypblc2"), m = 1)
}

\arguments{
  \item{K}{a state-by-problem indicator matrix representing the knowledge
    structure.  An element is one if the problem is contained in the state,
    and else zero.}
  \item{N.R}{a (named) vector of absolute frequencies of response patterns.}
  \item{method}{\code{MD} for minimum discrepancy estimation, \code{ML} for
    maximum likelihood estimation, \code{MDML} for minimum discrepancy
    maximum likelihood estimation.}
  \item{R}{a person-by-problem indicator matrix of unique response patterns.
    Per default inferred from the names of \code{N.R}.}
  \item{P.K}{the vector of initial parameter values for probabilities of
    knowledge states.}
  \item{beta}{the vector of initial parameter values for probabilities of a
    careless error.}
  \item{eta}{the vector of initial parameter values for probabilities of a
    lucky guess.}
  \item{errtype}{type of response errors that can occur: \code{error} for
    careless errors only, \code{guessing} for lucky guesses only, and
    \code{both} for both error types.}
  \item{errequal}{logical, if \code{TRUE} then items have identical error
    rates.}
  \item{incradius}{include knowledge states of distance from the minimum
    discrepant states less than or equal to \code{incradius}.}
  \item{tol}{tolerance, stopping criterion for iteration.}
  \item{maxiter}{the maximum number of iterations.}
  \item{incrule}{inclusion rule for knowledge states.  (See Details.)}
  \item{m}{exponent for hyperbolic inclusion rules.}
}

\details{
  See Doignon and Falmagne (1999) for details on the basic local independence
  model (BLIM) for probabilistic knowledge structures.

  Minimum discrepancy (MD) minimizes the number of expected response errors
  (careless errors or lucky guesses). Maximum likelihood maximizes the
  likelihood, possibly at the expense of inflating the error and guessing
  parameters.  Minimum discrepancy maximum likelihood (MDML) maximizes the
  likelihood subject to the constraint of minimum response errors.  See Heller
  and Wickelmaier (2004) for details on the parameter estimation methods.

  \code{blimMD} uses minimum discrepancy estimation only.  Apart from the
  hyperbolic inclusion rules, all of its functionality is also provided by
  \code{blim}.  It may be removed in the future.
}

\value{
  An object of class \code{blim} having the following components:
  \item{discrepancy}{the mean discrepancy.}
  \item{P.K}{the vector of estimated parameter values for probabilities of
    knowledge states.}
  \item{beta}{the vector of estimated parameter values for probabilities of
    a careless error.}
  \item{eta}{the vector of estimated parameter values for probabilities of a
    lucky guess.}
  \item{disc.tab}{the matrix of symmetric set distances.}
  \item{K}{the knowledge structure.}
  \item{N.R}{the vector of frequencies of response patterns.}
  \item{nitems}{the number of items.}
  \item{nstates}{the number of knowledge states.}
  \item{npatterns}{the number of response patterns.}
  \item{ntotal}{the number of respondents.}
  \item{nerror}{the number of response errors.}
  \item{npar}{the number of parameters.}
  \item{errtype}{the type of response errors.}
  \item{method}{the parameter estimation method.}
  \item{iter}{the number of iterations needed.}
  \item{loglik}{the log-likelihood.}
  \item{fitted.values}{the fitted response frequencies.}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (G2), the degrees of
    freedom, and the p-value of the corresponding chi-square distribution.}
}

\references{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.

  Heller, J., & Wickelmaier, F. (2004).
  Parameter estimation in probabilistic knowledge structures: the
  minimum-discrepancy method.
  Presented at the
  \emph{35th European Mathematical Psychology Group (EMPG) Meeting}, September
  2-4, Ghent, Belgium. 
}

\seealso{
  \code{\link{simulate.blim}}, \code{\link{plot.blim}},
  \code{\link{residuals.blim}}, \code{\link{logLik.blim}}.
}

\examples{
data(DoignonFalmagne7)
K   <- DoignonFalmagne7$K    # knowledge structure
N.R <- DoignonFalmagne7$N.R  # frequencies of response patterns

## Fit basic local independence model (BLIM) by different methods
blim(K, N.R, method="MD")    # minimum discrepancy estimation
blim(K, N.R, method="ML")    # maximum likelihood estimation by EM algorithm
blim(K, N.R, method="MDML")  # MDML estimation
}
\keyword{models}
