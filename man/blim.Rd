\name{blim}
\alias{blim}
\alias{blimMD}
\alias{anova.blim}
\alias{coef.blim}
\alias{deviance.blim}
\alias{logLik.blim}
\alias{nobs.blim}

\title{Basic Local Independence Models (BLIMs)}
\description{Fits a basic local independence model (BLIM) for probabilistic
  knowledge structures by minimum discrepancy maximum likelihood estimation.
}
\usage{
blim(K, N.R, method = c("MD", "ML", "MDML"), R = as.binmat(N.R),
     P.K = rep(1/nstates, nstates),
     beta = rep(0.1, nitems), eta = rep(0.1, nitems),
     betafix = rep(NA, nitems), etafix = rep(NA, nitems),
     betaequal = NULL, etaequal = NULL,
     randinit = FALSE, incradius = 0,
     tol = 1e-07, maxiter = 10000, zeropad = 16)

blimMD(K, N.R, R = as.binmat(N.R),
       betafix = rep(NA, nitems), etafix = rep(NA, nitems),
       incrule = c("minimum", "hypblc1", "hypblc2"), m = 1)

\method{anova}{blim}(object, \dots, test = c("Chisq", "none"))
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
  \item{beta, eta}{vectors of initial parameter values for probabilities of a
    careless error and a lucky guess, respectively.}
  \item{betafix, etafix}{vectors of fixed error and guessing parameter values;
    \code{NA} indicates a free parameter.}
  \item{betaequal, etaequal}{lists of vectors of problem indices; each vector
    represents an equivalence class: it contains the indices of problems for
    which the error or guessing parameters are constrained to be equal.  (See
    Examples.)}
  \item{randinit}{logical, if \code{TRUE} then initial parameter values are
    sampled uniformly with constraints.  (See Details.)}
  \item{incradius}{include knowledge states of distance from the minimum
    discrepant states less than or equal to \code{incradius}.}
  \item{tol}{tolerance, stopping criterion for iteration.}
  \item{maxiter}{the maximum number of iterations.}
  \item{zeropad}{the maximum number of items for which an incomplete
    \code{N.R} vector is completed and padded with zeros.}
  \item{incrule}{inclusion rule for knowledge states.  (See Details.)}
  \item{m}{exponent for hyperbolic inclusion rules.}
  \item{object}{an object of class \code{blim}, typically the result of a
    call to \code{blim}.}
  \item{test}{should the p-values of the chi-square distributions be
    reported?}
  \item{\dots}{additional arguments passed to other methods.}
}

\details{
  See Doignon and Falmagne (1999) for details on the basic local independence
  model (BLIM) for probabilistic knowledge structures.

  Minimum discrepancy (MD) minimizes the number of expected response errors
  (careless errors or lucky guesses). Maximum likelihood maximizes the
  likelihood, possibly at the expense of inflating the error and guessing
  parameters.  Minimum discrepancy maximum likelihood (MDML) maximizes the
  likelihood subject to the constraint of minimum response errors.  See Heller
  and Wickelmaier (2013) for details on the parameter estimation methods.

  If \code{randinit} is \code{TRUE}, initial parameter values are sampled
  uniformly with the constraint \code{beta + eta < 1} (Weisstein, 2013) for
  the error parameters, and with \code{sum(P.K) == 1} (Rubin, 1981) for the
  probabilities of knowledge states. Setting \code{randinit} to \code{TRUE}
  overrides any values given in the \code{P.K}, \code{beta}, and \code{eta}
  arguments.

  The degrees of freedom in the goodness-of-fit test are calculated as number
  of possible response patterns minus one or number of respondents, whichever
  is smaller, minus number of parameters.

  \code{blimMD} uses minimum discrepancy estimation only.  Apart from the
  hyperbolic inclusion rules, all of its functionality is also provided by
  \code{blim}.  It may be removed in the future.
}

\value{
  An object of class \code{blim} having the following components:
  \item{discrepancy}{the mean minimum discrepancy between response patterns
    and knowledge states.}
  \item{P.K}{the vector of estimated parameter values for probabilities of
    knowledge states.}
  \item{beta}{the vector of estimated parameter values for probabilities of
    a careless error.}
  \item{eta}{the vector of estimated parameter values for probabilities of a
    lucky guess.}
  \item{disc.tab}{the minimum discrepancy distribution.}
  \item{K}{the knowledge structure.}
  \item{N.R}{the vector of frequencies of response patterns.}
  \item{nitems}{the number of items.}
  \item{nstates}{the number of knowledge states.}
  \item{npatterns}{the number of response patterns.}
  \item{ntotal}{the number of respondents.}
  \item{nerror}{the number of response errors.}
  \item{npar}{the number of parameters.}
  \item{method}{the parameter estimation method.}
  \item{iter}{the number of iterations needed.}
  \item{loglik}{the log-likelihood.}
  \item{fitted.values}{the fitted response frequencies.}
  \item{goodness.of.fit}{the goodness of fit statistic including the
    likelihood ratio fitted vs. saturated model (G2), the degrees of
    freedom, and the p-value of the corresponding chi-square distribution.
    (See Details.)}
}

\references{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.

  Heller, J., & Wickelmaier, F. (2013).
  Minimum discrepancy estimation in probabilistic knowledge structures.
  \emph{Electronic Notes in Discrete Mathematics}, \bold{42}, 49--56.
  \doi{10.1016/j.endm.2013.05.145}

  Rubin, D.B. (1981). The Bayesian bootstrap.
  \emph{The Annals of Statistics}, \bold{9}(1), 130--134.
  \doi{10.1214/aos/1176345338}

  Weisstein, E.W. (2013, August 29). Triangle point picking.
  In \emph{MathWorld -- A Wolfram Web Resource}.
  Retrieved from
  \url{https://mathworld.wolfram.com/TrianglePointPicking.html}.
}

\seealso{
  \code{\link{simulate.blim}}, \code{\link{plot.blim}},
  \code{\link{residuals.blim}}, \code{\link{logLik.blim}},
  \code{\link{delineate}}, \code{\link{jacobian}}, \code{\link{endm}},
  \code{\link{probability}}, \code{\link{chess}}.
}

\examples{
data(DoignonFalmagne7)
K   <- DoignonFalmagne7$K    # knowledge structure
N.R <- DoignonFalmagne7$N.R  # frequencies of response patterns

## Fit basic local independence model (BLIM) by different methods
blim(K, N.R, method = "MD")    # minimum discrepancy estimation
blim(K, N.R, method = "ML")    # maximum likelihood estimation by EM
blim(K, N.R, method = "MDML")  # MDML estimation

## Parameter restrictions: beta_a = beta_b = beta_d, beta_c = beta_e
##                          eta_a =  eta_b = 0.1
m1 <- blim(K, N.R, method = "ML",
           betaequal = list(c(1, 2, 4), c(3, 5)),
              etafix = c(0.1, 0.1, NA, NA, NA))
m2 <- blim(K, N.R, method = "ML")
anova(m1, m2)

## See ?endm, ?probability, and ?chess for further examples.
}
\keyword{models}
