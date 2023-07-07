\name{slm}
\alias{slm}
\alias{getSlmPK}
\alias{coef.slm}
\alias{print.slm}

\title{Simple Learning Models (SLMs)}
\description{Fits a simple learning model (SLM) for probabilistic knowledge
  structures by minimum discrepancy maximum likelihood estimation.
}
\usage{
slm(K, N.R, method = c("MD", "ML", "MDML"), R = as.binmat(N.R),
    beta = rep(0.1, nitems), eta = rep(0.1, nitems),
    g = rep(0.1, nitems),
    betafix = rep(NA, nitems), etafix = rep(NA, nitems),
    betaequal = NULL, etaequal = NULL,
    randinit = FALSE, incradius = 0,
    tol = 1e-07, maxiter = 10000, zeropad = 16,
    checkK = TRUE)

getSlmPK(g, K, Ko)

\method{print}{slm}(x, P.Kshow = FALSE, parshow = TRUE,
      digits=max(3, getOption("digits") - 2), \dots)
}

\arguments{
  \item{K}{a state-by-problem indicator matrix representing the knowledge
    space.  An element is one if the problem is contained in the state,
    and else zero.}
  \item{N.R}{a (named) vector of absolute frequencies of response patterns.}
  \item{method}{\code{MD} for minimum discrepancy estimation, \code{ML} for
    maximum likelihood estimation, \code{MDML} for minimum discrepancy
    maximum likelihood estimation.}
  \item{R}{a person-by-problem indicator matrix of unique response patterns.
    Per default inferred from the names of \code{N.R}.}
  \item{beta, eta, g}{vectors of initial values for the error, guessing, and
    solvability parameters.}
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
  \item{checkK}{logical, if \code{TRUE} K is checked for well-gradedness.}
  \item{Ko}{a state-by-problem indicator matrix representing the outer fringe
    for each knowledge state in \code{K}; typically the result of a call to
    \code{getKFringe}.}
  \item{x}{an object of class \code{slm}, typically the result of a call to
    \code{slm}.}
  \item{P.Kshow}{logical, should the estimated distribution of knowledge
    states be printed?}
  \item{parshow}{logical, should the estimates of error, guessing, and
    solvability parameters be printed?}
  \item{digits}{a non-null value for \code{digits} specifies the minimum
    number of significant digits to be printed in values.}
  \item{\dots}{additional arguments passed to other methods.}
}

\details{
  See Doignon and Falmagne (1999) for details on the simple learning model
  (SLM) for probabilistic knowledge structures. The model requires a
  well-graded knowledge space \code{K}.

  An \code{slm} object inherits from class \code{blim}.  See \code{blim} for
  details on the function arguments.  The helper function \code{getSlmPK}
  returns the distribution of knowledge states \code{P.K}.
}

\value{
  An object of class \code{slm} and \code{blim}.  It contains all components
  of a \code{blim} object.  In addition, it includes:
  \item{g}{the vector of estimates of the solvability parameters.}
}

\references{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.
}

\seealso{
  \code{\link{blim}}, \code{\link{simulate.blim}}, \code{\link{getKFringe}},
  \code{\link{is.downgradable}}
}

\examples{
data(DoignonFalmagne7)
K   <- DoignonFalmagne7$K     # well-graded knowledge space
N.R <- DoignonFalmagne7$N.R   # frequencies of response patterns

## Fit simple learning model (SLM) by different methods
slm(K, N.R, method = "MD")    # minimum discrepancy estimation
slm(K, N.R, method = "ML")    # maximum likelihood estimation by EM
slm(K, N.R, method = "MDML")  # MDML estimation

## Compare SLM and BLIM
m1 <-  slm(K, N.R, method = "ML")
m2 <- blim(K, N.R, method = "ML")
anova(m1, m2)
}
\keyword{models}
