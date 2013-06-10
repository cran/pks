\name{print.blim}
\alias{print.blim}
\title{Print a blim Object}
\description{
  Prints the output of a \code{blim} model object.
}
\usage{
\method{print}{blim}(x, P.Kshow = FALSE, errshow = TRUE,
      digits=max(3, getOption("digits") - 2), \dots)
}
\arguments{
  \item{x}{an object of class \code{blim}, typically the result of a call to
    \code{\link{blim}}.}
  \item{P.Kshow}{logical, should the estimated distribution of knowledge
    states be printed?}
  \item{errshow}{logical, should the estimates of careless error and lucky
    guess parameters be printed?}
  \item{digits}{a non-null value for \code{digits} specifies the minimum
    number of significant digits to be printed in values.}
  \item{\dots}{further arguments passed to or from other methods.
    None are used in this method.}
}
\value{
  Returns the \code{blim} object invisibly.
}
\seealso{
  \code{\link{blim}}.
}
\examples{
data(DoignonFalmagne7)
 
blim1 <- blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R)
print(blim1, showP.K = TRUE)
}
\keyword{models}
