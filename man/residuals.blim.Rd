\name{residuals.blim}
\alias{residuals.blim}
\title{Residuals for Basic Local Independence Models}
\description{
  Computes deviance and Pearson residuals for \code{blim} objects.
}
\usage{
\method{residuals}{blim}(object, type = c("deviance", "pearson"), \dots)
}
\arguments{
  \item{object}{an object of class \code{blim}, typically the result of a
    call to \code{\link{blim}}.}
  \item{type}{the type of residuals which should be returned; the
     alternatives are: \code{"deviance"} (default) and \code{"pearson"}.}
  \item{\dots}{further arguments passed to or from other methods.
    None are used in this method.}
}
\details{
  See \code{\link{residuals.glm}} for details.
}
\value{
  A named vector of residuals having as many elements as response patterns.
}
\seealso{\code{\link{blim}}, \code{\link{residuals.glm}},
  \code{\link{plot.blim}}.}
\examples{
data(DoignonFalmagne7)
blim1 <- blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R)

sum( resid(blim1)^2 )                # likelihood ratio G2
sum( resid(blim1, "pearson")^2 )     # Pearson X2
}
\keyword{models}
