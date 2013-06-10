\name{plot.blim}
\alias{plot.blim}
\title{Diagnostic Plot for Basic Local Independence Models}
\description{
  Plots BLIM residuals against fitted values.
}
\usage{
\method{plot}{blim}(x, xlab = "Predicted response probabilities",
     ylab = "Deviance residuals", \dots)
}
\arguments{
  \item{x}{an object of class \code{blim}, typically the result of a
    call to \code{\link{blim}}.}
  \item{xlab, ylab, \dots}{graphical parameters passed to plot.}
}
\details{
  The deviance residuals are plotted against the predicted response
  probabilities for each response pattern.
}
\seealso{\code{\link{blim}}, \code{\link{residuals.blim}}.}
\examples{
## Compare MD and MDML estimation

data(DoignonFalmagne7)
blim1 <- blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R, method="MD")
blim2 <- blim(DoignonFalmagne7$K, DoignonFalmagne7$N.R, method="MDML")

par(mfrow = 1:2)      # residuals versus fitted values
plot(blim1, main = "MD estimation",   ylim = c(-4, 4))
plot(blim2, main = "MDML estimation", ylim = c(-4, 4))
}
\keyword{models}
