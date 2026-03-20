\name{hsgeometry}
\alias{hsgeometry}
\alias{angles}
\alias{circles}
\title{Responses to Geometry Problems and Knowledge Structures}
\description{
  Lakshminarayan (1995, 1996) presented high school geometry problems to 959
  undergraduate students before and after a lesson. The problems deal with
  basic concepts in geometry, such as angles, lengths, circles, and polygons.
  Doignon and Falmagne (1999, chap. 8) refer to a subset of these data.
}
\usage{data(hsgeometry)}
\format{
  Two lists, each consisting of three components:
  \describe{
    \item{\code{angles}}{a list with components \code{K}, \code{N.R}, and
      \code{changescore} for five problems mostly about angles.}
    \item{\code{circles}}{a list with components \code{K}, \code{N.R}, and
      \code{changescore} for four problems mostly about circles.}
    \item{\code{K}}{a state-by-problem indicator matrix representing the
      knowledge structure in Lakshminarayan (1995, 1996).}
    \item{\code{N.R}}{a matrix. The row names denote response patterns, the
      column entries their frequencies before and after the lesson.}
    \item{\code{changescore}}{a 3d table. For each problem, it contains the
      number of students who solved or did not solve the problem at pre- and
      posttest.}
  }
}
\note{
  Only those 800 responses used to estimate the model parameters are reported
  (the held-out sample is missing).

  Thanks to Kamakshi Lakshminarayan, Cord Hockemeyer, and Mike Regenwetter for
  their help with finding the data.
}
\source{
  Lakshminarayan, K (1996).
  \emph{A hybrid latent trait and latent class model of learning - Theoretical
  details and empirical application}.
  Technical report No. MBS 96-07. University of Califomia, Irvine.
}
\references{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.

  Lakshminarayan, K. (1995).
  \emph{Theoretical and empirical aspects of some stochastic learning models}.
  Doctoral dissertation. University of California, Irvine.
}
\examples{
data(hsgeometry)
dotchart(angles$N.R[, "posttest"], pch = 4,
         main = "High school geometry problems",
         xlab = "Response frequency (pre- [o] and posttest [x])")
points(angles$N.R[, "pretest"], 1:32)
mtext("Lakshminarayan (1995)", side = 3, line = 0.5)
}
\keyword{datasets}
