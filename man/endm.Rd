\name{endm}
\alias{endm}
\title{Responses and Knowledge Structures from Heller and Wickelmaier (2013)}
\description{
  Knowledge structures and 200 artificial responses to four problems are used
  to illustrate parameter estimation in Heller and Wickelmaier (2013).
}
\usage{data(endm)}
\format{
  A list consisting of three components:
  \describe{
    \item{\code{K}}{a state-by-problem indicator matrix representing the
      true knowledge structure that underlies the model that generated the
      data.}
    \item{\code{K2}}{a slightly misspecified knowledge structure.}
    \item{\code{N.R}}{a named numeric vector.  The names denote response
      patterns, the values denote their frequencies.}
  }
}
\source{
  Heller, J., & Wickelmaier, F. (2013).
  Minimum discrepancy estimation in probabilistic knowledge structures.
  \emph{Electronic Notes in Discrete Mathematics}, \bold{42}, 49--56.
  \doi{10.1016/j.endm.2013.05.145}
}
\examples{
data(endm)
endm$K    # true knowledge structure
endm$K2   # misspecified knowledge structure
endm$N.R  # response patterns

## Generate data from BLIM based on K
blim0 <- list(
     P.K = setNames(c(.1, .15, .15, .2, .2, .1, .1), as.pattern(endm$K)),
    beta = rep(.1, 4),
     eta = rep(.1, 4),
       K = endm$K,
  ntotal = 200)
class(blim0) <- "blim"
simulate(blim0)

## Fit BLIM based on K2
blim1 <- blim(endm$K2, endm$N.R, "MD")
}
\keyword{datasets}
