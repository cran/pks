\name{Taagepera}
\alias{Taagepera}
\alias{density97}
\alias{matter97}
\title{Responses and Knowledge Structures from Taagepera et al. (1997)}
\description{
  Taagepera et al. (1997) applied knowledge space theory to specific science
  problems.  The density test was administered to 2060 students, the
  conservation of matter test to 1620 students.  A subtest of five items each
  is included here.  The response frequencies were reconstructed from
  histograms in the paper.
}
\usage{data(Taagepera)}
\format{
  Two lists, each consisting of two components:
  \describe{
    \item{\code{density97}}{a list with components \code{K} and \code{N.R} for
      the density test.}
    \item{\code{matter97}}{a list with components \code{K} and \code{N.R} for
      the conservation of matter test.}
    \item{\code{K}}{a state-by-problem indicator matrix representing the
      hypothetical knowledge structure.  An element is one if the problem
      is contained in the state, and else zero.}
    \item{\code{N.R}}{a named numeric vector.  The names denote response
      patterns, the values denote their frequencies.}
  }
}
\source{
  Taagepera, M., Potter, F., Miller, G.E., & Lakshminarayan, K. (1997).
  Mapping students' thinking patterns by the use of knowledge space theory.
  \emph{International Journal of Science Education}, \bold{19}, 283--302.
}
\examples{
data(Taagepera)
density97$K     # density test knowledge structure
density97$N.R   # density test response patterns
matter97$K      # conservation of matter knowledge structure
matter97$N.R    # conservation of matter response patterns
}
\keyword{datasets}
