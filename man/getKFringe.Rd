\name{getKFringe}
\alias{getKFringe}
\title{Outer and Inner Fringes of a Knowledge Structure}
\description{
  Returns the outer or inner fringe for each state in a knowledge structure.
}
\usage{
getKFringe(K, nstates = nrow(K), nitems = ncol(K), outer = TRUE)
}
\arguments{
  \item{K}{a state-by-problem indicator matrix representing the knowledge
    structure.  An element is one if the problem is contained in the state,
    and else zero.}
  \item{nstates}{the number of knowledge states in \code{K}.}
  \item{nitems}{the number of items in \code{K}.}
  \item{outer}{logical. If \code{TRUE} return outer fringe, else return inner
    fringe.}
}
\details{
  The outer fringe of a knowledge state is the set of all items that can be
  learned from that state, such that adding an outer-fringe item to the state
  results in another state in \code{K},

  \deqn{K^O = \{q \notin K | K \cup \{q\} \in \mathcal{K}\}.}

  The inner fringe of a knowledge state is the set of all items that have been
  learned most recently to reach that state, such that deleting an
  inner-fringe item from the state results in another state in \code{K},

  \deqn{K^I = \{q \in K | K - \{q\} \in \mathcal{K}\}.}
}
\value{
  A state-by-problem indicator matrix representing the outer or inner fringe
  for each knowledge state in \code{K}.
}
\seealso{
  \code{\link{slm}}, \code{\link{simulate.blim}}.
}
\examples{
data(DoignonFalmagne7)

## Which items can be learned from each state?
getKFringe(DoignonFalmagne7$K)

## Which items in each state have been recently learned?
getKFringe(DoignonFalmagne7$K, outer = FALSE)
}
\keyword{models}
