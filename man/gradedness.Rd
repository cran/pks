\name{gradedness}
\alias{gradedness}
\alias{is.forward.graded}
\alias{is.backward.graded}
\alias{is.downgradable}
\title{Forward-/Backward-Gradedness and Downgradability of a Knowledge
       Structure}
\description{
  Checks if a knowledge structure is
  \itemize{
    \item{forward- or backward-graded in any item;}
    \item{downgradable.}
  }
}
\usage{
is.forward.graded(K)

is.backward.graded(K)

is.downgradable(K)
}
\arguments{
  \item{K}{a state-by-problem indicator matrix representing the knowledge
    structure.  An element is one if the problem is contained in the state,
    and else zero.  \code{K} should have non-empty \code{colnames}.}
}

\details{
  A knowledge structure \eqn{K} is forward-graded in item \eqn{q}, if
  \eqn{S \cup \{q\}} is in \eqn{K} for every state \eqn{S \in K}.
  A knowledge structure \eqn{K} is backward-graded in item \eqn{q}, if
  \eqn{S - \{q\}} is in \eqn{K} for every state \eqn{S \in K}.
  See Spoto, Stefanutti, and Vidotto (2012).

  A knowledge structure \eqn{K} is downgradable, if its inner fringe is empty
  only for a single state (the empty set).  See Doignon and Falmagne (2015).
}

\value{
  For forward- and backward-gradedness, a named logical vector with as many
  elements as columns in \code{K}.

  For downgradability, a single logical value.
}

\references{
  Doignon, J.-P., & Falmagne, J.-C. (2015).
  Knowledge spaces and learning spaces.
  \emph{arXiv}.
  \doi{https://doi.org/10.48550/arXiv.1511.06757}

  Spoto, A., Stefanutti, L., & Vidotto, G. (2012).
  On the unidentifiability of a certain class of skill multi map based
  probabilistic knowledge structures.
  \emph{Journal of Mathematical Psychology}, \bold{56}(4), 248--255.
  \doi{10.1016/j.jmp.2012.05.001}
}

\seealso{
  \code{\link{blim}}, \code{\link{jacobian}}, \code{\link{getKFringe}}.
}
\examples{
K <- as.binmat(c("0000", "1000", "1100", "1010", "0110", "1110", "1111"))
is.forward.graded(K)                  # forward-graded in a
is.backward.graded(K)                 # not backward-graded in a
is.downgradable(K)                    # not downgradable
all(K[, "a"] | getKFringe(K)[, "a"])  # every K or outer fringe contains a
}
\keyword{models}
