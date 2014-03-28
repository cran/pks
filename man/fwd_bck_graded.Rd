\name{gradedness}
\alias{gradedness}
\alias{is.forward.graded}
\alias{is.backward.graded}
\title{Forward- or Backward-Gradedness of a Knowledge Structure}
\description{
  Checks if a knowledge structure is forward- or backward-graded in any
  item.
}
\usage{
is.forward.graded(K)

is.backward.graded(K)
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
}

\value{
  A named logical vector with as many elements as columns in \code{K}.
}

\references{
  Spoto, A., Stefanutti, L., & Vidotto, G. (2012).
  On the unidentifiability of a certain class of skill multi map based
  probabilistic knowledge structures.
  \emph{Journal of Mathematical Psychology}, \bold{56}, 248--255.
}

\seealso{
  \code{\link{blim}}.
}
\examples{
K <- as.binmat(c("0000", "1000", "1100", "1010", "0110", "1110", "1111"))
is.forward.graded(K)   # forward-graded in a
is.backward.graded(K)  # not backward-graded in a
}
\keyword{models}
