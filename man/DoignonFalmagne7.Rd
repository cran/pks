\name{DoignonFalmagne7}
\alias{DoignonFalmagne7}
\title{Artificial Responses from Doignon and Falmagne (1999)}
\description{
  Fictitious data set from Doignon and Falmagne (1999, chap. 7).
  Response patterns of 1000 respondents to five problems.  Each
  respondent is assumed to be in one of nine possible states of
  the knowledge structure \code{K}.
}
\usage{data(DoignonFalmagne7)}
\format{
  A list consisting of two components:
  \describe{
    \item{\code{K}}{a state-by-problem indicator matrix representing the
      hypothetical knowledge structure.  An element is one if the problem
      is contained in the state, and else zero.}
    \item{\code{N.R}}{a named numeric vector.  The names denote response
      patterns, the values denote their frequencies.}
  }
}
\source{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.
}
\examples{
data(DoignonFalmagne7)
DoignonFalmagne7$K    # knowledge structure
DoignonFalmagne7$N.R  # response patterns
}
\keyword{datasets}
