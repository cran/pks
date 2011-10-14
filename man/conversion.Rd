\name{conversion}
\alias{conversion}
\alias{as.pattern}
\alias{as.binmat}
\title{Conversion between Representations of Responses or States}
\description{
  Converts between binary matrix and pattern representations of response
  patterns or knowledge states.
}
\usage{
as.pattern(R, freq = FALSE)

as.binmat(N.R, uniq = TRUE, col.names = NULL)
}
\arguments{
  \item{R}{an indicator matrix of response patterns or knowledge states.}
  \item{N.R}{either a (named) vector of absolute frequencies of response
    patterns, or a character vector of response patterns or knowledge states.}
  \item{freq}{logical, should the frequency of response patterns be
    reported?}
  \item{uniq}{logical, if \code{TRUE}, only the unique response patterns are
    returned.}
  \item{col.names}{column names for the state or response matrix.}
}
\value{
  \code{as.pattern} returns a vector of integers named by the response
  patterns if \code{freq} is \code{TRUE}, else a character vector.

  \code{as.binmat} returns an indicator matrix.
}
\seealso{
  \code{\link{blim}}.
}
\examples{
data(DoignonFalmagne7)
as.pattern(DoignonFalmagne7$K)
as.pattern(DoignonFalmagne7$K, freq = TRUE)

dim(as.binmat(DoignonFalmagne7$N.R))
dim(as.binmat(DoignonFalmagne7$N.R, uniq = FALSE))

## Knowledge structure as binary matrix
as.binmat(c("000", "100", "101", "111"))
}
\keyword{models}
