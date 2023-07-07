\name{conversion}
\alias{conversion}
\alias{as.pattern}
\alias{as.binmat}
\alias{is.subset}
\title{Conversion between Representations of Responses or States}
\description{
  Converts between binary matrix and pattern representations of response
  patterns or knowledge states.
}
\usage{
as.pattern(R, freq = FALSE, useNames = FALSE, as.set = FALSE,
           sep = "", emptyset = "{}", as.letters = NULL)

as.binmat(N.R, uniq = TRUE, col.names = NULL, as.logical = FALSE)

is.subset(R)
}
\arguments{
  \item{R}{an indicator matrix of response patterns or knowledge states.}
  \item{N.R}{either a (named) vector of absolute frequencies of response
    patterns; or a character vector of response patterns or knowledge states;
    or a \code{set} of sets representing the knowledge structure.}
  \item{freq}{logical, should the frequencies of response patterns be
    reported?}
  \item{uniq}{logical, if \code{TRUE}, only the unique response patterns are
    returned.}
  \item{useNames}{logical, return response patterns as combinations of
    item names.}
  \item{as.set}{logical, return response patterns as set of sets.}
  \item{sep}{character to separate the item names.}
  \item{emptyset}{string representing the empty set if \code{useNames}
    is \code{TRUE}.}
  \item{as.letters}{deprecated, use \code{useNames} instead.}
  \item{col.names}{column names for the state or response matrix.}
  \item{as.logical}{logical, return logical matrix of states.}
}
\value{
  \code{as.pattern} returns a vector of integers named by the response
  patterns if \code{freq} is \code{TRUE}, else a character vector. If
  \code{as.set} is \code{TRUE}, the return value is of class \code{set}.

  \code{as.binmat} returns an indicator matrix. If \code{as.logical} is
  \code{TRUE}, it returns a logical matrix.

  \code{is.subset} returns a logical incidence matrix of the subset relation
  among states.
}
\seealso{
  \code{\link{blim}}, \code{set} in package \code{sets}.
}
\examples{
data(DoignonFalmagne7)
K <- DoignonFalmagne7$K
as.pattern(K, freq = TRUE)
as.pattern(K)
as.pattern(K, useNames = TRUE)
as.pattern(K, as.set = TRUE)

N.R <- DoignonFalmagne7$N.R
dim(as.binmat(N.R))
dim(as.binmat(N.R, uniq = FALSE))

## Knowledge structure as binary matrix
as.binmat(c("000", "100", "101", "111"))
as.binmat(set(set(), set("a"), set("a", "c"), set("a", "b", "c")))
as.binmat(c("000", "100", "101", "111"), as.logical = TRUE)

## Subset relation incidence matrix
is.subset(K)

## Plotting the knowledge structure
if(requireNamespace("relations") &&
   requireNamespace("Rgraphviz")) {
  rownames(K) <- as.pattern(K, useNames = TRUE)
  plot(relations::as.relation(is.subset(K)), main = "")
}
}
\keyword{models}
