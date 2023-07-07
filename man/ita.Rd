\name{ita}
\alias{ita}
\alias{print.ita}
\title{Item Tree Analysis (ITA)}
\description{
  Item tree analysis (ITA) on a set of binary responses.
}
\usage{
ita(R, L = NULL, makeK = FALSE, search = c("local", "global"))
}
\arguments{
  \item{R}{a subject-by-problem indicator matrix representing the responses.}
  \item{L}{the threshold of violations acceptable for the precedence relation.
    If \code{NULL} (default), an optimal threshold is searched for.}
  \item{makeK}{should the corresponding knowledge structure be returned?}
  \item{search}{local (default) or global threshold search.}
}

\details{
  ITA seeks to establish a precedence relation among a set of binary items.
  For each pair of items \eqn{(p, q)}, it counts how often \eqn{p} is not
  solved if \eqn{q} is solved, which constitutes a violation of the relation.
  ITA searches for a threshold \code{L} for the maximum number of violations
  consistent with a (transitive) precedence relation.  Its attempts to
  minimize the total discrepancy between \code{R} and \code{K}.

  See van Leeuwe (1974) and Schrepp (1999) for details.
}

\value{
  An object of class \code{ita} having the following components:
  \item{K}{the knowledge structure corresponding to the precedence relation.}
  \item{discrepancy}{the discrepancy between \code{R} and \code{K} (fit),
    between \code{K} and \code{R} (complexity), and their sum (total).}
  \item{transitiveL}{the vector of transitive thresholds.}
  \item{searchL}{either \code{NULL} or the method used for threshold search.}
  \item{L}{the selected or requested threshold.}
  \item{P}{the precedence matrix containing the number of violations.}
  \item{I}{the precedence relation as a logical incidence matrix at threshold
    \code{L}.}
}

\references{
  Schrepp, M. (1999).
  On the empirical construction of implications between bi-valued test items.
  \emph{Mathematical Social Sciences},
  \bold{38}(3), 361--375.
  \doi{10.1016/S0165-4896(99)00025-6}

  Van Leeuwe, J.F. (1974).
  Item tree analysis.
  \emph{Nederlands Tijdschrift voor de Psychologie en haar Grensgebieden},
  \bold{29}(6), 475--483.
}

\seealso{
  \code{\link{blim}}.
}

\examples{
data(chess)

ita(chess$R)  # find (locally) optimal threshold L

i <- ita(chess$R, L = 6, makeK = TRUE)
identical(sort(as.pattern(i$K)),
          sort(as.pattern(chess$dst1)))

## Plotting the precedence relation
if(requireNamespace("relations") &&
   requireNamespace("Rgraphviz")) {
  plot(relations::as.relation(i$I))
}
}
\keyword{models}
