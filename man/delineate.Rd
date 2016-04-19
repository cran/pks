\name{delineate}
\alias{delineate}
\title{Delineate a Knowledge Structure by a Skill Function}
\description{
  Computes the knowledge structure delineated by a skill function.
}
\usage{
delineate(skillfun, itemID = 1)
}
\arguments{
  \item{skillfun}{a data frame or a matrix representing the skill function.
    It consists of an item indicator and a problem-by-skill indicator matrix.}
  \item{itemID}{index of the column in \code{skillfun} that holds the item
    indicator.}
}

\details{
  The skill function \eqn{(Q, S, \mu)} indicates for each item in \eqn{Q}
  which subsets of skills in \eqn{S} are required to solve the item.  Thus,
  \eqn{\mu(q)} is a set containing sets of skills.  An item may have multiple
  entries in \code{skillfun}, each in a separate row identified by the same
  \code{itemID}.

  See Doignon and Falmagne (1999, Chap. 4).
}

\value{
  A list of two components:
  \item{K}{the knowledge structure delineated by the skill function.}
  \item{classes}{a list of equivalence classes of competence states; the
    members of these classes are mapped onto the same knowledge state by the
    problem function induced by the skill function \eqn{\mu}.}
}

\references{
  Doignon, J.-P., & Falmagne, J.-C. (1999).
  \emph{Knowledge spaces}. Berlin: Springer.
}

\seealso{
  \code{\link{blim}}.
}

\examples{
# Skill function
# mu(e) = {{s, t}, {s, u}},  mu(f) = {{u}}
# mu(g) = {{s}, {t}},        mu(h) = {{t}}
sf <- read.table(header = TRUE, text = "
item s t u
   e 1 1 0
   e 1 0 1
   f 0 0 1
   g 1 0 0
   g 0 1 0
   h 0 1 0
")
delineate(sf)

## See ?probability for further examples.
}
\keyword{models}
