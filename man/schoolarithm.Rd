\name{schoolarithm}
\alias{schoolarithm}
\alias{fraction17}
\alias{subtraction13}
\title{Arithmetic Problems for Elementary and Middle School Students}
\description{
  The 23 fraction problems were presented to 191 first-level middle school
  students (about 11 to 12 years old).  A subset of 13 problems is included in
  Stefanutti and de Chiusole (2017).

  The eight subtraction problems were presented to 294 elementary school
  students and are described in de Chiusole and Stefanutti (2013).
}
\usage{data(schoolarithm)}
\format{
  \describe{
    \item{\code{fraction17}}{a person-by-problem indicator matrix representing
      the responses of 191 persons to 23 problems.  The responses are
      classified as correct (0) or incorrect (1).}
  }
  The 23 problems were:
  \itemize{
    \item{\code{p01} \eqn{%
      \big(\frac{1}{3} + \frac{1}{12}\big) : \frac{2}{9} = ?}{%
      (1/3 + 1/12) : 2/9 = ?}}
    \item{\code{p02} \eqn{%
      \big(\frac{3}{2} + \frac{3}{4}\big) \times \frac{5}{3} - 2 = ?}{%
      (3/2 + 3/4) x 5/3 - 2 = ?}}
    \item{\code{p03} \eqn{%
      \big(\frac{5}{6} + \frac{3}{14}\big) \times
           \big(\frac{19}{8} - \frac{3}{2}\big) = ?}{%
      (5/6 + 3/14) x (19/8 - 3/2) = ?}}
    \item{\code{p04} \eqn{%
      \big(\frac{1}{6} + \frac{2}{9}\big) - \frac{7}{36} = ?}{%
      (1/6 + 2/9) - 7/36 = ?}}
    \item{\code{p05} \eqn{%
      \frac{7}{10} + \frac{9}{10} = ?}{%
      7/10 + 9/10 = ?}}
    \item{\code{p06} \eqn{%
      \frac{8}{13} + \frac{5}{2} = ?}{%
      8/13 + 5/2 = ?}}
    \item{\code{p07} \eqn{%
      \frac{8}{12} + \frac{4}{15} = ?}{%
      8/12 + 4/15 = ?}}
    \item{\code{p08} \eqn{%
      \frac{2}{9} + \frac{5}{6} = ?}{%
      2/9 + 5/6 = ?}}
    \item{\code{p09} \eqn{%
      \frac{7}{5} + \frac{1}{5} = ?}{%
      7/5 + 1/5 = ?}}
    \item{\code{p10} \eqn{%
      \frac{2}{7} + \frac{3}{14} = ?}{%
      2/7 + 3/14 = ?}}
    \item{\code{p11} \eqn{%
      \frac{5}{9} + \frac{1}{6} = ?}{%
      5/9 + 1/6 = ?}}
    \item{\code{p12} \eqn{%
      \big(\frac{1}{12} + \frac{1}{3}\big) \times \frac{24}{15} = ?}{%
      (1/12 + 1/3) x 24/15 = ?}}
    \item{\code{p13} \eqn{%
      2 - \frac{3}{4} = ?}{%
      2 - 3/4 = ?}}
    \item{\code{p14} \eqn{%
      \big(4 + \frac{3}{4} - \frac{1}{2}\big) \times \frac{8}{6} = ?}{%
      (4 + 3/4 - 1/2) x 8/6 = ?}}
    \item{\code{p15} \eqn{%
      \frac{4}{7} + \frac{3}{4} = \frac{?}{28}}{%
      4/7 + 3/4 = ?/28}}
    \item{\code{p16} \eqn{%
      \frac{5}{8} - \frac{3}{16} = \frac{? - ?}{16}}{%
      5/8 - 3/16 = (? - ?)/16}}
    \item{\code{p17} \eqn{%
      \frac{3}{8} + \frac{5}{12} = \frac{? \times 3 + ? \times 5}{24}}{%
      3/8 + 5/12 = (? x 3 + ? x 5)/24}}
    \item{\code{p18} \eqn{%
      \frac{2}{7} + \frac{3}{5} = \frac{5 \times ? + 7 \times ?}{35}}{%
      2/7 + 3/5 = (5 x ? + 7 x ?)/35}}
    \item{\code{p19} \eqn{%
      \frac{2}{3} + \frac{6}{9} = \frac{?}{9} = \frac{?}{?}}{%
      2/3 + 6/9 = ?/9 = ?/?}}
    \item{\code{p20} Least common multiple \eqn{lcm(6, 8) = ?}}
    \item{\code{p21} \eqn{%
      \frac{7}{11} \times \frac{2}{3} = ?}{%
      7/11 x 2/3 = ?}}
    \item{\code{p22} \eqn{%
      \frac{2}{5} \times \frac{15}{4} = ?}{%
      2/5 x 15/4 = ?}}
    \item{\code{p23} \eqn{%
      \frac{9}{7} : \frac{2}{3} = ?}{%
      9/7 : 2/3 = ?}}
  }
  \code{subtraction13} is a data frame consisting of the following components:
  \describe{
    \item{\code{School}}{factor; school id.}
    \item{\code{Classroom}}{factor; class room id.}
    \item{\code{Gender}}{factor; participant gender.}
    \item{\code{Age}}{participant age.}
    \item{\code{R}}{a person-by-problem indicator matrix representing the
      responses of 294 persons to eight problems.}
  }
  The eight problems were:
  \itemize{
    \item{\code{p1} \eqn{73 - 58}}
    \item{\code{p2} \eqn{317 - 94}}
    \item{\code{p3} \eqn{784 - 693}}
    \item{\code{p4} \eqn{507 - 49}}
    \item{\code{p5} \eqn{253 - 178}}
    \item{\code{p6} \eqn{2245 - 418}}
    \item{\code{p7} \eqn{156 - 68}}
    \item{\code{p8} \eqn{3642 - 753}}
  }
}
\source{
  The data were made available by Debora de Chiusole, Andrea Brancaccio, and
  Luca Stefanutti.
}
\references{
  de Chiusole, D., & Stefanutti, L. (2013).
  Modeling skill dependence in probabilistic competence structures.
  \emph{Electronic Notes in Discrete Mathematics},
  \bold{42}, 41--48.
  \doi{https://doi.org/10.1016/j.endm.2013.05.144}

  Stefanutti, L., & de Chiusole, D. (2017).
  On the assessment of learning in competence based knowledge space theory.
  \emph{Journal of Mathematical Psychology},
  \bold{80}, 22--32.
  \doi{10.1016/j.jmp.2017.08.003}
}
\examples{
data(schoolarithm)

## Fraction problems used in Sefanutti and de Chiusole (2017)
R <- fraction17[, c(4:8, 10:11, 15:20)]
colnames(R) <- 1:13
N.R <- as.pattern(R, freq = TRUE)

## Conjunctive skill function in Table 1
sf <- read.table(header = TRUE, text = "
  item  a  b  c  d  e  f  g  h
     1  1  1  1  0  1  1  0  0
     2  1  0  0  0  0  0  1  1
     3  1  1  0  1  1  0  0  0
     4  1  1  0  0  1  1  1  1
     5  1  1  0  0  1  1  0  0
     6  1  1  1  0  1  0  1  1
     7  1  1  0  0  1  1  0  0
     8  1  1  0  0  1  0  1  1
     9  0  1  0  0  1  0  0  0
    10  0  1  0  0  0  0  0  0
    11  0  0  0  0  1  0  0  0
    12  1  1  0  0  1  0  1  1
    13  0  0  0  0  0  1  0  0
")
K <- delineate(sf)$K  # delineated knowledge structure
blim(K, N.R)

## Subtraction problems used in de Chiusole and Stefanutti (2013)
N.R <- as.pattern(subtraction13$R, freq = TRUE)

# Skill function in Table 1
# (f) mastering tens and hundreds; (g) mastering thousands; (h1) one borrow;
# (h2) two borrows; (h3) three borrows; (i) mastering the proximity of
# borrows; (j) mastering the presence of the zero; (k) mental calculation
sf <- read.table(header = TRUE, text = "
  item  f  g h1 h2 h3  i  j  k
     1  0  0  1  0  0  0  0  0
     2  1  0  1  0  0  0  0  0
     3  1  0  1  0  0  1  0  0
     4  1  0  1  1  1  0  1  0
     4  0  0  0  0  0  0  0  1
     5  1  0  1  1  1  1  0  0
     6  1  1  1  1  0  0  0  0
     7  1  0  1  1  1  1  0  0
     8  1  1  1  1  1  0  0  0
")
K <- delineate(sf)$K
blim(K, N.R)
}
\keyword{datasets}
