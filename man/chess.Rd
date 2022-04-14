\name{chess}
\alias{chess}
\title{Responses to Chess Problems and Knowledge Structures}
\description{
  Held, Schrepp and Fries (1995) derive several knowledge structures for the
  representation of 92 responses to 16 chess problems.  See Schrepp, Held and
  Albert (1999) for a detailed description of these problems.
}
\usage{data(chess)}
\format{
  A list consisting of five components:
  \describe{
    \item{\code{dst1}}{a state-by-problem indicator matrix representing the
      knowledge structure DST1.}
    \item{\code{dst3}}{the knowledge structure DST3.}
    \item{\code{dst4}}{the knowledge structure DST4.}
    \item{\code{N.R}}{a named integer vector.  The names denote response
      patterns, the values denote their frequencies.}
    \item{\code{R}}{a state-by-problem indicator matrix representing the raw
      responses.  Column names \code{hdbgXX} and \code{grazYY} identify
      responses collected in Heidelberg and Graz, respectively.}
  }
}
\note{
  The graphs of the precedence relations for DST1 and DST4 in Held et. al
  (1995) contain mistakes that have been corrected.  See examples.
}
\source{
  Held, T., Schrepp, M., & Fries, S. (1995).
  Methoden zur Bestimmung von Wissensstrukturen -- eine Vergleichsstudie.
  \emph{Zeitschrift fuer Experimentelle Psychologie},
  \bold{42}(2), 205--236.
}
\references{
  Schrepp, M., Held, T., & Albert, D. (1999).
  Component-based construction of surmise relations for chess problems.
  In D. Albert & J. Lukas (Eds.),
  \emph{Knowledge spaces: Theories, empirical research, and applications}
  (pp. 41--66).  Mahwah, NJ: Erlbaum.
}
\examples{
data(chess)
chess$dst1  # knowledge structure DST1

## Precedence relation (Held et al., 1995, p. 215) and knowledge space
P <- as.binmat(c("1111011101111001",   # s
               # "0100000000000000",   # gs   mistake in Abb. 3
                 "0111010100111000",   # gs   correction
                 "0011010000011000",   # egs
                 "0011010000011000",   # eegs
                 "0000110000000000",   # cs
                 "0000010000000000",   # gcs
                 "0011011100111000",   # ts
                 "0011010100011000",   # ges
                 "1111111111111111",   # f
                 "0111010101111000",   # gf
                 "0011010000111000",   # gff
                 "0000000000010000",   # ggff
                 "0000000000001000",   # ggf
                 "0111011101111101",   # ff
                 "0111011101111011",   # tf
                 "0011010100111001"),  # tff
               as.logical = TRUE)
dimnames(P) <- list("<" = colnames(chess$R), ">" = colnames(chess$R))
K <- rbind(0L, binary_closure(t(P)))
identical(sort(as.pattern(K)),
          sort(as.pattern(chess$dst1)))

blim(chess$dst1, chess$N.R)  # Tab. 1
}
\keyword{datasets}
