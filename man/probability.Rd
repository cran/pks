\name{probability}
\alias{probability}
\title{Problems in Elementary Probability Theory}
\description{
  This data set contains responses to problems in elementary probability
  theory observed before and after some instructions (the so-called learning
  object) were given.  Data were collected both in the lab and via an online
  questionnaire.  Of the 1127 participants eligible in the online study, 649
  were excluded because they did not complete the first set of problems
  (p101, ..., p112) or they responded too quickly or too slowly.  Based on
  similar criteria, further participants were excluded for the second set of
  problems, indicated by missing values in the variables b201, ..., b212.
  Problems were presented in random order.

  Participants were randomized to two conditions: an enhanced learning object
  including instructions with examples and a basic learning object without
  examples.  Instructions were given on four concepts: how to calculate the
  classic probability of an event (pb), the probability of the complement of
  an event (cp), of the union of two disjoint events (un), and of two
  independent events (id).

  The questionnaire was organized as follows:
  \describe{
    \item{Page 1}{Welcome page.}
    \item{Page 2}{Demographic data.}
    \item{Page 3}{First set of problems.}
    \item{Page 4 to 8}{Instructions (learning object).}
    \item{Page 9}{Second set of problems.}
    \item{Page 10}{Feedback about number of correctly solved problems.}
  }
}
\usage{data(probability)}
\format{
  A data frame with 504 cases and 68 variables:
  \itemize{
    \item{\code{case} a factor giving the case id, a five-digits code
      the fist digit denoting lab or online branch of the study, the last
      four digits being the case number.}
    \item{\code{lastpage} Which page of the questionnaire was reached
      before quitting?  The questionnaire consisted of ten pages.}
    \item{\code{mode} a factor; \code{lab} or \code{online} branch of
      study.}
    \item{\code{started} a timestamp of class POSIXlt.  When did participant
      start working on the questionnaire?}
    \item{\code{sex} a factor coding sex of participant.}
    \item{\code{age} age of participant.}
    \item{\code{educat} education as a factor with three levels:
      \code{1} secondary school or below; \code{2} higher education entrance
      qualification; \code{3} university degree.}
    \item{\code{fos} field of study.  Factor with eight levels:
      \code{ecla} economics, business, law; \code{else} miscellaneous;
      \code{hipo} history, politics; \code{lang} languages; \code{mabi}
      mathematics, physics, biology; \code{medi} medical science;
      \code{phth} philosophy, theology; \code{psco} psychology,
      computer science, cognitive science.}
    \item{\code{semester} ordered factor.  What semester are you in?}
    \item{\code{learnobj} a factor with two levels:
      \code{enhan} learning object enhanced with examples;
      \code{basic} learning object without examples.}
  }
  The twelve problems of the first part (before the learning object):
  \itemize{
    \item{\code{p101} A box contains 30 marbles in the following colors:
      8 red, 10 black, 12 yellow.  What is the probability that a randomly
      drawn marble is yellow?  (Correct: 0.40)}
    \item{\code{p102} A bag contains 5-cent, 10-cent, and 20-cent coins.
      The probability of drawing a 5-cent coin is 0.35, that of drawing a
      10-cent coin is 0.25, and that of drawing a 20-cent coin is 0.40.  What
      is the probability that the coin randomly drawn is not a 5-cent
      coin?  (0.65)}
    \item{\code{p103} A bag contains 5-cent, 10-cent, and 20-cent coins.
      The probability of drawing a 5-cent coin is 0.20, that of drawing a
      10-cent coin is 0.45, and that of drawing a 20-cent coin is 0.35.  What
      is the probability that the coin randomly drawn is a 5-cent coin or
      a 20-cent coin?  (0.55)}
    \item{\code{p104} In a school, 40\% of the pupils are boys and 80\% of
      the pupils are right-handed.  Suppose that gender and handedness are
      independent.  What is the probability of randomly selecting a
      right-handed boy?  (0.32)}
    \item{\code{p105} Given a standard deck containing 32 different cards,
      what is the probability of not drawing a heart?  (0.75)}
    \item{\code{p106} A box contains 20 marbles in the following colors:
      4 white, 14 green, 2 red.  What is the probability that a randomly
      drawn marble is not white?  (0.80)}
    \item{\code{p107} A box contains 10 marbles in the following colors:
      2 yellow, 5 blue, 3 red.  What is the probability that a randomly
      drawn marble is yellow or blue?  (0.70)}
    \item{\code{p108} What is the probability of obtaining an even number
      by throwing a dice?  (0.50)}
    \item{\code{p109} Given a standard deck containing 32 different cards,
      what is the probability of drawing a 4 in a black suit?  (Responses
      that round to 0.06 were considered correct.)}
    \item{\code{p110} A box contains marbles that are red or yellow, small
      or large.  The probability of drawing a red marble is 0.70 (lab: 0.30),
      the probability of drawing a small marble is 0.40.  Suppose that the
      color of the marbles is independent of their size.  What is the
      probability of randomly drawing a small marble that is not red?  (0.12,
      lab: 0.28)}
    \item{\code{p111} In a garage there are 50 cars.  20 are black and 10 are
      diesel powered.  Suppose that the color of the cars is independent of
      the kind of fuel.  What is the probability that a randomly selected car
      is not black and it is diesel powered?  (0.12)}
    \item{\code{p112} A box contains 20 marbles.  10 marbles are red, 6 are
      yellow and 4 are black.  12 marbles are small and 8 are large.  Suppose
      that the color of the marbles is independent of their size.  What is the
      probability of randomly drawing a small marble that is yellow or
      red?  (0.48)}
  }
  The twelve problems of the second part (after the learning object):
  \itemize{
    \item{\code{p201} A box contains 30 marbles in the following colors:
      10 red, 14 yellow, 6 green.  What is the probability that a randomly
      drawn marble is green?  (0.20)}
    \item{\code{p202} A bag contains 5-cent, 10-cent, and 20-cent coins.
      The probability of drawing a 5-cent coin is 0.25, that of drawing a
      10-cent coin is 0.60, and that of drawing a 20-cent coin is 0.15.  What
      is the probability that the coin randomly drawn is not a 5-cent
      coin?  (0.75)}
    \item{\code{p203} A bag contains 5-cent, 10-cent, and 20-cent coins.
      The probability of drawing a 5-cent coin is 0.35, that of drawing a
      10-cent coin is 0.20, and that of drawing a 20-cent coin is 0.45.  What
      is the probability that the coin randomly drawn is a 5-cent coin or
      a 20-cent coin?  (0.80)}
    \item{\code{p204} In a school, 70\% of the pupils are girls and 10\% of
      the pupils are left-handed.  Suppose that gender and handedness are
      independent. What is the probability of randomly selecting a
      left-handed girl?  (0.07)}
    \item{\code{p205} Given a standard deck containing 32 different cards,
      what is the probability of not drawing a club?  (0.75)}
    \item{\code{p206} A box contains 20 marbles in the following colors:
      6 yellow, 10 red, 4 green.  What is the probability that a randomly
      drawn marble is not yellow?  (0.70)}
    \item{\code{p207} A box contains 10 marbles in the following colors:
      5 blue, 3 red, 2 green.  What is the probability that a randomly
      drawn marble is blue or red?  (0.80)}
    \item{\code{p208} What is the probability of obtaining an odd number
      by throwing a dice?  (0.50)}
    \item{\code{p209} Given a standard deck containing 32 different cards,
      what is the probability of drawing a 10 in a red suit?  (Responses
      that round to 0.06 were considered correct.)}
    \item{\code{p210} A box contains marbles that are green or red, large
      or small The probability of drawing a green marble is 0.40,
      the probability of drawing a large marble is 0.20.  Suppose that the
      color of the marbles is independent of their size.  What is the
      probability of randomly drawing a large marble that is not
      green?  (0.12)}
    \item{\code{p211} In a garage there are 50 cars.  15 are white and 20 are
      diesel powered. Suppose that the color of the cars is independent of
      the kind of fuel.  What is the probability that a randomly selected car
      is not white and it is diesel powered?  (0.28)}
    \item{\code{p212} A box contains 20 marbles.  8 marbles are white, 4 are
      green and 8 are red.  15 marbles are small and 5 are large.  Suppose
      that the color of the marbles is independent of their size.  What is the
      probability of randomly drawing a large marble that is white or
      green?  (0.15)}
  }
  Further variables:
  \itemize{
    \item{\code{time01}, ..., \code{time10} the time (in s) spent on each
      page of the questionnaire.  In the lab branch of the study, participants
      started directly on Page 2.}
    \item{\code{b101}, ..., \code{b112} the twelve problems of the first
      part coded as correct (1) or error (0).}
    \item{\code{b201}, ..., \code{b212} the twelve problems of the second
      part coded as correct (1) or error (0).}
  }
}
\source{
  Data were collected by Pasquale Anselmi and Florian Wickelmaier at the
  Department of Psychology, University of Tuebingen, in February and March
  2010.
}
\examples{
data(probability)

## "Completer" sample
pb <- probability[!is.na(probability$b201), ]

## Response frequencies for first and second part
N.R1 <- as.pattern(pb[, sprintf("b1\%.2i", 1:12)], freq = TRUE)
N.R2 <- as.pattern(pb[, sprintf("b2\%.2i", 1:12)], freq = TRUE)

## Conjunctive skill function, one-to-one problem function
sf1 <- read.table(header = TRUE, text = "
  item cp id pb un
     1  0  0  1  0
     2  1  0  0  0
     3  0  0  0  1
     4  0  1  0  0
     5  1  0  1  0
     6  1  0  1  0
     7  0  0  1  1
     8  0  0  1  1
     9  0  1  1  0
    10  1  1  0  0
    11  1  1  1  0
    12  0  1  1  1
")

## Extended skill function
sf2 <- rbind(sf1, read.table(header = TRUE, text = "
  item cp id pb un
     2  0  0  0  1
     3  1  0  0  0
     6  0  0  1  1
     7  1  0  1  0
    12  1  1  1  0
"))

## Delineated knowledge structures
K1 <- delineate(sf1)$K
K2 <- delineate(sf2)$K

## After instructions, fit of knowledge structures improves
sapply(list(N.R1, N.R2), function(n) blim(K1, n)$discrepancy)
sapply(list(N.R1, N.R2), function(n) blim(K2, n)$discrepancy)
}
\keyword{datasets}
