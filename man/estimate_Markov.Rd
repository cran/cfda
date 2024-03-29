% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{estimate_Markov}
\alias{estimate_Markov}
\title{Estimate transition matrix and spent time}
\usage{
estimate_Markov(data)
}
\arguments{
\item{data}{data.frame containing \code{id}, id of the trajectory, \code{time}, time at which a change occurs and
\code{state}, associated state.}
}
\value{
list of two elements: \code{Q}, the estimated transition matrix, and \code{lambda},
the estimated time spent in each state
}
\description{
Calculates crude initial values for transition intensities by assuming that the data represent the exact transition
times of the Markov process.
}
\examples{
# Simulate the Jukes-Cantor model of nucleotide replacement
K <- 4
PJK <- matrix(1 / 3, nrow = K, ncol = K) - diag(rep(1 / 3, K))
lambda_PJK <- c(1, 1, 1, 1)
d_JK <- generate_Markov(n = 100, K = K, P = PJK, lambda = lambda_PJK, Tmax = 10)

# estimation
mark <- estimate_Markov(d_JK)
mark$P
mark$lambda
}
\seealso{
\link{plot.Markov}
}
\author{
Cristian Preda
}
