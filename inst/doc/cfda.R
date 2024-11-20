## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.width = 5, fig.height = 5)

## ----load, echo=TRUE, message=FALSE-------------------------------------------
library(cfda)
set.seed(42)

## ----genData, echo=TRUE-------------------------------------------------------
K <- 4
PJK <- matrix(1 / 3, ncol = K, nrow = K) - diag(rep(1 / 3, K))
lambda_PJK <- c(1, 1, 1, 1)

Tmax <- 10
n <- 100

d_JK <- generate_Markov(n = n, K = 4, P = PJK, lambda = lambda_PJK, Tmax = Tmax, labels = as.factor(c("A", "C", "G", "T")))

head(d_JK, 10)

## ----suumary, echo=TRUE-------------------------------------------------------
summary_cfd(d_JK)

## ----echo=TRUE----------------------------------------------------------------
duration <- compute_duration(d_JK)
head(duration)

## ----cutT, echo=TRUE----------------------------------------------------------
d_JKT <- cut_data(d_JK, Tmax = Tmax)

## ----plot, echo=TRUE, fig.height = 7------------------------------------------
plotData(d_JKT)

## ----echo=TRUE----------------------------------------------------------------
timeSpent <- compute_time_spent(d_JKT)
timeSpent[1:10, ]

## ----echo=TRUE----------------------------------------------------------------
boxplot(timeSpent)

## ----echo=TRUE----------------------------------------------------------------
nJump <- compute_number_jumps(d_JK)
head(nJump)

## ----echo=TRUE----------------------------------------------------------------
hist(nJump)

## ----echo=TRUE----------------------------------------------------------------
pt_evol <- estimate_pt(d_JKT)
pt_evol$pt[1:K, 1:10]
head(pt_evol$t)

## ----echo=TRUE----------------------------------------------------------------
plot(pt_evol, ribbon = TRUE)

## ----echo=TRUE----------------------------------------------------------------
statetable(d_JK)

## ----echo=TRUE----------------------------------------------------------------
mark <- estimate_Markov(d_JK)
mark

## ----echo=TRUE----------------------------------------------------------------
plot(mark)

## ----echo=TRUE----------------------------------------------------------------
m <- 20
basis <- create.bspline.basis(c(0, Tmax), nbasis = m, norder = 4)

## ----echo=TRUE----------------------------------------------------------------
fmca <- compute_optimal_encoding(d_JKT, basis, verbose = FALSE)
summary(fmca)

## ----echo=TRUE----------------------------------------------------------------
plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)

## ----echo=TRUE----------------------------------------------------------------
print(fmca$alpha[[1]])

## ----echo=TRUE----------------------------------------------------------------
plot(fmca)

## ----echo=TRUE----------------------------------------------------------------
encoding <- get_encoding(fmca, fdObject = TRUE)

## ----echo=TRUE----------------------------------------------------------------
plot(fmca, addCI = TRUE, coeff = 2, states = "A")

## ----echo=TRUE----------------------------------------------------------------
plotComponent(fmca, comp = c(1, 2), addNames = TRUE)

## ----echo=TRUE----------------------------------------------------------------
plotData(d_JKT[d_JKT$id %in% c(3, 14, 67, 84), ])

## ----echo=TRUE----------------------------------------------------------------
indicators <- reconstructIndicators(fmca)
head(indicators)

## ----echo=TRUE----------------------------------------------------------------
iInd <- 3
plotData(d_JKT[d_JKT$id == iInd, ])
plotIndicatorsReconstruction(indicators, id = iInd)

