% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{care}
\alias{care}
\title{Care trajectories}
\format{
A data.frame containing three columns:
\itemize{
  \item \emph{id} id of individuals (2929 different ids)
  \item \emph{time} number of months since the diagnosis
  \item \emph{state} new state.
}
}
\source{
https://larmarange.github.io/analyse-R/data/care_trajectories.RData
https://larmarange.github.io/analyse-R/trajectoires-de-soins.html
}
\usage{
data(care)
}
\description{
Care trajectories of patients diagnosed with a serious and chronic condition
}
\details{
In this study, patients were followed from the time they were diagnosed with a serious and chronic condition
and their care trajectories were tracked monthly from the time of diagnosis.
The status variable contains the care status of each individual for each month of follow-up.
Trajectories have different lengths.

The four states are:
\itemize{
 \item{D: diagnosed, but not in care}
 \item{C: in care, but not on treatment}
 \item{T: on treatment, but infection not suppressed}
 \item{S: on treatment and suppressed infection}
}
}
\examples{
data(care)
head(care)

plotData(care)

# Individuals has not the same length. In order to compute the encoding,
# we keep individuals with at least 18 months of history and work
# with the 18 first months.
duration <- compute_duration(care)
idToKeep <- as.numeric(names(duration[duration >= 18]))
care2 <- cut_data(care[care$id \%in\% idToKeep, ], 18)
head(care2)
\donttest{
# It is recommended to increase the number of cores to reduce computation time
set.seed(42)
basis <- create.bspline.basis(c(0, 18), nbasis = 10, norder = 4)
fmca <- compute_optimal_encoding(care2, basis, nCores = 2)

plotEigenvalues(fmca, cumulative = TRUE, normalize = TRUE)
plot(fmca)
plot(fmca, addCI = TRUE)
plotComponent(fmca, addNames = FALSE)
}
}
\seealso{
Other datasets: 
\code{\link{biofam2}},
\code{\link{flours}}
}
\concept{datasets}
\keyword{data}
