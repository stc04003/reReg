% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcf.R
\name{mcf}
\alias{mcf}
\title{The \code{mcf} function is imported from \code{reda}.}
\description{
An S4 class generic function that returns the mean cumulative function (mcf) estimates.
The function is imported from the \code{reda} package; see \code{?reda::mcf} for more details.
}
\examples{
data(simDat)
m <- mcf(Recur(t.start \%to\% t.stop, id, event, status) ~ x1, data = simDat)
plot(m, conf.int = TRUE)
}
