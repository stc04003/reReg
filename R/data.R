#' Simulated dataset for demonstration
#'
#' @description A simulated data frame with the following variables
#' \describe{
#'   \item{id}{subjects identification}
#'   \item{t.start}{start of the interval}
#'   \item{t.stop}{endpoint of the interval; when time origin is 0 this variable also marks the recurrence or terminal/censoring time}
#'   \item{status}{terminal event indicator; 1 if a terminal event is recorded}
#'   \item{event}{recurrent event indicator; 1 if a recurrent event is recorded}
#'   \item{x1}{baseline covariate generated from a standard uniform distribution}
#'   \item{x2}{baseline covariate generated from a standard uniform distribution (independent from \code{z1}}
#' }
#'
#' @details
#' See \code{\link{simGSC}} for instruction on simulating recurrent event data from
#' scale-change models.
#'
#' @usage data(simDat)
#' @docType data
#' @name simDat
#' @rdname simDat
#' @format A data frame with 874 rows and 7 variables.
NULL

