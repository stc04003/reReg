#' @name reSurv
#' @rdname reSurv
#' @title Create an \code{reSurv} Object
#'
#' @description Create a recurrent event survival object, used as a response variable in \code{reReg}.
#' This function is being deprecated in Version 1.1.6.
#' A recurrent event object is now being created with \code{Recur()}.
#' See '?Recur()' for details.
#'
#' @param time1 when "\code{time2}" is provided, this vector is treated as the starting time for the gap time between two successive recurrent events.
#' In the absence of "\code{time2}", this is the observation time of recurrence on calendar time scale, in which, the time corresponds to the time since entry/inclusion in the study.
#' @param time2 an optional vector for ending time for the gap time between two successive recurrent events.
#' @param event a binary vector used as the recurrent event indicator. \code{event = 1} for recurrent times.
#' @param status a binary vector used as the status indicator for the terminal event. \code{status = 0} for censored times.
#' @param id subject's id.
#' @param origin a numerical vector indicating the time origin of subjects.
#' When \code{origin} is a scalar, \code{reSurv} assumes all subjects have the same origin.
#' Otherwise, \code{origin} needs to be a numerical vector, with length equals to the number of subjects.
#' In this case, each element corresponds to different origins for different subjects.
#' This argument is only needed when "\code{time2}" is missing.
#' 
#' @rdname reSurv
#' @export
#' @example inst/examples/ex_reSurv.R
reSurv <- function(time1, time2, id, event, status, origin = 0) {
    warning("'reSurv()' is being deprecated in Version 1.1.7. Output is prepared by 'Recur()'.\n See '?Recur()' for details.\n")
    nArg <- length(match.call()) - 1
    if (nArg >= 5) return(Recur(time1 %2% time2, id, event, status))
    if (nArg < 5) return(Recur(time1, time2, id, event, origin = origin))
}

is.reReg <- function(x) inherits(x, "reReg")
