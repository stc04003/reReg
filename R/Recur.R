#' The \code{Recur} function is imported from \code{reda}.
#'
#' This is a placeholder.
#' 
#' @importFrom reda Recur
#' @noRd
#' @aliases Recur
#' @export Recur
NULL

#' The \code{is.Recur} function is imported from \code{reda}.
#'
#' This is a placeholder.
#' 
#' @importFrom reda is.Recur
#' @noRd
#' @aliases is.Recur
#' @export is.Recur
NULL

##' Recurrent Episodes
##'
##' Specify time segements or recurrent episodes by endpoints.
##'
##' This function is intended to be used for specifying the argument \code{time}
##' in function \code{\link{Recur}}.
##'
##' @name Recur-to
##'
##' @param time1 The left end-points of the recurrent episodes.
##' @param time2 The right end-points of the recurrent episodes.
##'
##' @return A list that consists of two elements named
##'     \code{"time1"} and \code{"time2"}.
##'
##' @export
`%to%` <- function(time1, time2) {
    list(time1 = time1, time2 = time2)
}


##' @rdname Recur-to
##' @export
`%2%` <- `%to%`
