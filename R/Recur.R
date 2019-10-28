#' The \code{Recur} function is imported from \code{reda}.
#'
#' See \code{?reda::Recur} for more details.
#'
#' @rdname Recur
#' @name Recur
#' @importFrom reda Recur
#' @aliases Recur
#' @export Recur
NULL

#' The \code{is.Recur} function is imported from \code{reda}.
#'
#' See \code{?reda::is.Recur} for more details.
#'
#' \itemize{
#' \item \code{\link[reda]{is.Recur}}: imported from the \code{reda}
#' package.
#' }
#'
#' @rdname Recur
#' @name Recur
#' @importFrom reda is.Recur
#' @aliases is.Recur
#' @export is.Recur
NULL

##' Recurrent Episodes
##'
##' Specify time segments or recurrent episodes by endpoints.
##'
##' This function is intended to be used for specifying the argument \code{time}
##' in function \code{\link{Recur}}.
##'
##' @name %to%
##' @rdname Recur-to
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
##' @name %2%
##' @export
`%2%` <- `%to%`
