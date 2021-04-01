#' The \code{Recur} function is imported from \code{reda}.
#'
#' Create a recurrent event survival object, used as a response variable in \code{reReg}.
#' This function is replacing the original \code{reSurv()} in version 1.1.6.
#' See \code{?reda::Recur} for more details.
#'
#' @rdname Recur
#' @name Recur
#' @aliases is.Recur
#' @importFrom reda Recur is.Recur check_Recur
#' @seealso \code{\link{\%2\%}}
#' @importMethodsFrom reda show summary
#' @exportMethod show summary
#' @export Recur is.Recur
#' @example inst/examples/ex_Recur.R
NULL

#' The \code{Recur} class is imported from \code{reda}.
#'
#' The class \code{Recur} is an S4 that represents a formula response for
#' recurrent event data model. See \code{reda} for details.
#'
#' @rdname Recur-class
#' @name Recur-class
#' @importClassesFrom reda Recur 
#' @exportClass Recur
NULL

#' The \code{summary.Recur} class is imported from \code{reda}.
#'
#' The class \code{summary.Recur} is an S4 that represents the summary of a \code{Recur} object.
#' See \code{reda} for details.
#'
#' @rdname summary.Recur-class
#' @name summary.Recur-class
#' @importClassesFrom reda summary.Recur 
#' @exportClass summary.Recur 
NULL

#' The \code{\%to\%} function is imported from \code{reda}
#'
#' This pipe operator specifies the time segments or recurrent episodes by endpoints.
#' See \code{reda} for more details.
#'
#' @name Recur-pipe
#' @rdname Recur-pipe
#' @aliases %to% %2%
#' @importFrom reda %to% %2%
#' @export %to% %2%
#' @example inst/examples/ex_Recur.R
NULL
