#' Check for whole number
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#' Calculate b %*% solve(A) %*% b
#' @noRd
bAib <- function(A, b) {
    crossprod(solve(t(chol(A)), b))
}

#' Calculate solve(A) %*% B %*% t(solve(A))
#' @noRd
AiBAi <- function(A, B) {
    tcrossprod(solve(A, t(chol(B))))
}

#' Calculate solve(t(A) %*% A) %*% t(A) %*% b
#' @noRd
Axb <- function(A, b) {
    solve(crossprod(A), crossprod(A, b))
}

#' Generalized printCoefmat for printing in summary.reReg
#' @noRd
printCoefmat2 <- function(tab) 
    printCoefmat(as.data.frame(tab), P.values = TRUE,
                 has.Pvalue = TRUE, signif.legend = FALSE)
