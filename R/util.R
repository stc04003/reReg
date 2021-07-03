#' Check for whole number
#' @noRd
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol


#' Calculate b %*% solve(A) %*% b
#' @noRd
bAib <- function(A, b) {
    tryCatch(crossprod(solve(t(chol(A)), b)),
             error = function(e) b %*% svd2ginv(A) %*% b)
}

#' Calculate solve(A) %*% B %*% t(solve(A))
#' @noRd
AiBAi <- function(A, B) {
    tryCatch(tcrossprod(solve(A, t(chol(B)))),
             error = function(e) svd2ginv(A) %*% B %*% t(svd2ginv(A)))
}

#' Calculate solve(t(A) %*% A) %*% t(A) %*% b
#' @noRd
Axb <- function(A, b) {
    tryCatch(solve(crossprod(A), crossprod(A, b)),
             error = function(e) svd2ginv(crossprod(A)) %*% crossprod(A, b))
}

#' Generalized printCoefmat for printing in summary.reReg
#' @noRd
printCoefmat2 <- function(tab) 
    printCoefmat(as.data.frame(tab), P.values = TRUE,
                 has.Pvalue = TRUE, signif.legend = FALSE)

#' Computes (generalized) inverse of a matrix from svd
#' This guard against error in applying solve(A) on a matrix whose inverse doesn't exist
#' @noRd
svd2ginv <- function(A) {
    tmp <- svd(A)
    p <- tmp$d > 1e-10
    tmp$id <- diag(1 / tmp$d[p], sum(p))
    tmp$tu <- t(tmp$u[,p])
    tmp$vv <- tmp$v[,p]
    Reduce("%*%", tmp[c("vv", "id", "tu")])
}
