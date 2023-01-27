#' reReg: Recurrent Event Regression
#'
#' The package offers a comprehensive collection of practical and easy-to-use tools for analyzing
#' recurrent event data, with or without the presence of a (possibly) correlated terminal event.
#' The modeling framework is based on a joint frailty scale-change model,
#' that encompasses many existing models, including the popular Cox-type models,
#' as special cases and accommodates informative censoring through a subject-specific frailty.
#' The implemented estimating procedure does not require any parametric assumption on the frailty
#' distribution.
#' The package allows the users to specify different model forms for both the recurrent event process
#' and the terminal event.
#' The package also includes tools for visualization of recurrent events and
#' simulation from the regression models.
#'
#' @aliases reReg-packages
#' @references Chiou, S.H., Xu, G., Yan, J. and Huang, C.-Y. (2023). Regression Modeling for Recurrent Events Possibly with an Informative Terminal Event Using {R} Package {reReg}.
#' \emph{Journal of Statistical Software}, \bold{105}(5): 1--34. 
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711--730.
#' @references Wang, M.-C., Qin, J., and Chiang, C.-T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association}, \bold{96}(455): 1057--1065.
#' @references Ghosh, D. and Lin, D.Y. (2002). Marginal Regression Models for Recurrent and Terminal Events. \emph{Statistica Sinica}: 663--688.
#' @references Ghosh, D. and Lin, D.Y. (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877--885.
#' @references Huang, C.-Y. and Wang, M.-C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association}, \bold{99}(468): 1153--1165.
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017).
#' Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association}, \bold{112}(518): 796--805.
#' @references Xu, G., Chiou, S.H., Yan, J., Marr, K., and Huang, C.-Y. (2019). Generalized Scale-Change Models for Recurrent Event
#' Processes under Informative Censoring. \emph{Statistica Sinica}, \bold{30}: 1773--1795.
#' 
#' @importFrom graphics lines par plot title boxplot
#' @importFrom stats approx ksmooth lm model.extract model.matrix pnorm printCoefmat quantile rbinom rexp rgamma rgeom rlnorm rnorm rpois
#' @importFrom stats runif rweibull var aggregate as.formula coef sd terms pchisq complete.cases stepfun
#' @importFrom methods getClass
#' @importFrom utils head
#' @importFrom MASS ginv
#' @importFrom nleqslv nleqslv
#' @importFrom BB BBsolve dfsane BBoptim
#' @importFrom SQUAREM squarem
#' @importFrom survival Surv basehaz coxph survfit
#' @importFrom ggplot2 ggplot geom_point aes geom_line geom_bar ggtitle scale_color_manual scale_shape_manual theme element_text labs facet_grid theme_bw coord_flip 
#' @importFrom ggplot2 element_blank element_line element_rect xlab ylab scale_x_continuous alpha guides guide_legend geom_step scale_color_discrete
#' @importFrom grDevices hcl
#' @importFrom parallel makeCluster clusterExport stopCluster detectCores parSapply
#' @importFrom Rcpp sourceCpp
#' 
#' @docType package
#' @useDynLib reReg
"_PACKAGE"
NULL


## #' @useDynLib reReg, .registration = TRUE
