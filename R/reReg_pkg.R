#' reReg: Recurrent Event Regression
#'
#' The package provides easy access to fit regression models to recurrent event data.
#' The available implementations allow users to explore recurrent data through event plot and cumulative sample mean function plot,
#' fit semiparametric regression models under different assumptions,
#' and simulate recurrent event data. 
#'
#' @aliases reReg-packages
#' @references Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.C. and Yan, J. (2017).
#' Joint Scale-change Models for Recurrent Events and Failure Time.
#' \emph{Journal of the American Statistical Association}, \bold{112}(518): 796--805.
#' @references Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric Regression for the Mean and Rate Functions of Recurrent Events.
#' \emph{Journal of the Royal Statistical Society: Series B (Methodological)}, \bold{62}: 711--730.
#' @references Wang, M.C., Qin, J., and Chiang, C.T. (2001). Analyzing Recurrent Event Data with Informative Censoring.
#' \emph{Journal of the American Statistical Association}, \bold{96}(455): 1057--1065.
#' @references Ghosh, D. and Lin, D.Y. (2003). Semiparametric Analysis of Recurrent Events Data in the Presence of Dependent Censoring.
#' \emph{Biometrics}, \bold{59}: 877--885.
#' @references Huang, C.-Y. and Wang, M.C. (2004). Joint Modeling and Estimation for Recurrent Event Processes and Failure Time Data.
#' \emph{Journal of the American Statistical Association}, \bold{99}(468): 1153--1165.
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
#' @importFrom tidyr unnest nest
#' @importFrom purrr map_lgl
#' @importFrom dplyr "%>%" group_by group_by_at filter summarize mutate full_join left_join arrange select tibble count bind_cols
#' @importFrom plyr ddply 
#' @importFrom grDevices hcl
#' @importFrom parallel makeCluster clusterExport stopCluster detectCores parSapply
#' 
#' @docType package
#' @useDynLib reReg, .registration = TRUE
"_PACKAGE"
NULL
