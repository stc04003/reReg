#' Rehospitalization Colorectal Cancer
#'
#' This contains rehospitalization times after surgery in patients diagnosed with colorectal cancer.
#'
#' @format A \code{data.frame} contains the following columns:
#' \describe{
#'   \item{id}{identificator of each subject (repeated for each recurrence).}
#'   \item{enum}{observation number within patient.}
#'   \item{t.start}{start of interval (0 or previous recurrence time).}
#'   \item{t.stop}{recurrence or censoring time.}
#'   \item{time}{interocurrence or censoring time.}
#'   \item{event}{rehospitalization status. All event are 1 for each subject excepting last one that it is 0.}
#'   \item{chemo}{Treatment (chemotherapy) indicator.}
#'   \item{sex}{Gender: Male or Female.}
#'   \item{dukes}{Tumour stage under Dukes's classification: A-B, C, or D.}
#'   \item{charlson}{Comorbidity Charlson's index, modelled as a time dependent covariate and classified into three groups: Index 0, Index 1-2, and Index >= 3.}
#'   \item{death}{death indicator: Dead = 1; alive = 0.}
#' }
#'
#' @docType data
#' @usage data(readmission)
#' @name readmission
#'
#' @source Gonzalez, J.R., Fernandez, E., Moreno, V., Ribes, J., Peris, M., Navarro, M., Cambray, M. and Borras, JM (2005). Sex differences in hospital readmission among colorectal cancer patients.
#' \emph{Journal of Epidemiology and Community Health}, \bold{59}, 6, 506-511.
NULL
