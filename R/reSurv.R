#' @name reSurv
#' @rdname reSurv
#' @title Create an \code{reSurv} Object
#'
#' @description Create a recurrent event survival object, used as a response variable in reReg model formula.
#'
#' @param time1 when "time2" is provided, this vector is treated as the starting time for the gap time between two successive recurrent events.
#' In the absence of "time2", this is the observation time of recurrence on calendar time scale, in which, the time corresponds to the time since entry/inclusion in the study.
#' @param time2 an optional vector for ending time for the gap time between two successive recurrent events.
#' @param event a binary vector used as the recurrent event indicator.
#' @param status a binary vector used as the status indicator for the terminal event.
#' @param id observation subject's id
#' @param x an \code{reSurv} object.
NULL

#' @rdname reSurv
#' @export
#' @examples
#' data(readmission)
#' with(readmission, reSurv(t.stop, event, death, id))
#' with(readmission, reSurv(t.start, t.stop, event, death, id))
reSurv <- function(time1, time2, event, status, id) {
    if (missing(time1)) stop("Must have a time argument.")
    if (!is.numeric(time1)) stop("Time argument (time1) must be numeric.")
    if (any(time1 < 0)) stop("Time argument (time1) must be positive.")
    if (any(missing(time2), missing(event), missing(status), missing(id))) {
        if (missing(event)) stop("Recurrent event indicator (event) is missing.")
        if (missing(status)) stop("Censoring indicator (status) is missing.")
        if (missing(id)) {
            id <- status
            status <- event
            event <- time2
            time2 <- NULL
        }
    } 
    if (!is.numeric(time2) & !is.null(time2)) stop("Time argument (time2) must be numeric.")
    if (any(time2 < 0) & !is.null(time2)) stop("Time argument (time2) must be positive.")
    if (!is.null(time2) & any(time1 > time2)) stop("Stop time (time2) must be > start time (time1).")
    if (length(event) != length(time1)) stop("Time argument and recurrent event indicator (event) are different lengths.")
    if (length(status) != length(time1)) stop("Time argument and status are different lengths.")
    ## Multiple event types?
    event2 <- NULL
    if (is.logical(event)) event <- as.numeric(event)
    else if (is.numeric(event)) {
        event2 <- event
        event <- ifelse(event == 0, 0, 1)
    }
    else stop("Event must be logical or numeric")
    ## Checking status
    if (is.logical(status)) status <- as.numeric(status)
    else if (is.numeric(status)) {
        temp <- (status == 0 | status == 1)
        status <- ifelse(temp, status, NA)
        if (any(is.na(status))) stop("Status must be 0 or 1)")
    }
    else stop("Status must be logical or numeric")
    ## Prepare DF
    if (is.null(time2))
        tab <- tibble(id = id, Time = time1, event = event, status = status, recType = event2)
    else tab <- tibble(id = id, Time = unlist(lapply(split(time2 - time1, id), cumsum)),
                           event = event, status = status, recType = event2)
    rownames(tab) <- NULL
    ## construct tibby
    tmp <- tab %>% mutate(status = status * ifelse(!event, 1, NA),
                          recType = ifelse(recType, recType, NA),
                          tij = Time * ifelse(event, 1, NA),
                          Yi = Time * ifelse(!event, 1, NA)) 
    x <- tmp %>% group_by(id) %>% filter(!is.na(tij)) %>% summarize(tij = list(tij), recType = list(recType))
    y <- tmp %>% group_by(id) %>% filter(!is.na(Yi)) %>% summarize(Yi = Yi, status = status)
    tab2 <- x %>% full_join(y, by = "id") %>% arrange(id)
    rc <- list(reDF = tab, reTb = tab2)
    class(rc) <- "reSurv"
    rc
}

#' @rdname reSurv
#' @export
is.reSurv <- function(x) inherits(x, "reSurv")

is.reReg <- function(x) inherits(x, "reReg")
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
