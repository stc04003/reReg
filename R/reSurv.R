#' @name reSurv
#' @rdname reSurv
#' @title Create an \code{reSurv} Object
#'
#' @description Create a recurrent event survival object, used as a response variable in \code{reReg}.
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
#' @param x an \code{reSurv} object.
NULL

#' @rdname reSurv
#' @export
#' @examples
#' data(readmission, package = "frailtypack")
#' attach(readmission)
#' reSurv(t.stop, id, event, death)
#' reSurv(t.start, t.stop, id, event, death)
#' detach(readmission)
reSurv <- function(time1, time2, id, event, status, origin = 0) {
    if (missing(time1) & missing(time2)) stop("Must have a time argument.")
    if (!is.numeric(time1)) stop("Time argument (time1) must be numeric.")
    if (any(time1 < 0)) stop("Time argument (time1) must be positive.")
    msArg <- sum(missing(time2), missing(event), missing(status), missing(id))
    n <- length(time1)
    if (msArg == 4) {
        eval(parse(text = default.reSurv(c("time2", "id", "event", "status"))))
    }
    if (msArg == 3) {
        if (missing(time2) & !missing(id)) 
            eval(parse(text = default.reSurv(c("time2", "event", "status"))))
        if (missing(time2) & !missing(event)) 
            eval(parse(text = default.reSurv(c("time2", "id", "status"))))
        if (missing(time2) & !missing(status)) 
            eval(parse(text = default.reSurv(c("time2", "id", "status"))))
        if (!missing(time2) & all(time2 >= time1)) 
            eval(parse(text = default.reSurv(c("event", "id", "status"))))
        if (!missing(time2) & !all(time2 >= time1)) {
            id <- time2
            eval(parse(text = default.reSurv(c("time2", "event", "status"))))
        }
    }
    if (msArg == 2) {
        if (missing(time2) & missing(id)) 
            eval(parse(text = default.reSurv(c("time2", "id"))))
        if (missing(time2) & missing(event)) 
            eval(parse(text = default.reSurv(c("time2", "event"))))
        if (missing(time2) & missing(status)) 
            eval(parse(text = default.reSurv(c("time2", "status"))))
        if (!missing(time2) & all(time2 >= time1)) {
            if (!missing(id))
                eval(parse(text = default.reSurv(c("event", "status"))))
            else if (!missing(event))
                eval(parse(text = default.reSurv(c("event", "id"))))
            else if (!missing(status))
                eval(parse(text = default.reSurv(c("id", "status"))))
        }
        if (!missing(time2) & !all(time2 >= time1)) {
            event <- id
            id <- time2
            eval(parse(text = default.reSurv(c("time2", "status"))))
        }
    }
    if (msArg == 1) {
        if (missing(time2)) time2 <- NULL
        if (missing(status) & !is.null(time2) & all(time2 >= time1)) {
            if (missing(id))
                eval(parse(text = default.reSurv(c("id"))))
            else if (missing(event))
                eval(parse(text = default.reSurv(c("event"))))
            else if (missing(status))
                eval(parse(text = default.reSurv(c("status"))))
        }
        if (missing(status) & !is.null(time2) & !all(time2 >= time1)) {
            status <- event
            event <- id
            id <- time2
            time2 <- NULL
        }
        if (missing(event)) stop("Recurrent event indicator (event) is missing.")
        if (missing(status)) stop("Censoring indicator (status) is missing.")
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
        tab <- tibble(id = id, Time = time1 - origin, event = event, status = status, recType = event2)
    else tab <- tibble(id = id, Time = unlist(lapply(split(time2 - time1, id), cumsum)) - origin,
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

#' @noRd
#' @keywords internal
#' Used to set default values
default.reSurv <- function(x) {
    text <- NULL
    if ("time2" %in% x) text <- c(text, "time2 <- NULL")
    if ("id" %in% x) text <- c(text, "id <- 1:n")
    if ("event" %in% x) text <- c(text, "event <- rep(1, n)")
    if ("status" %in% x) text <- c(text, "status <- rep(0, n)")
    text
}
