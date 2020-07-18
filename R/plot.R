globalVariables(c("time2", "Y", "Y.upper", "Y.lower", "id", "event", "MCF"))

#' Produce Event Plot or Mean Cumulative Function Plot
#'
#' Plot the event plot or the mean cumulative function (MCF) from an \code{Recur} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject" for event plot and "Cumulative mean" for MCF plot.}
#'   \item{main}{customizable title, the default value is "Recurrent event plot" when \code{mcf = FALSE} and
#' "Sample cumulative mean function plot" when \code{mcf = TRUE}.}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable legend title for recurrent event, default value is "Recurrent events".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is \code{NULL}.}
#'   \item{alpha}{between 0 and 1, controls the transparency of points.}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can also be passed down without specifying a \code{control} list. See \bold{Examples}.
#' 
#' @param x an object of class \code{Recur} returned by the \code{Recur()} function. See \code{?Recur} for creating \code{Recur} objects.
#' @param event.result an optional character string that is passed to the \code{plotEvents()} function as the \code{result} argument. See \code{\link{plotEvents}}.
#' This argument is used to specify whether the event plot is sorted by the subjects' terminal time.
#' The available options are
#' \describe{
#'   \item{\code{increasing}}{sort the terminal time from in ascending order (default). This places longer terminal times on top. }
#'   \item{\code{decreasing}}{sort the terminal time from in descending order (default). This places shorter terminal times on top. }
#'   \item{\code{asis}}{present the as is, without sorting.}
#' }
#' @param control a list of control parameters. See \bold{Details}.
#' @param mcf.smooth an optional logical value that is passed to the \code{plotMCF()} function as the \code{smooth} argument. See \code{\link{plotMCF}}.
#' This argument indicates whether to add a smooth curve obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param mcf.adjrisk an optional logical value that is passed to the \code{plotMCF()} function as the \code{adjrisk} argument. See \code{\link{plotMCF}}.
#' This argument indicates whether risk set size will be adjusted. If \code{mcf.adjrisk = TRUE}, subjects leave the risk set after terminal times as in the Nelson-Aalen estimator.
#' If \code{mcf.adjrisk = FALSE}, subjects remain in the risk set after terminal time. 
#' @param mcf an optional logical value indicating whether the mean cumulative function (MCF) will
#' be plotted instead of the event plot (default).
#' @param ... additional graphical parameters to be passed to methods.
#'  
#' @seealso \code{\link{Recur}}, \code{\link{plotEvents}}, \code{\link{plotMCF}}
#'
#' @references Nelson, W. B. (1995) Confidence Limits for Recurrence Data-Applied to Cost or Number of Product Repairs. \emph{Technometrics}, \bold{37}(2): 147--157.
#' 
#' @keywords Plots
#' @export
#'
#' @return A \code{ggplot} object.
#' @example inst/examples/ex_plot_reSurv.R
plot.Recur <- function(x, mcf = FALSE, event.result = c("increasing", "decreasing", "asis"),
                       mcf.adjrisk = TRUE, mcf.smooth = FALSE,
                       control = list(), ...) {
    result <- match.arg(event.result)
    if (!is.Recur(x)) stop("Response must be a `Recur` object.")
    if (!mcf) ctrl <- plotEvents.control()
    if (mcf) ctrl <- plotMCF.control()
    call <- match.call()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) call[[x]])
    }
    if (!mcf) {
        return(plotEvents(x, result = event.result, control = ctrl))
    }
    if (mcf)
        return(plotMCF(x, adjrisk = mcf.adjrisk, smooth = mcf.smooth, control = ctrl))
}

#' Produce Event Plots
#'
#' Plot the event plot for an \code{Recur} object.
#' The usage of the function is similar to that of \code{plot.Recur()} but with more flexible options.
#'
#' The argument \code{control} consists of options with argument defaults to a
#' list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{main}{customizable title, default value is "Recurrent event plot".}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable legend title for recurrent event, default value is "Recurrent events".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is \code{NULL}.}
#'   \item{alpha}{between 0 and 1, controls the transparency of points.}
#'   \item{legend}{a character string specifying the position of the legend (if any).
#'   The available options are "right", "left", "top", "bottom", and "none". The default value is "right".}
#' }
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator,
#' and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{Recur()}.
#' @param data an optional data frame in which to interpret the variables occurring in the "\code{formula}".
#' @param result an optional character string specifying whether the event plot is sorted by the subjects' terminal time. The available options are
#' \describe{
#'   \item{\code{increasing}}{sort the terminal time from in ascending order (default). This places longer terminal times on top. }
#'   \item{\code{decreasing}}{sort the terminal time from in descending order (default). This places shorter terminal times on top. }
#'   \item{\code{asis}}{present the as is, without sorting.}
#' }
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#'
#' @seealso \code{\link{Recur}}, \code{\link{plot.Recur}}
#' 
#' @keywords Plots
#' @export
#' 
#' @return A \code{ggplot} object.
#' 
#' @example inst/examples/ex_plot_event.R
plotEvents <- function(formula, data, result = c("increasing", "decreasing", "asis"),
                       control = list(), ...) {
    result <- match.arg(result)
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) eval(call[[x]]))
    }
    nX <- 0
    if (is.Recur(formula)) {
        DF <- as.data.frame(formula@.Data)
        vNames <- NULL
    } else {
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        if (!is.Recur(obj)) stop("Response must be a `Recur` object.")
        nX <- length(formula[[3]])
        if (formula[[3]] == 1) DF <- as.data.frame(obj@.Data)
        if (formula[[3]] != 1 && nX == 1) {
            if (missing(data)) DF <- cbind(obj@.Data, eval(formula[[3]], parent.frame()))
            if (!missing(data)) DF <- cbind(obj@.Data, eval(formula[[3]], data))
            colnames(DF) <- c(colnames(obj@.Data), paste0(formula[[3]], collapse = ""))
            DF <- as.data.frame(DF)
        }
        if (formula[[3]] != 1 && nX > 1) {
            DF <- as.data.frame(obj@.Data)
            if (missing(data)) {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], parent.frame()))
                }
            } else {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], data))
                }
            }
            vNames <- attr(terms(formula), "term.labels")
            colnames(DF) <- c(colnames(obj@.Data), vNames)
        }
        vNames <- attr(terms(formula), "term.labels")
        if (length(vNames) == 0) vNames <- NULL
    }    
    ## dat$status <- ifelse(is.na(dat$status), 0, dat$status)
    ## dat$Yi <- ifelse(is.na(dat$Yi), unlist(lapply(dat$tij, max)), dat$Yi)
    newIDtime2 <- function(dat, result = "increasing") {
        if (result == "asis") {
            tmp <- table(dat$id)
            dat$id <- rep(1:length(tmp), tmp[match(unique(dat$id), names(tmp))])
            return(dat)
        } else {
            tmp <- rank(dat$time2[dat$event == 0], ties.method = "first")
        }
        if (result == "decreasing") tmp <- length(tmp) - tmp + 1
        dat$id <- rep(tmp, table(dat$id))
        return(dat)
    }
    if (nX == 0 || formula[[3]] == 1) {
        DF <- newIDtime2(DF, result = result)
    } else {
        DF <- do.call(rbind,
                      lapply(split(DF, DF[, 7:ncol(DF)], drop = TRUE), newIDtime2, result = result))
        rownames(DF) <- NULL
    }
    if (ctrl$cex == "Default") sz <- 1 + 8 / (1 + exp(length(unique(DF$id)) / 30)) / max(1, nX)
    else sz <- ctrl$cex
    k <- length(unique(DF$event)) - 1 ## exclude event = 0
    shp.val <- c(17, rep(19, k))
    clr.val <- c(alpha("red", ctrl$alpha), hcl(h = seq(120, 360, length.out = k),
                                               l = 60, alpha = ctrl$alpha))
    rec.lab <- paste("r", 1:k, sep = "")
    if (k == 1) {
        shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.name)
    } 
    if (k > 1 & is.null(ctrl$recurrent.type)) {
        shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))
    }
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.type)
        } else {
            cat('The length of "recurrent.type" mismatched, default names are used.\n')
            shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))            
        }
    }
    if (nX > 0) {
        for (i in vNames) {
            DF[,i] <- factor(DF[,i], labels = paste(i, "=", unique(DF[,i])))
        }}
    names(shp.val) <- names(clr.val) <- c("terminal", rec.lab)
    gg <- ggplot(DF[DF$event == 0,], aes(id, time2)) +
        geom_bar(stat = "identity", fill = "gray75") +
        coord_flip() + 
        theme(axis.line.y = element_blank(),
              axis.title.y = element_text(vjust = 0),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (any(table(DF$id) > 0))
        gg <- gg + geom_point(data = DF[DF$event > 0,],
                              aes(id, time2,
                                  shape = factor(event, labels = rec.lab),
                                  color = factor(event, labels = rec.lab)),
                              size = sz)    
    if (sum(DF$terminal, na.rm = TRUE) > 0)
        gg <- gg + geom_point(data = DF[DF$terminal > 0,], 
                              aes(id, time2, shape = "terminal", color = "terminal"),
                              size = sz)    
    if (nX > 0 && formula[[3]] != 1)        
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "both")
    gg <- gg + scale_shape_manual(name = "", values = shp.val,
                                  labels = shp.lab, breaks = c("terminal", rec.lab)) +
        scale_color_manual(name = "", values = clr.val,
                           labels = shp.lab, breaks = c("terminal", rec.lab))
    gg + theme(panel.background = element_blank(),
               axis.line = element_line(color = "black"), legend.position = ctrl$legend, 
               legend.key = element_rect(fill = "white", color = "white")) +
        scale_x_continuous(expand = c(0, 1)) +
        ggtitle(ctrl$main) + labs(x = ctrl$ylab, y = ctrl$xlab) +
        guides(shape = guide_legend(override.aes = list(size = 2.7)))
}

#' Produce Cumulative Sample Mean Function Plots
#'
#' Plot the mean cumulative function (MCF) for an \code{Recur} object.
#' The usage of the function is similar to that of \code{plot.Recur()} but with more flexible options.
#'
#' When \code{adjrisk = TRUE}, the \code{plotMCF()} is equivalent to
#' the Nelson-Aalen estimator for the intensity function of the recurrent event process.
#' When \code{adjrisk = FALSE}, the \code{plotMCF()} does not adjust for the risk set and
#' assumes all subjects remain at risk after the last observed recurrent event.
#' This is known as the survivor rate function.
#' The argument \code{control} consists of options with argument defaults
#' to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Cumulative mean".}
#'   \item{main}{customizable title, default value is "Sample cumulative mean function plot".}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can also be passed down without specifying a \code{control} list.
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator, and the predictors on the right.
#' The response must be a recurrent event survival object returned by the \code{Recur()} function.
#' @param data an optional data frame in which to interpret the variables occurring in the "\code{formula}".
#' @param adjrisk an optional logical value that is passed to the \code{plotMCF()} function as the \code{adjrisk} argument. See \code{\link{plotMCF}}.
#' This argument indicates whether risk set size will be adjusted. If \code{mcf.adjrisk = TRUE}, subjects leave the risk set after terminal times as in the Nelson-Aalen estimator.
#' If \code{mcf.adjrisk = FALSE}, subjects remain in the risk set after terminal time. 
#' @param smooth an optional logical value indicating whether to add a smooth curve obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' This feature only works for data with one recurrent event type.
#' @param onePanel an optional logical value indicating whether the mean cumulative functions will be plotted in the same panel.
#' This is only useful when there are multiple recurrent event types or in the presence of (discrete) covariates.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#' 
#' @seealso \code{\link{Recur}}, \code{\link{plot.Recur}}
#' @keywords Plots
#' @export
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom ggplot2 guides guide_legend
#' @importFrom scam scam
#' 
#' @example inst/examples/ex_plot_MCF.R
plotMCF <- function(formula, data, onePanel = FALSE, adjrisk = TRUE,
                     smooth = FALSE, control = list(), ...) {
    call <- match.call()
    ctrl <- plotMCF.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) eval(call[[x]]))
    }
    nX <- 0
    if(is.Recur(formula)) {
        DF <- as.data.frame(formula@.Data)
        vNames <- NULL
    } else {
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        if (!is.Recur(obj)) stop("Response must be a `Recur` object.")
        DF <- as.data.frame(obj@.Data)
        nX <- length(formula[[3]])
        if (formula[[3]] != 1 && nX == 1) {
            if (missing(data)) DF <- cbind(obj@.Data, tmp = eval(formula[[3]], parent.frame()))
            if (!missing(data)) DF <- cbind(obj@.Data, tmp = eval(formula[[3]], data))
            colnames(DF) <-  c(colnames(obj@.Data), paste0(formula[[3]], collapse = ""))
            DF <- as.data.frame(DF)
        }
        if (formula[[3]] != 1 && nX > 1) {
            DF <- as.data.frame(obj@.Data)
            if (missing(data)) {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], parent.frame()))
                }
            } else {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], data))
                }
            }
            colnames(DF) <- c(colnames(obj@.Data), attr(terms(formula), "term.labels"))
        }
        vNames <- attr(terms(formula), "term.labels")
        if (length(vNames) == 0) vNames <- NULL
    }   
    ## dd <- subset(DF, select = c("event", vNames, "time2"))
    dd <- DF[,c("event", vNames, "time2")]
    dd <- dd[do.call(order, dd),]
    nn <- table(apply(dd, 1, paste, collapse = ""))
    dd <- unique(dd)
    dd$n <- as.integer(nn) ## tmp1 in the 1st version
    rownames(dd) <- NULL
    k <- length(unique(dd$event)) - 1    
    if (!is.null(vNames)) { ## any covariates/stratifications?
        dd2 <- DF[DF$event == 0, vNames, drop = FALSE]
        dd2 <- dd2[do.call(order, dd2),, drop = FALSE]
        nn <- table(apply(dd2, 1, paste, collapse = ""))
        dd2 <- unique(dd2)
        dd2$n <- as.integer(nn) ## tmp2 in the 1st version
        rownames(dd2) <- NULL
        dd$GrpInd <- match(apply(dd[,vNames, drop = FALSE], 1, paste, collapse = ""),
                           apply(dd2[,vNames, drop = FALSE], 1, paste, collapse = ""))
        tmp1 <- merge(dd, dd2, by = vNames)
        ## as.integer(eval(parse(text = paste(attr(terms(formula), "term.labels"), collapse = ":"))))
        rec0 <- tmp1[tmp1$event == 0,]
        dat0 <- do.call(rbind, lapply(split(tmp1, tmp1$GrpInd), function(x) {
            x$adjrisk = apply(x, 1, function(y)
                as.numeric(y['n.y']) - sum(rec0$n.x[as.numeric(y['time2'])> rec0$time2 & rec0$GrpInd == as.numeric(y['GrpInd'])]))
            return(x)}))
        dat0$n.x <- dat0$n.x * (dat0$event > 0)
        rec0$time2 <- rec0$n.x <- 0
        rec0$adjrisk <- 1
        dat0 <- unique(rbind(dat0, rec0))
        rownames(dat0) <- NULL
        ## Number of recurrent types
        if (k > 1) {
            tmp <- dat0[dat0$event == 0,]
            tmp <- tmp[rep(1:NROW(tmp), k),]
            tmp$event <- rep(1:k, each = sum(dat0$event == 0))
            dat0 <- rbind(dat0[dat0$event > 0,], tmp)            
        } else {
            dat0$event <- dat0$event[1]
        }
        dat0 <- dat0[order(dat0$event, dat0$GrpInd, dat0$time2),]
        if (adjrisk) {
            dat0 <- do.call(rbind, 
                            lapply(split(dat0, list(dat0$event, dat0$GrpInd)), function(x) {
                                x$mu <- x$n.x / x$adjrisk
                                x$MCF <- cumsum(x$mu)
                                return(x)}))
        } else {
            dat0 <- do.call(rbind, 
                            lapply(split(dat0, list(dat0$event, dat0$GrpInd)), function(x) {
                                x$mu <- x$n.x / x$n.y
                                x$MCF <- cumsum(x$mu) 
                                return(x)}))
        }
    } else { ## no covariates
        dd$n.y <- length(unique(DF$id))
        rec0 <- dd[dd$event == 0, ]
        dd$adjrisk <- apply(dd, 1, function(x) x[4] - sum(rec0$n[x[2] > rec0$time2]))
        dd$n <- dd$n * (dd$event > 0)
        rec0$time2 <- rec0$n <- 0
        rec0$adjrisk <- 1
        dat0 <- unique(rbind(dd, rec0))       
        dat0 <- dat0[order(dat0$time2),]
        dat0$mu <- dat0$n / (adjrisk * dat0$adjrisk + (!adjrisk) * dat0$n.y)
        dat0$MCF <- cumsum(dat0$mu)
    }
    if (k == 1) {
        dat0$event <- dat0$event[1]
        dat0$event <- factor(dat0$event, labels = ctrl$recurrent.name)
    }
    if (k > 1) {
        dat00 <- dat0[dat0$event == 0,]
        dat0 <- dat0[dat0$event > 0,]
        if (nrow(dat00) > 0) {
            dat00 <- dat00[rep(1:nrow(dat00), k),]
            dat00$event <- rep(1:k, each = nrow(dat00) / k)
            dat0 <- rbind(dat0, dat00)
        }
        if (is.null(ctrl$recurrent.type))
            dat0$event <- factor(dat0$event, labels = paste(ctrl$recurrent.name, 1:k))
        if (!is.null(ctrl$recurrent.type)) {
            if (length(ctrl$recurrent.type) == k) {
                dat0$event <- factor(dat0$event, labels = ctrl$recurrent.type)
            } else {
                cat('The length of "recurrent.type" mismatched, default names are used.\n')
                dat0$event <- factor(dat0$event, labels = paste(ctrl$recurrent.name, 1:k))
            }
        }
    }
    if (nX > 0) {
        for (i in vNames) {
            dat0[,i] <- factor(dat0[,i], labels = paste(i, "=", unique(dat0[,i])))
        }}
    gg <- ggplot(data = dat0, aes(x = time2, y = MCF))
    if (is.null(vNames) & k == 1) {
        gg <- gg + geom_step(size = ctrl$lwd)
    } else {
        if (!onePanel & k == 1) gg <- gg + geom_step(size = ctrl$lwd)
        if (!onePanel & k > 1) 
            gg <- gg + geom_step(aes(color = event), direction = "hv", size = ctrl$lwd) +
                guides(color = guide_legend(title = ctrl$recurrent.name))
        if (onePanel) {
            dat0$GrpInd <- factor(dat0$GrpInd)
            levels(dat0$GrpInd) <- apply(unique(dat0[,vNames, drop = FALSE]), 1, paste, collapse = ", ")
            gg <- gg + geom_step(aes(color = dat0$GrpInd), direction = "hv", size = ctrl$lwd) +
                guides(color = guide_legend(title = ""))
        }
    }
    if (!onePanel && nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "x")
    if (!onePanel & k == 1)
        gg <- gg + theme(legend.position="none")
    ## if (onePanel & k == 1)
    ##     gg <- gg + scale_color_discrete(name = "", labels = levels(interaction(vNames)))
    ## if (onePanel & k > 1) gg <- gg + scale_color_discrete(name = "")
    if (smooth & k == 1 & !onePanel) {
        if (is.null(dat0$GrpInd)) dat0$GrpInd <- 1
        dat0 <- do.call(rbind, lapply(split(dat0, dat0$GrpInd), function(x){
            x$bs <- scam(x$MCF ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values
            return(x)}))       
        gg <- gg + geom_line(aes(time2, y = dat0$bs), color = 4, size = ctrl$lwd)
        ## geom_smooth(method = "scam", formula = y ~ s(x, k = 10, bs = "mpi"), size = ctrl$lwd, se = FALSE)
    }
    ## gg <- gg + geom_smooth(method = "loess", size = ctrl$lwd, se = FALSE)
    if (smooth & k > 1) cat('Smoothing only works for data with one recurrent event type.\n')
    gg + theme(axis.line = element_line(color = "black"),
                legend.key = element_rect(fill = "white", color = "white")) +
        ggtitle(ctrl$main) + labs(y = ctrl$ylab, x = ctrl$xlab)
}

#' Plot the Baseline Cumulative Rate Function and the Baseline Cumulative Hazard Function
#'
#' Plot the baseline cumulative rate function and the baseline cumulative hazard function
#' (if applicable) for an \code{reReg} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty.}
#'   \item{main}{customizable title, default value are "Baseline cumulative rate and hazard function" when \code{baseline = "both"},
#' "Baseline cumulative rate function" when \code{baseline = "rate"}, and "Baseline cumulative hazard function" when \code{baseline = "hazard"}.}
#' }
## #' These arguments can also be passed down without specifying a \code{control} list.
#' 
#' @param x an object of class \code{reReg}, returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether to add a smooth curve obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param baseline a character string specifying which baseline function to plot.
#' \describe{
#'   \item{\code{baseline = "both"}}{plot both the baseline cumulative rate and the baseline cumulative hazard function (if applicable) in separate panels within the same display (default).}
#'   \item{\code{baseline = "rate"}}{plot the baseline cumulative rate function.}
#'   \item{\code{baseline = "hazard"}}{plot the baseline cumulative hazard function.}
#' }
#' @param control a list of control parameters. See \bold{Details}.
#' @param ... additional graphical parameters to be passed to methods.
#' 
#' @seealso \code{\link{reReg}}
#' @export
#' @keywords Plots
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom ggplot2 geom_smooth geom_step
#' @example inst/examples/ex_plot_reReg.R
plot.reReg <- function(x, baseline = c("both", "rate", "hazard"),
                       smooth = FALSE, control = list(), ...) {
    baseline <- match.arg(baseline)
    if (x$recType %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline functions not available for this method.")
    if (baseline == "both") {
        ctrl <- plot.reReg.control(main = "Baseline cumulative rate and cumulative hazard functions")
        if (x$method == "cox.LWYY")
            ctrl <- plot.reReg.control(main = "Baseline cumulative rate function")
        smooth  <- FALSE
    }
    if (baseline == "rate") ctrl <- plot.reReg.control(main = "Baseline cumulative rate function")
    if (baseline == "hazard") ctrl <- plot.reReg.control(main = "Baseline cumulative hazard function")
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) call[[x]])
    }
    if (!is.reReg(x)) stop("Response must be a `reReg` class")
    if (baseline == "rate")
        return(plotRate(x, smooth = smooth, control = ctrl))
    if (baseline == "hazard")
        return(plotHaz(x, smooth = smooth, control = ctrl))
    if (baseline == "both" & x$temType == ".") {
        cat(paste("Baseline cumulative hazard function is not available."))
        cat("\nOnly the baseline cumulative rate function is plotted.\n")
        return(plotRate(x, smooth = smooth, control = ctrl))
    }
    dat1 <- dat2 <- x$DF[, "time2", drop = FALSE]
    dat1$Y <- x$Lam0(dat1$time2)
    if (!is.null(x$Lam0.upper)) {
        dat1$Y.upper <- x$Lam0.upper(dat1$time2)
        dat1$Y.lower <- x$Lam0.lower(dat1$time2)
    }
    if (x$temType != ".") {
        dat2$Y <- x$Haz0(dat1$time2)
        if (!is.null(x$Haz0.upper)) {    
            dat2$Y.upper <- x$Haz0.upper(dat2$time2)
            dat2$Y.lower <- x$Haz0.lower(dat2$time2)
        }
        dat <- rbind(dat1, dat2)
        dat$group <- c(rep(1, nrow(dat1)), rep(2, nrow(dat2)))
        dat$group <- factor(dat$group, levels = 1:2,
                            labels = c("Baseline cumulative rate", "Baseline cumulative hazard"))
    }
    gg <- ggplot(data = dat, aes(x = time2, y = Y)) +
        facet_grid(group ~ ., scales = "free") +
        theme(axis.line = element_line(color = "black"),
              strip.text = element_text(face = "bold", size = 12))   
    if (smooth) {
        dat <- do.call(rbind, lapply(split(dat, dat$group), function(x){
            x$bs <- scam(x$Y ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values
            return(x)}))
        gg <- gg + geom_line(aes(time2, y = dat$bs), color = 4, size = ctrl$lwd)
        if (!is.null(x$Lam0.upper)) {
            dat <- do.call(rbind, lapply(split(dat, dat$group), function(x){
                x$bs.upper <- scam(x$Y.upper ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values
                return(x)}))
            gg <- gg + geom_line(aes(time2, y = dat$bs.upper), color = 4, size = ctrl$lwd, lty = 2)
        }
        if (!is.null(x$Lam0.lower)) {
            dat <- do.call(rbind, lapply(split(dat, dat$group), function(x) {
                x$bs.lower <- scam(x$Y.lower ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values
                return(x)}))
            gg <- gg + geom_line(aes(time2, y = dat$bs.lower), color = 4, size = ctrl$lwd, lty = 2)
        }
    } else {
        gg <- gg + geom_step()
        if (!is.null(x$Lam0.upper))
            gg <- gg + geom_step(aes(x = time2, y = Y.upper), lty = 2)+ 
                geom_step(aes(x = time2, y = Y.lower), lty = 2)
    }
    gg + ggtitle(ctrl$main) + labs(y = ctrl$ylab, x = ctrl$xlab)
}

#' Plotting the Baseline Cumulative Rate Function for the Recurrent Event Process
#'
#' Plot the baseline rate function for an \code{reReg} object.
#'
#' 
#' The \code{plotRate()} plots the estimated baseline cumulative rate function 
#' depending on the identifiability assumption.
#' When \code{type = "unrestricted"} (default), the baseline cumulative rate function
#' is plotted under the assumption \eqn{E(Z) = 1}.
#' When \code{type = "scaled"}, the baseline cumulative rate function is plotted
#' under the assumption \eqn{\Lambda(\min(Y^\ast, \tau)) = 1}.
#' When \code{type = "raw"}, the baseline cumulative rate function is plotted
#' under the assumption \eqn{\Lambda(\tau) = 1}.
#' See \code{?reReg} for the specification of the notations and underlying models.
#' 
#' The argument \code{control} consists of options with argument defaults
#' to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty.}
#'   \item{main}{customizable title, default value is "Baseline cumulative rate function".}
#' }
#' These arguments can also be passed down without specifying a \code{control} list. See \bold{Examples}.
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param type a character string specifying the type of rate function to be plotted.
#' Options are "unrestricted", "scaled", "raw". See \bold{Details}.
#' @param smooth an optional logical value indicating whether to add a smooth curve obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#'
#' @seealso \code{\link{reReg}} \code{\link{plot.reReg}}
#' @export
#'
#' @return A \code{ggplot} object.
#'
#' @keywords Plots
#' 
#' @example inst/examples/ex_plot_rate.R
plotRate <- function(x, type = c("unrestricted", "scaled", "raw"),
                     smooth = FALSE, control = list(), ...) {
    if (x$recType %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline cumulative rate function is not available.")
    ctrl <- plot.reReg.control(main = "Baseline cumulative rate function")
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) call[[x]])
    }
    type <- match.arg(type)
    if (!is.reReg(x)) stop("Response must be a `reReg` class")    
    dat <- x$DF[,"time2",drop = FALSE]
    if (type == "unrestricted") dat$Y <- x$Lam0(dat$time2) * exp(x$log.muZ)
    if (type == "scaled") dat$Y <- x$Lam0(dat$time2) / x$Lam0(max(dat$time2))
    if (type == "raw") dat$Y <- x$Lam0(dat$time2)
    if (!is.null(x$Lam0.upper)) {
        if (type == "raw") {
            dat$Y.upper <- x$Lam0.upper(dat$time2)
            dat$Y.lower <- x$Lam0.lower(dat$time2)
        }
        if (type == "unrestricted") {
            dat$Y.upper <- x$Lam0.upper(dat$time2) * exp(x$log.muZ)
            dat$Y.lower <- x$Lam0.lower(dat$time2) * exp(x$log.muZ)
        }
        if (type == "scaled") {
            dat$Y.upper <- x$Lam0.upper(dat$time2) / x$Lam0(max(dat$time2))
            dat$Y.lower <- x$Lam0.lower(dat$time2) / x$Lam0(max(dat$time2))
        }
    }    
    gg <- ggplot(data = dat, aes(x = time2, y = Y)) +
        theme(axis.line = element_line(color = "black"))
    if (smooth) {
        dat$bs <- scam(dat$Y ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
        gg <- gg + geom_line(aes(time2, y = dat$bs), color = 4)
        if (!is.null(x$Lam0.upper)) {
            dat$bs.upper <- scam(dat$Y.upper ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
            dat$bs.lower <- scam(dat$Y.lower ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
            gg <- gg + geom_line(aes(time2, y = dat$bs.upper), color = 4, lty = 2) + 
                geom_line(aes(time2, y = dat$bs.lower), color = 4, lty = 2)
        }
    } else {
        gg <- gg + geom_step()
        if (!is.null(x$Lam0.upper))
            gg <- gg + geom_step(aes(x = time2, y = Y.upper), lty = 2) +
                geom_step(aes(x = time2,  y = Y.lower), lty = 2)
    }
    gg + ggtitle(ctrl$main) + labs(x = ctrl$xlab, y = ctrl$ylab)
}

#' Plot the Baseline Cumulative Hazard Function for the Terminal Time
#'
#' Plot the baseline cumulative hazard function for an \code{reReg} object.
#' The 95\% confidence interval on the baseline cumulative rate function
#'
#' The argument \code{control} consists of options with argument
#' defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty.}
#'   \item{main}{customizable title, default value is "Baseline cumulative hazard function".}
#' }
#' These arguments can also be passed down without specifying a \code{control} list.
#' See \bold{Examples}.
#'
#' @param x an object of class \code{reReg}, returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether to add a smooth curve obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#' 
#' @seealso \code{\link{reReg}} \code{\link{plot.reReg}}
#' @export
#'
#' @return A \code{ggplot} object.
#' @keywords Plots
#' 
#' @example inst/examples/ex_plot_Haz.R
plotHaz <- function(x, smooth = FALSE, control = list(), ...) {
    if (x$recType %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline cumulative hazard function is not available.")
    if (x$temType == ".") {
        stop("Baseline cumulative hazard function is not available.")
    }
    ctrl <- plot.reReg.control(main = "Baseline cumulative hazard function")
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    namp <- names(match.call())
    if (any(namp %in% names(ctrl))) {
        namp <- namp[namp %in% names(ctrl)]
        ctrl[namp] <- lapply(namp, function(x) call[[x]])
    }
    if (!is.reReg(x)) stop("Response must be a `reReg` class.")
    dat <- x$DF[, "time2", drop = FALSE]
    dat$Y <- x$Haz0(dat$time2)
    if (!is.null(x$Haz0.upper)) {
        dat$Y.upper <- x$Haz0.upper(dat$time2)
        dat$Y.lower <- x$Haz0.lower(dat$time2)
    }
    gg <- ggplot(data = dat, aes(x = time2, y = Y)) +
        theme(axis.line = element_line(color = "black"))
    if (smooth) {
        dat$bs <- scam(dat$Y ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
        gg <- gg + geom_line(aes(time2, y = dat$bs), color = 4)
        if (!is.null(x$Lam0.upper)) {
            dat$bs.upper <- scam(dat$Y.upper ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
            dat$bs.lower <- scam(dat$Y.lower ~ s(dat$time2, k = 10, bs = "mpi"))$fitted.values
            gg <- gg + geom_line(aes(time2, y = dat$bs.upper), color = 4, lty = 2) +
                geom_line(aes(time2, y = dat$bs.lower), color = 4, lty = 2)
        }
    } else {
        gg <- gg + geom_step()
        if (!is.null(x$Lam0.upper))
            gg <- gg + geom_step(aes(x = time2,  y = Y.upper), lty = 2) +
                geom_step(aes(x = time2,  y = Y.lower), lty = 2)
    }
    gg + ggtitle(ctrl$main) + labs(x = ctrl$xlab, y = ctrl$ylab)
}

plotEvents.control <- function(xlab = "Time", ylab = "Subject",
                               main = "Recurrent event plot",
                               cex = "Default", 
                               terminal.name = "Terminal event",
                               recurrent.name = "Recurrent events",
                               recurrent.type = NULL, alpha = .7,
                               legend = "right") {
    list(xlab = xlab, ylab = ylab, main = main, cex = cex,
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type, alpha = alpha, legend = legend)
}

plotMCF.control <- function(xlab = "Time", ylab = "Cumulative mean",
                            main = "Sample cumulative mean function plot",
                            lwd = 1,
                            terminal.name = "Terminal event",
                            recurrent.name = "Recurrent events",
                            recurrent.type = NULL,
                            legend = "right") {
    list(xlab = xlab, ylab = ylab, main = main, lwd = lwd, 
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type, legend = legend)
}

plot.reReg.control <- function(xlab = "Time", ylab = "", main = "") {
    list(xlab = xlab, ylab = ylab, main = main)
}

