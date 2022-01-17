globalVariables(c("time1", "time2", "group", "Y", "Y.upper", "Y.lower", "id", "event", "origin", "MCF", "bs"))

#' Produce Event Plot or Mean Cumulative Function Plot
#'
#' Plot the event plot or the mean cumulative function (MCF) from an \code{Recur} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with
#' the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject" for event plot and
#' "Cumulative mean" for MCF plot.}
#'   \item{main}{customizable title, the default value is "Recurrent event plot"
#' when \code{mcf = FALSE} and
#' "Sample cumulative mean function plot" when \code{mcf = TRUE}.}
#'   \item{terminal.name}{customizable label for terminal event,
#' the default value is "Terminal event".}
#'   \item{recurrent.name}{customizable legend title for recurrent event,
#' the default value is "Recurrent events".}
#'   \item{recurrent.types}{customizable label for recurrent event type,
#' the default value is \code{NULL}.}
#'   \item{alpha}{between 0 and 1, controls the transparency of points.}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can be specified
#' outside of the \code{control} list. 
#' 
#' @param x an object of class \code{Recur} returned by the \code{Recur()} function.
#' See \code{?Recur} for creating \code{Recur} objects.
#' @param event.result an optional character string that is passed to the
#' \code{plotEvents()} function as the \code{result} argument. See \code{\link{plotEvents}}.
#' This argument is used to specify whether the event plot is sorted by the subjects' terminal time.
#' The available options are
#' \describe{
#'   \item{\code{increasing}}{sort the terminal time from in ascending order (default).
#' This places longer terminal times on top. }
#'   \item{\code{decreasing}}{sort the terminal time from in descending order.
#' This places shorter terminal times on top. }
#'   \item{\code{none}}{present the event plots as is, without sorting by the terminal times.}
#' }
#' @param event.calendarTime an optional logical value indicating whether to plot in calendar time.
#' When \code{event.calendarTime = FALSE} (default),
#' the event plot will have patient time on the x-axis.
#' @param control a list of control parameters. See \bold{Details}.
#' @param mcf.adjustRiskset an optional logical value that is passed to
#' the \code{mcf()} function as the \code{adjustRiskset} argument. 
#' This argument indicates whether risk set size will be adjusted.
#' If \code{mcf.adjustRiskset = TRUE}, subjects leave the risk set after terminal times
#' as in the Nelson-Aalen estimator.
#' If \code{mcf.adjustRiskset = FALSE}, subjects remain in the risk set after terminal time.
#' @param mcf.conf.int an optional logical value that is passed to
#' the \code{mcf()} function as the \code{conf.int} argument. See \code{\link{mcf}} for details.
#' @param mcf an optional logical value indicating whether the mean cumulative function (MCF) will
#' be plotted instead of the event plot. When \code{mcf = TRUE},
#' the \code{mcf} is internally called. See \code{\link{mcf}} for details.
#' @param ... additional graphical parameters to be passed to methods.
#'  
#' @seealso \code{\link{Recur}}, \code{\link{plotEvents}}, \code{\link{mcf}}
#'
#' @references Nelson, W. B. (1995) Confidence Limits for Recurrence Data-Applied to Cost
#' or Number of Product Repairs. \emph{Technometrics}, \bold{37}(2): 147--157.
#' 
#' @keywords Plots
#' @export
#'
#' @return A \code{ggplot} object.
#' @example inst/examples/ex_plot_Recur.R
#' @method plot Recur
plot.Recur <- function(x, mcf = FALSE,
                       event.result = c("increasing", "decreasing", "asis"),
                       event.calendarTime = FALSE, 
                       mcf.adjustRiskset = TRUE,
                       mcf.conf.int = FALSE,
                       control = list(), ...) {
    event.result <- match.arg(event.result)
    if (!is.Recur(x)) stop("Response must be a `Recur` object.")
    ctrl <- plotEvents.control()
    if (mcf & is.null(ctrl$ylab)) ctrl <- ctrl$ylab <- "Cumulative mean"
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
        return(plotEvents(x, result = event.result, calendarTime = event.calendarTime, control = ctrl))
    }
    if (mcf) {
        return(plot(mcf(x ~ 1, adjustRiskset = mcf.adjustRiskset), conf.int = mcf.conf.int))
        ## return(plotMCF(x, adjustRiskset = mcf.adjustRiskset, smooth = mcf.smooth, control = ctrl))
    }
}

#' Produce Event Plots
#'
#' Plot the event plot for an \code{Recur} object.
#' The usage of the function is similar to that of \code{plot.Recur()} but with more flexible options.
#'
#' The argument \code{control} consists of options with argument defaults to a list with
#' the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject" for event plot and
#' "Cumulative mean" for MCF plot.}
#'   \item{main}{customizable title, the default value is "Recurrent event plot"
#' when \code{mcf = FALSE} and
#' "Sample cumulative mean function plot" when \code{mcf = TRUE}.}
#'   \item{terminal.name}{customizable label for terminal event,
#' the default value is "Terminal event".}
#'   \item{recurrent.name}{customizable legend title for recurrent event,
#' the default value is "Recurrent events".}
#'   \item{recurrent.types}{customizable label for recurrent event type,
#' the default value is \code{NULL}.}
#'   \item{alpha}{between 0 and 1, controls the transparency of points.}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can be specified
#' outside of the \code{control} list. 
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator,
#' and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{Recur()}.
#' @param data an optional data frame in which to interpret the variables occurring in
#' the "\code{formula}".
#' @param result an optional character string specifying whether the event plot is
#' sorted by the subjects' terminal time. The available options are
#' \describe{
#'   \item{\code{increasing}}{sort the terminal time from in ascending order (default).
#' This places longer terminal times on top. }
#'   \item{\code{decreasing}}{sort the terminal time from in descending order.
#' This places shorter terminal times on top. }
#'   \item{\code{none}}{present the event plots as is, without sorting by the terminal times.}
#' }
#' @param calendarTime an optional logical value indicating whether to plot in calendar time.
#' When \code{calendarTime = FALSE} (default), the event plot will have patient time on the x-axis.
#' @param control a list of control parameters. See \bold{Details}.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#'
#' @seealso \code{\link{Recur}}, \code{\link{plot.Recur}}
#' 
#' @keywords Plots
#' @export
#' 
#' @return A \code{ggplot} object.
#' @importFrom ggplot2 geom_rect ggplot_build scale_y_continuous unit
#' @example inst/examples/ex_plot_event.R
plotEvents <- function(formula, data, result = c("increasing", "decreasing", "none"),
                       calendarTime = FALSE, control = list(), ...) {
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
        isDate <- "Date" %in% formula@time_class
        vNames <- NULL
    } else {
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        if (!is.Recur(obj)) stop("Response must be a `Recur` object.")
        nX <- length(formula[[3]])
        isDate <- "Date" %in% obj@time_class
        if (formula[[3]] == 1) DF <- as.data.frame(obj@.Data)
        if (formula[[3]] != 1 && nX == 1) {
            if (missing(data)) DF <- data.frame(obj@.Data, eval(formula[[3]], parent.frame()))
            if (!missing(data)) DF <- data.frame(obj@.Data, eval(formula[[3]], data))
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
        if (result == "none") {
            tmp <- table(dat$id)
            dat$id <- rep(1:length(tmp), tmp[match(unique(dat$id), names(tmp))])
            return(dat)
        } else {
            dat <- dat[order(dat$id),]
            ## if (all(dat$origin == 0))
            if (!calendarTime)
                tmp <- rank((dat$time2 - dat$origin)[dat$event == 0], ties.method = "first")
            else
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
    if (is.null(ctrl$cex)) sz <- 1 + 8 / (1 + exp(length(unique(DF$id)) / 30)) / max(1, nX)
    else sz <- ctrl$cex
    k <- length(unique(DF$event)) - 1 ## exclude event = 0
    shp.val <- c(17, rep(19, k))
    clr.val <- c(alpha("red", ctrl$alpha), hcl(h = seq(120, 360, length.out = k),
                                               l = 60, alpha = ctrl$alpha))
    if (k == 0) rec.lab <- NULL
    else rec.lab <- paste("r", 1:k, sep = "")
    if (k == 0) shp.lab <- ctrl$terminal.name
    if (k == 1) shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.name)
    if (k > 1 & is.null(ctrl$recurrent.type)) {
        shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))
    }
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            shp.lab <- c(ctrl$terminal.name, ctrl$recurrent.type)
        } else {
            message('The length of "recurrent.type" mismatched, default names are used.')
            shp.lab <- c(ctrl$terminal.name, paste(ctrl$recurrent.name, 1:k))            
        }
    }
    if (nX > 0) {
        for (i in vNames) {
            DF[,i] <- factor(DF[,i], labels = paste(i, "=", unique(DF[,i])))
        }}
    names(shp.val) <- names(clr.val) <- c("terminal", rec.lab)
    ## Bars
    if (calendarTime)
        gg <- ggplot(DF, aes(xmin = id - .45, xmax = id + .45, ymin = time1, ymax = time2)) +
            geom_rect(fill = "gray75") + coord_flip()
    else gg <- ggplot(DF[DF$event == 0,], aes(id, time2 - origin)) +
             geom_bar(stat = "identity", fill = "gray75") +
             coord_flip()
    ## event dots
    if (!calendarTime) DF$time2 <- DF$time2 - DF$origin
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
    ## Add theme and final touch ups
    if (ctrl$main != "") gg <- gg + ggtitle(ctrl$main) 
    gg <- gg + scale_shape_manual(name = "", values = shp.val,
                                  labels = shp.lab, breaks = c("terminal", rec.lab)) +
        scale_color_manual(name = "", values = clr.val,
                           labels = shp.lab, breaks = c("terminal", rec.lab)) +
        theme(panel.background = element_blank(),
              axis.line = element_line(color = "black"),
              legend.position = ctrl$legend.position, 
              legend.key = element_rect(fill = "white", color = "white"),
              axis.line.y = element_blank(),
              axis.title.y = element_text(vjust = 0),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              plot.title = element_text(size = 2 * ctrl$base_size),
              strip.text = element_text(size = ctrl$base_size),
              legend.text = element_text(size = 1.5 * ctrl$base_size),
              legend.title = element_text(size = 1.5 * ctrl$base_size),
              axis.text = element_text(size = ctrl$base_size),
              axis.title = element_text(size = 1.5 * ctrl$base_size)) +
        scale_x_continuous(expand = c(0, 1)) +
        labs(x = ctrl$ylab, y = ctrl$xlab) +
        guides(shape = guide_legend(override.aes = list(size = 2.7)))
    if (isDate & calendarTime) {
        xl <- ggplot_build(gg)$layout$panel_params[[1]]$x$breaks
        xl <- xl[!is.na(xl)]
        gg <- gg + scale_y_continuous(breaks = xl, labels = as.Date(xl, origin = "1970-01-01"))
    }
    gg
}

#' Produce Cumulative Sample Mean Function Plots
#'
#' Plot the mean cumulative function (MCF) for an \code{Recur} object.
#' The usage of the function is similar to that of \code{plot.Recur()} 
#' but with more flexible options.
#'
#' When \code{adjustRiskset = TRUE}, the \code{plotMCF()} is equivalent to
#' the Nelson-Aalen estimator for the intensity function of the recurrent event process.
#' When \code{adjustRiskset = FALSE}, the \code{plotMCF()} does not adjust for the risk set and
#' assumes all subjects remain at risk after the last observed recurrent event.
#' This is known as the survivor rate function.
#' The argument \code{control} consists of options with argument defaults
#' to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Cumulative mean".}
#'   \item{main}{customizable title, default value is "Sample cumulative mean function plot".}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can also be
#' specified outside of the \code{control} list.
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator,
#' and the predictors on the right.
#' The response must be a recurrent event survival object returned by the \code{Recur()} function.
#' @param data an optional data frame in which to interpret the variables occurring in
#' the "\code{formula}".
#' @param adjustRiskset an optional logical value that is passed to the \code{plotMCF()} function
#' as the \code{adjustRiskset} argument. See \code{\link{plotMCF}}.
#' This argument indicates whether risk set size will be adjusted. If \code{mcf.adjustRiskset = TRUE},
#' subjects leave the risk set after terminal times as in the Nelson-Aalen estimator.
#' If \code{mcf.adjustRiskset = FALSE}, subjects remain in the risk set after terminal time. 
#' @param smooth an optional logical value indicating whether to add a smooth curve
#' obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' This feature only works for data with one recurrent event type.
#' @param onePanel an optional logical value indicating whether the mean cumulative functions
#' will be plotted in the same panel.
#' This is only useful when there are multiple recurrent event types or
#' in the presence of (discrete) covariates.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab}, \code{main}, and more. See \bold{Details}.
#' 
#' @seealso \code{\link{Recur}}, \code{\link{plot.Recur}}
#' @keywords Plots
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom ggplot2 guides guide_legend
#' @importFrom scam scam
#'
#' @noRd
#' 
#' @example inst/examples/ex_plot_MCF.R
plotMCF <- function(formula, data, adjustRiskset = TRUE, onePanel = FALSE, 
                     smooth = FALSE, control = list(), ...) {
    call <- match.call()
    ctrl <- plotEvents.control()
    if (is.null(ctrl$ylab)) ctrl$ylab <- "Cumulative mean"
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
    nn <- unlist(aggregate(time2~., data = dd, table)$time2)
    dd <- unique(dd)
    dd$n <- as.integer(nn)
    rownames(dd) <- NULL
    k <- length(unique(dd$event)) - 1    
    if (!is.null(vNames)) { ## any covariates/stratifications?
        dd2 <- DF[DF$event == 0, vNames, drop = FALSE]
        dd2 <- dd2[do.call(order, dd2),,drop = FALSE]
        nn <- c(t(table(dd2)))
        dd2 <- unique(dd2)
        dd2$n <- as.integer(nn) ## tmp2 in the 1st version
        rownames(dd2) <- NULL
        dd$GrpInd <- match(apply(dd[,vNames, drop = FALSE], 1, paste, collapse = ""),
                           apply(dd2[,vNames, drop = FALSE], 1, paste, collapse = ""))
        tmp1 <- merge(dd, dd2, by = vNames)
        ## as.integer(eval(parse(text = paste(attr(terms(formula), "term.labels"), collapse = ":"))))
        rec0 <- tmp1[tmp1$event == 0,]
        dat0 <- do.call(rbind, lapply(split(tmp1, tmp1$GrpInd), function(x) {
            x$adjustRiskset = apply(x, 1, function(y)
                as.numeric(y['n.y']) - sum(rec0$n.x[as.numeric(y['time2'])> rec0$time2 &
                                                    rec0$GrpInd == as.numeric(y['GrpInd'])]))
            return(x)}))
        dat0$n.x <- dat0$n.x * (dat0$event > 0)
        rec0$time2 <- rec0$n.x <- 0
        rec0$adjustRiskset <- 1
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
        if (adjustRiskset) {
            dat0 <- do.call(rbind, 
                            lapply(split(dat0, list(dat0$event, dat0$GrpInd)), function(x) {
                                x$mu <- x$n.x / x$adjustRiskset
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
        dd$adjustRiskset <- apply(dd, 1, function(x) x[4] - sum(rec0$n[x[2] > rec0$time2]))
        dd$n <- dd$n * (dd$event > 0)
        rec0$time2 <- rec0$n <- 0
        rec0$adjustRiskset <- 1
        dat0 <- unique(rbind(dd, rec0))       
        dat0 <- dat0[order(dat0$time2),]
        dat0$mu <- dat0$n / (adjustRiskset * dat0$adjustRiskset + (!adjustRiskset) * dat0$n.y)
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
                message('The length of "recurrent.type" mismatched, default names are used.')
                dat0$event <- factor(dat0$event, labels = paste(ctrl$recurrent.name, 1:k))
            }
        }
    }
    if (nX > 0) {
        for (i in vNames) {
            dat0[,i] <- factor(dat0[,i], labels = paste(i, "=", unique(dat0[,i])))
        }}
    dat0 <- dat0[complete.cases(dat0),]
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
        gg <- gg + theme(legend.position = "none")
    ## if (onePanel & k == 1)
    ##     gg <- gg + scale_color_discrete(name = "", labels = levels(interaction(vNames)))
    ## if (onePanel & k > 1) gg <- gg + scale_color_discrete(name = "")
    if (smooth & k == 1 & !onePanel) {
        if (is.null(dat0$GrpInd)) dat0$GrpInd <- 1
        dat0 <- do.call(rbind, lapply(split(dat0, dat0$GrpInd), function(x){
            x$bs <- scam(x$MCF ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values
            return(x)}))       
        gg <- gg + geom_line(data = dat0, aes(time2, y = bs), color = 4, size = ctrl$lwd)
        ## geom_smooth(method = "scam", formula = y ~ s(x, k = 10, bs = "mpi"), size = ctrl$lwd, se = FALSE)
    }
    ## gg <- gg + geom_smooth(method = "loess", size = ctrl$lwd, se = FALSE)
    if (smooth & k > 1) message('Smoothing only works for data with one recurrent event type.')
    if (ctrl$main != "") gg <- gg + ggtitle(ctrl$main) 
    gg + theme(axis.line = element_line(color = "black"),
               legend.position = ctrl$legend.position,
               legend.key = element_rect(fill = "white", color = "white"),
               plot.title = element_text(size = 2 * ctrl$base_size),
               strip.text = element_text(size = ctrl$base_size),
               legend.text = element_text(size = 1.5 * ctrl$base_size),
               legend.title = element_text(size = 1.5 * ctrl$base_size),
               axis.text = element_text(size = ctrl$base_size),
               axis.title = element_text(size = 1.5 * ctrl$base_size)) +
        labs(y = ctrl$ylab, x = ctrl$xlab)
}

#' Plot the Baseline Cumulative Rate Function and the Baseline Cumulative Hazard Function
#'
#' Plot the baseline cumulative rate function and the baseline cumulative hazard function
#' (if applicable) for an \code{reReg} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list
#' with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty.}
#'   \item{main}{customizable title, default value are "Baseline cumulative rate and
#' hazard function" when \code{baseline = "both"},
#' "Baseline cumulative rate function" when \code{baseline = "rate"},
#' and "Baseline cumulative hazard function" when \code{baseline = "hazard"}.}
#' }
## #' These arguments can also be passed down without specifying a \code{control} list.
#' 
#' @param x an object of class \code{reReg}, returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether to add a smooth curve
#' obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param baseline a character string specifying which baseline function to plot.
#' \describe{
#'   \item{\code{baseline = "both"}}{plot both the baseline cumulative rate and
#' the baseline cumulative hazard function (if applicable) in separate panels
#' within the same display (default).}
#'   \item{\code{baseline = "rate"}}{plot the baseline cumulative rate function.}
#'   \item{\code{baseline = "hazard"}}{plot the baseline cumulative hazard function.}
#' }
## #' @param rateType a character string specifying the type of rate function to be plotted.
## #' This argument is passed to the \code{plotRate()} function as the \code{type} argument.
## #' See \code{\link{plotRate}}.
## #' Options are "unrestricted", "scaled", "bounded". 
#' @param newdata an optional data frame contains variables to include in the calculation
#' of the cumulative rate function.
#' If omitted, the baseline rate function will be plotted.
#' @param frailty an optional vector to specify the shared frailty for \code{newdata}.
#' If \code{newdata} is given and \code{frailty} is not specified, the
#' @param showName an optional logical value indicating whether to label the curves
#' when \code{newdata} is specified.
#' @param control a list of control parameters. See \bold{Details}.
#' @param ... additional graphical parameters to be passed to methods.
#' 
#' @seealso \code{\link{reReg}}
#' @export
#' @keywords Plots
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom ggplot2 geom_smooth geom_step ggplotGrob
#' @importFrom grid grid.draw
#' @example inst/examples/ex_plot_reReg.R
#' @method plot reReg
plot.reReg <- function(x,
                       baseline = c("both", "rate", "hazard"),
                       ## type = c("unrestricted", "bounded", "scaled"),
                       smooth = FALSE, newdata = NULL, frailty = NULL, showName = FALSE,
                       control = list(), ...) {
    baseline <- match.arg(baseline)
    ## type <- match.arg(type)
    type <- "unrestricted"
    if (x$typeRec %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline functions not yet available for this method.")
    if (x$typeRec %in% c("cox.LWYY", "cox.HH"))
        ctrl <- plotEvents.control(ylab = "Rate")
    if (baseline %in% c("both", "rate"))
        ctrl <- plotEvents.control(ylab = "Rate")
    if (baseline == "hazard")
        ctrl <- plotEvents.control(ylab = "Hazard")
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
        return(plotRate(x, smooth = smooth, type = type,
                        newdata = newdata, frailty = frailty, showName = showName, control = ctrl))
    if (baseline == "hazard")
        return(plotHaz(x, smooth = smooth,
                       newdata = newdata, frailty = frailty, showName = showName, control = ctrl))
    if (x$typeTem == ".") {
        ## cat(paste("Baseline cumulative hazard function is not available."))
        ## cat("\nOnly the baseline cumulative rate function is plotted.\n")
        return(plotRate(x, smooth = smooth, type = type,
                        newdata = newdata, frailty = frailty, showName = showName, control = ctrl))
    }
    g1 <- plotRate(x, smooth = smooth, type = type,
                   newdata = newdata, frailty = frailty, showName = showName, control = ctrl)
    g2 <- plotHaz(x, smooth = smooth, type = type,
                  newdata = newdata, frailty = frailty, showName = showName, control = ctrl)
    if (ctrl$main != "") g1 <- g1 + ggtitle(ctrl$main) 
    g1 <- g1 + ylab("Rate") + xlab("") +
        ## facet_grid(factor(rep(1, nrow(x$DF)), levels = 1, labels = "Baseline cumulative rate")) +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              plot.margin = unit(c(0, 0, 0, 0), "cm"))
    g2 <- g2 + ylab("Hazard") + xlab("Time") + 
        ## facet_grid(factor(rep(1, nrow(x$DF)), levels = 1, labels = "Baseline cumulative hazard")) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))        
    g1 <- ggplotGrob(g1)
    g2 <- ggplotGrob(g2)
    grid.draw(rbind(g1, g2))
    ## gg + ggtitle(ctrl$main)
}

#' Plotting the Baseline Cumulative Rate Function for the Recurrent Event Process
#'
#' Plot the baseline cumulative rate function for an \code{reReg} object.
#'
#' 
#' The \code{plotRate()} plots the estimated baseline cumulative rate function 
#' depending on the identifiability assumption.
#' When \code{type = "unrestricted"} (default), the baseline cumulative rate function
#' is plotted under the assumption \eqn{E(Z) = 1}.
#' When \code{type = "scaled"}, the baseline cumulative rate function is plotted
#' under the assumption \eqn{\Lambda(\min(Y^\ast, \tau)) = 1}.
#' When \code{type = "bounded"}, the baseline cumulative rate function is plotted
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
#' These arguments can also be specified outside of the \code{control} list.
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param type a character string specifying the type of rate function to be plotted.
#' Options are "unrestricted", "scaled", "bounded". See \bold{Details}.
#' @param smooth an optional logical value indicating whether to add a smooth curve
#' obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param newdata an optional data frame contains variables to include in the calculation
#' of the cumulative rate function.
#' If omitted, the baseline rate function will be plotted.
#' @param frailty an optional vector to specify the shared frailty for \code{newdata}.
#' If \code{newdata} is given and \code{frailty} is not specified, the
#' @param showName an optional logical value indicating whether to label the curves
#' when \code{newdata} is specified.
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
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom directlabels geom_dl
#' 
#' @example inst/examples/ex_plot_rate.R
plotRate <- function(x, newdata = NULL, frailty = NULL, showName = FALSE, 
                     type = c("unrestricted", "bounded", "scaled"),
                     smooth = FALSE, control = list(), ...) {
    if (x$typeRec %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline cumulative rate function is not available")
    if (is.null(frailty)) frailty <- exp(x$log.muZ)
    if (length(frailty) > 1 & !is.null(newdata) && length(frailty) != nrow(newdata))
        stop("newdata and frailty are different lengths")
    ## ctrl <- plot.reReg.control(main = "Baseline cumulative rate function")
    ctrl <- plotEvents.control(ylab = "Rate")
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
    if (x$typeRec == "nonparametric") {
        type <- "bounded"
        ctrl$ylab = "MCF Estimates"
    }
    if (!is.reReg(x)) stop("Response must be a `reReg` class")    
    dat <- x$DF[,"time2",drop = FALSE]
    if (is.null(newdata)) {    
        if (type == "unrestricted") dat$Y <- x$Lam0(dat$time2) * exp(x$log.muZ)
        if (type == "scaled") dat$Y <- x$Lam0(dat$time2) / x$Lam0(max(dat$time2))
        if (type == "bounded") dat$Y <- x$Lam0(dat$time2)
        if (!is.null(x$Lam0.upper)) {
            if (type == "bounded") {
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
            gg <- gg + geom_line(aes(x = time2, y = dat$bs), color = 4)
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
                    geom_step(aes(x = time2, y = Y.lower), lty = 2)
        }        
    }
    if (!is.null(newdata)) {
        if (!is.null(x$xlevels)) {
            for(i in which(names(newdata) %in% names(x$xlevels)))
                newdata[,i] <- factor(newdata[,i], levels = x$xlevels[[i]])
            newdata <- model.matrix(~., newdata)
        }
        X <- as.matrix(unique(newdata[,match(x$varNames, colnames(newdata)), drop = FALSE]))
        if (ncol(X) != length(x$varNames))
            stop(paste0("Variables ",
                        paste(setdiff(x$varNames, colnames(newdata)), collapse = ", "),
                        " are missing"))
        p <- ncol(X)
        exa1 <- exa2 <- 1
        if (x$typeRec == "cox") exa2 <- exp(X %*% x$par1)
        if (x$typeRec == "ar") exa1 <- exp(X %*% x$par1)
        if (x$typeRec == "am") exa1 <- exa2 <- exp(X %*% x$par1)
        if (x$typeRec == "sc") {
            exa1 <- exp(X %*% x$par1)
            exa2 <- exp(X %*% x$par2)
        }
        exa1 <- rep(drop(exa1), each = nrow(dat))
        exa2 <- rep(drop(exa2), each = nrow(dat))        
        Y <- frailty * x$Lam0(dat$time2 * exa1) * exa2 / exa1
        if (!is.null(x$Lam0.upper)) {
            Y.upper <- frailty * x$Lam0.upper(dat$time2 * exa1) * exa2 / exa1
            Y.lower <- frailty * x$Lam0.lower(dat$time2 * exa1) * exa2 / exa1
            dat <- data.frame(time2 = dat$time2,
                              id = rep(rownames(X), each = nrow(dat)),
                              Y = Y, Y.upper = Y.upper, Y.lower = Y.lower)
        } else
            dat <- data.frame(time2 = dat$time2, id = rep(rownames(X), each = nrow(dat)), Y = Y)
        gg <- ggplot(data = dat, aes(x = time2, y = Y, group = id)) +
            theme(axis.line = element_line(color = "black"))
        if (smooth) {
            dat$bs <- unlist(lapply(split(dat, dat$id), function(x)
                scam(x$Y ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
            gg <- gg + geom_line(aes(x = time2, y = dat$bs, group = id), color = 4)
            if (!is.null(x$Lam0.upper)) {
                dat$bs.upper <- unlist(lapply(split(dat, dat$id), function(x)
                    scam(x$Y.upper ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
                dat$bs.lower <- unlist(lapply(split(dat, dat$id), function(x)
                    scam(x$Y.lower ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
                gg <- gg +
                    geom_line(aes(x = time2, y = dat$bs.upper, group = id), color = 4, lty = 2) + 
                    geom_line(aes(x = time2, y = dat$bs.lower, group = id), color = 4, lty = 2)
            }
        } else {
            gg <- gg + geom_step()
            if (!is.null(x$Lam0.upper))
                gg <- gg +
                    geom_step(aes(x = time2, y = Y.upper, group = id), lty = 2) +
                    geom_step(aes(x = time2, y = Y.lower, group = id), lty = 2)
        }
        if (showName) 
            gg <- gg + geom_dl(aes(label = paste(" Obs. =", id)), method = "last.bumpup") +
                scale_x_continuous(limits = c(0, max(dat$time2) * 1.1))
    }
    if (ctrl$main != "") gg <- gg + ggtitle(ctrl$main) 
    gg <- gg + labs(x = ctrl$xlab, y = ctrl$ylab) +
        theme(plot.title = element_text(size = 2 * ctrl$base_size),
              strip.text = element_text(size = ctrl$base_size),
              legend.text = element_text(size = 1.5 * ctrl$base_size),
              legend.title = element_text(size = 1.5 * ctrl$base_size),
              axis.line = element_blank(),
              axis.text = element_text(size = ctrl$base_size),
              axis.title = element_text(size = 1.5 * ctrl$base_size))
    attr(gg, "from") <- "reReg"
    return(gg)
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
#'
#' @param x an object of class \code{reReg}, returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether to add a smooth curve
#' obtained from a monotone increasing P-splines implemented in package \code{scam}.
#' @param newdata an optional data frame contains variables to include in the calculation
#' of the cumulative rate function.
#' If omitted, the baseline rate function will be plotted.
#' @param frailty an optional vector to specify the shared frailty for \code{newdata}.
#' If \code{newdata} is given and \code{frailty} is not specified, the
#' @param showName an optional logical value indicating whether to label the curves
#' when \code{newdata} is specified.
#' @param type a character string specifying the type of rate function to be plotted.
#' Options are "unrestricted", "scaled", "bounded". See \bold{Details}.
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
plotHaz <- function(x, newdata = NULL, frailty = NULL, showName = FALSE,
                    type = c("unrestricted", "bounded", "scaled"),
                    smooth = FALSE, control = list(), ...) {
    if (x$typeRec %in% c("cox.GL", "cox.LWYY", "am.GL"))
        stop("Baseline cumulative hazard function is not available.")
    if (x$typeTem == ".") {
        stop("Baseline cumulative hazard function is not available.")
    }
    if (is.null(frailty)) frailty <- exp(x$log.muZ)
    if (length(frailty) > 1 & !is.null(newdata) && length(frailty) != nrow(newdata))
        stop("newdata and frailty are different lengths")
    ## ctrl <- plot.reReg.control(main = "Baseline cumulative hazard function")
    ctrl <- plotEvents.control(ylab = "Hazard")
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
    if (!is.reReg(x)) stop("Response must be a `reReg` class.")
    dat <- x$DF[, "time2", drop = FALSE]
    if (is.null(newdata)) {
        if (type == "unrestricted") dat$Y <- x$Haz0(dat$time2) * exp(x$log.muZ)
        if (type == "scaled") dat$Y <- x$Haz0(dat$time2) / x$Haz0(max(dat$time2))
        if (type == "bounded") dat$Y <- x$Haz0(dat$time2)
        ## dat$Y <- x$Haz0(dat$time2)
        if (!is.null(x$Haz0.upper)) {
            if (type == "bounded") {
                dat$Y.upper <- x$Haz0.upper(dat$time2)
                dat$Y.lower <- x$Haz0.lower(dat$time2)
            }
            if (type == "unrestricted") {
                dat$Y.upper <- x$Haz0.upper(dat$time2) * exp(x$log.muZ)
                dat$Y.lower <- x$Haz0.lower(dat$time2) * exp(x$log.muZ)
            }
            if (type == "scaled") {
                dat$Y.upper <- x$Haz0.upper(dat$time2) / x$Haz0(max(dat$time2))
                dat$Y.lower <- x$Haz0.lower(dat$time2) / x$Haz0(max(dat$time2))
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
                gg <- gg + geom_step(aes(x = time2,  y = Y.upper), lty = 2) +
                    geom_step(aes(x = time2,  y = Y.lower), lty = 2)
        }
    }
    if (!is.null(newdata)) { 
        X <- as.matrix(unique(newdata[,match(x$varNames, names(newdata)), drop = FALSE]))
        if (ncol(X) != length(x$varNames))
            stop(paste0("Variables ",
                        paste(setdiff(x$varNames, names(newdata)), collapse = ", "),
                        " are missing"))
        p <- ncol(X)
        exb1 <- exb2 <- 1
        if (x$typeTem == "cox") exb2 <- exp(X %*% x$par3)
        if (x$typeTem == "ar") exb1 <- exp(X %*% x$par3)
        if (x$typeTem == "am") exb1 <- exb2 <- exp(X %*% x$par3)
        if (x$typeTem == "sc") {
            exb1 <- exp(X %*% x$par3)
            exb2 <- exp(X %*% x$par4)
        }
        exb1 <- rep(drop(exb1), each = nrow(dat))
        exb2 <- rep(drop(exb2), each = nrow(dat))        
        Y <- frailty * x$Haz0(dat$time2 * exb1) * exb2 / exb1
        if (!is.null(x$Haz0.upper)) {
            Y.upper <- frailty * x$Haz0.upper(dat$time2 * exb1) * exb2 / exb1
            Y.lower <- frailty * x$Haz0.lower(dat$time2 * exb1) * exb2 / exb1
            dat <- data.frame(time2 = dat$time2, id = rep(rownames(X), each = nrow(dat)),
                              Y = Y, Y.upper = Y.upper, Y.lower = Y.lower)
        } else
            dat <- data.frame(time2 = dat$time2, id = rep(rownames(X), each = nrow(dat)), Y = Y)
        gg <- ggplot(data = dat, aes(x = time2, y = Y, group = id)) +
            theme(axis.line = element_line(color = "black"))
        if (smooth) {
            dat$bs <- unlist(lapply(split(dat, dat$id), function(x)
                scam(x$Y ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
            gg <- gg + geom_line(aes(x = time2, y = dat$bs, group = id), color = 4)
            if (!is.null(x$Haz0.upper)) {
                dat$bs.upper <- unlist(lapply(split(dat, dat$id), function(x)
                    scam(x$Y.upper ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
                dat$bs.lower <- unlist(lapply(split(dat, dat$id), function(x)
                    scam(x$Y.lower ~ s(x$time2, k = 10, bs = "mpi"))$fitted.values))
                gg <- gg +
                    geom_line(aes(x = time2, y = dat$bs.upper, group = id), color = 4, lty = 2) + 
                    geom_line(aes(x = time2, y = dat$bs.lower, group = id), color = 4, lty = 2)
            }
        } else {
            gg <- gg + geom_step()
            if (!is.null(x$Haz0.upper))
                gg <- gg +
                    geom_step(aes(x = time2, y = Y.upper, group = id), lty = 2) +
                    geom_step(aes(x = time2, y = Y.lower, group = id), lty = 2)
        }
        if (showName) 
            gg <- gg + geom_dl(aes(label = paste(" Obs. =", id)), method = "last.bumpup") +
                scale_x_continuous(limits = c(0, max(dat$time2) * 1.1))
    }
    if (ctrl$main != "") gg <- gg + ggtitle(ctrl$main) 
    gg <- gg + labs(x = ctrl$xlab, y = ctrl$ylab) +
        theme(plot.title = element_text(size = 2 * ctrl$base_size),
              strip.text = element_text(size = ctrl$base_size),
              legend.text = element_text(size = 1.5 * ctrl$base_size),
              legend.title = element_text(size = 1.5 * ctrl$base_size),
              axis.line = element_blank(),
              axis.text = element_text(size = ctrl$base_size),
              axis.title = element_text(size = 1.5 * ctrl$base_size))
    attr(gg, "from") <- "reReg"
    return(gg)
}

#' Plot options for plotEvents
#'
#' This function provides the plotting options for the \code{plotEvents()} function.
#'
#' @param xlab a character string indicating the label for the x axis.
#' The default value is "Time".
#' @param ylab a character string indicating the label for the y axis.
#' The default value is "Subject".
#' @param main a character string indicating the title of the plot.
#' @param terminal.name a character string indicating the label for the terminal event
#' displayed in the legend. The default value is "Terminal event".
#' @param recurrent.name a character string indicating the label for the recurrent event
#' displayed in the legend. The default value is "Recurrent events".
#' @param recurrent.type a factor indicating the labels for the different recurrent event types.
#' This option is only available when there are more than one types of recurrent events.
#' The default value is "Recurrent events 1", "Recurrent events 2", ....
#' @param legend.position a character string specifies the position of the legend.
#' The available options are "none", "left", "right", "bottom", "top",
#' or a two-element numeric vector specifies the coordinate of the legend.
#' This argument is passed to the \code{ggplot} theme environment.
#' The default value is "top".
#' @param base_size a numerical value to specify the base font size, given in pts.
#' This argument is passed to the \code{ggplot} theme environment.
#' The default value is 12.
#' @param cex a numerical value specifies the size of the points. 
#' @param alpha a numerical value specifies the transparency of the points. 
#' 
#' @seealso \code{\link{plotEvents}}
#' @export
plotEvents.control <- function(xlab = NULL, ylab = NULL,
                               main = NULL, 
                               terminal.name = NULL,
                               recurrent.name = NULL, 
                               recurrent.type = NULL, 
                               legend.position = "top", base_size = 12,
                               cex = NULL, alpha = .7) {
    if (is.null(ylab)) ylab <- "Subject"
    if (is.null(xlab)) xlab <- "Time"
    if (is.null(main)) main <- ""
        ## main <- "Recurrent event plot"
    if (is.null(terminal.name)) terminal.name <-  "Terminal event"
    if (is.null(recurrent.name)) recurrent.name <- "Recurrent events"
    list(xlab = xlab, ylab = ylab, main = main, cex = cex,
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type, alpha = alpha,
         legend.position = legend.position, base_size = base_size)
}

#' Function used to combine baseline functions in one plot
#'
#' Combine different plots into one.
#'
#' @param ... \code{ggplot} objects created by plotting \code{reReg} objects.
#' @param legend.title an optional character string to specify the legend title.
#' @param legend.labels an optional character string to specify the legend labels.
#' @param control a list of control parameters. 
#' 
#' @export
#' @keywords Plots
#' 
#' @example inst/examples/ex_basebind.R
basebind <- function(..., legend.title, legend.labels, control = list()) {
    gglst <- list(...)
    if (any(sapply(gglst, function(x) attr(x, "from")) != "reReg"))
        stop("Plots must be created from reReg objects")
    if (missing(legend.title)) legend.title <- ""
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (is.null(ctrl$ylab)) ctrl$ylab <- ""
    nargs <- length(gglst)
    if (missing(legend.labels)) legend.labels <- 1:nargs
    if (length(legend.labels) != nargs) {
        message('The length of "name" mismatched, default names are used.')
        legend.labels <- 1:nargs
    }
    d <- do.call(rbind, lapply(gglst, function(x) x$data))
    nobs <- sapply(gglst, function(x) nrow(x$data))
    d$group <- factor(rep(1:nargs, nobs), labels = legend.labels)
    gg <- ggplot(data = d, aes(x = time2, y = Y, color = group)) + geom_step()
    if (!is.null(d$Y.upper)) 
        gg <- gg + geom_step(aes(x = time2, y = Y.upper), lty = 2) +
            geom_step(aes(x = time2, y = Y.lower), lty = 2)
    gg + labs(x = gglst[[1]]$labels$x, y = gglst[[1]]$labels$y, color = legend.title) +
        theme(plot.title = element_text(size = 2 * ctrl$base_size),
              strip.text = element_text(size = ctrl$base_size),
              legend.position = ctrl$legend.position,
              legend.text = element_text(size = 1.5 * ctrl$base_size),
              legend.title = element_text(size = 1.5 * ctrl$base_size),
              axis.line = element_blank(),
              axis.text = element_text(size = ctrl$base_size),
              axis.title = element_text(size = 1.5 * ctrl$base_size))
}
