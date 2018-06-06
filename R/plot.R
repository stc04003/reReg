globalVariables(c("Time", "Yi", "id", "recType", "status", "tij"))

#' Produce Event Plots
#'
#' Plot the event plot for an \code{reSurv} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{title}{customizable title, default value is "Recurrent event plot".}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable label for recurrent event, default value is "Recurrent event plot".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is NULL.}
#'   \item{alpha}{controls the transparency of points.}
#' }
#' 
#' @param x an object of class \code{reSurv}, usually returned by the \code{reSurv} function.
#' @param data an optional data frame in which to interpret the variables occurring in the "formula".
#' @param order an optional logical value. If "TRUE", the event plot is sorted by the terminal times; the default value is TRUE.
#' @param return.grob an optional logical value. If "TRUE", a \code{ggplot2} plot grob will be returned.
#' @param onePanel an optinoal logical value. If "TRUE", only one graphical panel will be displayed.
#' @param control a list of control parameters.
#' @param CMF an optional logical value. If "TRUE", the cumulative sample mean function will be plotted.
#' If "FALSE", the event plot will be plotted. The default is \code{CMF = FALSE}.
#' @param ... for future developments.
#' @seealso \code{\link{reSurv}}
#' @keywords plot.reSurv
#' @export
#' 
#' @examples
#' data(readmission)
#' reObj <- with(subset(readmission, id <= 10), reSurv(t.stop, event, death, id))
#' ## Default labels
#' plot(reObj)
#' plot(reObj, order = FALSE)
#' ## User specified labels
#' plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))
#'
#' ## With multiple hypothetical event types
#' set.seed(1)
#' reObj2 <- with(readmission, reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id))
#' plot(reObj2)
#'
#' ## With multiple hypothetical event types
#' set.seed(1)
#' reObj2 <- with(readmission, reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id))
#' plot(reObj2)
plot.reSurv <- function(x, data, CMF = FALSE, order = TRUE, onePanel = FALSE, return.grob = FALSE, control = list(), ...) {
    if (!CMF)
        return(plotEvents(x, data, order = order, return.grob = return.grob, control = control))
    if (CMF)
        return(plotCMF(x, data, onePanel = onePanel, return.grob = return.grob, control = control))
}

#' Produce Event Plots
#'
#' Plot the event plot for an \code{reSurv} object, with more flexible options.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{title}{customizable title, default value is "Recurrent event plot".}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable label for recurrent event, default value is "Recurrent event plot".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is NULL.}
#'   \item{alpha}{controls the transparency of points.}
#' }
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator, and the terms on the right.
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data an optional data frame in which to interpret the variables occurring in the "formula".
#' @param order an optional logical value. If "TRUE", the event plot is sorted by the terminal times; the default value is TRUE.
#' @param return.grob an optional logical value. If "TRUE", a \code{ggplot2} plot grob will be returned.
#' @param control a list of control parameters.
#' @param ... for future developments.
#' @seealso \code{\link{reSurv}}, \code{\link{plot.reSurv}}
#' @keywords plot.reSurv
#' @export
#' 
#' @examples
#' data(readmission)
#' plotEvents(reSurv(t.stop, event, death, id) ~ 1, data = readmission)
#'
#' ## Separate plots by gender
#' plotEvents(reSurv(t.stop, event, death, id) ~ sex, data = readmission)
#'
#' ## Separate plots by gender and chemo type
#' plotEvents(reSurv(t.stop, event, death, id) ~ sex + chemo, data = readmission)
#'
#' ## With multiple hypothetical event types
#' plotEvents(reSurv(t.stop, event * sample(1:3, 861, TRUE), death, id) ~
#'   sex + chemo, data = readmission)
plotEvents <- function(formula, data, order = TRUE, return.grob = FALSE, control = list(), ...) {
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    call <- match.call()
    nX <- 0
    if (is.reSurv(formula)) {dat <- formula$reTb
    } else {
        ## if (length(attr(terms(formula), "term.labels")) > 1)
        ##     stop("The current vision can only handle one covaraite.")
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        dat <- obj$reTb
        nX <- length(formula[[3]])
        if (formula[[3]] != 1 && nX == 1) {
            if (missing(data)) DF <- cbind(obj$reDF, eval(formula[[3]], parent.frame()))
            if (!missing(data)) DF <- cbind(obj$reDF, eval(formula[[3]], data))
            colnames(DF) <- c(names(obj$reDF), paste0(formula[[3]], collapse = ""))
            suppressMessages(dat <- left_join(obj$reTb, unique(select(DF, id, paste(formula[[3]])))))
        }
        if (formula[[3]] != 1 && nX > 1) {
            DF <- obj$reDF
            if (missing(data)) {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], parent.frame()))
                }
            } else {
                for (i in 2:nX) {
                    DF <- cbind(DF, eval(formula[[3]][[i]], data))
                }
            }
            colnames(DF) <- c(names(obj$reDF), sapply(2:nX, function(x) paste0(formula[[3]][[x]], collapse = "")))
            suppressMessages(dat <- left_join(obj$reTb, unique(select(DF, id, paste(formula[[3]][-1])))))
        }
    }
    if (order) {
        if (nX == 0 || formula[[3]] == 1) dat <- dat %>% mutate(id = rank(Yi, ties.method = "first"))
        else dat <- dat %>% group_by_at(6:ncol(dat)) %>% mutate(id = rank(Yi, ties.method = "first")) 
    }
    sz <- 1 + 8 / (1 + exp(length(unique(dat$id)) / 30))
    k <- length(unique(unlist(dat$recType)))
    shp.val <- c(17, rep(19, k))
    clr.val <- c(alpha("red", ctrl$alpha), hcl(h = seq(120, 360, length.out = k), l = 60, alpha = ctrl$alpha))
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
    names(shp.val) <- names(clr.val) <- c("Yi", rec.lab)
    ## ggplot starts here
    gg <- ggplot(dat, aes(id, Yi)) +
        geom_bar(stat = "identity", fill = "gray75") +
        geom_point(data = dat %>% filter(status > 0),
                   aes(id, Yi, shape = "Yi", color = "Yi"), size = sz) +
        geom_point(data = dat %>% filter(!map_lgl(tij, is.null)) %>%
                       unnest(tij, recType), # %>% select(id, tij, recType),
                   aes(id, tij, shape = factor(recType, labels = rec.lab), 
                       color = factor(recType, labels = rec.lab)),
                   size = sz) + 
        ## position = position_jitter(w = 0.1, h = 0)) +
        scale_shape_manual(
            name = "", values = shp.val,
            labels = shp.lab,
            breaks = c("Yi", rec.lab)) +
        scale_color_manual(
            name = "", values = clr.val,
            labels = shp.lab,
            breaks = c("Yi", rec.lab)) +
        coord_flip() + 
        theme(axis.line.y = element_blank(),
              axis.title.y = element_text(vjust = 0),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "both")
    if (!return.grob) {
        gg + theme(panel.background = element_blank(),
                   axis.line = element_line(colour = "black"),
                   legend.key = element_rect(fill = "white", color = "white")) +
            scale_x_continuous(expand = c(0, 1)) +
            ggtitle(ctrl$title) + labs(x = ctrl$ylab, y = ctrl$xlab) +
            guides(shape = guide_legend(override.aes = list(size = 2.7)))
    } else {return(ggplotGrob(gg))}
}

#' Produce Cumulative Sample Mean Function Plot
#'
#' Plot the cumulative sample mean function for an \code{reSurv} object, with more flexible options.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{title}{customizable title, default value is "Recurrent event plot".}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable label for recurrent event, default value is "Recurrent event plot".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is NULL.}
#'   \item{alpha}{controls the transparency of points.}
#' }
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator, and the terms on the right.
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data an optional data frame in which to interpret the variables occurring in the "formula".
#' @param return.grob an optional logical value. If "TRUE", a \code{ggplot2} plot grob will be returned.
#' @param adjrisk an optional logical value. If "TRUE", event times with \code{event == 0} will be treated as
#' (independent) censoring time, and risk set sizes are adjusted accordingly. 
#' @param onePanel an optinoal logical value. If "TRUE", only one graphical panel will be displayed.
#' @param control a list of control parameters.
#' @param ... for future developments.
#' @seealso \code{\link{reSurv}}, \code{\link{plot.reSurv}}
#' @keywords plot.reSurv
#' @export
#'
#' @importFrom dplyr summarise rowwise
#' 
#' @examples
#' data(readmission)
#' plotCMF(reSurv(t.stop, event, death, id) ~ 1, data = readmission)
plotCMF <- function(formula, data, onePanel = FALSE, return.grob = FALSE, adjrisk = TRUE, control = list(), ...) {
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (!("title" %in% namc)) ctrl$title = "Sample cumulative mean function plot"
    if (!("xlab" %in% namc)) ctrl$xlab = "Time"
    if (!("ylab" %in% namc)) ctrl$ylab = "Cumulative mean"
    call <- match.call()
    nX <- 0
    if(is.reSurv(formula)) {
        dat1 <- formula$reDF
    } else {
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
        dat1 <- obj$reDF
        nX <- length(formula[[3]])
        if (formula[[3]] != 1 && nX == 1) {
            if (missing(data)) dat1 <- bind_cols(dat1, tmp = eval(formula[[3]], parent.frame()))
            if (!missing(data)) dat1 <- bind_cols(dat1, tmp = eval(formula[[3]], data))
            names(dat1) <- c(names(dat1)[1:5], paste0(formula[[3]], collapse = ""))
        }
        if (formula[[3]] != 1 && nX > 1) {
            if (missing(data)) {
                for (i in 2:nX) {
                    dat1 <- cbind(dat1, eval(formula[[3]][[i]], parent.frame()))
                }
            } else {
                for (i in 2:nX) {
                    dat1 <- cbind(dat1, eval(formula[[3]][[i]], data))
                }
            }
            names(dat1) <- c(names(dat1)[1:5], sapply(2:nX, function(x) paste0(formula[[3]][[x]], collapse = "")))
        }
    }
    rText1 <- paste("dat1 %>% count(", paste(names(dat1)[5:ncol(dat1)], collapse = ","),", Time)")
    tmp1 <- eval(parse(text = rText1))
    if (ncol(dat1) > 5) { ## any covariates?
        rText2 <- paste("dat1 %>% group_by(", paste(names(dat1)[6:ncol(dat1)], collapse = ","),")",
                        "%>% summarise(n = length(unique(id)))")
        tmp2 <- eval(parse(text = rText2))
        tmp1 <- left_join(tmp1, tmp2, by = paste(names(dat1)[6:ncol(dat1)])) %>%
            mutate(GrpInd = as.integer(eval(parse(text = paste(attr(terms(formula), "term.labels"), collapse = ":")))))
        rec0 <- tmp1 %>% filter(recType == 0)
        dat0 <- tmp1 %>% filter(recType > 0) %>% rowwise() %>%
            mutate(adjrisk = n.y - sum(rec0$n.x[Time > rec0$Time & rec0$GrpInd == GrpInd])) %>% unrowwise
        if (adjrisk) {
            dat0 <- eval(parse(text = paste("dat0 %>% group_by(", paste(names(dat1)[5:ncol(dat1)], collapse = ","), ")"))) %>%
                mutate(mu = n.x / adjrisk, CMF = cumsum(mu))
        } else {
            dat0 <- eval(parse(text = paste("dat0 %>% group_by(", paste(names(dat1)[5:ncol(dat1)], collapse = ","), ")"))) %>%
                mutate(mu = n.x / n.y, CMF = cumsum(mu))
        }
        ## dat0 <- eval(parse(text = paste("dat0 %>% group_by(", paste(names(dat1)[5:ncol(dat1)], collapse = ","), ")")))
        ## dat0 <- dat0 %>% filter(recType > 0) %>% mutate(CMF = cumsum(mu))
    } else {
        rec0 <- tmp1 %>% filter(recType == 0)
        dat0 <- tmp1 %>% filter(recType > 0) %>% rowwise() %>%
            mutate(n.y = length(unique(dat1$id)), adjrisk = n.y - sum(rec0$n[Time > rec0$Time]))  %>%
            unrowwise
        if (adjrisk) {
            dat0 <- dat0 %>% mutate(mu = n / adjrisk, CMF = cumsum(mu))
        } else {
            dat0 <- dat0 %>% mutate(mu = n / n.y, CMF = cumsum(mu))
        }
    }
    k <- length(unique(unlist(dat0$recType)))
    if (k == 1) dat0$recType <- factor(dat0$recType, label = ctrl$recurrent.name)
    if (k > 1 & is.null(ctrl$recurrent.type))
        dat0$recType <- factor(dat0$recType, label = paste(ctrl$recurrent.name, 1:k))
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            dat0$recType <- factor(dat0$recType, label = ctrl$recurrent.type)
        } else {
            cat('The length of "recurrent.type" mismatched, default names are used.\n')
            dat0$recType <- factor(recType, label = paste(ctrl$recurrent.name, 1:k))
        }
    }
    if (!onePanel) {
    gg <- ggplot(data = dat0, aes(x = Time, y = CMF)) +
        geom_step(aes(color = recType), direction = "hv")
    } else {
        rText <- paste("geom_step(aes(color = interaction(",
                       paste(names(dat0)[5:ncol(dat1) - 4], collapse = ","),
                       ")), direction = \"hv\")")
        gg <- ggplot(data = dat0, aes(x = Time, y = CMF)) +
            eval(parse(text = rText))
    }
    if (!onePanel && nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "x")
    if (!onePanel & k == 1)
        gg <- gg + theme(legend.position="none")
    if (onePanel & k == 1)
        gg <- gg + scale_color_discrete(name = "",
                                        labels = levels(interaction(dat0[,(5 + 1):ncol(dat1) - 4])))
    if (onePanel & k > 1)
        gg <- gg + scale_color_discrete(name = "")
    if (!return.grob) {
        gg +  theme(axis.line = element_line(colour = "black"),
                    legend.key = element_rect(fill = "white", color = "white")) +
            ## scale_color_discrete(name = "") +
            ## scale_y_continuous(expand = c(0, 1)) +
            ggtitle(ctrl$title) + labs(y = ctrl$ylab, x = ctrl$xlab) 
    } else {
        return(ggplotGrob(gg))
    }
}

plotEvents.control <- function(xlab = "Time", ylab = "Subject", title = "Recurrent event plot",
                               terminal.name = "Terminal event",
                               recurrent.name = "Recurrent event",
                               recurrent.type = NULL, alpha = .7) {
    list(xlab = xlab, ylab = ylab, title = title, 
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type, alpha = alpha)
}

#' Plotting Estimated Baseline Cumulative Rate Functions
#'
#' Plot the estimated baseline cumulative rate function and harzard function if available.
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param ... for future methods.
#'
#' @seealso \code{\link{reReg}}
#' @export
#' @keywords plot.reReg
#' @examples
#' data(readmission)
#' fit <- reReg(reSurv(t.stop, event, death, id) ~ sex + chemo,
#'              data = subset(readmission, id < 50),
#'              method = "am.XCHWY", se = "resampling", B = 20)
#' plot(fit)
plot.reReg <- function(x, ...) {
    options(warn = -1)
    if (!is.reReg(x)) stop("Response must be a reReg class")
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH") 
    par(mfrow = c(2, 1))
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.ly <- max(ly, lyU, lyL, na.rm = TRUE)
    win.hy <- max(hy, hyU, hyL, na.rm = TRUE)
    plot(t0, ly, type = "s",  xlab = "", ylab = "", ylim = c(0, win.ly),
         main = "Baseline Cumulative Rate Function")
    if (any(!is.na(x$lamU))) {
        lines(t0, lyU, col = 2)
        lines(t0, lyL, col = 2)
    }
    title(ylab = expression(hat(Lambda)[0](t)), xlab = "Time", line = 2.2)
    ##     title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.ly, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    ## op <- par(ask=TRUE)
    ## plot(t, hy, type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, win.hy),
    ##      main = "Baseline Cumulative Hazard Function")
    plot(t0, hy, type = "s",  xlab = "", ylab = "", ylim = c(0, win.hy),
         main = "Baseline Cumulative Hazard Function")
    if (any(!is.na(x$hazU))) {
        lines(t0, hyU, col = 2, lty = 2)
        lines(t0, hyL, col = 2, lty = 2)
    }
    title(ylab = expression(hat(H)[0](t)), xlab = "Time", line = 2.2)
    ## title(xlab = "time", line = 1.6, cex.lab = 0.8)
    ## axis(1, at = seq(0, round(max(t), 1), length.out = 11), 
    ##      cex.axis = 0.8, tck = -0.015, mgp = c(3, .3, 0))
    ## axis(2, at = seq(0, round(win.hy, 2), length.out = 11),
    ##      las = 2, cex.axis = 0.8, tck = -0.015, mgp = c(3, .4, 0))
    ## par(op)
    par(mfrow = c(1, 1))
    out <- c(x, list(t0 = t0, lam = ly, lamU = lyU, lamL = lyL, haz = hy, hazU = hyU, hazL = hyL))
    options(warn = 0)
    invisible(out)
}

#' Plotting the baseline rate function
#'
#' Plot the baseline rate function after fitting \code{reReg}.
#' 
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{title}{customizable title, default value is "Recurrent event plot".}
#' }
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param control a list of control parameters.
#' @param ... for future developments.
#'
#' @seealso \code{\link{reReg}}
#' @export
#'
#' @examples
#' ## readmission data
#' data(readmission)
#' set.seed(123)
#' fit <- reReg(reSurv(t.stop, event, death, id) ~ sex + chemo,
#'                    data = subset(readmission, id < 50),
#'                    method = "am.XCHWY", se = "resampling", B = 20)
#' ## Plot both the baseline cumulative rate and hazard function
#' plot(fit)
#' ## Plot baseline cumulative rate function
#' plotRate(fit)
#' ## Plot with user-specified labels
#' plotRate(fit, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title")) 
plotRate <- function(x, control = list(), ...) {
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH") 
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$title == "Recurrent event plot") # default value
        ctrl$title <- "Baseline Cumulative Rate Function"
    if (!is.reReg(x)) stop("Response must be a reReg class")
    options(warn = -1)
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.ly <- max(ly, lyU, lyL, na.rm = TRUE)
    plot(t0, ly, type = "s",  xlab = "", ylab = "", ylim = c(0, win.ly),
         main = ctrl$title)
    if (any(!is.na(x$lamU))) {
        lines(t0, lyU, col = 2, lty = 2, "s")
        lines(t0, lyL, col = 2, lty = 2, "s")
    }
    if (ctrl$ylab == "Subject") #default value 
        title(ylab = expression(hat(Lambda)[0](t)), xlab = ctrl$xlab, line = 2.2)
    else title(ylab = ctrl$ylab, xlab = ctrl$xlab, line = 2.2)
    options(warn = 0)
}

#' Plotting the baseline hazard function
#'
#' Plot the baseline hazard function after fitting \code{reReg}.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject".}
#'   \item{title}{customizable title, default value is "Recurrent event plot".}
#' }
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param control a list of control parameters.
#' @param ... for future developments.
#' 
#' @seealso \code{\link{reReg}}
#' @export
#'
#' @examples
#' ## readmission data
#' data(readmission)
#' set.seed(123)
#' fit <- reReg(reSurv(t.stop, event, death, id) ~ sex + chemo,
#'              data = subset(readmission, id < 50),
#'              method = "am.XCHWY", se = "resampling", B = 20)
#' ## Plot both the baseline cumulative rate and hazard function
#' plot(fit)
#' ## Plot baseline cumulative hazard function
#' plotHaz(fit)
#' ## Plot with user-specified labels
#' plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", title = "User title"))  
plotHaz <- function(x, control = list(), ...) {
    if (x$method == "sc.XCYH") stop("Rate function plot is not yet available for method = sc.XCYH")
    ctrl <- plotEvents.control()
    namc <- names(control)
    if (!all(namc %in% names(ctrl))) 
        stop("unknown names in control: ", namc[!(namc %in% names(ctrl))])
    ctrl[namc] <- control
    if (ctrl$title == "Recurrent event plot") # default value
        ctrl$title <- "Baseline Cumulative Hazard Function"
    if (!is.reReg(x)) stop("Response must be a reReg class")
    options(warn = -1)
    t0 <- x$t0
    ly <- x$lam
    lyU <- x$lamU
    lyL <- x$lamL
    hy <- x$haz
    hyU <- x$hazU
    hyL <- x$hazL
    win.hy <- max(hy, hyU, hyL, na.rm = TRUE)
    plot(t0, hy, type = "s",  xlab = "", ylab = "", ylim = c(0, win.hy),
         main = ctrl$title)
    if (any(!is.na(x$hazU))) {
        lines(t0, hyU, col = 2, lty = 2, "s")
        lines(t0, hyL, col = 2, lty = 2, "s")
    }
    if (ctrl$ylab == "Subject") #default value 
        title(ylab = expression(hat(Lambda)[0](t)), xlab = ctrl$xlab, line = 2.2)
    else title(ylab = ctrl$ylab, xlab = ctrl$xlab, line = 2.2)
    options(warn = 0)
}

unrowwise <- function(x) {
  class(x) <- c( "tbl_df", "data.frame")
  x
}
