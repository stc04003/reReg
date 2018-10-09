globalVariables(c("Time", "Yi", "id", "recType", "status", "tij", "n.y", "GrpInd", "n.x", "mu", "n", "CSM"))
globalVariables(c("Y", "Y.upper", "Y.lower", "group")) ## global variables for plot.reReg

#' Produce Event Plot or Cumulative Sample Mean Function Plot
#'
#' Plot the event plot or the cumulative sample mean (CSM) function for an \code{reSurv} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is "Subject" for event plot and "Cumulative mean" for CSM plot.}
#'   \item{main}{customizable title, the default value is "Recurrent event plot" when \code{CSM = FALSE} and
#' "Sample cumulative mean function plot" when \code{CSM = TRUE}.}
#'   \item{terminal.name}{customizable label for terminal event, default value is "Terminal event".}
#'   \item{recurrent.name}{customizable legend title for recurrent event, default value is "Recurrent events".}
#'   \item{recurrent.types}{customizable label for recurrent event type, default value is \code{NULL}.}
#'   \item{alpha}{between 0 and 1, controls the transparency of points.}
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can also be passed down without specifying a \code{control} list.
#' 
#' @param x an object of class \code{reSurv}, usually returned by the \code{reSurv} function.
## #' @param data an optional data frame in which to interpret the variables occurring in the "formula".
#' @param order an optional logical value indicating whether the event plot (when \code{CSM = FALSE})
#' will be sorted by the terminal times.
## #' @param return.grob an optional logical value indicating whether a \code{ggplot2} plot grob will be returned.
#' @param control a list of control parameters. See \bold{Details}.
#' @param CSM an optional logical value indicating whether the cumulative sample mean (CSM) function will
#' be plotted instead of the event plot (default).
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#' 
#' @seealso \code{\link{reSurv}}
#' 
#' @keywords Plots
#' @export
#'
#' @return A \code{ggplot} object.
#' @examples
#' data(readmission, package = "frailtypack")
#' reObj <- with(subset(readmission, id <= 10), reSurv(t.stop, id, event, death))
#'
#' ## Event plots:
#' ## Default labels
#' plot(reObj)
#' plot(reObj, order = FALSE)
#' ## User specified labels
#' plot(reObj, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))
#'
#' ## With multiple hypothetical event types
#' set.seed(1)
#' reObj2 <- with(readmission, reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death))
#' plot(reObj2)
#'
#' ## CSM plots
#' plot(reObj, CSM = TRUE)
plot.reSurv <- function(x, CSM = FALSE, order = TRUE, control = list(), ...) {
    if (!is.reSurv(x)) stop("Response must be a reSurv class")
    if (!CSM) ctrl <- plotEvents.control()
    if (CSM) ctrl <- plotCSM.control()
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
    if (!CSM) {
        return(plotEvents(x, order = order, control = ctrl))
    }
    if (CSM)
        return(plotCSM(x, onePanel = TRUE, control = ctrl))
}

#' Produce Event Plots
#'
#' Plot the event plot for an \code{reSurv} object.
#' The function is similar to \code{plot.reSurve} but with more flexible options.
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
#' }
#' The \code{xlab}, \code{ylab} and \code{main} parameters can also be passed down without specifying a \code{control} list. See \bold{Examples}.
#' 
#' @param formula  a formula object, with the response on the left of a "~" operator,
#' and the predictors on the right.
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data an optional data frame in which to interpret the variables occurring in the "\code{formula}".
#' @param order an optional logical value indicating whether the event plot will be sorted by the terminal times.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#'
#' @seealso \code{\link{reSurv}}, \code{\link{plot.reSurv}}
#' 
#' @keywords Plots
#' @export
#' 
#' @return A \code{ggplot} object.
#' 
#' @examples
#' data(readmission, package = "frailtypack")
#' plotEvents(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
#'
#' ## Separate plots by gender
#' plotEvents(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
#'
#' ## Separate plots by gender and chemo type
#' plotEvents(reSurv(t.stop, id, event, death) ~ sex + chemo, data = readmission)
#'
#' ## With multiple hypothetical event types
#' plotEvents(reSurv(t.stop, id, event * sample(1:3, 861, TRUE), death) ~
#'   sex + chemo, data = readmission)
plotEvents <- function(formula, data, order = TRUE, control = list(), ...) {
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
    if (is.reSurv(formula)) {dat <- formula$reTb
    } else {
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
    if (ctrl$cex == "Default") sz <- 1 + 8 / (1 + exp(length(unique(dat$id)) / 30))
    else sz <- ctrl$cex
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
    gg <- ggplot(dat, aes(id, Yi)) +
        geom_bar(stat = "identity", fill = "gray75") +
        coord_flip() + 
        theme(axis.line.y = element_blank(),
              axis.title.y = element_text(vjust = 0),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
    if (sum(!is.na(dat$tij)) > 0) 
        gg <- gg + geom_point(data = dat %>% filter(!map_lgl(tij, is.null)) %>%
                                  unnest(tij, recType), # %>% select(id, tij, recType),
                              aes(id, tij, shape = factor(recType, labels = rec.lab), 
                                  color = factor(recType, labels = rec.lab)), size = sz)  
            ## ## position = position_jitter(w = 0.1, h = 0)) +
            ## scale_shape_manual(name = "", values = shp.val,
            ##     labels = shp.lab, breaks = c("Yi", rec.lab)) +
            ## scale_color_manual(name = "", values = clr.val,
            ##     labels = shp.lab, breaks = c("Yi", rec.lab))
    if (sum(dat$status) > 0)
        gg <- gg + geom_point(data = dat %>% filter(status > 0),
                              aes(id, Yi, shape = "Yi", color = "Yi"), size = sz) 
    if (nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "both")
    gg <- gg + scale_shape_manual(name = "", values = shp.val,
                                  labels = shp.lab, breaks = c("Yi", rec.lab)) +
        scale_color_manual(name = "", values = clr.val,
                           labels = shp.lab, breaks = c("Yi", rec.lab))
    gg + theme(panel.background = element_blank(),
               axis.line = element_line(color = "black"),
               legend.key = element_rect(fill = "white", color = "white")) +
        scale_x_continuous(expand = c(0, 1)) +
        ggtitle(ctrl$main) + labs(x = ctrl$ylab, y = ctrl$xlab) +
        guides(shape = guide_legend(override.aes = list(size = 2.7)))
}

#' Produce Cumulative Sample Mean Function Plots
#'
#' Plot the cumulative sample mean function (CSM) for an \code{reSurv} object.
#' The function is similar to \code{plot.reSurv} but with more flexible options.
#'
#' When \code{adjrisk = TRUE}, the \code{plotCSM} is equivalent to
#' the Nelson-Aalen estimator for the intensity function of the recurrent event process.
#' When \code{adjrisk = FALSE}, the \code{plotCSM} does not adjust for the risk set and
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
#' The response must be a recurrent event survival object as returned by function \code{reSurv}.
#' @param data an optional data frame in which to interpret the variables occurring in the "\code{formula}".
#' @param adjrisk an optional logical value indicating whether risk set will be adjusted. See \bold{Details}.
#' @param onePanel an optional logical value indicating whether cumulative sample means (CSM) will be plotted in the same panel.
#' This is useful when comparing CSM from different groups.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#' 
#' @seealso \code{\link{reSurv}}, \code{\link{plot.reSurv}}
#' @keywords Plots
#' @export
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom dplyr summarise rowwise
#' @importFrom ggplot2 guides guide_legend
#' 
#' @examples
#' data(readmission, package = "frailtypack")
#' plotCSM(reSurv(t.stop, id, event, death) ~ 1, data = readmission)
#' plotCSM(reSurv(t.stop, id, event, death) ~ sex, data = readmission)
#' plotCSM(reSurv(t.stop, id, event, death) ~ sex, data = readmission, onePanel = TRUE)
plotCSM <- function(formula, data, onePanel = FALSE, adjrisk = TRUE, control = list(), ...) {
    call <- match.call()
    ctrl <- plotCSM.control()
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
                mutate(mu = n.x / adjrisk, CSM = cumsum(mu))
        } else {
            dat0 <- eval(parse(text = paste("dat0 %>% group_by(", paste(names(dat1)[5:ncol(dat1)], collapse = ","), ")"))) %>%
                mutate(mu = n.x / n.y, CSM = cumsum(mu))
        }
    } else {
        rec0 <- tmp1 %>% filter(recType == 0)
        dat0 <- tmp1 %>% filter(recType > 0) %>% rowwise() %>%
            mutate(n.y = length(unique(dat1$id)), adjrisk = n.y - sum(rec0$n[Time > rec0$Time]))  %>%
            unrowwise
        if (adjrisk) {
            dat0 <- dat0 %>% mutate(mu = n / adjrisk, CSM = cumsum(mu))
        } else {
            dat0 <- dat0 %>% mutate(mu = n / n.y, CSM = cumsum(mu))
        }
    }
    k <- length(unique(unlist(dat0$recType)))
    if (k == 1) dat0$recType <- factor(dat0$recType, labels = ctrl$recurrent.name)
    if (k > 1 & is.null(ctrl$recurrent.type))
        dat0$recType <- factor(dat0$recType, labels = paste(ctrl$recurrent.name, 1:k))
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            dat0$recType <- factor(dat0$recType, labels = ctrl$recurrent.type)
        } else {
            cat('The length of "recurrent.type" mismatched, default names are used.\n')
            dat0$recType <- factor(recType, labels = paste(ctrl$recurrent.name, 1:k))
        }
    }
    gg <- ggplot(data = dat0, aes(x = Time, y = CSM))
    if (ncol(dat1) == 5 & length(unique(dat1$recType)) == 2) {
        gg <- gg + geom_step(size = ctrl$lwd)
    } else {
        if (!onePanel & length(unique(dat1$recType)) == 2) gg <- gg + geom_step(size = ctrl$lwd)
        if (!onePanel & length(unique(dat1$recType)) > 2) 
            gg <- gg + geom_step(aes(color = recType), direction = "hv", size = ctrl$lwd) +
                guides(color = guide_legend(title = ctrl$recurrent.name))
        if (onePanel) {
            rText <- paste("geom_step(aes(color = interaction(",
                           paste(names(dat0)[5:ncol(dat1) - 4], collapse = ","),
                           ")), direction = \"hv\")")
            gg <- gg + eval(parse(text = rText))
        }
    }
    if (!onePanel && nX > 0 && formula[[3]] != 1) 
        gg <- gg + facet_grid(as.formula(paste(formula[3], "~.", collapse = "")),
                              scales = "free", space = "free", switch = "x")
    if (!onePanel & k == 1)
        gg <- gg + theme(legend.position="none")
    if (onePanel & k == 1)
        gg <- gg + scale_color_discrete(name = "",
                                        labels = levels(interaction(dat0[,(5 + 1):ncol(dat1) - 4])))
    if (onePanel & k > 1) gg <- gg + scale_color_discrete(name = "")
    
    gg + theme(axis.line = element_line(color = "black"),
                legend.key = element_rect(fill = "white", color = "white")) +
        ggtitle(ctrl$main) + labs(y = ctrl$ylab, x = ctrl$xlab)
    ## scale_color_discrete(name = "") +
    ## scale_y_continuous(expand = c(0, 1)) +
}

#' Plotting Baseline Cumulative Rate Function and Baseline Cumulative Hazard Function
#'
#' Plot the baseline cumulative rate function and the baseline cumulative hazard function
#' (if applicable) for an \code{reReg} object.
#'
#' The argument \code{control} consists of options with argument defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty, e.g., "".}
#'   \item{main}{customizable title, default value are "Baseline cumulative rate and hazard function" when \code{baseline = "both"},
#' "Baseline cumulative rate function" when \code{baseline = "rate"}, and "Baseline cumulative hazard function" when \code{baseline = "hazard"}.}
#' }
## #' These arguments can also be passed down without specifying a \code{control} list.
#' 
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether to add a smooth curve (\emph{loess} smooth).
#' @param baseline a character string specifying which baseline function to plot.
#' If \code{baseline = "both"} (default), both the baseline cumulative rate and baseline cumulative hazard function will be plotted in separate panels within the same display;
#' if \code{baseline = "rate"}, only the baseline cumulative rate function will be plotted;
#' if \code{baseline = "hazard"}, only the baseline cumulative hazard function will be plotted.
#' @param control a list of control parameters. See \bold{Details}.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#' 
#' @seealso \code{\link{reReg}}
#' @export
#' @keywords Plots
#'
#' @return A \code{ggplot} object.
#' 
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 geom_smooth geom_step
#' @examples
#' data(readmission, package = "frailtypack")
#' fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
#'              data = subset(readmission, id < 50))
#' plot(fit)
#' plot(fit, baseline = "rate")
#' plot(fit, baseline = "rate", xlab = "Time (days)")
plot.reReg <- function(x, baseline = c("both", "rate", "hazard"),
                       smooth = FALSE, control = list(), ...) {
    baseline <- match.arg(baseline)
    if (baseline == "both") ctrl <- plotRe.control()
    if (baseline == "rate") ctrl <- plotRate.control()
    if (baseline == "hazard") ctrl <- plotHaz.control()
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
    if (!is.reReg(x)) stop("Response must be a reReg class")
    if (baseline == "rate")
        return(plotRate(x, smooth = smooth, control = ctrl))
    if (baseline == "hazard")
        return(plotHaz(x, smooth = smooth, control = ctrl))
    if (baseline == "both" & x$method %in% c("cox.LWYY", "sc.XCYH")) {
        cat(paste("Baseline cumulative hazard function is not available for method = ", x$method, ".", sep = ""))
        cat("\nOnly the baseline cumulative rate function is plotted.\n")
        return(plotRate(x, smooth = smooth, control = ctrl))
    }
    if (baseline == "both" & !(x$method %in% c("cox.LWYY", "sc.XCYH"))) {
        if (is.null(x$rate0.upper)) {
            dat1 <- as_tibble(x$DF) %>% mutate(Y = x$rate0(Time)) %>% select(Time, Y)
            dat2 <- as_tibble(x$DF) %>% mutate(Y = x$haz0(Time)) %>% select(Time, Y)
        } else {
            dat1 <- as_tibble(x$DF) %>%
                mutate(Y = x$rate0(Time), Y.upper = x$rate0.upper(Time), Y.lower = x$rate0.lower(Time)) %>%
                select(Time, Y, Y.upper, Y.lower)
            dat2 <- as_tibble(x$DF) %>%
                mutate(Y = x$haz0(Time), Y.upper = x$haz0.upper(Time), Y.lower = x$haz0.lower(Time)) %>%
                select(Time, Y, Y.upper, Y.lower)
        }
        dat <- bind_rows(dat1, dat2, .id = "group") %>%
            mutate(group = factor(group, levels = 1:2,
                                  labels = c("Baseline cumulative rate",
                                             "Baseline cumulative hazard")))
        gg <- ggplot(data = dat, aes(x = Time, y = Y)) +
            facet_grid(group ~ ., scales = "free") +
            theme(axis.line = element_line(color = "black"),
                  strip.text = element_text(face = "bold", size = 12))   
        if (smooth) {
            gg <- gg + geom_smooth(se = FALSE, method = "loess", col = 1)
            if (!is.null(x$rate0.upper))
                gg <- gg + geom_smooth(aes(x = Time, y = Y.upper),
                                       col = 1, se = FALSE, method = "loess", lty = 2) +
                    geom_smooth(aes(x = Time, y = Y.lower),
                                col = 1, se = FALSE, method = "loess", lty = 2)
        } else {
            gg <- gg + geom_step()
            if (!is.null(x$rate0.upper))
                gg <- gg + geom_step(aes(x = Time, y = Y.upper), lty = 2)+ 
                    geom_step(aes(x = Time, y = Y.lower), lty = 2)
        }
        gg + ggtitle(ctrl$main) + labs(y = ctrl$ylab, x = ctrl$xlab)
    }
}

#' Plotting the Baseline Cumulative Rate Function for the Recurrent Event Process
#'
#' Plot the baseline rate function for an \code{reReg} object.
#' 
#' The argument \code{control} consists of options with argument defaults
#' to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty, e.g., "".}
#'   \item{main}{customizable title, default value is "Baseline cumulative rate function".}
#' }
#' These arguments can also be passed down without specifying a \code{control} list. See \bold{Examples}.
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether the \emph{loess} smoothing will be applied.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#'
#' @seealso \code{\link{reReg}} \code{\link{plot.reReg}}
#' @export
#'
#' @return A \code{ggplot} object.
#'
#' @keywords Plots
#' 
#' @examples
#' ## readmission data
#' data(readmission, package = "frailtypack")
#' set.seed(123)
#' fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
#'                    data = subset(readmission, id < 50),
#'                    method = "am.XCHWY", se = "resampling", B = 20)
#' ## Plot both the baseline cumulative rate and hazard function
#' plot(fit)
#' ## Plot baseline cumulative rate function
#' plotRate(fit)
#' ## Plot with user-specified labels
#' plotRate(fit, xlab = "User xlab", ylab = "User ylab", main = "User title")
#' plotRate(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))
plotRate <- function(x, smooth = FALSE, control = list(), ...) {
    ctrl <- plotRate.control()
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
    if (!is.reReg(x)) stop("Response must be a reReg class")
    if (is.null(x$rate0.upper)) {
        dat <- as_tibble(x$DF) %>% mutate(Y = x$rate0(Time)) %>% select(Time, Y)
    } else {
        dat <- as_tibble(x$DF) %>%
            mutate(Y = x$rate0(Time), Y.upper = x$rate0.upper(Time), Y.lower = x$rate0.lower(Time)) %>%
            select(Time, Y, Y.upper, Y.lower)
    }
    gg <- ggplot(data = dat, aes(x = Time, y = Y)) +
        theme(axis.line = element_line(color = "black"))
    if (smooth) {
        gg <- gg + geom_smooth(se = FALSE, method = "loess", col = 1)
        if (!is.null(x$rate0.upper))
            gg <- gg + geom_smooth(aes(x = Time,  y = Y.upper),
                                   col = 1, se = FALSE, method = "loess", lty = 2) +
                geom_smooth(aes(x = Time,  y = Y.lower),
                            col = 1, se = FALSE, method = "loess", lty = 2) 
    } else {
        gg <- gg + geom_step()
        if (!is.null(x$rate0.upper))
            gg <- gg + geom_step(aes(x = Time,  y = Y.upper), lty = 2) +
                geom_step(aes(x = Time,  y = Y.lower), lty = 2)
    }
    gg + ggtitle(ctrl$main) + labs(x = ctrl$xlab, y = ctrl$ylab)
}

#' Produce the Baseline Cumulative Hazard Function for the Censoring Time
#'
#' Plot the baseline cumulative hazard function for an \code{reReg} object.
#'
#' The argument \code{control} consists of options with argument
#' defaults to a list with the following values:
#' \describe{
#'   \item{xlab}{customizable x-label, default value is "Time".}
#'   \item{ylab}{customizable y-label, default value is empty, e.g., "".}
#'   \item{main}{customizable title, default value is "Baseline cumulative hazard function".}
#' }
#' These arguments can also be passed down without specifying a \code{control} list.
#' See \bold{Examples}.
#'
#' @param x an object of class \code{reReg}, usually returned by the \code{reReg} function.
#' @param smooth an optional logical value indicating whether the \emph{loess} smoothing will be applied.
#' @param control a list of control parameters.
#' @param ... graphical parameters to be passed to methods.
#' These include \code{xlab}, \code{ylab} and \code{main}.
#' 
#' @seealso \code{\link{reReg}} \code{\link{plot.reReg}}
#' @export
#'
#' @return A \code{ggplot} object.
#' @keywords Plots
#' 
#' @examples
#' ## readmission data
#' data(readmission, package = "frailtypack")
#' set.seed(123)
#' fit <- reReg(reSurv(t.stop, id, event, death) ~ sex + chemo,
#'              data = subset(readmission, id < 50),
#'              method = "am.XCHWY", se = "resampling", B = 20)
#' ## Plot both the baseline cumulative rate and hazard function
#' plot(fit)
#' ## Plot baseline cumulative hazard function
#' plotHaz(fit)
#' ## Plot with user-specified labels
#' plotHaz(fit, control = list(xlab = "User xlab", ylab = "User ylab", main = "User title"))  
plotHaz <- function(x, smooth = FALSE, control = list(), ...) {
    if (x$method %in% c("cox.LWYY", "sc.XCYH")) {
        cat(paste("Baseline cumulative hazard function is not available for method = ", x$method, ".", sep = ""))
        return()
    }
    ctrl <- plotHaz.control()
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
    if (!is.reReg(x)) stop("Response must be a reReg class")
    if (is.null(x$haz0.upper)) {
        dat <- as_tibble(x$DF) %>% mutate(Y = x$haz0(Time)) %>% select(Time, Y)
    } else {
        dat <- as_tibble(x$DF) %>%
            mutate(Y = x$haz0(Time), Y.upper = x$haz0.upper(Time), Y.lower = x$haz0.lower(Time)) %>%
            select(Time, Y, Y.upper, Y.lower)
    }
    gg <- ggplot(data = dat, aes(x = Time, y = Y)) +
        theme(axis.line = element_line(color = "black"))
    if (smooth) {
        gg <- gg + geom_smooth(se = FALSE, method = "loess", col = 1)
        if (!is.null(x$rate0.upper))
            gg <- gg + geom_smooth(aes(x = Time,  y = Y.upper),
                                   col = 1, se = FALSE, method = "loess", lty = 2) +
                geom_smooth(aes(x = Time,  y = Y.lower),
                            col = 1, se = FALSE, method = "loess", lty = 2) 
    } else {
        gg <- gg + geom_step()
        if (!is.null(x$rate0.upper))
            gg <- gg + geom_step(aes(x = Time,  y = Y.upper), lty = 2) +
                geom_step(aes(x = Time,  y = Y.lower), lty = 2)
    }
    gg + ggtitle(ctrl$main) + labs(x = ctrl$xlab, y = ctrl$ylab)
}

unrowwise <- function(x) {
  class(x) <- c( "tbl_df", "data.frame")
  x
}

plotEvents.control <- function(xlab = "Time", ylab = "Subject",
                               main = "Recurrent event plot",
                               cex = "Default", 
                               terminal.name = "Terminal event",
                               recurrent.name = "Recurrent events",
                               recurrent.type = NULL, alpha = .7) {
    list(xlab = xlab, ylab = ylab, main = main, cex = cex,
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type, alpha = alpha)
}

plotCSM.control <- function(xlab = "Time", ylab = "Cumulative mean",
                            main = "Sample cumulative mean function plot",
                            lwd = 1,
                            terminal.name = "Terminal event",
                            recurrent.name = "Recurrent events",
                            recurrent.type = NULL) {
    list(xlab = xlab, ylab = ylab, main = main, lwd = lwd, 
         terminal.name = terminal.name, recurrent.name = recurrent.name,
         recurrent.type = recurrent.type)         
}

plotHaz.control <- function(xlab = "Time", ylab = "", main = "Baseline cumulative hazard function") {
    list(xlab = xlab, ylab = ylab, main = main)
}

plotRate.control <- function(xlab = "Time", ylab = "", main = "Baseline cumulative rate function") {
    list(xlab = xlab, ylab = ylab, main = main)
}

plotRe.control <- function(xlab = "Time", ylab = "", main = "Baseline cumulative rate and hazard functions") {
    list(xlab = xlab, ylab = ylab, main = main)
}

