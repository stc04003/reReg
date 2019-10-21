library(reReg)

args(Recur)
?Recur
Recur(3:5)
Recur(6:1, id = rep(1:2, 3))
is.Recur(Recur(3:5))
    
time1 <- c(1, 5, 7)
time2 <- c(3, 7, 9)
ex3 <- Recur(time1 %to% time2, id = c("A1", "A1", "A2"))
ex3
head(ex3)

## #########################################################
## Testing readmission data
## Some subjects in this data have multiple terminal events 
## #########################################################
data(readmission, package = "frailtypack")
head(readmission)
attach(readmission)

Recur(t.start %2% t.stop, id, event, death)
Recur(t.stop, id, event, death)

reSurv(t.start, t.stop, id, event, death)
reSurv(t.stop, id, event, death)

detach(readmission)

testDF <- function(dat) {
    ## `dat` comes from simSC
    ## `fit1` is from reSurv
    ## `fit2` is from Recur 
    fit1 <- with(dat, reSurv(Time, id, event, status))
    fit2 <- with(dat, Recur(Time, id, event, status))
    reDF2 <- as.data.frame(fit2@.Data)
    all(identical(fit1$reDF$Time, reDF2$time2),
        ## all.equal(fit1$reDF$id, reDF2$id),
        all.equal(fit1$reDF$id, as.numeric(levels(fit2@ID))[fit2@ID]),
        identical(fit1$reDF$recType, reDF2$event),
        identical(fit1$reDF$status, reDF2$terminal))
}

testDF(simSC(200, c(-1, 1), c(-1, 1)))
table(replicate(1000, testDF(simSC(200, c(-1, 1), c(-1, 1)))))

## ID in different orders
test1 <- function() {
    dat <- simSC(200, c(-1, 1), c(-1, 1))
    dat <- subset(dat, id %in% sample(1:200, 50))
    testDF(dat)
}

table(replicate(1000, test1()))

test2 <- function() {
    dat <- simSC(200, c(-1, 1), c(-1, 1))
    dat$id <- rep(sample(1:200, 200), table(dat$id))
    testDF(dat)
}

table(replicate(1000, test2()))

## Different event types

dat <- simSC(200, c(-1, 1), c(-1, 1))
dat$event <- dat$event * sample(1:3, nrow(dat), TRUE)
testDF(dat)

head(dat)

debug(testDF)


## #########################################################
## Testing plots
## #########################################################
set.seed(1)
n <- 20
dat <- simSC(n, c(-1, 1), c(-1, 1))
dat$x3 <- rep(sample(0:1, n, TRUE), table(dat$id))
dat

fit1 <- with(dat, reSurv(Time, id, event, status))
fit2 <- with(dat, Recur(Time, id, event, status))
str(fit2)

plotCSM(with(dat, reSurv(Time, id, event, status)))
plotCSM2(with(dat, Recur(Time, id, event, status)))

plotCSM(reSurv(Time, id, event, status) ~ 1, data = dat)
plotCSM2(Recur(Time, id, event, status) ~ 1, data = dat)

plotCSM(reSurv(Time, id, event, status) ~ x1, data = dat)
plotCSM2(Recur(Time, id, event, status) ~ x1, data = dat)

plotCSM(reSurv(Time, id, event, status) ~ x1 + x3, data = dat)
plotCSM2(Recur(Time, id, event, status) ~ x1 + x3, data = dat)

fit2 <- with(dat, Recur(Time, id, event, status))
debug(plotCSM2)

plotCSM2 <- function(formula, data, onePanel = FALSE, adjrisk = TRUE,
                     smooth = FALSE, control = list(), ...) {
    call <- match.call()
    ctrl <- reReg:::plotCSM.control()
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
    } else {
        if (missing(data)) obj <- eval(formula[[2]], parent.frame())
        else obj <- eval(formula[[2]], data)
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
    }
    vNames <- attr(terms(formula), "term.labels")
    dd <- subset(DF, select = c("event", vNames, "time2"))
    dd <- dd[do.call(order, dd),]
    nn <- table(apply(dd, 1, paste, collapse = ""))
    dd <- unique(dd)
    dd$n <- as.integer(nn) ## tmp1 in the 1st version
    rownames(dd) <- NULL
    k <- length(unique(dd$event)) - 1    
    if (ncol(DF) > 6) { ## any covariates/stratifications?
        dd2 <- subset(DF, event == 0, select = vNames)
        dd2 <- dd2[do.call(order, dd2),]
        nn <- table(apply(dd2, 1, paste, collapse = ""))
        dd2 <- unique(dd2)
        dd2$n <- as.integer(nn) ## tmp2 in the 1st version
        rownames(dd2) <- NULL
        tmp1 <- merge(dd, dd2, by = vNames)
        tmp1$GrpInd <- match(apply(dd[,vNames], 1, paste, collapse = ""),
                             apply(dd2[,vNames], 1, paste, collapse = ""))
        ## as.integer(eval(parse(text = paste(attr(terms(formula), "term.labels"), collapse = ":"))))
        rec0 <- subset(tmp1, event == 0)
        dat0 <- do.call(rbind, lapply(split(tmp1, tmp1$GrpInd), function(x) {
            x$adjrisk = apply(x, 1, function(y)
                y[6] - sum(rec0$n.x[y[4] > rec0$time2 & rec0$GrpInd == y[7]]))
            return(x)}))
        dat0$n.x <- dat0$n.x * (dat0$event > 0)
        rec0$time2 <- rec0$n.x <- 0
        rec0$adjrisk <- 1
        dat0 <- unique(rbind(dat0, rec0))
        ## ############################################
        
        if (k > 1) {
            tmp <- dat0 %>% filter(recType == 0)
            dat0 <- bind_rows(dat0 %>% filter(recType > 0),
                              do.call(rbind, lapply(split(tmp, tmp$GrpInd), function(x) x[rep(1:NROW(x), k),] %>% mutate(recType = rep(1:k, each = NROW(x))))))
        } else {
            dat0$recType <- dat0$recType[1]
        }
        if (adjrisk) {
            dat0 <- dat0 %>% arrange(recType, GrpInd, Time) %>% group_by(recType, GrpInd) %>% mutate(mu = n.x / adjrisk, CSM = cumsum(mu))
        } else {
            dat0 <- dat0 %>% arrange(recType, GrpInd, Time) %>% group_by(recType, GrpInd) %>% mutate(mu = n.x / n.y, CSM = cumsum(mu))
        }
    } else { ## no covariates
        rec0 <- tmp1 %>% filter(recType == 0)
        dat0 <- tmp1 %>% rowwise() %>%
            mutate(n.y = length(unique(dat1$id)), adjrisk = n.y - sum(rec0$n[Time > rec0$Time])) %>% unrowwise %>%
            mutate(n = n * (recType > 0))
        dat0 <- bind_rows(dat0, rec0 %>% mutate(Time = 0, n = 0, adjrisk = 1)) %>% distinct()
        dat0$recType <- dat0$recType[1]
        if (adjrisk) {
            dat0 <- dat0 %>% arrange(Time) %>% mutate(mu = n / adjrisk, CSM = cumsum(mu))
        } else {
            dat0 <- dat0 %>% arrange(Time) %>% mutate(mu = n / n.y, CSM = cumsum(mu))
        }
    }
    if (k == 1) dat0$recType <- factor(dat0$recType, labels = ctrl$recurrent.name)
    if (k > 1 & is.null(ctrl$recurrent.type))
        dat0$recType <- factor(dat0$recType, labels = paste(ctrl$recurrent.name, 1:k))
    if (k > 1 & !is.null(ctrl$recurrent.type)) {
        if (length(ctrl$recurrent.type) == k) {
            dat0$recType <- factor(dat0$recType, labels = ctrl$recurrent.type)
        } else {
            cat('The length of "recurrent.type" mismatched, default names are used.\n')
            dat0$recType <- factor(dat0$recType, labels = paste(ctrl$recurrent.name, 1:k))
        }
    }
    gg <- ggplot(data = dat0, aes(x = Time, y = CSM))
    if (ncol(dat1) == 5 & k == 1) {
        gg <- gg + geom_step(size = ctrl$lwd)
    } else {
        if (!onePanel & k == 1) gg <- gg + geom_step(size = ctrl$lwd)
        if (!onePanel & k > 1) 
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
    if (smooth & k == 1) gg <- gg + geom_smooth(method = "loess", size = ctrl$lwd, se = FALSE)
    if (smooth & k > 1) cat('Smoothing only works for data with one recurrent event type.\n')
    gg + theme(axis.line = element_line(color = "black"),
                legend.key = element_rect(fill = "white", color = "white")) +
        ggtitle(ctrl$main) + labs(y = ctrl$ylab, x = ctrl$xlab)
    ## scale_color_discrete(name = "") +
    ## scale_y_continuous(expand = c(0, 1)) +
}

