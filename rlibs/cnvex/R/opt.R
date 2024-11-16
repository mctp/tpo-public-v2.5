CNVEX_EMPTY_GRID <- structure(list(p = 0.99, P = 1.98,
                                   L = NA_real_, P1 = NA_real_, cP1 = NA_real_),
                                   row.names = c(NA, 1L), class = c("data.table", "data.frame"))
CNVEX_EMPTY_CAND <- structure(list(cand = "00", p = 0.99, P = 1.98,
                                   L = NA_real_, P1 = NA_real_, cP1 = NA_real_),
                                   row.names = c(NA, 1L), class = c("data.table", "data.frame"))
CNVEX_EMPTY_FINE <- structure(list(cand = "00", data = NA, converged = TRUE, p = 0.99, P = 1.98),
                                   row.names = c(NA, 1L), class = c("data.table", "data.frame"))
CNVEX_EMPTY_EVAL <- structure(list(cand = "00", data = NA, converged = TRUE,
                                   p = 0.99, P = 1.98, nlr = NA_integer_, naf = NA_integer_,
                                   tL = NA_real_, hL = NA_real_, aL = NA_real_, bL = NA_real_,
                                   P1 = NA_real_, cP = NA_real_, sP = NA_real_, mC = NA_integer_, aD = NA_real_,
                                   CN = NA_real_, C0 = NA_real_, C1 = NA_real_, C2 = NA_real_, C3 = NA_real_,
                                   C4 = NA_real_, C5 = NA_real_, C6 = NA_real_, C7 = NA_real_, C8 = NA_real_,
                                   pA = NA_real_, aE = NA_real_,
                                   ploh = NA_real_, podd = NA_real_, plow = NA_real_, pdip = NA_real_, pdel = NA_real_,
                                   rank = 1
                                   ), row.names = c(NA, 1L), class = c("data.table", "data.frame"))

.opt.ploidy <- function(C, sC, nC, len, max.C) {
    weighted.mean(ifelse(sC < max.C + 0.5, C * (nC/2), sC * (nC/2)), as.double(len))
}

.opt.clonal.ploidy <- function(C, nC, len) {
    weighted.mean(C * (nC/2), as.double(len))
}

.opt.subclonal.ploidy <- function(sC, nC, len) {
    weighted.mean(sC * (nC/2), as.double(len))
}

.opt.grid.localmax <- function(grid, res, x, y, var) {
    log_info("Searching for local maxima at resolution={res} and likelihood='{var}'")
    mx <- as.matrix(dcast.data.table(grid, reformulate(x, response=y), value.var=var)[,-1])
    r <- raster(mx)
    ## nearest odd integer >= to 3
    rres <- max(2*floor((nrow(r)*res)/2)+1, 3)
    cres <- max(2*floor((ncol(r)*res)/2)+1, 3)
    wind <- matrix(1, nrow=rres, ncol=cres)
    log_info("Local search window: {rres}x{cres}")
    localmax <- focal(r, fun = max_na_rm, w = wind, pad=TRUE, padValue=NA)
    cand <- grid[Which(localmax==r, cells=TRUE)]
    cand <- cbind(cand=sprintf("%03d", seq_len(nrow(cand))), cand)
    log_info("Found {nrow(cand)} candidates")
    cand <- rbind(cand, CNVEX_EMPTY_CAND[0])
    return(cand)
}

.opt.llik.lr.C <- function(lrC, seg, p, P, monotonic, opts) {
    ## optimize C
    C.pick <- .llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$opt.max.C)[seg,nomatch=0]
    if (monotonic) {
        P1 <- C.pick[(len / nlr < opts$opt.max.len.per.probe) & nC != 0, .opt.ploidy(C, sC, nC, len, opts$opt.max.C)]
        C.pick <- .llik.lrC.outer(lrC, p, P1, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$opt.max.C)[seg,nomatch=0]
    }
    cand <- C.pick[(len / nlr < opts$opt.max.len.per.probe) & nC != 0,.(
        p = ..p,
        P = ..P,
        tL = sum(llik), ## total per-tile likelihood
        sL = sum(dnorm(aD, sd=lrC$sd[1]/sqrt(nlr))), ## total per-segment likelihood
        hL = sum(pmin(2 * aD, 1) * nlr * log(0.5)), ## heuristic clonality penalty
        aL = NA_real_,
        P1 = .opt.ploidy(C, sC, nC, len, opts$opt.max.C),
        cP1 = .opt.clonal.ploidy(C, nC, len),
        sP1 = .opt.subclonal.ploidy(sC, nC, len)
    )]
    return(cand)
}

.opt.llik.af.K <- function(lrC, af, seg, p, P, opts) {
    ## optimize K (given C)
    C.pick <- .llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$opt.max.C)[seg,nomatch=0]
    afC <- .af.grid.afC(C.pick[,.(seg, C)], af, p, P, opts$opt.max.C)
    K.pick <- .llik.afC.outer(afC, opts$opt.p.af.anom, opts$opt.dp.af.max, TRUE)[seg,nomatch=0]
    cand <- C.pick[(len / nlr < opts$opt.max.len.per.probe) & nC != 0,.(
        p = ..p,
        P = ..P,
        tL = sum(llik), ## total per-tile likelihood
        sL = sum(dnorm(aD, sd=lrC$sd[1]/sqrt(nlr))), ## total per-segment likelihood
        hL = sum(pmin(2 * aD, 1) * nlr * log(0.5)), ## heuristic clonality penalty
        aL = NA_real_,
        P1 = .opt.ploidy(C, sC, nC, len, opts$opt.max.C),
        cP1 = .opt.clonal.ploidy(C, nC, len),
        sP1 = .opt.subclonal.ploidy(sC, nC, len)
    )]
    cand$aL <- sum(K.pick$llik)
    return(cand)
}

.opt.llik.lr.p <- function(lrC, seg, pi, Pi, first, opts) {
    ## optimize p
    p.lox <- max(pi-opts$opt.fine.p.off, opts$opt.p.lo)
    p.hix <- min(pi+opts$opt.fine.p.off, opts$opt.p.hi)
    pP <- cbind(p=seq(p.lox, p.hix, opts$opt.fine.p.res), P=Pi)
    tmp <- rbindlist(lapply(seq_len(nrow(pP)), function(i) {
        .opt.llik.lr.C(lrC, seg, pP[[i,1]], pP[[i,2]], monotonic=TRUE, opts)
    }))
    cand <- tmp[order(-(tL+sL+hL)), .SD[1]]
    ## log_debug("Optimizing candidate {format(pi, digits=3)}->{format(cand$p, digits=3)} {format(cand$P, digits=3)}")
    return(cand)
}

.opt.llik.af.p <- function(lrC, af, seg, pi, Pi, opts) {
    ## optimize p
    p.lox <- max(pi-opts$opt.fine.p.off, opts$opt.p.lo)
    p.hix <- min(pi+opts$opt.fine.p.off, opts$opt.p.hi)
    pP <- cbind(p=seq(p.lox, p.hix, opts$opt.fine.p.res), P=Pi)
    tmp <- rbindlist(lapply(seq_len(nrow(pP)), function(i) {
        .opt.llik.af.K(lrC, af, seg, pP[[i,1]], pP[[i,2]], opts)
    }))
    cand <- tmp[order(-aL), .SD[1]]
    ## log_debug("Optimizing candidate {format(pi, digits=3)}->{format(cand$p, digits=3)} {format(cand$P, digits=3)}")
    return(cand)
}

.get.opt.heuristic.llik <- function(opts) {
    if (opts$opt.grid.llik=="heuristic-v1") {
        ll.fun <- function(llik, aD, nlr, grid.subc.exp=opts$opt.grid.subc.exp) {
            sum((1 - (2 * pmin(abs(aD), 0.5)))**grid.subc.exp * llik) ## heuristic v1
        }
    } else if (opts$opt.grid.llik=="heuristic-v2") {
        ll.fun <- function(llik, aD, nlr, grid.subc.prob=opts$opt.grid.subc.prob) {
            sum(llik) + sum(pmin(2 * aD, 1) * nlr * log(grid.subc.prob)) ## heuristic v2
        }
    } else {
        ll.fun <- function(llik, aD, nlr) {
            sum(llik) ## default
        }
    }
    return(ll.fun)
}

.opt.grid <- function(data, opts) {
    ## grid over all p and P combinations
    pP <- .lr.grid.pP(opts$opt.grid.p.res, opts$opt.p.lo, opts$opt.p.hi, opts$opt.P.lo, opts$opt.P.hi)
    ## grid over all C's for each segment
    ## only.hq = TRUE, want to fit to hq data
    lrC <- .lr.grid.lrC(data$lr, data$stats$sd.lr, opts$opt.grid.max.C, opts$opt.only.hq)
    log_info("Grid search purity:{opts$opt.p.lo}-{opts$opt.p.hi} ploidy:{opts$opt.P.lo}-{opts$opt.P.hi} resolution:{opts$opt.grid.p.res} max.c:{opts$opt.grid.max.C}")
    .opt.heuristic.llik <- .get.opt.heuristic.llik(opts)
    ## calculate likelihood
    grid <- rbindlist(foreach(i=seq_len(nrow(pP))) %dopar% {
        ## nomatch=0, because some seg may not have hq
        x <- .llik.lrC.outer(lrC, pP[[i,1]], pP[[i,2]], opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$opt.max.C)[data$seg,nomatch=0]
        ## pick segments with enough probe density, ignore Y in males
        y <- x[(len / nlr < opts$opt.max.len.per.probe) & nC != 0,.(
                p = pP[[i,1]], # input purity
                P = pP[[i,2]], # input ploidy
                L = .opt.heuristic.llik(llik, aD, nlr),
                P1 = .opt.ploidy(C, sC, nC, len, opts$opt.max.C),
                cP1 = .opt.clonal.ploidy(C, nC, len)
                )]
        log_debug("grid-point purity:{format(pP[[i,1]], digits=3)} ploidy:{format(pP[[i,2]], digits=3)} L:{format(y$L, digits=2)} P1:{format(y$P1, digits=2)}")
        y
    })
    grid <- rbind(grid, CNVEX_EMPTY_GRID[0])
    return(grid)
}

.opt.find.cand <- function(grid, opts) {
    cand <- .opt.grid.localmax(grid, res=opts$opt.cand.res, x="p", y="P", var="L")
    return(cand)
}

.opt.cand.lr <- function(lrC, af, seg, cand, opts, first=FALSE) {
    trace <- NULL
    iter <- 0
    pi <- cand[,p]
    Pi <- cand[,P1]
    cPi <- cand[,cP1]
    i <- cand[,cand]
    log_debug("Optimizing candidate {i} <lr> starting purity={format(pi, digits=3)} ploidy={format(Pi, digits=3)} first={first}")
    while(TRUE) {
        iter <- iter + 1
        ## log_debug("Optimizing candidate lr={i} iteration={iter} input p={format(pi, digits=3)} P={format(Pi, digits=3)} cP={format(cPi, digits=3)}")
        tracej <- .opt.llik.lr.p(lrC, seg, pi, Pi, first, opts)
        pj <- tracej[,p]
        Pj <- tracej[,P1]
        cPj <- tracej[,cP1]
        tracej$converged <- (cPi==cPj)
        trace <- rbind(trace, tracej)
        if ((iter==opts$opt.cand.max.iter) || tracej$converged || first) {
            break
        }
        pi <- pj
        Pi <- Pj
        cPi <- cPj
    }
    trace[,":="(cand=i, iter=.I)]
    cand <- trace[order(!converged, -(tL+hL+sL))][1]
    log_debug("Optimizing candidate {i} <lr> best purity={format(cand$p, digits=3)} ploidy={format(cand$P1, digits=3)}")
    return(cand)
}

.opt.cand.af <- function(lrC, af, seg, cand, opts, first=FALSE) {
    trace <- NULL
    iter <- 0
    pi <- cand[,p]
    Pi <- cand[,P1]
    i <- cand[,cand]
    log_debug("Optimizing candidate {i} <af> starting purity={format(pi, digits=3)} ploidy={format(Pi, digits=3)} first={first}")
    while(TRUE) {
        iter <- iter + 1
        ## log_debug("Candidate af candidate={i} iteration={iter} pi={format(pi, digits=3)}")
        tracej <- .opt.llik.af.p(lrC, af, seg, pi, Pi, opts)
        pj <- tracej[,p]
        tracej$converged <- (pi==pj)
        trace <- rbind(trace, tracej)
        if ((iter==opts$opt.cand.max.iter) || tracej$converged || first) {
            break
        }
        pi <- pj
    }
    trace[,":="(cand=i, iter=.I)]
    cand <- trace[order(!converged, -aL)][1]
    log_debug("Optimizing candidate {i} <af> best purity={format(cand$p, digits=3)} ploidy={format(cand$P1, digits=3)}")
    return(cand)
}

.opt.cand <- function(i, lrC, af, seg, cand, opts) {
    log_info("Optimizing candidate {i}")
    cand.lr0 <- cand
    ## first-iteration
    log_debug("Optimizing candidate {i} first-iteration <lr>")
    cand.lr1 <- .opt.cand.lr(lrC, af, seg, cand.lr0, opts, first=TRUE)
    log_debug("Optimizing candidate {i} first-iteration <af>")
    cand.af1 <- .opt.cand.af(lrC, af, seg, cand.lr0, opts, first=TRUE)
    ## full optimization
    log_debug("Optimizing candidate {i} full-optimization <lr>")
    cand.lr2 <- .opt.cand.lr(lrC, af, seg, cand.lr1, opts, first=FALSE)
    log_debug("Optimizing candidate {i} full-optimization <af>")
    cand.af2 <- .opt.cand.af(lrC, af, seg, cand.af1, opts, first=FALSE)
    ## last-iteration
    log_debug("Optimizing candidate {i} last-iteration <lr-af>")
    cand.lr3 <- .opt.cand.af(lrC, af, seg, cand.lr2, opts, first=TRUE)
    log_debug("Optimizing candidate {i} last-iteration <af-lr>")
    cand.af3 <- .opt.cand.lr(lrC, af, seg, cand.af2, opts, first=TRUE)
    ##
    if (cand.lr2$converged) {
        cand.lr <- data.table(
            cand=sprintf("%s-lr", i),
            data="lr",
            converged=TRUE,
            p=cand.lr3$p,
            P=cand.lr3$P
        )
    } else {
        cand.lr <- data.table(
            cand=sprintf("%s-lr", i),
            data="lr",
            converged=FALSE,
            p=cand.lr1$p,
            P=cand.lr1$P
        )
    }
    if (cand.af2$converged) {
        cand.af <- data.table(
            cand=sprintf("%s-af", i),
            data="af",
            converged=TRUE,
            p=cand.af3$p,
            P=cand.af3$P
        )
    } else {
        cand.af <- data.table(
            cand=sprintf("%s-af", i),
            data="af",
            converged=FALSE,
            p=cand.af1$p,
            P=cand.af1$P
        )
    }
    fine <- rbind(cand.lr, cand.af)
    fine <- rbind(fine, CNVEX_EMPTY_FINE[0])
    return(fine)
}

.opt.fine <- function(data, cand, opts) {
    log_info("Optimizing {nrow(cand)} candidates")
    lrC <- .lr.grid.lrC(data$lr, data$stats$sd.lr, opts$opt.max.C, opts$opt.only.hq)
    fine <- rbindlist(foreach(i=seq_len(nrow(cand))) %dopar% {
        .opt.cand(cand[i]$cand, lrC, data$af, data$seg, cand[i], opts)
    })
    return(fine)
}

.opt <- function(data, nogrid, nofine, opts) {
    if (!nogrid) {
        ## grid search
        grid <- .opt.grid(data, opts)
        ## find candidates
        cand <- .opt.find.cand(grid, opts)
    } else {
        grid <- CNVEX_EMPTY_GRID
        cand <- CNVEX_EMPTY_CAND
    }
    if (!nogrid && !nofine) {
        ## fine-tuning
        fine <- .opt.fine(data, cand, opts)
        ## evaluation
        eval <- .eval(data, fine, opts)
    } else {
        fine <- CNVEX_EMPTY_FINE
        eval <- CNVEX_EMPTY_EVAL
    }
    models <- list(grid=grid, cand=cand, fine=fine, eval=eval)
    return(models)
}
