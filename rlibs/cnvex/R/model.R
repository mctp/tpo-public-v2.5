.mod.fix.XY <- function(mod) {
    mod[, ":="(
        C=round(C * nC/2),
        sC=sC * nC/2,
        raw.sC=raw.sC * nC/2,
        aD=aD * nC/2
    )]
    return(mod)
}

.mod.fix.LOH <- function(mod) {
    ## prune unreliable LOH calls
    col.count <- ncol(mod)
    mod[nC %in% c(1, 0), K:=0]
    mod[C %in% c(1, 0), K:=0]
    mod[,":="(nextC = c(NA_real_,C[-.N]), prevC = c(C[2:.N],NA_real_))]
    mod[,":="(nextK = c(NA_real_,K[-.N]), prevK = c(K[2:.N],NA_real_))]
    mod[,newK := fcase(
        (C == nextC & C == prevC & nextK <= prevK), nextK,
        (C == nextC & C == prevC & nextK > prevK), prevK,
        (C == nextC & C != prevC), nextK,
        (C != nextC & C == prevC), prevK,
        (C != nextC & C != prevC), NA_real_
    )]
    mod[(naf<5 | is.na(naf)) & (K==0), K := newK] # TODO: any case where naf=NA and K==0?
    ## correct for C=1 (possible to be NA after prune)
    mod[nC %in% c(1, 0), K:=0]
    mod[C %in% c(1, 0), K:=0]
    mod <- mod[,1:col.count]
    return(mod)
}

.mod.fix.anom <- function(mod) {
    ## ignore K/aL for small amplicons
    ## attempts to fix #58
    mod[anom > 0.5 & sC > C+2 & naf < 50,
        ":="(K=NA, aL=NA)
    ]
}

.mod.fix.pP <- function(mod, p, P) {
    if (is.na(p) || is.na(P)) {
        mod[,
            ":="(C=NA, K=NA)
            ]
    }
    return(mod)
}

.fix.model <- function(mod, p, P, opts) {
    mod <- .mod.fix.XY(mod)
    mod <- .mod.fix.LOH(mod)
    mod <- .mod.fix.anom(mod)
    mod <- .mod.fix.pP(mod, p, P)
    return(mod)
}

.get.model.lrC <- function(lrC, af, seg, p, P, opts) {
    C.pick <- .llik.lrC.outer(lrC, p, P, opts$opt.max.sC, opts$opt.p.lr.anom, TRUE, opts$opt.max.C)
    setnames(C.pick, "llik", "tL")
    if (!is.null(af) && nrow(af) > 0) {
        afC <- .af.grid.afC(C.pick, af, p, P, opts$opt.max.C)
        K.pick <- .llik.afC.outer(afC, opts$opt.p.af.anom, opts$opt.dp.af.max, TRUE)
        setnames(K.pick, "llik", "aL")
        setkey(K.pick, seg, C)
        setkey(C.pick, seg, C)
        CK.pick <- merge(K.pick, C.pick, all = TRUE)
    } else {
        setkey(C.pick, seg, C)
        CK.pick <- C.pick
        CK.pick[,":="(K=NA_integer_, aL=NA_real_, anom=NA_real_, mse=NA_real_)]
    }
    setkey(seg, seg)
    mod <- seg[CK.pick]
    return(mod)
}

.get.model <- function(data, p, P, only.hq, opts) {
    lrC <- .lr.grid.lrC(data$lr, data$stats$sd.lr, opts$opt.max.C, only.hq)
    mod <- .get.model.lrC(lrC, data$af, data$seg, p, P, opts)
    mod <- .fix.model(mod, p, P, opts)
    return(mod)
}

.get.models <- function(data, fine, only.hq, opts) {
    lrC <- .lr.grid.lrC(data$lr, data$stats$sd.lr, opts$opt.max.C, only.hq)
    mods <- foreach(i=seq_len(nrow(fine))) %dopar% {
        mod <- .get.model.lrC(lrC, data$af, data$seg, p=fine[i,p], P=fine[i,P], opts)
        mod <- .fix.model(mod, p=fine[i,p], P=fine[i,P], opts)
    }
    return(mods)
}

.eval.model <- function(mod, opts) {
    tmp <- mod[(len / nlr < opts$opt.max.len.per.probe)]
    .opt.heuristic.llik <- .get.opt.heuristic.llik(opts)
    eval <- tmp[,.(
        nlr = sum(nlr), ## total number of tiles
        naf = sum(naf), ## total number of SNPs
        ## likelihoods
        tL = sum(tL), ## total segment log-likelihood
        hL = .opt.heuristic.llik(tL, aD, nlr),
        aL = sum(aL, na.rm=TRUE),
        ## BIC inspired LL correction
        bL = (2*sum(tL, na.rm=TRUE) - log(sum(nlr, na.rm=TRUE)) * nrow(tmp) * weighted.mean(C, len)) +
             (2*sum(aL, na.rm=TRUE) - log(sum(naf, na.rm=TRUE)) * nrow(tmp) * weighted.mean(C, len)),
        ## ploidy
        P1 = weighted.mean(ifelse(sC < opts$opt.max.C + 0.5, C, sC), len), ## total heuristic ploidy
        cP = weighted.mean(C, len), ## total clonal ploidy
        sP = weighted.mean(sC, len), ## sub-clonal ploidy
        mC = weighted.median(C, as.numeric(len)), ## median C
        aD = weighted.mean(aD, len), ## average sub-clonal allele Deviation
        ## copy-number
        CN = weighted.mean(raw.sC < -0.5, len), ## proportion negative copy-number
        C0 = weighted.mean(C==0, len), C1 = weighted.mean(C==1, len), C2 = weighted.mean(C==2, len),
        C3 = weighted.mean(C==3, len), C4 = weighted.mean(C==4, len), C5 = weighted.mean(C==5, len),
        C6 = weighted.mean(C==6, len), C7 = weighted.mean(C==7, len), C8 = weighted.mean(C>=8, len),
        ## BAF
        pA = weighted.mean(anom, naf), ## proportion Anomalous
        aE = weighted.mean(mse, naf), ## mean-squered allele Error
        ## heuristics, compute means over autosomes
        ## - ignore X and Y in XY
        ## - ignore Y in XX
        ploh = weighted.mean(nC == 2 & K %in% 0, len, na.rm=TRUE), # proportion LOH, assume NA het
        podd = weighted.mean(nC == 2 & C %in% c(1,3,5,7,9), len), # proportion odd copy number
        plow = weighted.mean(nC == 2 & C %in% c(0,1,2,3), len), # proportion low copy number
        pdip = weighted.mean(nC == 2 & C == 2 & K %in% 1, len, na.rm = TRUE), # proportion diploid, assume NA homo
        pdel = weighted.mean(nC == 2 & C == 0, len) ## proportion deleted
    )]
    return(eval)
}

.eval.models.inner <- function(mods, opts) {
    eval <- rbindlist(foreach(mod=mods) %dopar% {
        .eval.model(mod, opts)
    })
    return(eval)
}

.rank.models.inner <- function(eval, opts) {
    eval.in <- eval
    ## top-3 models
    ## diploid model, make sure not too many homozygous deletions, rank by aL+tL
    eval.lo <- eval.in[P<2.7]
    eval.lo <- head(eval.lo[order(pdel>0.005, -(aL+tL))], 1)
    ## tetraploid model, make sure we have some odd-copy chromosomes, rank by aL+tL
    eval.mi <- eval.in[P>=2.7 & P <=4.7]
    eval.mi <- head(eval.mi[order(ploh < 0.025, pdel>0.005, podd < 0.025, -(aL+tL))], 1)
    ## high-ploidy model, make sure we have some LOH and some chromosomes C<=3, dislike models with most chromosomes odd
    eval.hi <- eval.in[P>4.7]
    eval.hi <- head(eval.hi[order(ploh < 0.025, podd > 0.5 | podd < 0.025, plow < 0.025, -(aL+hL))], 1)
    eval.t3 <- rbind(eval.lo, eval.mi, eval.hi)
    ## among the 3 top, pick based on BIC
    eval.t3[order(-bL), rank:=1:.N]
    ## bottom N-models
    setkey(eval.in, cand)
    setkey(eval.t3, cand)
    eval.bn <- eval.in[!eval.t3]
    eval.bn[order(-bL), rank:=(1:.N) + nrow(eval.t3)]
    eval <- rbind(eval.t3, eval.bn)[order(rank)]
    return(eval)
}

.eval.models <- function(mods, fine, opts) {
    eval <- .eval.models.inner(mods, opts)
    eval <- cbind(fine, eval)
    eval <- .rank.models.inner(eval, opts)
    eval <- rbind(eval, CNVEX_EMPTY_EVAL[0]) ## this is just a check to keep updated
    return(eval)
}

.eval <- function(data, fine, opts) {
    ## got model for each fine candidate
    ## only.hq == TRUE, because we want to rank models
    ## based on their hq llik
    mods <- .get.models(data, fine, opts$opt.only.hq, opts)
    ## evaluate and rank each model
    eval <- .eval.models(mods, fine, opts)
    return(eval)
}
