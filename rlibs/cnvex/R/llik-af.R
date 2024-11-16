.af.grid.MC <- function(p, P, max.C) {
    ## M - number of chromosome with variant
    ## C - total number of chromosomes
    ## Ef - expected frequency of variant
    MC <- as.data.table(expand.grid(M=0:max.C,C=0:max.C))[M<=C]
    MC[,":="(
        K=pmin(M,C-M),
        Ef=(p * M + 1 * (1-p)) / (p * C + 2 * (1-p)),
        prior.M=ifelse(M==C-M, log(1), log(1/2))
    )]
    return(MC[])
}

.af.grid.afC <- function(C.pick, af, p, P, max.C) {
    ## get M K Ef's for each segment C
    MC <- .af.grid.MC(p, P, max.C)
    setkey(C.pick, C)
    setkey(MC, C)
    MC.pick <- MC[C.pick[,.(seg, C)],allow.cartesian=TRUE]
    ## get M K Ef's for each variant
    setkey(MC.pick, seg)
    setkey(af, seg)
    ## some af's will be in segments without
    afC <- MC.pick[af,allow.cartesian=TRUE,nomatch=0]
    setkey(C.pick, seg)
    return(afC)
}

.llik.afC.inner <- function(afC, p.af.anom, dp.af.max) {
    afC[,":="(
            beta=dbeta(
                Ef,
                shape1=   AF  * pmin(DP, dp.af.max) + 1,
                shape2=(1-AF) * pmin(DP, dp.af.max) + 1,
                log=TRUE
            ),
            unif=1 # dbeta(x, 1, 1) or dunif(x, 0, 1)
    )]
    bllik <- afC[,beta + prior.M + log(1-p.af.anom)]
    ullik <- afC[,unif + prior.M + log(p.af.anom)]
    afC[,":="(
        llik=plog_sum_exp(bllik, ullik),
        anom=bllik < ullik,
        se=(Ef-AF)**2
    )]
    ## pick best M per K C
    ## in regions where sC >> C, this will result in incorrect C
    ## and unrealistic K (prior forces it to be balanced #57).
    x <- afC[order(-llik),.SD[1],.(idx, K, C)]
    ## total/average by K C
    x <- x[,.(
        llik=sum(llik),
        anom=mean(anom),
        mse=mean(se)
    ), .(seg, K, C)]
    return(x[])
}

.llik.afC.outer <- function(afC, p.af.anom, dp.af.max, collapse) {
    x <- .llik.afC.inner(afC, p.af.anom, dp.af.max)
    if (collapse) {
        x <- x[order(-llik),.SD[1],.(seg, C)]
    }
    setkey(x, seg)
    return(x[])
}

.llik.afC.inner.full <- function(af, MC, p.af.anom, dp.af.max) {
    ## af data grid generation
    tmp1 <- rbindlist(rep(list(af), nrow(MC)))[order(idx)]
    tmp2 <- rbindlist(rep(list(MC), nrow(af)))
    afC <- cbind(tmp1, tmp2)
    x <- .llik.afC.inner(afC, p.af.anom, dp.af.max)
    return(x[])
}

.llik.afC.outer.full <- function(af, p, P, max.C, p.af.anom, dp.af.max, collapse) {
    ## get af
    MC <- .af.grid.MC(p, P, max.C)
    afs <- split(af, af$seg)
    x <- rbindlist(lapply(afs, function(af) {
        .llik.afC.inner.full(af, MC, p.af.anom, dp.af.max)
    }))
    if (collapse) {
        x <- x[order(-llik),.SD[1],.(seg, C)]
    }
    return(x[])
}
