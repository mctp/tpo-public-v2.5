## pP grid generation
.lr.grid.pP <- function(grid.p.res, p.lo, p.hi, P.lo, P.hi) {
    p <- seq(p.lo, p.hi, grid.p.res)
    P <- seq(P.lo, P.hi, length.out=length(p))
    pP <- as.matrix(expand.grid(p=p, P=P))
    return(pP)
}

## lr data grid generation
.lr.grid.lrC <- function(lr, sd.lr, max.C, only.hq=FALSE) {
    ## only keep tiles with coverage for llik calculations
    lr <- lr[!is.na(lr)]
    ## only keep hq tiles for llik calculations
    ## this creates a problem because some segments will not have any hq tiles
    if (only.hq) {
        lr <- lr[(hq)]
    }
    ## The 2 makes sure that we generate sufficiently high C's for haploid
    ## chromosomes with nC==1
    lrC <- expand.grid(lr=lr$lr, C=0:(max.C * 2))
    lrC$seg <- lr$seg
    lrC$sd <- sd.lr
    lrC$nC <- lr$nC
    setDT(lrC)
    ## if nc==1 keep all C's (doubled)
    lrC <- lrC[nC==1 | C<=max.C]
    setkey(lrC, C, seg)
    return(lrC[])
}

## lr likelihood function
.llik.lrC.inner <- function(lrC, p, P, p.lr.anom) {
    lrC[,":="(
        ## normal distribution to model clonal variants
        norm = dnorm(lr, mean=log2((p * C + (1-p) * 2) / ((p * P) + (1-p) * 2)),
                     sd=sd, log=TRUE),
        ## uniform distribution to model anomalies
        unif = dunif(0, min=-6, max=6, log=TRUE)
    )]
    lrC[,":="(
        llik=plog_sum_exp(norm + log(1-p.lr.anom), unif + log(p.lr.anom))
        ## subc=unif>norm
    )]
    x <- lrC[,.( ## sum by C and seg
        lr = mean(lr),
        llik = sum(llik),
        ## subc = mean(subc),
        ## sd = sd[1],
        nC = median(nC)
    ), by=.(C, seg)]
    return(x)
}

.llik.lrC.outer <- function(lrC, p, P, max.sC, p.lr.anom, collapse, max.C) {
    ## compute likelihood for each segment and each C
    x <- .llik.lrC.inner(lrC, p, P, p.lr.anom)
    ## ML sub-clonal copy-number
    x[,raw.sC:=(2^(lr) * ((p * P) + (1-p) * 2) - ((1-p) * 2)) / p]
    x[,sC:=pmax(pmin(raw.sC, max.sC / (nC/2)), 0)]
    x[,aD:=abs(C - sC)]
    ## for each segment pick C with highest llik
    if (collapse) {
      x <- x[order(-llik, ifelse(lr>0, -C, C)),.SD[1],by=seg]
      ## replace C with max.C for high sC values
      ## for high C, sC is calculated by lr and is immune
      ## to possible model failures
      x[sC > max.C, C := max.C]
    }
    setkey(x, seg)
    return(x[])
}
