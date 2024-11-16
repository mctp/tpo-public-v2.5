.breakpoint.stat <- function(bpt, eff.frg) {
    ## dt selectors
    sl.sel <- quote(l3=="sl")
    ts.sel <- quote(l3=="ts")
    art.sel <- quote(art==TRUE)
    dup.sel <- quote(art.5=="nil" & art.3=="nil" & l2=="dist" & l3=="nd" & !d2a)
    spn.dist.nd.sel <- quote(l1=="spn" & l2=="dist" & l3=="nd")
    spn.prox.nd.sel <- quote(l1=="spn" & l2=="prox" & l3=="nd")
    dist.nd.sel <- quote(!art & l2=="dist" & l3=="nd")
    prox.nd.sel <- quote(!art & l2=="prox" & l3=="nd")
    dist.lig.nd.sel <- quote(l1=="spn" & l2=="dist" & l3=="nd" & !d2a)
    prox.lig.nd.sel <- quote(l1=="spn" & l2=="prox" & l3=="nd" & !d2a)
    spn.dist.sel <- quote(l1=="spn" & l2=="dist")
    enc.dist.sel <- quote(l1=="enc" & l2=="dist")
    spn.dist.gtag.sel <- quote(l1=="spn" & l2=="dist" & ((str_sub(mot.5,7,8)=="GT") & (str_sub(mot.3,7,8)=="AG")))
    enc.dist.gtag.sel <- quote(l1=="enc" & l2=="dist" & ((str_sub(mot.5,7,8)=="GT") & (str_sub(mot.3,7,8)=="AG")))
    spn.prox.d2a.dup.sel <- quote(l1=="spn" & l2=="prox" & d2a & topo=="dup")
    ## artifact rate
    art.rte <- bpt[eval(art.sel), sum(sum.jnc)] / eff.frg
    ## duplication rate
    dup.rte <- bpt[eval(dup.sel), 1-sum(unq.sum.jnc)/sum(sum.jnc)]
    ## ligation rate
    spn.dist.rte <- bpt[eval(spn.dist.nd.sel), sum(sum.jnc)] / bpt[eval(dist.nd.sel), sum(sum.jnc)]
    spn.prox.rte <- bpt[eval(spn.prox.nd.sel), sum(sum.jnc)] / bpt[eval(prox.nd.sel), sum(sum.jnc)]
    spn.dist.lig <- bpt[eval(dist.lig.nd.sel), sum(sum.jnc)]
    spn.prox.lig <- bpt[eval(prox.lig.nd.sel), sum(sum.jnc)]
    dist.lig <- spn.dist.lig / spn.dist.rte
    prox.lig <- spn.prox.lig / spn.prox.rte
    lig.rte <- (dist.lig+prox.lig) / eff.frg
    ## trans-splice log fold-change
    spn.dist <- bpt[eval(spn.dist.sel), sum(sum.jnc)]
    enc.dist <- bpt[eval(enc.dist.sel), sum(sum.jnc)]
    spn.dist.gtag <- bpt[eval(spn.dist.gtag.sel), sum(sum.jnc)]
    enc.dist.gtag <- bpt[eval(enc.dist.gtag.sel), sum(sum.jnc)]
    gtag.exp.rte <- (enc.dist.gtag/enc.dist)
    gtag.obs.rte <- (spn.dist.gtag/spn.dist)
    ts.lfc <- log2(gtag.obs.rte / gtag.exp.rte)
    ## trans-splice rate
    gtag.exp.cnt <- spn.dist * gtag.exp.rte
    gtag.obs.cnt <- spn.dist.gtag
    ts.rte.lo <- (bpt[eval(ts.sel), sum(sum.jnc)] / spn.dist.rte) / eff.frg
    ts.rte.hi <- (max(gtag.obs.cnt - gtag.exp.cnt, 0) / spn.dist.rte) / eff.frg
    ## back-splice rate
    spn.prox.d2a <- bpt[eval(spn.prox.d2a.dup.sel), sum(sum.jnc)]
    prox.d2a <- spn.prox.d2a / spn.prox.rte
    bs.rte <- prox.d2a / eff.frg
    ## stem-loop rate
    sl.rte <- bpt[eval(sl.sel), sum(sum.jnc)] / eff.frg
    bpt.surv <- bpt[,.(
        n=.N,
        hq.sum.jnc=mean(hq.sum.jnc), unq.sum.jnc=mean(unq.sum.jnc), sum.jnc=mean(sum.jnc),
        avg.err.5=mean(avg.err.5), avg.err.3=mean(avg.err.3),
        avg.low.5=mean(avg.low.5), avg.low.3=mean(avg.low.3),
        max.ovr.5=mean(max.ovr.5), max.ovr.3=mean(max.ovr.3),
        unq.ovr.5=mean(unq.ovr.5), unq.ovr.3=mean(unq.ovr.3),
        orf=mean(orf), spn=mean(type>-1),
        hq.bpt=mean(hq.bpt), hi.bpt=mean(hi.bpt)
    ), by=.(
           l1, l2, l3, topo, dst, d2a,
           gene.5=(cls.5 != "."), gene.3=(cls.3 != "."),
           ok.5=(art.5=="nil"), ok.3=(art.3=="nil")
       )
    ]
    res <- list(art.rte=art.rte, dup.rte=dup.rte, lig.rte=lig.rte,
                ts.rte.lo=ts.rte.lo, ts.rte.hi=ts.rte.hi, ts.lfc=ts.lfc,
                bs.rte=bs.rte, sl.rte=sl.rte,
                s2e.prox.rte=spn.prox.rte, s2e.dist.rte=spn.dist.rte,
                bpt.surv=bpt.surv)
    return(res)
}

.locus.stat <- function(bpt, ann) {
    tmp.x <- data.table(
        locus_id=factor(),
        enc_dist_lig=integer(),
        enc_prox_lig=integer(),
        spn_dist_lig=integer(),
        spn_dist_mot=integer(),
        spn_prox_lig=integer(),
        spn_prox_mot=integer()
    )
    tmp.5 <- bpt[,.(N=sum(sum.jnc)), by=.(locus_id=locus_id.5.1, l1, l2, l3, d2a)]
    tmp.5 <- dcast.data.table(tmp.5, locus_id~l1+l2+ifelse(d2a, "mot", "lig"), value.var = "N", fun.aggregate = sum)
    tmp.5 <- rbind(tmp.5, tmp.x, fill=TRUE)[,colnames(tmp.x),with=FALSE]
    setnames(tmp.5, c("locus_id", paste0(colnames(tmp.5)[-1], ".5")))
    tmp.3 <- bpt[,.(N=sum(sum.jnc)), by=.(locus_id=locus_id.3.1, l1, l2, l3, d2a)]
    tmp.3 <- dcast.data.table(tmp.3, locus_id~l1+l2+ifelse(d2a, "mot", "lig"), value.var = "N", fun.aggregate = sum)
    tmp.3 <- rbind(tmp.3, tmp.x, fill=TRUE)[,colnames(tmp.x),with=FALSE]
    setnames(tmp.3, c("locus_id", paste0(colnames(tmp.3)[-1], ".3")))
    setkey(tmp.5, locus_id)
    setkey(tmp.3, locus_id)
    tmp.53 <- merge(tmp.5, tmp.3, all=TRUE)
    tmp.53 <- tmp.53[ann$loci$locus_id]
    tmp.53 <- naTo0(tmp.53)
    setkey(tmp.53, locus_id)
    return(list(loc.surv=tmp.53))
}

#' @export
statReport <- function(bun, spl, bam, ann, par) {
    bpt.stat <- .breakpoint.stat(bun$bpt, bam$alig.reads / 2)
    loc.stat <- .locus.stat(bun$bpt, ann)
    stat <- list(bpt.stat=bpt.stat, loc.stat=loc.stat)
}
