.poolIcaDenoise <- function(t.cov, n.cov, sex, gpool, opts) {
    sspool <- .sexSpecificPool(gpool, sex)

    S.aut <- sspool$projection.aut[,seq(1,min(ncol(sspool$projection.aut),opts$pool.n.comp))]
    S.sex <- sspool$projection.sex[,seq(1,min(ncol(sspool$projection.sex),opts$pool.n.comp))]

    # adjust tumor values
    if (sex=="female") {
        t.cov <- t.cov * sspool$adjust.coef
    }

    lr <- log2(t.cov / sspool$cov.med)
    lr <- .correctBias(lr, gpool$target, gpool$hq, "gc", gpool$gc, opts)

    x.aut <- .fixLogRatioBounds(lr[(sspool$filter&(!sspool$sex.chr))])
    x.sex <- .fixLogRatioBounds(lr[(sspool$filter)])

    ## denoise step
    x.aut <- x.aut - as.vector(tcrossprod(x.aut %*% t(ginv(S.aut)), S.aut))
    x.sex <- x.sex - as.vector(tcrossprod(x.sex %*% t(ginv(S.sex)), S.sex))

    lr[(sspool$filter&(!sspool$sex.chr))] <- x.aut
    lr[(sspool$filter&(sspool$sex.chr))] <- x.sex[sspool$sex.chr[sspool$filter]]
    lr[!(sspool$filter)] <- NA_real_
    return(lr)
}

.poolPcaDenoise <- function(t.cov, n.cov, sex, gpool, opts) {
    sspool <- .sexSpecificPool(gpool, sex)

    ## adjust tumor values
    if (sex=="female") {
        t.cov <- t.cov * sspool$adjust.coef
    }

    lr <- log2(t.cov / sspool$cov.med)
    lr <- .correctBias(lr, gpool$target, gpool$hq, "gc", gpool$gc, opts)

    x.aut <- .fixLogRatioBounds(lr[(sspool$filter&(!sspool$sex.chr))])
    x.sex <- .fixLogRatioBounds(lr[(sspool$filter)])

    ## denoise step
    P.aut <- sspool$projection.aut[,seq(1,min(ncol(sspool$projection.aut),opts$pool.n.comp))]
    x.aut <- x.aut - as.vector(tcrossprod(x.aut %*% P.aut, P.aut))

    if (gpool$method!="pca-nosex") {
        P.sex <- sspool$projection.sex[,seq(1,min(ncol(sspool$projection.sex),opts$pool.n.comp))]
        x.sex <- x.sex - as.vector(tcrossprod(x.sex %*% P.sex, P.sex))
    }

    lr[(sspool$filter&(!sspool$sex.chr))] <- x.aut
    lr[(sspool$filter&(sspool$sex.chr))] <- x.sex[sspool$sex.chr[sspool$filter]]
    lr[!(sspool$filter)] <- NA_real_

    return(lr)
}

.poolLogRatio <- function(t.cov, n.cov, sex, pool, opts) {
    ## check for normal coverage
    if (all(is.na(n.cov))) {
        if (sex == "male") {
            n.cov <- pool$male$cov.med
        } else {
            n.cov <- pool$female$cov.med
        }
    } else {
        n.cov <- n.cov
    }
    ## pool denoised log-ratio
    if (pool$method %in% c("pca", "pca-nosex")) {
        lr.pool <- .poolPcaDenoise(t.cov, n.cov, sex, pool, opts)
    } else if (pool$method %in% c("ica", "ica-nosex")) {
        lr.pool <- .poolIcaDenoise(t.cov, n.cov, sex, pool, opts)
    }
    ## pool median log-ratio
    lr.pair <- .rawLogRatio(t.cov, n.cov, opts)
    ## hi.nvar threshold only on autosomes, same for both sex
    hi.nvar <- quantile(pool$male$nvar[!pool$sex.chr], opts$pool.hi.nvar)
    ## for tiles with high nvar (e.g. 95th percentile) use pool-median log ratio
    if (sex == "male") {
        lr <- ifelse(pool$male$nvar > hi.nvar, lr.pair, lr.pool)
    } else {
        lr <- ifelse(pool$female$nvar > hi.nvar, lr.pair, lr.pool)
    }
    ## TODO: check this carefully
    lr <- ifelse(is.na(lr), lr.pair, lr)
    return(lr)
}
