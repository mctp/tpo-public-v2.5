#' @export
robustImport <- function(fn, seqi, ...) {
    .robust.import(fn, seqi, ...)
}

#' @export
getGobj <- function(genome, fasta, refs=TRUE) {
    gobj <- .gobj(genome, fasta, refs)
    return(gobj)
}

#' @export
getOpts <- function(settings, opts=list()) {
    get.opts(settings, opts)
}

#' @export
importFeat <- function(gtf, feat, gobj) {
    ann <- .robust.import(gtf, gobj$seqi)
    feat <- ann[ann$type==feat]
    return(feat)
}

#' @export
importVcf <- function(vcf, tidx, nidx, gobj, opts) {
    func <- getVcfImporter(opts$caller)
    var <- unlist(GRangesList(unname(lapply(names(vcf), function(name) {
        fn <- vcf[[name]]
        v <- func(fn, tidx, nidx, gobj, opts)
        mcols(v)$SOURCE <- name
        return(v)
    }))))
    return(var)
}

#' @export
importBam <- function(capture, t.bam, n.bam, gobj, opts) {
    if (!is.null(capture)) {
        target.fun <- getTargetTiles
    } else {
        target.fun <- getGenomeTiles
    }
    tile <- target.fun(capture, gobj, opts)

    ## normalize by sequencing depth
    if (!is.null(t.bam)) {
        t.cov.raw <- .runMosdepthTile(t.bam, tile, gobj$fasta, opts$cores)
    } else {
        t.cov.raw <- NA_real_
    }
    if (!is.null(n.bam)) {
        n.cov.raw <- .runMosdepthTile(n.bam, tile, gobj$fasta, opts$cores)
    } else {
        n.cov.raw <- NA_real_
    }
    mcols(tile) <- cbind(mcols(tile), cbind(t.cov.raw, n.cov.raw))
    return(tile)
}

#' @export
createCnv <- function(vcf, capture, t.bam, n.bam, gobj, opts) {
    if (!is.null(t.bam) && !is.null(n.bam)) {
        tidx <- 1
        nidx <- 2
    } else if (!is.null(t.bam) && is.null(n.bam)) {
        tidx <- 1
        nidx <- NULL
    } else if (is.null(t.bam) && !is.null(n.bam)) {
        tidx <- NULL
        nidx <- 1
    }
    var <- importVcf(vcf, tidx, nidx, gobj, opts)
    tile <- importBam(capture, t.bam, n.bam, gobj, opts)
    cnv <- list(tile=tile, var=var)
    return(cnv)
}

#' @export
upgradeCnv <- function(cnv, opts) {
    if (is.null(cnv$var$population.af)) {
        ## for CNVX objects done prior to CNVEX_POPAF implementation
        ## set to default value if CNVEX_POPAF is not provided
        ## see: callers.R importSentieon
        cnv$var$population.af <- 0.5
    }
    return(cnv)
}

#' @export
addHqTile <- function(cnv, opts) {
    cnv$tile <- .annotateHQ(cnv$tile, opts)
    return(cnv)
}

#' @export
addPassVariant <- function(cnv, opts) {
    if ("t.GT" %in% names(mcols(cnv$var))) {
        cnv$var$t.PASS <- passTumor(cnv$tile, cnv$var, opts)
    }
    if ("n.GT" %in% names(mcols(cnv$var))) {
        cnv$var$n.PASS <- passNormal(cnv$tile, cnv$var, opts)
    }
    return(cnv)
}

#' @export
addSex <- function(cnv, gobj, opts) {
    if (is.null(opts$sex)) {
        cnv$sex <- .detect.sex(cnv$var, cnv$tile, gobj, opts)
    } else {
        cnv$sex <- opts$sex
    }
    cnv$tile$nC <- .get.sexcopy(cnv$sex, cnv$tile, gobj)
    return(cnv)
}

#' @export
addNormCoverage <- function(cnv, opts) {
    t.cov.raw <- cnv$tile$t.cov.raw
    t.cov <- .normCoverage(t.cov.raw, cnv$tile, opts)
    n.cov.raw <- cnv$tile$n.cov.raw
    if (!is.null(n.cov.raw)) {
        n.cov <- .normCoverage(n.cov.raw, cnv$tile, opts)
    } else {
        n.cov <- NA_real_
    }
    mcols(cnv$tile) <- cbind(mcols(cnv$tile), cbind(t.cov, n.cov))
    return(cnv)
}

#' @export
addLogRatio <- function(cnv, pool, opts) {
    nas <- rep(NA_real_, length(cnv$tile))
    cnv$tile$t.lr <- nas
    cnv$tile$n.lr <- nas
    cnv$tile$tn.lr <- nas
    cnv$tile$tp.lr <- nas
    cnv$tile$np.lr <- nas
    if (!all(is.na(cnv$tile$t.cov))) {
        cnv$tile$t.lr <- .polishLogRatio(.rawLogRatio(cnv$tile$t.cov, 1, opts), cnv$tile, opts)
    }
    if (!all(is.na(cnv$tile$n.cov))) {
        cnv$tile$n.lr <- .polishLogRatio(.rawLogRatio(cnv$tile$n.cov, 1, opts), cnv$tile, opts)
    }
    if (!all(is.na(cnv$tile$t.cov)) && !all(is.na(cnv$tile$n.cov))) {
        cnv$tile$tn.lr <- .polishLogRatio(.rawLogRatio(cnv$tile$t.cov, cnv$tile$n.cov, opts), cnv$tile, opts)
    }
    if (!all(is.na(cnv$tile$t.cov)) && !is.null(pool)) {
        if (!all(is.na(cnv$tile$n.cov))) {
            n.cov <- cnv$tile$n.cov
        } else {
            n.cov <- nas
        }
        cnv$tile$tp.lr <- .polishLogRatio(.poolLogRatio(cnv$tile$t.cov, n.cov, cnv$sex, pool, opts), cnv$tile, opts)
    }
    if (!all(is.na(cnv$tile$n.cov)) && !is.null(pool)) {
        cnv$tile$np.lr <- .polishLogRatio(.poolLogRatio(cnv$tile$n.cov, nas, cnv$sex, pool, opts), cnv$tile, opts)
    }
    return(cnv)
}

#' @export
modelCnv <- function(sample, cnv, pool, gobj, opts) {
    mcnv <- .modelCnv(sample, cnv, pool, gobj, opts)
    return(mcnv)
}

#' @export
modelSeg <- function(mcnv, pool, opts) {
    seg <- jointSegment(mcnv$tile, opts)
    seg <- pruneSegments(seg, mcnv$tile, mcnv$stats, opts)
    return(seg)
}

#' @export
modelOpt <- function(mcnv, seg, nogrid, nofine, opts) {
    log_info("Running modelOpt nogrid={nogrid} nofine={nofine}")
    log_debug("tiles:{length(mcnv$tile)} variants:{length(mcnv$var)} segments:{length(seg)}")
    ## get data
    data <- .opt.data(mcnv, seg, opts)
    ## optimize models
    opt <- .opt(data, nogrid, nofine, opts)
    return(opt)
}

#' @export
modelGetOne <- function(mcnv, seg, purity, ploidy, opts) {
    log_info("Running modelGet for purity={purity} and ploidy={ploidy}")
    log_debug("tiles:{length(mcnv$tile)} variants:{length(mcnv$var)} segments:{length(seg)}")
    data <- .opt.data(mcnv, seg, opts)
    ## only.hq == FALSE, because we want to fit all segments
    mod <- .get.model(data, purity, ploidy, FALSE, opts)
    return(mod)
}

#' @export
modelGetEval <- function(mcnv, seg, eval, n=99, opts) {
    n <- min(n, nrow(eval))
    log_info("Running modelGetEval for top {n} out of {nrow(eval)} models")
    log_debug("tiles:{length(mcnv$tile)} variants:{length(mcnv$var)} segments:{length(seg)}")
    data <- .opt.data(mcnv, seg, opts)
    mods <- foreach(i=seq_len(n)) %dopar% {
        ## only.hq == FALSE, because we want to fit all segments
        .get.model(data, eval$p[i], eval$P[i], FALSE, opts)
    }
    return(mods)
}

#' @export
getBlocks <- function(digest) {
    blk <- unlist(GRangesList(sapply(split(digest$seg, digest$fit$block), range)))
    return(blk)
}

#' @export
modelDigest <- function(mcnv, seg, mod, purity, ploidy, eval=list(), opts, log) {

    ## add block id to model
    blk <- .blockSegmentsModel(seg, mcnv$tile, mod, opts)
    mod$block <- subjectHits(findOverlaps(seg, blk, select="all"))

    ## this internal function is for backwards compatibility
    .format.fit <- function(mod) {
        fit <- mod[,.(seg, block, C, K, lr, tL, aL, d=aD, anom, mse, nlr, naf, len, sC, nC)]
        setkey(fit, seg)
        return(fit)
    }
    fit <- .format.fit(mod)

    ## merge data
    mcnv$purity <- eval$p
    mcnv$ploidy <- eval$P
    mcnv$seg <- seg
    mcnv$fit <- fit
    mcnv$eval <- eval
    mcnv$opts <- opts
    mcnv$log <- log
    ## add segment
    mcnv$var$seg <- findOverlaps(mcnv$var, seg, maxgap=opts$tile.shoulder-1, select="first")
    mcnv$tile$seg <- findOverlaps(mcnv$tile, seg, select="first")
    ## add fit tile
    mcnv$tile$sC <- fit[mcnv$tile$seg]$sC
    mcnv$tile$C <- fit[mcnv$tile$seg]$C
    mcnv$tile$K <- fit[mcnv$tile$seg]$K
    mcnv$tile$seg.lr <- fit[mcnv$tile$seg]$lr
    mcnv$tile$seg.len <- fit[mcnv$tile$seg]$len
    mcnv$tile$seg.nlr <- fit[mcnv$tile$seg]$nlr
    mcnv$tile$seg.naf <- fit[mcnv$tile$seg]$naf
    ## add fit var
    mcnv$var$sC <- fit[mcnv$var$seg]$sC
    mcnv$var$C <- fit[mcnv$var$seg]$C
    mcnv$var$K <- fit[mcnv$var$seg]$K
    mcnv$var$seg.lr <- fit[mcnv$var$seg]$lr
    mcnv$var$seg.len <- fit[mcnv$var$seg]$len
    mcnv$var$seg.nlr <- fit[mcnv$var$seg]$nlr
    mcnv$var$seg.naf <- fit[mcnv$var$seg]$naf
    return(mcnv)
}

#' @export
createPoolData <- function(cnv.fns, gobj, opts) {
    pd <- .importPoolData(cnv.fns, gobj, opts)
    return(pd)
}

#' @export
createPool <- function(pd, opts) {
    pool <- .createPool(pd, opts)
    return(pool)
}




#### CINs

#' @export
cinLST <- function(digest,gobj,opts) {
  lst <- .cinLST(digest,gobj,opts)
  return(lst)
}

#' @export
cinNtAI <- function(digest,gobj,opts) {
  ntai <- .cinNtAI(digest,gobj,opts)
  return(ntai)
}

#' @export
cinNLOH <- function(digest,gobj,opts) {
  nloh <- .cinNLOH(digest,gobj,opts)
  return(nloh)
}

#' @export
cinWgii <- function(digest) {
  LST <- .cinWGii(digest)
  return(LST)
}

#' @export
cinFG <- function(DIGEST, PLOIDY) {
  fg <- .cinFG(DIGEST,PLOIDY)
  return(fg)
}

#' @export
cinSI <- function(DIGEST, GOBJ, OPTS) {
  si <- .cinSI(DIGEST, GOBJ, OPTS)
  return(si)
}
