.modelCnv <- function(sample, cnv, pool, gobj, opts) {
    sel <- .sample.switch(sample, opts)
    tile <- cnv$tile[,c("arm", "target", "unmasked", "hq", "gap", "nC")]
    tile$lr <- mcols(cnv$tile)[,sel$lr]
    if (!is.null(pool)) {
        if (cnv$sex=="male") {
            tile$nvar <- pool$male$nvar
        } else {
            tile$nvar <- pool$female$nvar
        }
    }
    if (length(cnv$var)>0) {
        col.sel <- c("SOURCE", "QUAL")
        var <- cnv$var[,col.sel]
        var$AF <- mcols(cnv$var)[,sel$af]
        var$DP <- mcols(cnv$var)[,sel$dp]
        var <- var[mcols(cnv$var)[,sel$pass]]
    } else {
        var <- cnv$var
        var$SOURCE <- character(0)
        var$QUAL <- integer(0)
        var$AF <- numeric(0)
        var$DP <- integer(0)
    }
    var <- .drop.var(tile, var, opts)
    tile <- .add.baf(tile, var, opts)
    stats <- .modelStats(tile, var)
    cnv <- list(sample=sample, tile=tile, var=var, stats=stats, sex=cnv$sex)
    return(cnv)
}

.drop.var <- function(tile, var, opts) {
    max.per.tile <- opts$baf.max.per.tile
    if (!is.null(max.per.tile) && is.finite(max.per.tile)) {
        hits <- findOverlaps(var, tile)
        var.df <- mcols(var[queryHits(hits)])
        hits.dt <- as.data.table(cbind(var.df[,c("QUAL"),drop=FALSE], hits))
        ## pick variant with highest QUAL
        sel <- hits.dt[order(-QUAL),
                       .(queryHits=head(queryHits, max.per.tile)), subjectHits]
        var.sel <- var[sort(sel$queryHits)]
    } else {
        var.sel <- var
    }
    return(var.sel)
}

.add.baf <- function(tile, var, opts) {
    hits <- .getHits(tile, var, opts)
    if (length(hits)>0) {
        bad <- ifelse(var[queryHits(hits)]$AF < 0.5,
                   round(   var[queryHits(hits)]$AF  * var[queryHits(hits)]$DP),
                   round((1-var[queryHits(hits)]$AF) * var[queryHits(hits)]$DP))
        tmp <- data.table(
            idx=subjectHits(hits),
            bad=bad,
            depth=var[queryHits(hits)]$DP
        )
        setkey(tmp, idx)
        tmp <- tmp[J(1:length(tile))]
        tmp <- tmp[,.(
            baf=sum(bad)/sum(depth),
            depth=sum(depth),
            n=length(na.omit(bad))
        ), by=idx]
        tile$baf <- .smoothOutliers(tmp$baf, tile, opts)
        tile$baf.depth <- tmp$depth
        tile$baf.n <- tmp$n
    } else {
        tile$baf <- NA_real_
        tile$baf.depth <- NA_real_
        tile$baf.n <- NA_real_
    }
    return(tile)
}

.modelStats <- function(tile, var) {
    ## exclude low-quality regions
    tile.hq <- tile[tile$hq]
    hq.lr <- na.omit(tile.hq$lr)
    if (length(hq.lr)>3) {
        sd.lr <- estimateSd(hq.lr)
    } else {
        sd.lr <- 0
    }
    hq.baf <- na.omit(tile.hq$baf)
    if (length(hq.baf)>3) {
        sd.baf <- estimateSd(hq.baf)
    } else {
        sd.baf <- 0
    }
    stats <- list(sd.lr=sd.lr, sd.baf=sd.baf)
    return(stats)
}

.sample.switch <- function(sample, opts) {
    if (sample=="tumor") {
        af <- "t.AF"
        dp <- "t.DP"
        pass <- "t.PASS"
        if (opts$lr.tumor=="pool") {
            lr <- "tp.lr"
        } else if (opts$lr.tumor=="pair") {
            lr <- "tn.lr"
        } else {
            lr <- "t.lr"
        }
    } else if (sample=="normal") {
        af <- "n.AF"
        dp <- "n.DP"
        pass <- "n.PASS"
        if (opts$lr.normal=="pool") {
            lr <- "np.lr"
        } else {
            lr <- "n.lr"
        }
    }
    out <- list(sample=sample, af=af, dp=dp, pass=pass, lr=lr)
    return(out)
}
