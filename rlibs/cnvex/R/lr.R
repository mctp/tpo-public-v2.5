.smoothOutliers <- function(y, tile, opts) {
    ## remove gross outliers
    y.tar.n <- sum(!is.na(y[tile$target]))
    if (y.tar.n>1) {
        y[tile$target] <- .runSmooth(y[tile$target], tile[tile$target])
    }
    y.off.n <- sum(!is.na(y[!tile$target]))
    if (y.off.n>1) {
        y[!tile$target] <- .runSmooth(y[!tile$target], tile[!tile$target])
    }
    return(y)
}

.rawLogRatio <- function(t.cov, n.cov, opts) {
    lr.raw <- log2(t.cov / n.cov)
    lr.raw <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
    return(lr.raw)
}

.scaleLogRatio <- function(lr, gt, opts) {
    weight <- width(gt)
    lr.mean.tgt <- weighted.mean(lr[gt$target & gt$hq], weight[gt$target & gt$hq], na.rm=TRUE)
    lr.mean.off <- weighted.mean(lr[!gt$target & gt$hq], weight[!gt$target & gt$hq], na.rm=TRUE)
    lr[gt$target] <- lr[gt$target] - lr.mean.tgt
    lr[!gt$target] <- lr[!gt$target] - lr.mean.off
    return(lr)
}

.smoothLogRatioInner <- function(lr, tile, opts) {
    lr.smooth <- unname(unlist(lapply(split(lr, tile$arm), function(y) {
        if (length(y) >= opts$lr.smooth.window) {
            y.smooth <- suppressWarnings(
                hybrid.filter(y, opts$lr.smooth.window, method=c("MEAN", "MED"))$level[["MEAN"]]
            )
        } else {
            y.smooth <- y
        }
        return(y.smooth)
    })))
    return(lr.smooth)
}

.smoothLogRatio <- function(lr, tile, opts) {
    lr.smooth <- lr
    if (opts$lr.smooth=="hybrid") {
        ## hybrid smooth only on targeted
        lr.smooth[tile$target] <-
            .smoothLogRatioInner(lr.smooth[tile$target], tile[tile$target], opts)
    } else if (opts$lr.smooth=="outlier") {
        lr.smooth <- .smoothOutliers(lr.smooth, tile, opts)
    }
    return(lr.smooth)
}

.polishLogRatio <- function(lr, tile, opts) {
    lr <- .scaleLogRatio(lr, tile, opts)
    if (opts$bias.gc.correct) {
        lr <- .correctBias(lr, tile$target, tile$hq, "gc", tile$gc, opts)
    }
    lr <- .smoothLogRatio(lr, tile, opts)
    return(lr)
}
