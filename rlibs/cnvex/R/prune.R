pruneSegments <- function(seg, tile, stats, opts) {
    ## add seg
    tile$seg <- queryHits(findOverlaps(seg, tile, select="all"))
    ## seperate arms
    arms <- split(tile, tile$arm)
    ## run prune for each arm
    tmp <- foreach(arm=arms) %dopar% {
        .pruneArm(arm, stats$sd.lr, stats$sd.baf, opts)
    }
    ## compute genome-wide unique indexes
    tmp <- paste(tile$arm, unlist(tmp))
    tmp <- as.integer(factor(tmp, levels=unique(tmp)))
    ## create segmentation ranges
    new.seg <- unname(unlist(range(split(tile, tmp))))
    return(new.seg)
}

.pruneArm <- function(arm, sd.lr, sd.baf, opts) {
    ## find breakpoints
    seg1 <- which(diff(arm$seg) == 1)
    ## start pruning
    if (opts$prune.len) {
        repeat({
            ## prune segments below diff threshold
            best1 <- .segPruneArmLen(seg1, arm, sd.lr, sd.baf, opts$prune.lr.lo.threshold, opts$prune.baf.lo.threshold,
                                     opts$prune.lr.len.penalty, opts$prune.baf.len.penalty)
            seg1 <- setdiff(seg1, best1)
            if (length(best1) == 0) {break}
        })
    }
    if (opts$prune.nvar) {
        repeat({
            ## require higher diff for segments with higher variance in normals
            best1 <- .segPruneArmNvar(seg1, arm, sd.lr, sd.baf, opts$prune.lr.lo.threshold, opts$prune.baf.lo.threshold,
                                      opts$prune.lr.nvar.penalty, opts$prune.baf.nvar.penalty)
            seg1 <- setdiff(seg1, best1)
            if (length(best1) == 0) {break}
        })
    }
    if (opts$prune.hq) {
        repeat({
            ## require higher diff for segments with a smaller proportion of hq tiles
            best1 <- .segPruneArmHq(seg1, arm, sd.lr, sd.baf, opts$prune.lr.lo.threshold, opts$prune.baf.lo.threshold,
                                    opts$prune.lr.hq.penalty, opts$prune.baf.hq.penalty)
            seg1 <- setdiff(seg1, best1)
            if (length(best1) == 0) {break}
        })
    }
    ## append segmentation
    bpt1 <- c(1, seg1 + 1)
    len1 <- diff(c(bpt1, length(arm)+1))
    idx1 <- rep(seq_along(bpt1), len1)
    return(idx1)
}

.segPruneArmLen <- function(seg0, arm, sd.lr, sd.baf, lr.lo.threshold, baf.lo.threshold, lr.len.penalty, baf.len.penalty) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
        ## compute breakpoint stats
        beg0 <- c(1, seg0+1)
        end0 <- c(seg0, length(arm))
        len0 <- diff(c(0, end0))
        idx0 <- rep(seq_along(beg0), len0)
        lr0 <- split(arm$lr, idx0)
        baf0 <- split(arm$baf, idx0)
        baf.n0 <- split(arm$baf.n, idx0)
        stat0 <- data.table(
            seg = seg0,
            lr.diff = sapply(2:length(beg0), function(i) absMedDiff(lr0[[i]], lr0[[i-1]])),
            baf.diff = sapply(2:length(beg0), function(i) absMedDiff(baf0[[i]], baf0[[i-1]])),
            lr.minlen = sapply(2:length(beg0), function(i) min(len0[i], len0[i-1])),
            baf.minlen = sapply(2:length(beg0), function(i) minSum(baf.n0[[i]], baf.n0[[i-1]]))
        )
        stat0[is.na(lr.diff), lr.diff:=0]
        stat0[is.na(baf.diff), baf.diff:=0]
        ## pick weakest breakpoint to merge
        stat0 <- stat0[order(lr.diff, baf.diff)]
        seg1 <- stat0[(
            ## keep all breakpoints above sd.lr and sd.baf penalty
            lr.diff <= sd.lr * (lr.lo.threshold + pmin(lr.len.penalty / lr.minlen, lr.len.penalty)) &
            baf.diff <= sd.baf * (baf.lo.threshold + pmin(baf.len.penalty / baf.minlen, baf.len.penalty))
        ), seg]
        best1 <- head(seg1, 1)
    }
    return(best1)
}

.segPruneArmNvar <- function(seg0, arm, sd.lr, sd.baf, lr.lo.threshold, baf.lo.threshold, lr.nvar.penalty, baf.nvar.penalty) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
        ## compute breakpoint stats
        beg0 <- c(1, seg0+1)
        end0 <- c(seg0, length(arm))
        len0 <- diff(c(0, end0))
        idx0 <- rep(seq_along(beg0), len0)
        lr0 <- split(arm$lr, idx0)
        baf0 <- split(arm$baf, idx0)
        nvar0 <- split(arm$nvar, idx0)
        stat0 <- data.table(
            seg = seg0,
            lr.diff = sapply(2:length(beg0), function(i) absMedDiff(lr0[[i]], lr0[[i-1]])),
            baf.diff = sapply(2:length(beg0), function(i) absMedDiff(baf0[[i]], baf0[[i-1]])),
            max.nvar = sapply(2:length(beg0), function(i) max(mean(nvar0[[i]]), mean(nvar0[[i-1]])))
        )
        stat0[is.na(lr.diff), lr.diff:=0]
        stat0[is.na(baf.diff), baf.diff:=0]
        ## pick weakest breakpoint to merge
        stat0 <- stat0[order(lr.diff, baf.diff)]
        seg1 <- stat0[(
            ## keep all breakpoints above sd.lr and sd.baf penalty
            lr.diff <= sd.lr * (lr.lo.threshold + lr.nvar.penalty * sigmoid(max.nvar,max(quantile(max.nvar,0.95),median(max.nvar)+1.5*mad(max.nvar)),median(max.nvar)+1.5*mad(max.nvar))) &
            baf.diff <= sd.baf * (baf.lo.threshold + baf.nvar.penalty * sigmoid(max.nvar,max(quantile(max.nvar,0.95),median(max.nvar)+1.5*mad(max.nvar)),median(max.nvar)+1.5*mad(max.nvar)))
        ), seg]
        best1 <- head(seg1, 1)
    }
    return(best1)
}

.segPruneArmHq <- function(seg0, arm, sd.lr, sd.baf, lr.lo.threshold, baf.lo.threshold, lr.hq.penalty, baf.hq.penalty) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
        ## compute breakpoint stats
        beg0 <- c(1, seg0+1)
        end0 <- c(seg0, length(arm))
        len0 <- diff(c(0, end0))
        idx0 <- rep(seq_along(beg0), len0)
        lr0 <- split(arm$lr, idx0)
        baf0 <- split(arm$baf, idx0)
        hq0 <- split(arm$hq, idx0)
        stat0 <- data.table(
            seg = seg0,
            lr.diff = sapply(2:length(beg0), function(i) absMedDiff(lr0[[i]], lr0[[i-1]])),
            baf.diff = sapply(2:length(beg0), function(i) absMedDiff(baf0[[i]], baf0[[i-1]])),
            min.hq = sapply(2:length(beg0), function(i) min(sum(hq0[[i]])/length(hq0[[i]]), sum(hq0[[i-1]])/length(hq0[[i-1]])))
        )
        stat0[is.na(lr.diff), lr.diff:=0]
        stat0[is.na(baf.diff), baf.diff:=0]
        ## pick weakest breakpoint to merge
        stat0 <- stat0[order(lr.diff, baf.diff)]
        seg1 <- stat0[(
            ## keep all breakpoints above sd.lr and sd.baf penalty
            lr.diff <= sd.lr * (lr.lo.threshold + lr.hq.penalty * sigmoid(min.hq,0,1)) &
            baf.diff <= sd.baf * (baf.lo.threshold + baf.hq.penalty * sigmoid(min.hq,0,1))
        ), seg]
        best1 <- head(seg1, 1)
    }
    return(best1)
}
