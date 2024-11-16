jointSegment <- function(tile, opts) {
    ## segment each arm
    arms <- split(tile, tile$arm)
    tmp <- foreach(arm=arms) %dopar% {
        .jointSegArm(arm, opts)
    }
    ## compute genome-wide unique indexes
    tmp <- paste(tile$arm, unlist(tmp))
    tmp <- as.integer(factor(tmp, levels=unique(tmp)))
    ## create segmentation ranges
    seg <- unname(unlist(range(split(tile, tmp))))
    return(seg)
}

.jointSegArm <- function(arm, opts) {
    ## initial segmentation
    seg1 <- .runJointSeg(arm, opts$seg.method, opts$tile.width, opts$baf.max.eff.dp, opts$seg.cbs.lr,
                         opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.cbs.weighted, opts$seg.only.target)
    ## append segmentation
    bpt1 <- c(1, seg1 + 1)
    len1 <- diff(c(bpt1, length(arm)+1))
    idx1 <- rep(seq_along(bpt1), len1)
    return(idx1)
}

.runJointSeg <- function(arm, method, tile.width, baf.max.eff.dp, cbs.lr, cbs.baf,
                         rbs.selection, cbs.weighted, only.target) {
    arm.lr <- arm$lr
    arm.baf <- arm$baf
    if (only.target) {
        arm.lr[!arm$target] <- NA_real_
        arm.baf[!arm$target] <- NA_real_
    }
    ## make sure we have enough points to segment
    len.min <- 2
    seg0 <- integer()
    if (method %in% c("RBS", "DynamicProgramming")) {
        n.points.lr <- sum(!is.na(arm.lr))
        n.points.baf <- sum(!is.na(arm.baf))
        if (n.points.lr >= 2*len.min) {
            if (!is.na(tile.width)) {
                tile.per.mb <- 1e6 / tile.width
                opt.K <- ceiling(n.points.lr / tile.per.mb)
            } else {
                opt.K <- max(2, floor(n.points.lr / (2*len.min)))
            }
            if (n.points.baf >= 2 *len.min) {
                data <- cbind(arm.lr, arm.baf)
            } else {
                data <- arm.lr
            }
            seg0 <- suppressWarnings(sort(unique(jointSeg(
                data, method=method, modelSelectionMethod=rbs.selection, K=opt.K)$bestBkp)))
        }
    } else if (method=="CBS") {
        ## add data weight
        if (cbs.weighted) {
            arm.lr.weight <- arm$unmasked
            arm.baf.weight <- sqrt(pmin(arm$baf.depth, baf.max.eff.dp))
        } else {
            arm.lr.weight <- NULL
            arm.baf.weight <- NULL
        }
        n.points.lr <- sum(!is.na(arm.lr))
        if (n.points.lr >= 2*len.min) {
            seg0.lr <- .runCBS(arm.lr, arm.lr.weight, args=cbs.lr)
        } else {
            seg0.lr <- c()
        }
        n.points.baf <- sum(!is.na(arm.baf))
        if (n.points.baf >= 2*len.min) {
            seg0.baf <- .runCBS(arm.baf, arm.baf.weight, args=cbs.baf)
        } else {
            seg0.baf <- c()
        }
        seg0 <- unique(sort(c(seg0.lr, seg0.baf)))
    } else {
        stop("Joint segmentation method not supported.")
    }
    return(seg0)
}
