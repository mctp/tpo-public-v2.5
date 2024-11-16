.blockSegmentsModel <- function(seg, tile, mod, opts) {
    ## add seg, C, K
    hits <- findOverlaps(seg, tile, select="all")
    tile.seg <- queryHits(hits)
    tile$seg <- tile.seg
    tile$C <- mod[tile.seg]$C
    tile$K <- mod[tile.seg]$K
    ## seperate arms
    arms <- split(tile, tile$arm)
    ## run blocking for each arm
    tmp <- foreach(arm=arms) %dopar% {
        .blockArmModel(arm, opts)
    }
    ## compute genome-wide unique indexes
    tmp <- paste(tile$arm, unlist(tmp))
    tmp <- as.integer(factor(tmp, levels=unique(tmp)))
    ## create segmentation ranges
    new.seg <- unname(unlist(range(split(tile, tmp))))
    return(new.seg)
}

.blockArmModel <- function(arm, opts) {
    ## find breakpoints
    seg1 <- which(diff(arm$seg) == 1)
    ## start blocking
    if (TRUE) {
        repeat({
            ## block segments below diff threshold
            best1 <- .segBlockArmModel(seg1, arm)
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

.segBlockArmModel <- function(seg0, arm) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
        ## compute breakpoint stats
        beg0 <- c(1, seg0+1)
        end0 <- c(seg0, length(arm))
        len0 <- diff(c(0, end0))
        idx0 <- rep(seq_along(beg0), len0)
        C0 <- split(arm$C, idx0)
        K0 <- split(arm$K, idx0)
        stat0 <- data.table(
            seg = seg0,
            C.diff = sapply(2:length(beg0), function(i) C0[[i]][1]!=C0[[i-1]][1]),
            K.diff = sapply(2:length(beg0), function(i) K0[[i]][1]!=K0[[i-1]][1])
        )
        stat0[is.na(C.diff), C.diff:=FALSE]
        stat0[is.na(K.diff), K.diff:=FALSE]
        ## pick weakest breakpoint to merge
        stat0 <- stat0[order(C.diff, K.diff)]
        seg1 <- stat0[(
            ## keep all breakpoints with same C and K
            C.diff == FALSE &
            K.diff == FALSE
        ), seg]
        best1 <- head(seg1, 1)
    }
    return(best1)
}
