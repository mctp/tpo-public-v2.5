.getTargetHits <- function(tile, snp, opts) {
    hits <- findOverlaps(snp, tile, maxgap = opts$tile.shoulder-1)
    hits <- hits[tile[subjectHits(hits)]$target] # prefer assignment to target
    hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    return(hits)
}

.getGenomeHits <- function(tile, snp, opts) {
    hits <- findOverlaps(snp, tile, maxgap = opts$tile.shoulder-1)
    return(hits)
}

.getHits <- function(tile, snp, opts) {
    if (any(!tile$target)) {
        hits <- .getTargetHits(tile, snp, opts)
    } else {
        hits <- .getGenomeHits(tile, snp, opts)
    }
    return(hits)
}

.addBafTile <- function(tile, var, opts) {
    ## filter variants to SNPs
    snp <- var[var$t.PASS]
    hits <- .getHits(tile, snp, opts)
    if (length(hits)>0) {
        bad <- ifelse(snp[queryHits(hits)]$t.AF < 0.5,
                   round(   snp[queryHits(hits)]$t.AF  * snp[queryHits(hits)]$t.DP),
                   round((1-snp[queryHits(hits)]$t.AF) * snp[queryHits(hits)]$t.DP))
        tmp <- data.table(
            idx=subjectHits(hits),
            bad=bad,
            depth=snp[queryHits(hits)]$t.DP
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
