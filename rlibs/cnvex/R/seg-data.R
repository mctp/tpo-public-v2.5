segData <- function(sel, tile, var, opts) {
    out <- tile[,c("arm", "target", "unmasked", "hq")]
    out$lr <- mcols(tile)[,sel$lr]
    snp <- var[,c("QUAL")]
    snp$AF <- mcols(var)[,sel$af]
    snp$DP <- mcols(var)[,sel$dp]
    snp <- snp[mcols(var)[,sel$pass]]
    hits <- .getHits(out, snp, opts)
    if (length(hits)>0) {
        bad <- ifelse(snp[queryHits(hits)]$AF < 0.5,
                   round(snp[queryHits(hits)]$AF  * snp[queryHits(hits)]$DP),
                   round((1-snp[queryHits(hits)]$AF) * snp[queryHits(hits)]$DP))
        tmp <- data.table(
            idx=subjectHits(hits),
            bad=bad,
            depth=snp[queryHits(hits)]$DP
        )
        setkey(tmp, idx)
        tmp <- tmp[J(1:length(out))]
        tmp <- tmp[,.(
            baf=sum(bad)/sum(depth),
            depth=sum(depth),
            n=length(na.omit(bad))
        ), by=idx]
        out$baf <- .smoothOutliers(tmp$baf, out, opts)
        out$baf.depth <- tmp$depth
        out$baf.n <- tmp$n
    } else {
        out$baf <- NA_real_
        out$baf.depth <- NA_real_
        out$baf.n <- NA_real_
    }
    return(out)
}
