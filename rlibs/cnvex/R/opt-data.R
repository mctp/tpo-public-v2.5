.opt.data <- function(mcnv, seg, opts) {
    var.seg <- findOverlaps(mcnv$var, seg, select="first", maxgap=opts$tile.shoulder-1)
    tile.seg <- findOverlaps(mcnv$tile, seg, select="first")
    if (length(mcnv$var)>0) {
        ## TODO: it is possible that this will return a 0-row table, fix.
        af <- .af.opt.data(mcnv$var, var.seg)
    } else {
        af <- data.table(seg=integer(0), idx=integer(0), SOURCE = character(0),
                         QUAL = numeric(0), AF = numeric(0), DP = integer(0))
    }
    lr <- .lr.opt.data(mcnv$tile, tile.seg)
    seg <- .seg.opt.data(seg, mcnv$tile, tile.seg)
    data <- list(lr=lr, af=af, seg=seg, stats=mcnv$stats)
    return(data)
}

.lr.opt.data <- function(tile, tile.seg) {
    lrt <- data.table(
        seg=tile.seg,
        lr=tile$lr,
        target=tile$target,
        nC=tile$nC,
        hq=tile$hq
    )
    return(lrt)
}

.af.opt.data <- function(var, var.seg) {
    snpt <- data.table(
        seg=var.seg,
        idx=seq_len(length(var)),
        as.data.table(mcols(var))
    )
    return(snpt)
}

.seg.opt.data <- function(seg, tile, tile.seg) {
    tmp <- as.data.table(mcols(tile))
    tmp$seg <- tile.seg
    tmp$len <- width(tile)
    seg.len <- width(seg)
    segt <- tmp[,.(
        len=seg.len[seg],
        gap=weighted.mean(gap, len),
        ntot=.N,
        nlr=sum(!is.na(lr)),
        naf=sum(baf.n),
        nhq=sum(hq)
    ),by=seg]
    setkey(segt, seg)
    return(segt)
}
