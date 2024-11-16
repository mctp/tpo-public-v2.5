.annotateTiles <- function(genome.tile, gobj) {

    ## blacklist
    bl <- GenomicRanges::reduce(granges(gobj$refs$blacklist))

    tmp <- findOverlaps(genome.tile, bl)
    tmp <- data.table(
        tile=queryHits(tmp),
        blacklist=width(pintersect(genome.tile[queryHits(tmp)], bl[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(blacklist), blacklist:=0]
    tmp <- tmp[,.(blacklist=sum(blacklist)),by=tile]
    genome.tile$blacklist <- tmp$blacklist / width(genome.tile)

    ## masking
    tmp <- findOverlaps(genome.tile, gobj$refs$unmask.strict)
    tmp <- data.table(
        tile=queryHits(tmp),
        unmasked=width(pintersect(genome.tile[queryHits(tmp)], gobj$refs$unmask.strict[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(unmasked), unmasked:=0]
    tmp <- tmp[,.(unmasked=sum(unmasked)),by=tile]
    genome.tile$unmasked <- tmp$unmasked / width(genome.tile)

    ## giab
    giab <- GenomicRanges::reduce(granges(gobj$refs$giab.isdifficult))

    tmp <- findOverlaps(genome.tile, giab)
    tmp <- data.table(
        tile=queryHits(tmp),
        giab=width(pintersect(genome.tile[queryHits(tmp)], giab[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(giab), giab:=0]
    tmp <- tmp[,.(giab=sum(giab)),by=tile]
    genome.tile$giab.difficults <- tmp$giab / width(genome.tile)

    ## GC content
    tmp <- getSeq(gobj$seq, genome.tile)
    genome.tile$gc <- letterFrequency(tmp, "GC", as.prob=TRUE)[,1]

    ## has assembly gaps
    genome.tile$gap <- letterFrequency(tmp, "N", as.prob=TRUE)[,1]

    ## cytobands
    tmp <- findOverlaps(genome.tile, gobj$refs$cytoband, select="first")
    genome.tile$cytoband <- gobj$refs$cytoband[tmp]$name

    ## order arms
    genome.tile <- sort(genome.tile)
    genome.tile$arm <- paste0(seqnames(genome.tile), str_sub(genome.tile$cytoband, 1, 1))
    genome.tile$arm <- factor(genome.tile$arm, unique(genome.tile$arm), ordered=TRUE)

    ## TODO: BiasCorrection
    ## ## methylation
    ## tmp <- findOverlaps(genome.tile, meth)
    ## tmp <- data.table(
    ##     tile=queryHits(tmp),
    ##     meth=meth$score[subjectHits(tmp)]/100
    ## )
    ## setkey(tmp, tile)
    ## tmp <- tmp[J(seq_along(genome.tile))]
    ## tmp <- tmp[,mval:=log2((meth+1)/(2-meth))]
    ## tmp <- tmp[,.(meth=mean(meth)),by=tile]
    ## genome.tile$meth <- (tmp$meth/100)

    ## TODO: BiasCorrection
    ## ## replication time
    ## tmp <- findOverlaps(genome.tile, gobj$refs$repliseq)
    ## tmp <- data.table(
    ##     tile=queryHits(tmp),
    ##     repliseqcord=subjectHits(tmp),
    ##     rept=repliseq$score[subjectHits(tmp)]
    ## )
    ## setkey(tmp, tile)
    ## tmp <- tmp[J(seq_along(cnv$tile))]
    ## tmp <- tmp[,.(rept=mean(rept)),by=tile]
    ## tmp <- tmp[rept<=0,rept:=NA_real_]
    ## genome.tile$reptime <- (tmp$rept/100)

    return(genome.tile)
}

.annotateHQ <- function(tile, opts) {
    ## define HQ tile
    t.cov.raw <- tile$t.cov.raw
    t.cov.raw[is.na(t.cov.raw)] <- 0
    n.cov.raw <- tile$n.cov.raw
    n.cov.raw[is.na(n.cov.raw)] <- 0
    tile.totalcov <- (t.cov.raw + n.cov.raw) * width(tile)
    tile$hq <-
        tile$gap <= opts$tile.hq.max.gap &
        tile$unmasked >= opts$tile.hq.min.unmasked &
        tile$blacklist <= opts$tile.hq.max.blacklist &
        tile$giab.difficults <= opts$tile.hq.max.giab.difficults & ## TODO: this threshold needs tuning
        tile.totalcov >= opts$tile.hq.min.totalcov
    return(tile)
}

.tilegenome <- function(seqi, tile.width) {
    tile <- tileGenome(seqi, tilewidth=tile.width, cut.last.tile.in.chrom=TRUE)
    return(tile)
}

getGenomeTiles <- function(tgt.fn, gobj, opts) {
    ## tile genome
    genome.tile <- .tilegenome(gobj$seqi, opts$tile.width)
    genome.tile$target <- TRUE
    genome.tile <- .annotateTiles(genome.tile, gobj)
    return(genome.tile)
}

getTargetTiles <- function(tgt.fn, gobj, opts) {
    ## tile targets
    tgt <- .robust.import(tgt.fn, gobj$seqi)
    tgt$name <- NULL
    tgt$score <- NULL
    tgt$target <- TRUE
    gap <- keepSeqlevels(gaps(tgt), seqlevels(gobj$seqi), pruning.mode="coarse")
    gap <- gap[strand(gap)=="*"]
    gap <- gap[width(gap)>2*opts$tile.shoulder+opts$tile.min.gap]
    gap <- gap[!(gap %over% gobj$refs$centromere)]
    gap <- gap-opts$tile.shoulder
    gap$target <- FALSE
    target.tile <- sort(c(tgt, gap))
    target.tile <- .annotateTiles(target.tile, gobj)
    return(target.tile)
}
