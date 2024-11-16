kpGetOpts <- function(settings, opts=list()) {
    kp.opts <- get.opts(settings, opts)
}

kpInit <- function(genome, kp.opts, zoom=NULL) {
    ## PARAMS
    pp <- getDefaultPlotParams(plot.type = 4)
    pp <- modifyList(pp, kp.opts)
    ## CYTOBANDS
    if (pp$cyto.type=="arm") {
        color.table <- getCytobandColors()
        color.table["gpos"] <- pp$cyto.pcol
        color.table["gneg"] <- pp$cyto.qcol
        cyto <- .kpArm(genome)
        cyto$gieStain <- ifelse(str_sub(cyto$name, -1)=="p", "gpos", "gneg")
    } else if (pp$cyto.type=="band") {
        cyto <- NULL
    }
    ## PLOT
    if (is.character(zoom) && (length(zoom)>1 || !grepl(":", zoom))) {
        kp <- plotKaryotype(genome=genome, chromosomes=zoom, plot.type=4,
                            ideogram.plotter=NULL, labels.plotter=NULL, cytobands=cyto, plot.params=pp)
    } else {
        kp <- plotKaryotype(genome=genome, zoom=zoom, plot.type=4,
                            ideogram.plotter=NULL, labels.plotter=NULL, cytobands=cyto, plot.params=pp)
    }
    ## names
    kpAddChromosomeNames(kp, str_replace(kp$chromosomes, "chr", ""), cex=pp$chromosome.cex)
    ## cytobands
    kpAddCytobands(kp, color.table=color.table, lwd=pp$cyto.lwd, lty=pp$cyto.lty)
    return(kp)
}

.kpArm <- function(genome) {
    tmp <- karyoploteR::getCytobands(genome)
    tmp <- keepStandardChromosomes(tmp, pruning.mode="coarse")
    tmp <- tmp[seqnames(tmp) != "chrM"]
    tmp <- sort(unlist(GenomicRanges::reduce(split(tmp, str_sub(tmp$name, 1, 1)))))
    tmp$name <- paste0(seqnames(tmp), names(tmp))
    tmp <- sort(unname(tmp))
    return(tmp)
}

.kpAddMainTitleCustom <- function(karyoplot, main = NULL, ...) {
    if (!is.null(main)) {
        karyoplot$beginKpPlot()
        on.exit(karyoplot$endKpPlot())
        bb <- getMainTitleBoundingBox(karyoplot)
        x <- bb$x0
        y <- (bb$y0 + bb$y1)/2
        graphics::text(x = x, y = y, adj=c(0, 0.5), labels = main, ...)
    }
    invisible(karyoplot)
}

.kpTitlePlot <- function(dig, title=NULL, kp, kp.opts) {
    title <- paste(
        title,
        sprintf("purity:%.0f%% ploidy:%.1f sex:%s sample:%s", dig$purity*100, dig$ploidy, dig$sex, dig$sample)
    )
    .kpAddMainTitleCustom(kp, title)
}

.kpLrPlot <- function(r0, r1, tile, purity, ploidy, kp, kp.opts) {
    ## PARAMS
    pp <- getDefaultPlotParams(plot.type = 4)
    pp <- modifyList(pp, kp.opts)
    ## boxes
    if (pp$box.odd || pp$box.even) {
        box.all <- GRanges(seqinfo(tile))
        if (pp$box.odd) {
            box.odd <- box.all[seq(1, length(box.all), 2)]
            kpRect(kp, box.odd, r0=r0, r1=r1, y0=0.0, y1=1.0, col=pp$box.oddcol, lwd=0, lty=0)
        }
        if (pp$box.even) {
            box.even <- box.all[seq(2, length(box.all), 2)]
            kpRect(kp, box.even, r0=r0, r1=r1, y0=0.0, y1=1.0, col=pp$box.evencol, lwd=0, lty=0)
        }
    }
    col <- colTypeSwitch(tile, pp$point.lr.col.type, "black")
    y  <-  tile$lr
    tick <- pp$axis.lr.tick
    if (pp$axis.lr.scale=="abs") {
        tick.pos <- copyToRatio(tick, purity, ploidy)
    } else if (pp$axis.lr.scale=="log") {
        tick.pos <- tick
    }
    if (pp$axis.lr.free) {
        ymin <- min(c(tile$lr, tick.pos-0.33), na.rm=TRUE)
        ymax <- max(c(tile$lr, tick.pos+0.33), na.rm=TRUE)
    } else {
        ymin <- -4
        ymax <- +4
        y[y < -4] <- -4
        y[y > +4] <- +4
    }
    kpPoints(kp, data=tile, y=y, r0=r0, r1=r1, ymin=ymin, ymax=ymax, cex=pp$point.lr.cex, col=col)
    ## axes
    tick.len <- pp$tick.len * sum(width(kp$plot.region))
    kpAxis(kp, r0=r0, r1=r1, ymin=ymin, ymax=ymax, tick.len=tick.len, tick.pos=tick.pos, cex=pp$axis.lr.cex, labels=tick)
}

.kpBafPlot <- function(r0, r1, var, kp, kp.opts) {
    ## PARAMS
    pp <- getDefaultPlotParams(plot.type = 4)
    pp <- modifyList(pp, kp.opts)
    ## box
    if (pp$box.odd || pp$box.even) {
        box.all <- GRanges(seqinfo(var))
        if (pp$box.odd) {
            box.odd <- box.all[seq(1, length(box.all), 2)]
            kpRect(kp, box.odd, r0=r0, r1=r1, y0=0.0, y1=1.0, col=pp$box.oddcol, lwd=0, lty=0)
        }
        if (pp$box.even) {
            box.even <- box.all[seq(2, length(box.all), 2)]
            kpRect(kp, box.even, r0=r0, r1=r1, y0=0.0, y1=1.0, col=pp$box.evencol, lwd=0, lty=0)
        }
    }
    ## points
    if (pp$axis.baf.free) {
        dbaf <- max(abs(var$AF-0.5))
        ymax <- 0.5 + dbaf
        ymin <- 0.5 - dbaf
        tick.pos <- c(ymin+0.1, 0, ymax-0.1)
        tick.lab <- sprintf("%.0f%%", tick.pos*100)
    } else {
        ymin <- 0.25
        ymax <- 0.75
        tick.pos <- c(0.25, 0.5, 0.75)
        tick.lab <- c("25%", "50%", "75%")
    }
    if (length(var)!=0) {
        col <- colTypeSwitch(var, pp$point.baf.col.type, "black")
        kpPoints(kp, data=var, y=var$AF, r0=r0, r1=r1, ymin=0, ymax=1, cex=pp$point.baf.cex, col=col)
    }
    ## axes
    tick.len <- pp$tick.len * sum(width(kp$plot.region))
    kpAxis(kp, r0=r0, r1=r1, ymin=0, ymax=1, tick.len=tick.len, tick.pos=tick.pos, labels=tick.lab, cex=pp$axis.baf.cex)
}

kpLrPlot <- function(dig, kp, kp.opts) {
    seqlevelsStyle(dig$tile) <- "UCSC"
    .kpLrPlot(0, 1, dig$tile[dig$tile$target], dig$purity, dig$ploidy, kp, kp.opts)
}

kpBafPlot <- function(dig, kp, kp.opts) {
    seqlevelsStyle(dig$var) <- "UCSC"
    .kpBafPlot(0, 1, dig$var, kp, kp.opts)
}

kpLrBafPlot <- function(dig, title=NULL, kp, kp.opts) {
    seqlevelsStyle(dig$var) <- "UCSC"
    seqlevelsStyle(dig$tile) <- "UCSC"
    .kpTitlePlot(dig, title, kp, kp.opts)
    .kpBafPlot(0, 0.485, dig$var, kp, kp.opts)
    .kpLrPlot(0.515, 1, dig$tile[dig$tile$target], dig$purity, dig$ploidy, kp, kp.opts)
}

kpSegStats <- function(seg, meta, kp, kp.opts, zoom=NULL,
                       tracks=c("gain", "loss", "loh", "diff"),
                       gain.col=c("rel.gain", "abs.gain", "rel.amp", "abs.amp"),
                       loss.col=c("rel.loss", "abs.loss", "rel.del", "abs.del"),
                       diff.col=c("rel.diff", "rel.dpct"),
                       chr.col="seqnames", sample.col="sample",
                       ymax.loss=1.0, ymax.gain=1.0) {
    ## gain type
    gain.col <- match.arg(gain.col)
    loss.col <- match.arg(loss.col)
    diff.col <- match.arg(diff.col)

    ## prepare segment table and ranges
    seg.dt <- as.data.table(seg)
    seg.dt <- seg.dt[,.(
        seqnames=get(chr.col),
        start=start,
        end=end
    )]
    seg.gr <- toGRanges(seg.dt)
    ## integrate all segments across samples
    seg.dj <- disjoin(seg.gr)
    hit.dj <- as.data.table(findOverlaps(seg.dj, seg.gr))
    hit.dj.anno <- cbind(hit.dj, seg[hit.dj$subjectHits])
    ## add ploidy information
    setkeyv(meta, sample.col)
    setkeyv(hit.dj.anno, sample.col)
    hit.dj.anno[meta, ":="(ploidy=ploidy)]
    ##
    rel.thr <- ifelse(is.null(kp.opts$SegStats.rel.thr), 0.95, kp.opts$SegStats.rel.thr)

    dj.stats <- hit.dj.anno[,.(
        n=.N,
        loh=mean(K==0,na.rm=TRUE),
        abs.gain=mean(C>nC,na.rm=TRUE),
        abs.loss=mean(C<nC,na.rm=TRUE),
        abs.amp=mean(C>=7 & width<10e6,na.rm=TRUE),
        abs.del=mean(C==0 & width<10e6,na.rm=TRUE),
        rel.gain=mean(C-(ploidy * nC/2)>  rel.thr,na.rm=TRUE),
        rel.loss=mean(C-(ploidy * nC/2)< -rel.thr,na.rm=TRUE),
        rel.amp=mean(C-(ploidy * nC/2)>  3 * rel.thr,na.rm=TRUE),
        rel.del=mean(C-(ploidy * nC/2)< -2 * rel.thr,na.rm=TRUE),
        rel.diff=pmax(pmin(mean(C-(ploidy * nC/2), na.rm=TRUE), 4), -2),
        rel.dpct=mean(C-(ploidy * nC/2)>  rel.thr,na.rm=TRUE)-mean(C-(ploidy * nC/2)< -rel.thr,na.rm=TRUE)
    ),queryHits]
    
    dj.stats[n<10,
    ":="(loh=0, abs.gain=0, abs.loss=0, abs.amp=0, abs.del=0, rel.gain=0, rel.loss=0, rel.diff=0, rel.dpct=0)
    ]
    dj.stats[is.nan(loh), loh:=0]
    dj.stats[is.nan(abs.gain), abs.gain:=0]
    dj.stats[is.nan(abs.loss), abs.loss:=0]
    dj.stats[is.nan(abs.amp), abs.amp:=0]
    dj.stats[is.nan(abs.del), abs.del:=0]
    dj.stats[is.nan(rel.gain), rel.gain:=0]
    dj.stats[is.nan(rel.loss), rel.loss:=0]
    dj.stats[is.nan(rel.diff), rel.diff:=0]
    dj.stats[is.nan(rel.dpct), rel.dpct:=0]
    dj.stats <- dj.stats[order(queryHits)]

    mcols(seg.dj) <- cbind(mcols(seg.dj), dj.stats)
    ##
    n.track <- length(tracks)
    for (i in 1:n.track) {
        if (tracks[i]=="loh") {
            at <- autotrack(i, n.track, r0=0.0, r1=1.0, margin=0.1)
            kpArea(kp, seg.dj, y=seg.dj$loh, r0=at$r0, r1=at$r1, ymin=0, ymax=1, col=CNVEX_COPY_COL["K0"])
            kpAxis(kp, ymin = 0, ymax = 1, r0=at$r0, r1=at$r1, lwd=2.5, side=1, cex=1, numticks = 2, labels=c("0%", "100%"))
            kpAddLabels(kp, "loh", r0=at$r0, r1=at$r1, label.margin=0.025, cex=1.5)
        }
        if (tracks[i]=="loss") {
            at <- autotrack(i, n.track, r0=0.0, r1=1.0, margin=0.1)
            if (ymax.loss=="free") {
                ymax <- min(max(mcols(seg.dj)[[loss.col]]) + 0.1, 1.0)
            } else {
                ymax <- ymax.loss
            }
            ylab <- sprintf("%d%%", round(ymax * 100))
            kpArea(kp, seg.dj, y=mcols(seg.dj)[[loss.col]], r0=at$r0, r1=at$r1, ymin=0, ymax=ymax, col=CNVEX_COPY_COL["C1"])
            kpAxis(kp, ymin = 0, ymax = ymax, r0=at$r0, r1=at$r1, lwd=2.5, side=1, cex=1, numticks = 2, labels=c("0%", ylab))
            kpAddLabels(kp, "loss", r0=at$r0, r1=at$r1, label.margin=0.025, cex=1.5)
        }
        if (tracks[i]=="gain") {
            at <- autotrack(i, n.track, r0=0.0, r1=1.0, margin=0.1)
            if (ymax.gain=="free") {
                ymax <- min(max(mcols(seg.dj)[[gain.col]]) + 0.1, 1.0)
            } else {
                ymax <- ymax.gain
            }
            ylab <- sprintf("%d%%", round(ymax * 100))
            kpArea(kp, seg.dj, y=mcols(seg.dj)[[gain.col]], r0=at$r0, r1=at$r1, ymin=0, ymax=ymax, col=CNVEX_COPY_COL["C3"])
            kpAxis(kp, ymin = 0, ymax = ymax, r0=at$r0, r1=at$r1, lwd=2.5, side=1, cex=1, numticks = 2, labels=c("0%", ylab))
            kpAddLabels(kp, "gain", r0=at$r0, r1=at$r1, label.margin=0.025, cex=1.5)
        }
        if (tracks[i]=="diff") {
            data <- mcols(seg.dj)[[diff.col]]
            sel.up <- data > 0
            data.up <- data
            data.up[!sel.up] <- 0
            data.dn <- data
            data.dn[sel.up] <- 0
            if (diff.col=="rel.diff") {
                at <- autotrack(i, n.track, r0=0.0, r1=1.0, margin=0.1)
                kpArea(kp, seg.dj, y=data.up, r0=at$r0, r1=at$r1, ymin=-2, ymax=+4, col=CNVEX_COPY_COL["C3"])
                kpArea(kp, seg.dj, y=data.dn, r0=at$r0, r1=at$r1, ymin=-2, ymax=+4, col=CNVEX_COPY_COL["C1"])
                kpAxis(kp, ymin = -2, ymax = +4, r0=at$r0, r1=at$r1, lwd=2.5, side=1, cex=1, numticks = 3, tick.pos=c(-2,0,4), labels=c("-2", "0", "+4"))
                kpAddLabels(kp, "rel-diff", r0=at$r0, r1=at$r1, label.margin=0.025, cex=1.5)
            } else {
                at <- autotrack(i, n.track, r0=0.0, r1=1.0, margin=0.1)
                kpArea(kp, seg.dj, y=data.up, r0=at$r0, r1=at$r1, ymin=-1, ymax=+1, col=CNVEX_COPY_COL["C3"])
                kpArea(kp, seg.dj, y=data.dn, r0=at$r0, r1=at$r1, ymin=-1, ymax=+1, col=CNVEX_COPY_COL["C1"])
                kpAxis(kp, ymin = -1, ymax = +1, r0=at$r0, r1=at$r1, lwd=2.5, side=1, cex=1, numticks = 3, tick.pos=c(-1,0,1), labels=c("-1", "0", "+1"))
                kpAddLabels(kp, "rel-dpct", r0=at$r0, r1=at$r1, label.margin=0.025, cex=1.5)
            }
        }
    }
}



kpSegHeatmap <- function(seg, kp, kp.opts, zoom=NULL, sample.col="sample", chr.col="seqnames", data.col="C", data.arm=NULL,
                        col.type="C", wt.col="black", row.order = c("by-name", "by-data-sort", "by-data-cluster"), row.names=FALSE
                        ) {
    ## input validation
    row.order <- match.arg(row.order)
    stopifnot(sample.col %in% colnames(seg))
    stopifnot(chr.col %in% colnames(seg))
    stopifnot("C" %in% colnames(seg))
    stopifnot("K" %in% colnames(seg))
    ## create segment data.table and GRanges
    seg.dt <- as.data.table(seg)
    seg.dt <- seg.dt[,.(
        seqnames=get(chr.col),
        start=start,
        end=end,
        data=get(data.col),
        sample=get(sample.col),
        nlr,
        col=colTypeSwitch(seg.dt, col.type, wt.col, TRUE)
    )]
    seg.gr <- toGRanges(seg.dt)
    
    ## add arm information
    arm <- .kpArm(unique(genome(kp$genome)))
    seg.dt$arm <- arm[findOverlaps(seg.gr, arm, select="first")]$name

    ## row order
    if (row.order=="by-name") {
        samp <- unique(seg.dt$sample)
    } else {
        ## D-matrix
        seg.dt[,D:=data]
        D.dt <- seg.dt[,.(D=weighted.mean(as.numeric(data), as.numeric(nlr), na.rm=TRUE)),.(sample, arm)]
        D.dt <- dcast.data.table(D.dt, sample ~ arm, value.var="D")
        D.mx <- as.matrix(D.dt[,-1])
        D.mx[is.na(D.mx)] <- 0
        rownames(D.mx) <- D.dt$sample
        if (!is.null(data.arm)) {
            D.mx <- D.mx[, data.arm, drop=FALSE]
        }
        if (row.order=="by-data-sort") {
            ord <- order(rowMeans(D.mx))
        }
        if (row.order=="by-data-cluster") {
            D.hc <- hclust(dist(D.mx))
            ord <- D.hc$order
        }
        samp <- D.dt$sample[ord]
    }

    ## plot
    n.samp <- length(samp)
    for (i in 1:n.samp) {
        s <- samp[i]
        dt <- seg.dt[(sample==s)]
        at <- autotrack(i, n.samp, margin=0.0, r1=1.0)
        kpPlotRegions(kp, dt[,.(seqnames, start, end)], r0=at$r0, r1=at$r1, col=dt$col)
        if (row.names) {
            kpAddLabels(kp, labels=s, r0=at$r0, r1=at$r1)
        }
        
    }
    kp$..samp <- samp
    return(kp)
}
