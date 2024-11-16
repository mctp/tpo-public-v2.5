#' @export
getGgOpts <- function(settings, opts=list()) {
    gg.opts <- get.opts(settings, opts)
}

.ggPlotDataCov <- function(dig, genes, off, gg.opts) {
    tile <- dig$tile
    all.chr <- seqlevels(dig$tile)
    cov.dt <- data.table(
        chr=as.character(seqnames(tile)),
        start=floor((start(tile)+end(tile))/2),
        end=floor((start(tile)+end(tile))/2),
        type="COV",
        as.data.table(mcols(tile))
    )
    cov.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    cov.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]
    if (!gg.opts$lr.off.target) {
      cov.dt <- cov.dt[(target)]
    }
    ## add gene
    cov.gr <- with(cov.dt, GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
    ol <- findOverlaps(cov.gr, genes)
    cov.dt$gene <- NA
    cov.dt$gene[queryHits(ol)] <- genes[subjectHits(ol)]$gene_name
    return(cov.dt)
}

.ggPlotDataSnp <- function(dig, genes, off, gg.opts) {
    snp <- dig$var
    if (length(snp) == 0) {
        snp.dt <- data.table(
            chr=character(0),
            start=integer(0),
            end=integer(0),
            type=character(0),
            idx=integer(0),
            SOURCE=character(0),
            QUAL=integer(0),
            AF=numeric(0),
            DP=integer(0),
            seg=integer(0),
            sC=numeric(0),
            C=numeric(0), # why?
            K=integer(0),
            seg.lr=numeric(0),
            seg.len=integer(0),
            seg.nlr=integer(0),
            seg.naf=numeric(0), # why?
            start.off=integer(0),
            end.off=integer(0),
            gene=character(0),
            lr=numeric(0)
        )
        return(snp.dt)
    }
    snp.dt <- data.table(
        chr=as.character(seqnames(snp)),
        start=start(snp),
        end=start(snp),
        type="BAF",
        idx=(1:length(snp)),
        as.data.table(mcols(snp))
    )
    all.chr <- seqlevels(dig$tile)
    snp.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    snp.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]
    ## genes
    snp.gr <- with(snp.dt,GRanges(seqnames=chr,ranges=IRanges(start=start,end=end)))
    ol <- findOverlaps(snp.gr,genes)
    snp.dt$gene <- NA
    snp.dt$gene[queryHits(ol)] <- genes[subjectHits(ol)]$gene_name
    ## lr
    snp.dt.gr <- GRanges(snp.dt)
    overlaps <- nearest(x = snp.dt.gr, subject = dig$tile, select = "all")
    snp.dt$lr[overlaps@from] <- dig$tile$lr[overlaps@to]
    return(snp.dt)
}

.ggPlotDataSeg <- function(dig, genes, off, gg.opts) {
    seg <- dig$seg
    fit <- dig$fit
    if (is.null(fit)) {
        fit <- data.table(
            seg=1:length(seg),
            C=NA_integer_,
            K=NA_integer_,
            lr=NA_real_,
            tL=NA_real_,
            aL=NA_real_,
            d=NA_real_,
            anom=NA_real_,
            mse=NA_real_,
            nlr=NA_real_,
            naf=NA_real_,
            len=NA_real_,
            sC=NA_real_
        )
    }
    all.chr <- seqlevels(dig$tile)
    seg.dt <- data.table(
        chr=as.character(seqnames(seg)),
        start=start(seg),
        end=end(seg),
        type="SEG",
        seg = 1:length(seg))
    seg.dt <- merge(seg.dt, fit, by="seg", all.x = TRUE)
    seg.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    seg.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]
    setkey(seg.dt, seg)
    return(seg.dt)
}

.ggGetChrOff <- function(seg) {
    tmp <- cumsum(as.numeric(seqlengths(seqinfo(seg))))
    off <- c(0, head(tmp, -1))
    names(off) <- seqnames(seqinfo(seg))
    chr.off <- data.table(chr=names(off), chr.start=off+1, chr.end=tmp, chr.col=rep(c("odd","even"), length.out=length(tmp)))
    return(list(off=off, chr.off=chr.off))
}

#' @export
ggPlotData <- function(dig, genes, gg.opts) {
    mgrid <- modelGrid(dig$purity, dig$ploidy, gg.opts$abs.range)
    mlines <- modelLines(dig$purity, gg.opts$baf.lines$M, gg.opts$baf.lines$C)
    ## chromosome offsets
    off_chr.off <- .ggGetChrOff(dig$seg)
    chr.off <- off_chr.off$chr.off
    off <- off_chr.off$off
    ## plotting tables
    cov.dt <- .ggPlotDataCov(dig, genes, off, gg.opts)
    snp.dt <- .ggPlotDataSnp(dig, genes, off, gg.opts)
    seg.dt <- .ggPlotDataSeg(dig, genes, off, gg.opts)
    ## final object
    pd <- list(
        cov=cov.dt,
        snp=snp.dt,
        seg=seg.dt,
        off=chr.off, mgrid=mgrid, mlines=mlines)
    return(pd)
}

.ggZoomXlim <- function(off, zoom) {
    if (is.null(zoom)) {
        xlim <- NULL
    } else {
        zoom <- GRanges(zoom)
        shift <- off$chr.start-1
        names(shift) <- off$chr
        zoom.off <- shift(zoom, shift=shift[head(as.character(seqnames(zoom)),1)])
        xlim <- c(start(zoom.off), end(zoom.off))
    }
    return(xlim)
}

.baseLrPlot <- function(off, mgrid, zoom, gg.opts) {
    if (gg.opts$lr.range.dynamic) {
        ymin <- head(mgrid, 1)$lr - 0.25
        ymax <- tail(mgrid, 1)$lr + 0.25
    } else {
        ymin <- gg.opts$lr.range.limit[1]
        ymax <- gg.opts$lr.range.limit[2]
    }
    mgrid.y <- mgrid[C %in% gg.opts$lr.lines]
    lr.breaks <- (gg.opts$lr.range.limit[1]+1):(gg.opts$lr.range.limit[2]-1)
    plt <- ggplot(aes=aes(text="placeholder")) +
        coord_cartesian(ylim=c(ymin,ymax)) +
        ## box
        scale_x_continuous(expand=c(0,0), sec.axis=dup_axis()) +
        geom_rect(aes(xmin=chr.start, xmax=chr.end, ymin=-Inf, ymax=Inf, fill=chr.col), data=off) +
        scale_fill_manual(values=c(gg.opts$box.evencol, gg.opts$box.oddcol), guide=FALSE) +
        ## y-axis
        geom_hline(aes(yintercept=lr), mgrid.y, size=gg.opts$lr.lines.size, col=gg.opts$lr.lines.col, alpha=gg.opts$lr.lines.alpha) +
        scale_y_continuous("Copy Number", expand=c(0,0), breaks=mgrid.y$lr, labels=mgrid.y$C,
                           sec.axis=dup_axis(name="Log2 Cov Ratio", breaks=lr.breaks, labels=lr.breaks)) +
        ## theme
        theme_pubr() +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.line.y.left = element_line(size = 0.5),
            axis.line.y.right = element_line(size = 0.5),
            axis.line.x.top = element_line(size = 0.5),
            axis.line.x.bottom = element_line(size = 0.5),
            legend.position="none"
        )
    if (gg.opts$lr.chr.names) {
        plt <- plt + geom_text(aes(
                         x=(chr.start+chr.end)/2,
                         y=ymax-(ymax-ymin)*ifelse(1:nrow(off) %% 2, gg.opts$lr.chr.yoff[1], gg.opts$lr.chr.yoff[2]),
                         label=str_replace(chr, "chr", gg.opts$lr.chr.prefix)
                     ), data=off, color=gg.opts$lr.chr.col)
    }
    if(is.null(zoom)) {
      plt <- plt + theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
      )
    }
    return(plt)
}

.ggLrPlot <- function(pd, zoom, gg.opts) {
    ## colors
    pd.cov <- pd$cov[!is.na(pd$cov$lr)]
    pk <- colTypeSwitch(pd.cov, gg.opts$lr.col.type, gg.opts$lr.col.wt, decode=FALSE)
    ## make plot
    plt <- .baseLrPlot(pd$off, pd$mgrid, zoom, gg.opts)
    ## make sure not points are lost
    pd.cov[pd.cov$lr < plt$coordinates$limits$y[1]]$lr <- plt$coordinates$limits$y[1]
    pd.cov[pd.cov$lr > plt$coordinates$limits$y[2]]$lr <- plt$coordinates$limits$y[2]
    #subset the data and plot by coordinate if zooming
    if(!is.null(zoom)) {
      zoom <- GRanges(zoom)
      chr <- as.character(seqnames(zoom)[1])
      w <- which(pd.cov$start>=start(zoom) & pd.cov$end<=end(zoom) & pd.cov$chr==chr)
      pd.cov <- pd.cov[w,]
      pk$keys <- pk$keys[w]
      plt <- plt +
          geom_point(
            aes(x=start, y=lr, color=pk$keys),
            pd.cov, size=0.75, alpha=1
          ) +
          scale_color_manual(values=pk$palette, guide=FALSE) +
          scale_x_continuous(limits=c(start(zoom),end(zoom)))
    } else {
      plt <- plt +
          geom_point(aes(x=start.off, y=lr, color=pk$keys), pd.cov, size=0.75, alpha=1) +
          scale_color_manual(values=pk$palette, guide=FALSE)
    }
    return(plt)
}

.baseBafPlot <- function(off, mgrid, mlines, zoom, gg.opts) {
    if (gg.opts$baf.range.dynamic) {
        ymin <- max(0, 0+min(mgrid$af))
        ymax <- 1-ymin
    } else {
        ymin <- gg.opts$baf.range.limit[1]
        ymax <- gg.opts$baf.range.limit[2]
    }
    plt <- ggplot() +
        coord_cartesian(ylim=c(ymin,ymax)) +
        ## box
        geom_rect(aes(xmin=chr.start, xmax=chr.end, ymin=-Inf, ymax=Inf, fill=chr.col), data=off) +
        scale_fill_manual(values=c(gg.opts$box.evencol, gg.opts$box.oddcol), guide=FALSE) +
        ## y-axis
        geom_hline(aes(yintercept=baf), mlines, size=gg.opts$baf.lines.size, color=gg.opts$baf.lines.col,alpha=gg.opts$baf.lines.alpha) +
        scale_y_continuous("Copies", expand=c(0,0), breaks=mlines$baf, labels=mlines[[gg.opts$baf.lines.lab]],
                           sec.axis=dup_axis(name="BAF", breaks=c(0.3,0.5,0.7), labels=c("0.3", "0.5", "0.8"))
                           ) +
        ## x-axis
        scale_x_continuous(expand=c(0,0), sec.axis=dup_axis()) +
        theme_pubr() +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.line.y.left = element_line(size = 0.5),
            axis.line.y.right = element_line(size = 0.5),
            axis.line.x.top = element_line(size = 0.5),
            axis.line.x.bottom = element_line(size = 0.5),
            legend.position="none"
        )
    if (gg.opts$baf.chr.names) {
        plt <- plt + geom_text(aes(
                         x=(chr.start+chr.end)/2,
                         y=ymax-ifelse(1:nrow(off) %% 2, gg.opts$baf.chr.yoff[1], gg.opts$baf.chr.yoff[2]),
                         label=str_replace(chr, "chr", gg.opts$baf.chr.prefix),
                         ), data=off, color=gg.opts$baf.chr.col)
    }
    if(is.null(zoom)) {
      plt <- plt + theme(
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
      )
    }
    return(plt)
}

.ggBafPlot <- function(pd, zoom, gg.opts) {
    ## subsample points
    snp <- pd$snp
    snp <- snp[sample(nrow(snp))]
    snp <- snp[,head(.SD,gg.opts$point.per.segment),by=seg]
    ## colors
    pk <- colTypeSwitch(snp, gg.opts$baf.col.type, gg.opts$baf.col.wt, decode=FALSE)
    ## make plot
    plt <- .baseBafPlot(pd$off, pd$mgrid, pd$mlines, zoom, gg.opts)
    #subset the data and plot by coordinate if zooming
    if(!is.null(zoom)) {
      zoom <- GRanges(zoom)
      chr <- as.character(seqnames(zoom)[1])
      w <- which(snp$start.off>=start(zoom) & snp$end<=end(zoom) & snp$chr==chr)
      snp <- snp[w,]
      pk$keys <- pk$keys[w]
      plt <- plt +
          geom_point(
            aes(x=start, y=AF, color=pk$keys,gene=gene,seg=seg,dp=DP,af=AF,C=C,K=K),
            snp, size=0.75, alpha=1
          ) +
          scale_color_manual(values=pk$palette, guide=FALSE) +
          scale_x_continuous(limits=c(start(zoom),end(zoom)))
          #coord_cartesian(xlim=c(start(zoom),end(zoom)))
    } else {
    plt <- plt +
        geom_point(aes(x=start.off, y=AF, color=pk$keys), snp, size=0.75, alpha=1) +
        scale_color_manual(values=pk$palette, guide=FALSE)
    }
    return(plt)
}

.ggSegPlot <- function(pd, gg.opts) {

  # params
  seg <- pd$seg

  # base plot (no data)
  plt <- .baseLrPlot(pd$off, pd$mgrid, zoom=NULL, gg.opts)

  # segment plot
  plt <- plt +
    geom_segment(aes(x=start.off, xend=end.off, y=lr, yend=lr, col=factor(seg %% 3)), seg, size=2) +
    scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none")

  return(plt)
}

#' @export
ggLrPlot <- function(dig, genes, zoom=NULL, gg.opts) {
    pd <- ggPlotData(dig, genes, gg.opts)
    lr.plt <- .ggLrPlot(pd, zoom, gg.opts)
    return(lr.plt)
}

#' @export
ggBafPlot <- function(dig, genes, zoom=NULL, gg.opts) {
    pd <- ggPlotData(dig, genes, gg.opts)
    baf.plt <- .ggBafPlot(pd, zoom, gg.opts)
    return(baf.plt)
}

#' @export
ggLrBafPlot <- function(dig, genes, zoom=NULL, gg.opts) {
    pd <- ggPlotData(dig, genes, gg.opts)
    lr.plt <- .ggLrPlot(pd, zoom, gg.opts)
    baf.plt <- .ggBafPlot(pd, zoom, gg.opts)
    return(ggarrange(lr.plt, baf.plt, ncol=1, newpage = FALSE, draw = FALSE))
}

#' @export
ggCopyBafGrid <- function(dig, genes, gg.opts) {
  pd <- ggPlotData(dig, genes, gg.opts)
  grid <- unique(pd$mgrid[,c("C","lr")])
  plt <- ggplot(pd$snp, aes(x=AF,y=lr)) +
    stat_binhex(aes(fill = log(..count..)),bins=100, alpha = 0.75) +
    geom_point(aes(x=AF,y=seg.lr), alpha=0.15,col="blue",size=0.5) +
    geom_point(aes(x=af,y=lr),color="red",size=2,data=pd$mgrid) +
    scale_y_continuous(labels=grid$C,breaks=grid$lr, limits = c(min(grid$lr)-.25, max(grid$lr)+0.25)) +
    scale_fill_gradient(low = "white", high = "black", guide = "none") +
    xlab("BAF") + ylab("Copy Number") +
    theme_pubr(base_size=14) +
    theme(legend.position = "none")
  return(plt)
}

#' @export
ggSegPlot <- function(dig, genes, gg.opts) {
  pd <- ggPlotData(dig, genes, gg.opts)
  seg.plt <- .ggSegPlot(pd, gg.opts)
}

#' Plots summary plot for model picking
#'
#' @param dig a digest exported from cnvex
#' @param genes a granges object
#' @param gg.opts obtained from getGgOpts()
#' 
#' @export
ggArrangedPlot <- function(dig, genes, gg.opts) {
  # plots
  seg.plt <- ggSegPlot(dig, genes, gg.opts)
  bafgrid.plt <- ggCopyBafGrid(dig, genes, gg.opts)
  lrbaf.plt <- ggLrBafPlot(dig, genes, NULL, gg.opts)
  # text
  text <- sprintf("Purity: %s \nPloidy: %s", round(dig$purity,3), round(dig$ploidy,3))
  text.plt <- ggparagraph(text = text, face = "italic", size = 24, color = "black")
  # arranged
  plt <- cowplot::plot_grid(seg.plt, cowplot::plot_grid(bafgrid.plt, text.plt), lrbaf.plt, nrow =3)
  return(plt)
}

#' @export
ggGridPlot <- function(digest, gg.opts) {
    ## prepare grid
    grid <- modelGrid(digest$purity, digest$ploidy, gg.opts$abs.range)
    gridy  <- unique(grid[,.(C,lr)])
    tmp <- copyCol(grid$C)
    grid$col <- factor(names(tmp), levels=names(CNVEX_COPY_COL), ordered=TRUE)
    ## prepare variant data table
    vdt <- as.data.table(digest$var)
    tile.tgt <- digest$tile[digest$tile$target==TRUE]
    hits <- findOverlaps(digest$var, tile.tgt + digest$opts$tile.shoulder, select="first")
    vdt$lr <- tile.tgt[hits]$lr
    vdt$seg.lr <- tile.tgt[hits]$seg.lr
    ## limit maximum number of points-per-segment
    if (point.per.segment < Inf) {
        vdt <- vdt[,.SD[1:min(gg.opts$point.per.segment, nrow(.SD))], seg]
    }
    ## limit maximum number of segment colors:
    vdt[,seg:=factor(seg %% gg.opts$grid.max.seg.cols)]
    ## plot
    plt <- ggplot(vdt) +
        aes(x=AF) +
        geom_point(aes(y=lr), size=2, color="gray", fill="gray", alpha=0.5, shape=21) +
        geom_point(aes(y=seg.lr, fill=seg), size=1, color="black", alpha=0.50, shape=21) +
        geom_point(aes(color=col, y=lr), data=grid, size=4) +
        scale_color_manual(values=CNVEX_COPY_COL, guide=FALSE) +
        scale_fill_brewer(palette="Paired", guide=FALSE) +
        scale_y_continuous(breaks = gridy$lr, labels = gridy$C) +
        ylab("Copy Number") +
        xlab("B-allele frequency") +
        theme_pubr(base_size=14)
    return(plt)
}
