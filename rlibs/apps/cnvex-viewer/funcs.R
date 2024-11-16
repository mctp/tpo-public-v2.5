STRING_COL <- c(
    "green", "#4500ACFF", "#6B58EEFF", "black", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15",
    "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D",
    "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D", "#67000D"
)
names(STRING_COL) <- as.character(-1:(length(STRING_COL)-2))

seg_goi_format <- function(seg,goi=NULL,v4=NULL){
  if(length(goi)>0){
    seg$notable_genes <- sapply(seg$gene_names, 
                            function(x){
                              foo <- unlist(x)
                              g <- foo[foo %in% goi]
                              if(length(g)==0){'None'}else{paste(g,collapse=',')}
                            }
                          )

  }
  seg$notable_genes <- as.character(seg$notable_genes)
  seg$gene_names <- sapply(seg$gene_names,
                      function(x){
                        foo <- unlist(x)
                        if(!is.null(v4)){
                          foo <- foo[foo %in% v4]
                        }
                        if(length(foo)>5){
                          paste0(c(foo[1:5],sprintf('(%d others)',length(foo)-5)),collapse=',')
                        }else{
                          paste0(foo,collapse=',')
                        }
                      }
                    )
  return(seg)
}

.pd.filter <- function(pd, cnv.filters, min.cov.cov, min.cov.snp) {
    fcb <- "blacklisted" %in% cnv.filters & pd$cov$blacklist
    fcm <- "cov-masked" %in% cnv.filters & pd$cov$unmasked < 0.75
    fcg <- "gapped" %in% cnv.filters & pd$cov$gap > 0
    fct <- "off-target" %in% cnv.filters & !pd$cov$target
    if (all(is.na(pd$cov$n.cov.raw))) {
        n.cov.raw <- Inf
    } else {
        n.cov.raw <- pd$cov$n.cov.raw
    }    
    fcc <- (pd$cov$t.cov.raw < min.cov.cov | n.cov.raw < min.cov.cov) & pd$cov$target
    fsm <- "cov-masked" %in% cnv.filters & pd$snp$mask.strict
    if (is.null(pd$snp$n.DP)) {
        n.DP <- Inf
    } else {
        n.DP <- pd$snp$n.DP
    }
    fsc <- pd$snp$t.DP < min.cov.snp | n.DP < min.cov.snp    
    filter <- list(
        cov=!(fcb | fcm | fcg | fct | fcc),
        snp=!(fsm | fsc)
    )
    filter <- list(cov=TRUE, snp=TRUE)
    return(filter)
}

.baseLogRatioPlot <- function(off, Clr, ymin, ymax) {
    plt <- ggplot(aes=aes(text=rep("zztop", nrow(off)))) + 
        geom_rect(aes(xmin=chr.start, xmax=chr.end, ymin=-Inf, ymax=Inf, fill=chr.col), data=off) +
        scale_fill_manual(values=c("white", "#EEEEEE"), guide="none") +
        geom_text(aes(x=(chr.start+chr.end)/2, y=ymax-(ymax-ymin)*ifelse(1:nrow(off) %% 2, 0.04, 0.07), label=chr), data = off) +
        coord_cartesian(ylim=c(ymin,ymax)) +
        geom_hline(aes(yintercept=lr), Clr, linewidth=0.25) +
        scale_y_continuous("log2(tumor/normal)", expand = c(0, 0), breaks=Clr$lr, labels=Clr$C) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_pubr() +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank()
        )
    return(plt)
}

.baseBafSimplePlot <- function(off, MCbaf) {
    plt <- ggplot() + 
        geom_rect(aes(xmin=chr.start, xmax=chr.end, ymin=-Inf, ymax=Inf, fill=chr.col), off, alpha = 1) +
        scale_fill_manual(values=c("white", "#EEEEEE"), guide="none") +
        geom_hline(aes(yintercept=  baf), MCbaf, linewidth=0.25, color="red") +
        scale_y_continuous("BAF", expand = c(0, 0), breaks=MCbaf$baf, labels=MCbaf$lab) +
        scale_x_continuous(expand = c(0, 0)) +
        theme_pubr() +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.line.x = element_blank()
        )
        return(plt)
}

.baseLogRatioVars <- function(purity, ploidy, range, Msel=0, Csel=1) {
    p <- purity
    D <- (ploidy * p) + 2 * (1-p)
    C <- seq(range[1], range[2], by=1)
    Clr <- data.table(C=factor(C), lr=log2((p * C + (1-p) * 2)  / D))
    ymin <- log2((p * min(C) + (1-p) * 2) / D) - 0.25
    ymax <- log2((p * max(C) + (1-p) * 2) / D) + 0.25
    baf <- (p * Msel + 1 * (1-p)) / (p * Csel + 2 * (1-p))
    MCbaf <- data.table(
        baf=c(baf,1-baf),
        lab=rep(paste(Msel,Csel,sep="/"), 2)
    )
    vars <- list(p=p, D=D, C=C, Clr=Clr, ymin=ymin, ymax=ymax, MCbaf=MCbaf)
    return(vars)
}

.feat2gene <- function(feat, gene) {
    tmp <- as.data.table(findOverlaps(feat, gene))
    setkey(tmp, queryHits)
    tmp <- cbind(tmp, as.data.table(mcols(gene[tmp$subjectHits])[,c("gene_id", "gene_name")]))
    tmp <- tmp[,.(gene_ids=list(gene_id), gene_names=list(gene_name)), queryHits]
    tmp <- tmp[J(seq_along(feat)), .(gene_ids, gene_names)]
    return(tmp)
}

plotData <- function(tile, var, seg, fit, gene, gobj, opts) {
    snp <- var
    all.chr <- seqlevels(gobj$seqi)
    tmp <- cumsum(as.numeric(seqlengths(seqinfo(seg))))
    off <- c(0, head(tmp, -1))
    names(off) <- seqnames(seqinfo(seg))
    chr.off <- data.table(chr=names(off), chr.start=off+1, chr.end=tmp, chr.col=rep(c("A","B"), length.out=length(tmp)))
    cov.dt <- data.table(
        chr=as.character(seqnames(tile)),
        start=floor((start(tile)+end(tile))/2),
        end=floor((start(tile)+end(tile))/2),
        type="COV",
        seg=GenomicRanges::findOverlaps(tile, seg, select = "first"),
        gene_names=.feat2gene(tile, gene)$gene_names,
        as.data.table(mcols(tile))
    )
    cov.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    cov.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]
    
    if (length(snp)>0) {
        snp.dt <- data.table(
            chr=as.character(seqnames(snp)),
            start=start(snp),
            end=start(snp),
            type="BAF",
            idx=(1:length(var)),
            tile=findOverlaps(snp, tile, maxgap=opts$tile.shoulder-1, select="first"),
            seg=findOverlaps(snp, seg, maxgap=opts$tile.shoulder-1, select = "first"),
            gene_names=.feat2gene(snp, gene)$gene_names,
            as.data.table(mcols(snp)[,!(colnames(mcols(snp)) %in% c("REF", "ALT"))])
        )
    } else {
        snp.dt <- data.table(
            chr=character(0),
            start=integer(0),
            end=integer(0),
            type=character(0),
            idx=integer(0),
            tile=integer(0),
            seg=integer(0),
            gene_names=character(0),
            t.AF=numeric(0)
        )
    }
    snp.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    snp.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]

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
    seg.dt <- data.table(
        chr=as.character(seqnames(seg)),
        start=start(seg),
        end=end(seg),
        type="SEG",
        gene_names=.feat2gene(seg, gene)$gene_names,
        fit
    )
    seg.dt[,chr:=factor(chr, all.chr, ordered=TRUE)]
    seg.dt[,":="(
        start.off=start+off[chr],
        end.off=end+off[chr]
    )]
    setkey(seg.dt, seg)
    cov.dt$C <- seg.dt[J(cov.dt$seg)]$C
    snp.dt$C <- seg.dt[J(snp.dt$seg)]$C
    cov.dt$K <- seg.dt[J(cov.dt$seg)]$K
    snp.dt$K <- seg.dt[J(snp.dt$seg)]$K
    pd <- list(cov=cov.dt, snp=snp.dt, seg=seg.dt, off=chr.off)
    return(pd)
}

plotLogRatio <- function(pd, purity=NULL, ploidy=NULL, sel.chr=NULL, sel.data="tile", sel.col="segment", lr.range=c(-4,4), C.range=c(0, 8), max.point=1000) {
    ## select data from chromosomes
    if (!is.null(sel.chr)) {
        cov <- pd$cov[chr %in% sel.chr]
        off <- pd$off[chr %in% sel.chr]
        seg <- pd$seg[chr %in% sel.chr]
    } else {
        cov <- pd$cov
        off <- pd$off
        seg <- pd$seg
    }
    ## subsample points
    cov <- cov[sample(nrow(cov))]
    cov <- cov[,head(.SD,max.point),by=seg]
    ## compute plot limits

    if (is.null(purity) || is.null(ploidy)) {
        Clr <- data.table(C=seq(lr.range[1], lr.range[2]), lr=seq(lr.range[1], lr.range[2]))
        ymin <- lr.range[1]
        ymax <- lr.range[2]
    } else {
        vars <- .baseLogRatioVars(purity, ploidy, C.range)
        Clr <- vars$Clr
        ymin <- vars$ymin
        ymax <- vars$ymax
    }
    ## base plot (no data)
    plt <- .baseLogRatioPlot(off, Clr, ymin, ymax)
    ## data is tile
    if (sel.data=="tile") {
        if (sel.col=="segment") {
            plt <- plt +
                geom_point(aes(x=start.off, y=lr, col=factor(seg %% 3)), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
                scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none") +
                geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
        } else if (sel.col=="sC") {
            plt <- plt +
                geom_point(aes(x=start.off, y=lr, col=sC), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
                scale_color_gradient2(low="blue", mid="black", high="red", midpoint=2, limits=c(0,6), oob = scales::squish) +
                geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
        } else if (sel.col=="C") {
            plt <- plt +
                geom_point(aes(x=start.off, y=lr, color=as.character(C)), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
                scale_color_manual(values=STRING_COL, guide="none") +
                geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
        } else if (sel.col=="CK") {
            cov[C>1 & K==0, C:=-1]
            cov[,C:=as.character(C)]
            plt <- plt +
                geom_point(aes(x=start.off, y=lr, color=C, text=gene_names), cov, size=0.75, alpha=1, na.rm=TRUE) +
                scale_color_manual(values=STRING_COL, guide="none")
        } else {
            plt <- plt + 
                geom_point(aes_string(x="start.off", y="lr", col=sel.col), cov, size=0.75, na.rm=TRUE) +
                scale_color_gradient(low="black", high="red", guide="none")
        }
    }
    ## data is segment
    if (sel.data=="segment") {
        if (sel.col=="segment") {
            plt <- plt + 
                geom_segment(aes(x=start.off, xend=end.off, y=lr, yend=lr, col=factor(seg %% 3)), seg, size=2) + 
                scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none")
        } else {
            plt <- plt +
                geom_segment(aes_string(x="start.off", xend="end.off", y="lr", yend="lr", col=lr-Clr[C==2]$lr), seg, size=2) +
                scale_color_gradient(low="black", high="red", guide="none")
        }
    }
    return(plt)    
}

plotBaf <- function(pd, purity=NULL, ploidy=NULL, sel.col="segment", sel.chr=NULL, baf.range=c(0,1), C.range=c(0,8), Msel=0, Csel=1, max.point=250) {
    if (!is.null(sel.chr)) {
        snp <- pd$snp[chr %in% sel.chr]
        off <- pd$off[chr %in% sel.chr]
    } else {
        snp <- pd$snp
        off <- pd$off
    }

    #backwards compatible with old pre-germline data
    if('t.AF' %in% colnames(snp)){
      snp %>% dplyr::rename(AF=t.AF) -> snp
    }
      
    ## subsample points
    snp <- snp[sample(nrow(snp))]
    snp <- snp[,head(.SD,max.point),by=seg]
    ## compute plot grid lines
    if (is.null(purity) || is.null(ploidy)) {
        MCbaf <- data.table(baf=seq(baf.range[1], baf.range[2], 0.1), lab=seq(baf.range[1], baf.range[2], 0.1))
    } else {
        vars <- .baseLogRatioVars(purity, ploidy, C.range, Msel, Csel)
        MCbaf <- vars$MCbaf
    }    
    plt <- .baseBafSimplePlot(off, MCbaf)
    if (sel.col=="segment") {
        plt <- plt + 
            geom_point(aes(x=start.off, y=AF, col=factor(seg %% 3)), snp, size=0.75, alpha=1, na.rm=TRUE) +
            scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none")
    }
    if (sel.col=="C") {
        plt <- plt +
            geom_point(aes(x=start.off, y=AF, col=as.character(C)), snp, size=0.75, alpha=1, na.rm=TRUE) +
            scale_color_manual(values=STRING_COL, guide="none")
    }
    if (sel.col=="CK") {
        snp[C>1 & K==0, C:=-1]
        snp[,C:=as.character(C)]
        plt <- plt +
            geom_point(aes(x=start.off, y=AF, col=C, text=gene_names), snp, size=0.75, alpha=1, na.rm=TRUE) +
            scale_color_manual(values=STRING_COL, guide="none")
    }
    return(plt)
}

plotGC <- function(pd, n=Inf) {
    tmp <- pd$cov[1:min(nrow(pd$cov), n)]
    if (!all(is.na(tmp$n.cov))) {
        tmp[,lr.raw:=cnvex:::.rawLogRatio(tmp$t.cov, tmp$n.cov)]
        tmp[,lr:=tmp$tn.lr]
    } else {
        tmp[,lr.raw:=cnvex:::.rawLogRatio(tmp$t.cov, 1)]
        tmp[,lr:=tmp$t.lr]
    }
    tmp <- melt(tmp[,.(gc, blacklist, target, lr.raw, lr, lr.off=lr.raw-lr)], id.vars=c("gc", "blacklist", "target"))
    tmp[,delta:=FALSE]
    tmp[variable=="lr.off", delta:=TRUE]
    tmp[variable=="lr.off", variable:="lr.raw"]
    tmp[,target.lab:=ifelse(target, "on-target", "off-target")]
    tmp[,lr.lab:=ifelse(variable=="lr.raw", "raw", "adjusted")]
    gc.plt <- ggplot(tmp) + aes(x=gc, y=value, color=delta) +
        facet_grid(lr.lab~target.lab) +
        geom_point(alpha=0.05, size=0.1, na.rm=TRUE) +
        scale_color_manual(values=c("black", "red"), guide="none") +
        coord_cartesian(xlim=c(0.25, 0.75), ylim=c(-3, 3)) +
        geom_smooth(color = "blue", size = 1, method="gam", formula = y ~ s(x, bs = "cs"), na.rm=TRUE) +
        ylab("log2(tumor/normal)") +
        xlab("GC [%]") +
        theme_pubr(base_size=14) +
        scale_x_continuous(labels=scales::percent_format(accuracy = 1))
    return(gc.plt)
}

plotGrid <- function(grid, cand, opts, var="L", better="higher") {
    plt <- ggplot(grid) +
        aes(x=p, y=P) +
        geom_tile(aes(fill=L)) +
        scale_x_continuous(breaks=seq(0, 1, 0.1)) +
        scale_y_continuous(breaks=seq(opts$opt.P.lo, opts$opt.P.hi, 1)) +
        scale_fill_gradient2(low="blue", mid="white", high="red", name="llik", midpoint=quantile(grid$L, 0.75)) +
        scale_color_manual(values=c("black", "red")) +
        geom_point(aes(color=converged), data=cand, size=3, na.rm=TRUE) +
        xlab("purity") + 
        ylab("ploidy") +
        theme_pubr(legend="right")
    return(plt)
}

plotNvar <- function(pool, opts) {
    ## sex sepecific nvar
    nvar.f <- pool$female$nvar
    nvar.m <- pool$male$nvar
    
    plt.f.ont <- ggplot() +
        geom_histogram(aes(x = nvar.f[pool$target]), binwidth = 0.00025) +
        coord_cartesian(xlim = c(0,0.1)) +
        scale_x_continuous(breaks=seq(0, 1, 0.01)) +
        theme_linedraw() +
        theme(text = element_text(size=7.5)) +
        xlab("Variance") +
        ggtitle("Female nvar distribution", subtitle = "On Target")
    
    plt.f.oft <- ggplot() +
        geom_histogram(aes(x = nvar.f[!pool$target]), binwidth = 0.00025) +
        coord_cartesian(xlim = c(0,0.1)) +
        scale_x_continuous(breaks=seq(0, 1, 0.01)) +
        theme_linedraw() +
        theme(text = element_text(size=7.5)) +
        xlab("Variance") +
        ggtitle("Female nvar distribution", subtitle = "Off Target")
    
    plt.m.ont <- ggplot() +
        geom_histogram(aes(x = nvar.m[pool$target]), binwidth = 0.00025) +
        coord_cartesian(xlim = c(0,0.1)) +
        scale_x_continuous(breaks=seq(0, 1, 0.01)) +
        theme_linedraw() +
        theme(text = element_text(size=7.5)) +
        xlab("Variance") +
        ggtitle("Male nvar distribution", subtitle = "On Target")
    
    plt.m.oft <- ggplot() +
        geom_histogram(aes(x = nvar.m[!pool$target]), binwidth = 0.00025) +
        coord_cartesian(xlim = c(0,0.1)) +
        scale_x_continuous(breaks=seq(0, 1, 0.01)) +
        theme_linedraw() +
        theme(text = element_text(size=7.5)) +
        xlab("Variance") +
        ggtitle("Male nvar distribution", subtitle = "Off Target")
    
    plt <- arrangeGrob(plt.f.ont,plt.f.oft, plt.m.ont, plt.m.oft, nrow = 2)
    return(plt)
}


## temporary function
plotGCncovLr<- function(cnv, pooldata, opts) {
    tmp = data.table()
    if (cnv$sex == "female") {
        cov.med = rowMedians(pooldata$female$cov,na.rm = T)
    } else {
        cov.med = rowMedians(pooldata$male$cov,na.rm = T)
    }
    tmp[,cov.noGCnorm := log2(cnv$tile$n.cov/cov.med)]
    tmp[,cov.GCnorm :=.correctBias(log2(cnv$tile$n.cov/cov.med),pooldata$target,pooldata$hq,"gc",pooldata$gc,opts)]
    tmp[,cov.off := cov.noGCnorm-cov.GCnorm]
    tmp[,gc := cnv$tile$gc]
    tmp = melt(tmp,id.vars = "gc")
    
    gc.plt <- ggplot(tmp) + aes(x=gc, y=value) +
        facet_grid(.~variable) +
        geom_point(alpha=0.05, size=0.1, na.rm=TRUE) +
        coord_cartesian(xlim=c(0.25, 0.75), ylim=c(-3, 3)) +
        ylab("log2(tumor/normal)") +
        xlab("GC [%]") +
        theme_pubr() +
        scale_x_continuous(labels=scales::percent_format(accuracy = 1))
    # geom_smooth(color = "blue",size = 1)
    
    return(gc.plt)
}

lr_compare = function(sample_number,cnv.fns, pool, chrs, opts) {
    cnv = read_rds(cnv.fns[sample_number])
    cnv$tile$filter = pool$female$filter
    cnv$tile$reptime = pool$reptime
    cnv <- addPoolLogRatio(cnv, pool, opts)
    cnv <- addScaleLogRatio(cnv, opts)
    cnv <- addBiasLogRatio(cnv, opts)
    cnv <- addSmoothLogRatio(cnv, opts)
    cnv <- addJointSegment(cnv, opts)
    
    cnv2 <- addPairLogRatio(cnv, opts)
    cnv2 <- addScaleLogRatio(cnv2, opts)
    cnv2 <- addBiasLogRatio(cnv2, opts)
    cnv2 <- addSmoothLogRatio(cnv2, opts)
    cnv2 <- addJointSegment(cnv2, opts)
    
    p1 = lr_plot(cnv,chr=chrs)+ggtitle(paste("pool_sample",sample_number))
    p2 = lr_plot(cnv2,chr=chrs,start = )+ggtitle(paste("pair_sample",sample_number))
    
    grid.arrange(p1,p2)
}

addSeg = function(cnv){
    seg_tiles <- queryHits(findOverlaps(cnv$tile,cnv$seg))
    seg_ids <- subjectHits(findOverlaps(cnv$tile,cnv$seg))
    cnv$tile$tad = NA_real_
    cnv$tile$seg <- findOverlaps(cnv$tile,cnv$seg,select = "first")
    # cnv$tile$seg <- rep(1,length(cnv$tile))
    cnv$tile$seg = seg_ids
    return(cnv)
}

#plot lr vs coordinates:
lr_plot = function(file,chr = "all", value = "lr"){
    file = addSeg(file)
    if(chr == "all") {
        chr = paste0("chr",c(1:22,"X","Y"))
    }
    if(is.null(file$tile$filter)) {
        file$tile$filter = file$tile$hq
    }
    data = as.data.frame(file$tile) %>% dplyr::select(seqnames,start,end,t.cov,n.cov,lr,seg,arm,blacklist, filter)
    data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% dplyr::select(-strand)
    data = data %>% left_join(data_seg,by = "seg")
    data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
    data = data %>% group_by(seg) %>% mutate(seg_lr_med = median(lr,na.rm = T)) %>% ungroup()
    data$blacklist = ifelse(data$blacklist > 0, 1,0)
    if(value == "lr") {value = data$lr} else {
        value = data$n.cov
        data$seg_lr = 0}
    data = data %>% mutate(val = value)
    p = data %>% dplyr::filter(seqnames.x %in% chr) %>% ggplot()+
        facet_grid(.~seqnames.x, scales="free_x")+
        geom_point(size = 0.3,aes(x = start.x,y = val, color = as.factor(filter)), na.rm=TRUE)+
        geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
        theme_linedraw()+
        coord_cartesian(ylim = c(-3,3))+
        geom_vline(xintercept = 1)
    
    return(p)
}

## testing segmentation changes
testSegChanges <- function(seg.old,seg.new,tile, cnv, chr = "all") {
    tile$seg.old = queryHits(findOverlaps(seg.old, tile, select="all"))
    tile$seg.new = queryHits(findOverlaps(seg.new, tile, select="all"))
    bp.old <- which(diff(tile$seg.old) == 1)
    bp.new <- which(diff(tile$seg.new) == 1)
    bp.diff <- bp.old[!(bp.old %in% bp.new)]
    tile$removed <- FALSE
    tile$removed[bp.diff] <- TRUE
    segments.diff <- tile$seg.old[c(bp.diff,(bp.diff+1))]
    
    #plot lr vs coordinates:
    cnv$tile$filter = cnv$tile$hq
    file = cnv
    file = addSeg(file)
    if(chr == "all") {
        chr = paste0("chr",c(1:22,"X","Y"))
    }
    data = as.data.frame(file$tile)
    data_seg = as.data.frame(file$seg) %>% mutate(seg = row_number()) %>% dplyr::select(-strand)
    data = data %>% left_join(data_seg,by = "seg")
    data = data %>% group_by(seg) %>% mutate(seg_lr = mean(lr,na.rm = T)) %>% ungroup()
    data = data %>% group_by(seg) %>% mutate(seg_lr_med = median(lr,na.rm = T)) %>% ungroup()
    data$blacklist = ifelse(data$blacklist > 0, 1,0)
    if(value == "lr") {val = data$lr} else {
        val = data$n.cov
        data$seg_lr = 0}
    data = data %>% mutate(val = val)
    p = data %>% dplyr::filter(seqnames.x %in% chr) %>% ggplot()+
        facet_grid(.~seqnames.x, scales="free_x")+
        geom_point(size = 0.3,aes(x = start.x,y = val, color = as.factor(filter)), na.rm=TRUE)+
        geom_segment(aes(x = start.y,xend = end.y,y = seg_lr,yend = seg_lr),colour = 'red')+
        theme_linedraw()+
        coord_cartesian(ylim = c(-3,3))+
        geom_vline(xintercept = 1)
    vline.data = data.frame(z = data$start.x[bp.diff], seqnames.x = data$seqnames.x[bp.diff])
    vline.data = vline.data %>% filter(seqnames.x %in% chr)
    p2 = p+geom_vline(aes(xintercept = z),vline.data, color = "green", alpha = 0.4, linewidth = 0.5) +
        labs(caption = paste(as.character(nrow(vline.data)),"BP has been merged"))
    return(p2)
}


plotLogRatioQuick <- function(tile, seg, sample.type, purity=NULL, ploidy=NULL, sel.chr=NULL, lr.range=c(-4,4), C.range=c(0, 8), opts) {
  # seqlevelsStyle(seg)  <- "NCBI"
  # seqlevelsStyle(tile) <- "NCBI"
  seg.plot <- as.data.table(seg)
  sel <- cnvex:::.sample.switch(sample.type, opts) # TODO: remove te cnvex:::
  ## average lr per segment
  hits <- as.data.table(findOverlaps(seg, tile))
  hits$lr.tile <- mcols(tile)[, sel$lr]
  seg.lr <- hits[, .(lr = mean(lr.tile, na.rm=TRUE)), by = queryHits]
  ## fix
  seg.plot$lr            <- seg.lr$lr
  seg.plot[, backgroundCol := ifelse(rleid(seqnames)%%2 == 1, "white", "#eaeaea")] 
  seg.plot[, segCol := (rleid(.I)%%3)] 
  
  ## pick Clr
  if (is.null(purity) || is.null(ploidy)) {
    Clr <- data.table(C=seq(lr.range[1], lr.range[2]), lr=seq(lr.range[1], lr.range[2]))
    ymin <- lr.range[1]
    ymax <- lr.range[2]
  } else {
    vars <- .baseLogRatioVars(purity, as.numeric(ploidy), C.range)
    Clr <- vars$Clr
    ymin <- vars$ymin
    ymax <- vars$ymax
  }
  
  
  
  if(!is.null(sel.chr)) {
    if(nchar(sel.chr) > 2) {
      sel.chr <- substr(sel.chr, 4,nchar(sel.chr))
    }
    seg.plot <- seg.plot[seqnames %in% sel.chr]
  }
  
  p <- seg.plot %>% ggplot()+
    geom_rect(aes(xmin=-Inf, xmax=+Inf, ymin=-Inf, ymax=Inf, fill=backgroundCol)) +
    scale_fill_manual(values=c("#EEEEEE","white"), guide="none")+
    geom_segment(aes(x = start, xend = end, y = lr, yend = lr, color = as.factor(segCol)), lwd = 2)+
    facet_grid(. ~ seqnames, scales = "free_x", space = "free_x")+
    theme_pubclean()+
    theme(panel.spacing = unit(0, "lines"), legend.position = "none", axis.text.x = element_blank(), 
          axis.title.x = element_blank(), axis.ticks.x = element_blank())+
    geom_hline(yintercept = Clr$lr, color = "grey", linetype = "dashed")+
    scale_y_continuous(breaks = Clr$lr, labels = Clr$C)+
    coord_cartesian(ylim = c(min(seg.plot$lr), max(seg.plot$lr)))+
    scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none")+
    coord_cartesian(ylim = c(min(Clr$lr - 0.1), max(Clr$lr + 0.1)))
  
  return(p)
}


library(scattermore)
plotLogRatioScatterMore <- function(pd, purity=NULL, ploidy=NULL, sel.chr=NULL, sel.data="tile", sel.col="segment", lr.range=c(-4,4), C.range=c(0, 8), max.point=1000) {
  ## select data from chromosomes
  if (!is.null(sel.chr)) {
    cov <- pd$cov[chr %in% sel.chr]
    off <- pd$off[chr %in% sel.chr]
    seg <- pd$seg[chr %in% sel.chr]
  } else {
    cov <- pd$cov
    off <- pd$off
    seg <- pd$seg
  }
  ## subsample points
  cov <- cov[sample(nrow(cov))]
  cov <- cov[,head(.SD,max.point),by=seg]
  ## compute plot limits
  
  if (is.null(purity) || is.null(ploidy)) {
    Clr <- data.table(C=seq(lr.range[1], lr.range[2]), lr=seq(lr.range[1], lr.range[2]))
    ymin <- lr.range[1]
    ymax <- lr.range[2]
  } else {
    vars <- .baseLogRatioVars(purity, ploidy, C.range)
    Clr <- vars$Clr
    ymin <- vars$ymin
    ymax <- vars$ymax
  }
  ## base plot (no data)
  plt <- .baseLogRatioPlot(off, Clr, ymin, ymax)
  ## data is tile
  if (sel.data=="tile") {
    if (sel.col=="segment") {
      plt <- plt +
        geom_scattermore(aes(x=start.off, y=lr, col=factor(seg %% 3)), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
        scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none") +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
    } else if (sel.col=="sC") {
      plt <- plt +
        geom_scattermore(aes(x=start.off, y=lr, col=sC), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
        scale_color_gradient2(low="blue", mid="black", high="red", midpoint=2, limits=c(0,6), oob = scales::squish) +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
    } else if (sel.col=="C") {
      plt <- plt +
        geom_scattermore(aes(x=start.off, y=lr, color=as.character(C)), cov[(target)], size=0.75, alpha=1, na.rm=TRUE) +
        scale_color_manual(values=STRING_COL, guide="none") +
        geom_point(aes(x=start.off, y=lr), cov[(!target)], size=0.75, alpha=1, na.rm=TRUE)
    } else if (sel.col=="CK") {
      cov[C>1 & K==0, C:=-1]
      cov[,C:=as.character(C)]
      plt <- plt +
        geom_scattermore(aes(x=start.off, y=lr, color=C, text=gene_names), cov, size=0.75, alpha=1, na.rm=TRUE) +
        scale_color_manual(values=STRING_COL, guide="none")
    } else {
      plt <- plt + 
        geom_scattermore(aes_string(x="start.off", y="lr", col=sel.col), cov, size=0.75, na.rm=TRUE) +
        scale_color_gradient(low="black", high="red", guide="none")
    }
  }
  ## data is segment
  if (sel.data=="segment") {
    if (sel.col=="segment") {
      plt <- plt + 
        geom_segment(aes(x=start.off, xend=end.off, y=lr, yend=lr, col=factor(seg %% 3)), seg, size=2) + 
        scale_color_manual(values=c("#6495ED", "#DD8080", "#CDE2B8"), guide="none")
    } else {
      plt <- plt +
        geom_segment(aes_string(x="start.off", xend="end.off", y="lr", yend="lr", col=lr-Clr[C==2]$lr), seg, size=2) +
        scale_color_gradient(low="black", high="red", guide="none")
    }
  }
  return(plt)    
}
