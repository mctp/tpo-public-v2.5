#' @export
.setCNVColors <- function(segs,sex){
  #set CNV colors
  CNVEX_COPY_COL <- c(
      K0="#78b09c", KN="#000000", CN="#690033",
      C0="#4500AC", C1="#6B58EE", C2="#000000", C3="#FC9272", C4="#FB6A4A",
      C5="#EF3B2C", C6="#CB181D", C7="#A50F15", C8="#67000D"
  )
  segs$col <- NA
  segs$col[segs$C==3] <- CNVEX_COPY_COL['C3']
  segs$col[segs$C==4] <- CNVEX_COPY_COL['C4']
  segs$col[segs$C==5] <- CNVEX_COPY_COL['C5']
  segs$col[segs$C==6] <- CNVEX_COPY_COL['C6']
  segs$col[segs$C==7] <- CNVEX_COPY_COL['C7']
  segs$col[segs$C>=8] <- CNVEX_COPY_COL['C8']
  segs$col[segs$C<4 & segs$C>1 & segs$K==0 & !segs$chr %in% c('chrX','chrY')] <- CNVEX_COPY_COL['K0']
  segs$col[segs$C==2 & (segs$K==1 | is.na(segs$K))] <- CNVEX_COPY_COL['KN']
  segs$col[segs$C==1] <- CNVEX_COPY_COL['C1']
  segs$col[segs$C==0] <- CNVEX_COPY_COL['C0']
  segs$C[segs$C>6] <- 6

  #color for sex chromosomes
  if(sex=='male'){
    segs$col[segs$C==1 & segs$chr %in% c('chrX','chrY')] <- CNVEX_COPY_COL['KN']
    segs$col[segs$C==2 & segs$chr %in% c('chrX','chrY')] <- CNVEX_COPY_COL['C3']
  }

  return(segs)

}


#' make circos plots of the CNV, SV and Fusion data, as available
#'
#' @param mv a TPO dtVault
#' @param savef [boolean] Save the pdf
#' @param gnm genome [hg38,mm10]
#' @return if save==FALSE, return the circos plot to the default graphics device, if save==TRUE, return nothing
#' 
#' @export
circosPlot <- function(v,savef=TRUE,gnm='hg38'){
  require(circlize)

  #get cases
  cases <- unique(v$meta$id)

  for(case in cases){
    # Sex
    sex <- v$meta %>% filter(id==!!case) %>% .$sex
    # Get the CNV
    segs <- v$tables$segment %>% 
              filter(id %in% case) %>% 
              select(chr,start,end,sC, C, K)
    segs <- .setCNVColors(segs, sex)

    #get the SV
    #remove any which map to an unknown or alternate contig
    species <- ifelse(gnm %in% c('hg38'), 'Homo sapiens','Mus musculus')
    chrs <- GenomeInfoDb::extractSeqlevels(species=species, style='UCSC')
    chrs <- chrs[chrs!='chrM']
    if(sex=='female'){
      chrs <- chrs[chrs!='chrY']
    }
    sv <- v$tables$structural %>% filter(id %in% case) %>% 
            filter(chr1 %in% chrs & chr2 %in% chrs)

    #get the fusions
    fus <- v$tables$fusion %>% filter(id %in% case) %>% 
              select(chr.3,pos.3,chr.5,pos.5) %>% 
              filter(chr.3 %in% chrs & chr.5 %in% chrs)

    # save or write to device
    if(savef){
      pdf(file=sprintf('%s-circos.pdf',case))
    }
      
    # Initialize
    circos.par("track.height" = 0.06)
    circos.initializeWithIdeogram(species=gnm, chromosome.index=chrs)

    if(nrow(segs)>0){
      circos.genomicTrack(segs, panel.fun = function(region, value, ...) {
              circos.genomicLines(region, value, numeric.column='sC', type='segment', col = value$col, lwd=2.5, ...)
      })
    }

    # Add SV
    if(nrow(sv)>0){
      #translocation
      bed1 <- sv  %>% filter(topo=='translocation') %>% select(chr1,pos1)
      bed2 <- sv  %>% filter(topo=='translocation') %>% select(chr2,pos2)
      circos.genomicLink(bed1, bed2, col = 'black', border = NA)
      #deletion
      bed1 <- sv %>% filter(topo=='deletion') %>% select(chr1,pos1) 
      bed2 <- sv %>% filter(topo=='deletion') %>% select(chr2,pos2) 
      circos.genomicLink(bed1, bed2, col = 'dodgerblue', border = NA, h=0.05)
      #duplication
      bed1 <- sv %>% filter(topo=='duplication') %>% select(chr1,pos1) 
      bed2 <- sv %>% filter(topo=='duplication') %>% select(chr2,pos2) 
      circos.genomicLink(bed1, bed2, col = 'firebrick', border = NA, h=0.05)
    }

    # Add Fusions
    if(nrow(fus)>0){
      fustmp <- fus %>% filter(chr.3!=chr.5)
      bed1 <- fustmp %>% select(chr.3,pos.3) 
      bed2 <- fustmp %>% select(chr.5,pos.5) 
      circos.genomicLink(bed1, bed2, col = 'darkgoldenrod', border = NA, h=.4, lty='dashed')
    }

    title(case)
    circos.clear()
    if(savef){
      dev.off()
    }
  }
}


#' @export
.cnvSegToBin <- function(x,bins){
   #get a GR of all bins split up if they contain a segment break
   d <- disjoin(c(x,bins))
   d2 <- subsetByOverlaps(d,x)
   d2$id <- unique(x$id)
   # annotate split bins with tile name 
   ol <- findOverlaps(d2,bins)
   d2$name <- NA
   d2$name[queryHits(ol)] <- bins$name[subjectHits(ol)]
   # annotate split bins with C,K
   ol <- findOverlaps(d2,x)
   d2$seg[queryHits(ol)] <- x$seg[subjectHits(ol)]
   d2$C[queryHits(ol)] <- x$C[subjectHits(ol)]
   d2$K[queryHits(ol)] <- x$C[subjectHits(ol)]
   d2$width <- width(d2)
   # Quantify average C per bin
   out <- as.data.frame(d2) %>% 
            select(chr=seqnames,id,name,C,width) %>% 
            group_by(chr,id,name) %>%
            summarize(avgC=sum(C*width)/sum(width)) %>%
            ungroup()
  #output a copy of bins, annotated with average C per bin
  ot <- bins
  ot$chr <- seqnames(ot)
  elementMetadata(ot) <- merge(elementMetadata(ot), out, by=c('chr','name'), all.x=T)
  return(ot)
}

#' @export
.prepareCohortCNV <- function(mv,bins){
  #quantify CNV by tile or cytoband
  cnvgr <- with(mv$tables$segment,GRanges(seqnames=chr, ranges=IRanges(start=start, end=end)))
  cnvgr$C <- mv$tables$segment$C
  cnvgr$K <- mv$tables$segment$C
  cnvgr$id <- mv$tables$segment$id
  cnvgr$seg <- mv$tables$segment$seg

  #convert segment space to bin space
  cnvgr_s <- split(cnvgr, cnvgr$id)
  binseg <- GRangesList(sapply(cnvgr_s,function(x){.cnvSegToBin(x,bins)}))

  #unlist and munge
  allbins <- unlist(binseg)
  names(allbins) <- paste0(allbins$id,seqnames(allbins),allbins$name)
  allbins <- as.data.frame(allbins) %>% 
              merge(mv$meta %>% select(id,ploidy,sex), by='id') %>%
              filter(!(chr=='chrY' & sex=='female')) %>%
              select(chr=seqnames,start,end,id,name,avgC,ploidy,sex)

  # to make it easier, double C on X,Y for males
  w <- with(allbins, which(sex=='male' & chr %in% c('chrX','chrY')))
  allbins$avgC[w] <- 2*allbins$avgC[w]

  #add colors
  cnv_col_simple <- circlize::colorRamp2(breaks=c(0,2,4,8), colors=c("#4500AC","black","#FB6A4A","#67000D"))
  allbins$col <- cnv_col_simple(allbins$avgC)

  return(allbins)
}

#' @export
.prepareCohortSV <- function(mv,bins){

  sv <- mv$tables$structural %>%
    filter(chr1 %in% seqnames(bins) & chr2 %in% seqnames(bins)) %>%
    mutate(chr1=as.character(chr1), chr2=as.character(chr2)) 
  svgr1 <- with(sv, GRanges(seqnames=chr1, ranges=IRanges(start=pos1, width=1)))
  svgr2 <- with(sv, GRanges(seqnames=chr2, ranges=IRanges(start=pos2, width=1)))

  #annotate with the bin
  .bincenter <- function(x){(end(x)+start(x))/2}
  ol <- findOverlaps(svgr1, bins)
  sv$bin1[queryHits(ol)] <- bins$name[subjectHits(ol)]
  sv$bincenter1[queryHits(ol)] <- .bincenter(bins[subjectHits(ol)])

  ol <- findOverlaps(svgr2, bins)
  sv$bin2[queryHits(ol)] <- bins$name[subjectHits(ol)]
  sv$bincenter2[queryHits(ol)] <- .bincenter(bins[subjectHits(ol)])

  #tabulate by bins 
  sv <- sv %>% tidyr::unite(SV,bin1,bin2,topo,sep='|',remove=F)
  recur <- table(sv$SV)
  sv <- merge(sv, as.data.frame(recur), by.x='SV', by.y='Var1')
  #
  return(sv)
}

#' @export
.prepareCohortFusion <- function(mv,bins){
  fus <- mv$tables$fusion %>%
    filter(chr.5 %in% seqnames(bins) & chr.3 %in% seqnames(bins)) %>%
    mutate(chr.5=as.character(chr.5), chr.3=as.character(chr.3)) 
  fusgr3 <- with(fus, GRanges(seqnames=chr.5, ranges=IRanges(start=pos.5, width=1)))
  fusgr5 <- with(fus, GRanges(seqnames=chr.3, ranges=IRanges(start=pos.3, width=1)))

  #annotate with the bin
  .bincenter <- function(x){(end(x)+start(x))/2}
  ol <- findOverlaps(fusgr5, bins)
  fus$bin5[queryHits(ol)] <- bins$name[subjectHits(ol)]
  fus$bincenter5[queryHits(ol)] <- .bincenter(bins[subjectHits(ol)])

  ol <- findOverlaps(fusgr3, bins)
  fus$bin3[queryHits(ol)] <- bins$name[subjectHits(ol)]
  fus$bincenter3[queryHits(ol)] <- .bincenter(bins[subjectHits(ol)])

  #tabulate by bins
  fus <- fus %>% tidyr::unite(fus,bin5,gene_names.5.1,gene_names.5.2,bin3,gene_names.3.1,gene_names.3.2,sep='|',remove=F)
  recur <- table(fus$fus)
  fus <- merge(fus, as.data.frame(recur), by.x='fus', by.y='Var1')

  return(fus)
}

#' make circos plots of the CNV, SV and Fusion data, summarizing recurrent events
#'
#' @param mv a TPO dtVault 
#' @param bin_size currently only supports 'cytoband'
#' @param file [NULL] pdf to save, if NULL will use default graphics device
#' @param cnv,sv,fusion booleans for which data types to display [Default: TRUE]
#' @return if save==FALSE, return the circos plot to the default graphics device, if save==TRUE, return nothing
#' 
#' @export
cohortCicosPlot <- function(mv, bin_size='cytoband', file=NULL, title=NULL, cnv=TRUE,sv=TRUE,fusion=TRUE, gnm='hg38'){
  #setup the bins
  if(bin_size=='cytoband'){
    bins = rtracklayer::import.bed(system.file("extdata/hg38/cytoband.bed.gz", package="cnvex"))
    chrs <- GenomeInfoDb::extractSeqlevels(species='Homo sapiens', style='UCSC')
    chrs <- chrs[chrs!='chrM']
    bins <- bins[seqnames(bins) %in% chrs]
    bins$chr <- seqnames(bins)
    elementMetadata(bins) <- elementMetadata(bins) %>% as.data.frame() %>% tidyr::unite(name,chr,name,sep='-')
  }else{
    #todo: bin size
    stop("only 'cytoband' is supported")
  }
  


  #save pdf
  if(!is.null(file)){
    pdf(file=file)
  }
    
  # Initialize circos
  circos.par("track.height" = 0.1, "gap.degree"=1.5, "track.margin"=c(0.005,0.005))
  circos.initializeWithIdeogram(
    species=gnm, 
    chromosome.index=chrs,
    plotType=c('ideogram','labels')
  )

  # Add CNV
  if(cnv){
    # Get interesting CNV
    plotcnv <- .prepareCohortCNV(mv,bins)
    plotcnv$col <- colorspace::adjust_transparency(plotcnv$col, 0.6)
    plotcnv$avgC[plotcnv$avgC>8] <- 8
    if(nrow(plotcnv)>0){
      circos.genomicTrack(plotcnv, ylim=c(0,8), panel.fun = function(region, value, ...) {
              circos.genomicLines(region, value, 
                numeric.column='avgC', type='segment', col = value$col, lwd=1, ...
                )
              circos.yaxis(
                at=c(0,2,4,8),
                labels.cex=0.01,
                tick=T,side='left',
                lwd=0.5,
                tick.length= convert_x(0.5, "mm")
              )

      })
      circos.yaxis(
        at=c(0,2,4,8),
        labels.cex=0.3,
        labels.col=c("#67000D","black","#FB6A4A","#67000D"),
        tick=T,
        side='left',
        lwd=0.5,
        sector.index='chr5',
        tick.length= convert_x(0.125, "mm")
      )
    }
  }else{plotcnv <- NULL}

  # Add SV
  if(sv){
    # Get interesting SV
    plotsv <- .prepareCohortSV(mv,bins) %>%
      filter(AFT>0.05 & ADT>7) %>%
      filter(triage=='PASS')

    if(nrow(plotsv)>0){
      #deletion
      del <- plotsv %>% 
        filter(topo=='deletion') %>% 
        select(chr=chr1,start=bincenter1,end=bincenter1,SV,Freq) 

      #duplication
      amp <- plotsv %>% 
        filter(topo=='duplication') %>% 
        select(chr=chr1,start=bincenter1,end=bincenter1,SV,Freq)

      #inversion
      inv <- plotsv %>% 
        filter(topo=='inversion') %>% 
        select(chr=chr1,start=bincenter1,end=bincenter1,SV,Freq)

      tloc1 <- plotsv %>% 
        filter(topo=='translocation') %>% 
        select(chr=chr1,start=bincenter1,end=bincenter1,SV)
      tloc2 <- plotsv %>% 
        filter(topo=='translocation') %>% 
        select(chr=chr2,start=bincenter2,end=bincenter2,SV)
      tloc <- rbind(tloc1,tloc2)


      #amp/del density
      to_plot <- list()
      col <- c()
      if(nrow(del)>0){
        to_plot <- c(to_plot, list(del))
        col <- c(col,colorspace::adjust_transparency('dodgerblue',0.5))
      }
      if(nrow(amp)>0){
        to_plot <- c(to_plot, list(amp))
        col <- c(col,colorspace::adjust_transparency('firebrick',0.5))
      }
      if(nrow(inv)>0){
        to_plot <- c(to_plot, list(inv))
        col <- c(col,colorspace::adjust_transparency('darkgreen',0.5))
      }
      circos.genomicDensity(
        data=to_plot, 
        col = col,
        track.height=0.08,
        count_by='number'
      )
      circos.genomicDensity(
        data=tloc,
        col=colorspace::adjust_transparency('purple4',0.5),
        track.height=0.08,
        count_by='number'
      )


      #translocation
      tl <- plotsv %>% 
        filter(topo=='translocation') %>% 
        select(chr1,bincenter1,chr2,bincenter2,SV,Freq) %>% 
        unique()
      bed1 <- tl %>% select(chr1,pos=bincenter1)
      bed2 <- tl %>% select(chr2,pos=bincenter2)
      circos.genomicLink(bed1, bed2,
        col = colorspace::adjust_transparency('black',0.8),
        lwd=2*sqrt(tl$Freq/max(tl$Freq)), border = NA
      )
    }
  }else{plotsv <- NULL}

  # Add Fusions
  if(fusion){
    # Get interesting Fusions
    plotfus <- .prepareCohortFusion(mv,bins) %>%
      filter(gmap.valid | mm2.valid) %>%
      filter(tot.jnc+tot.enc>7) %>%
      filter(hq.bpt)
    if(nrow(plotfus)>0){
      # same chr
      fustmp <- plotfus  %>% 
        select(chr.3,bincenter3,chr.5,bincenter5,Freq) %>% 
        unique() %>% 
        mutate(lwd=2*sqrt(Freq/max(Freq)))
      fustmp1 <- fustmp %>% filter(chr.5==chr.3)
      bed1 <- fustmp1 %>% select(chr.3,bincenter3) 
      bed2 <- fustmp1 %>% select(chr.5,bincenter5) 
      circos.genomicLink(bed1, bed2, 
        col = colorspace::adjust_transparency('darkgoldenrod',0.8), 
        border = NA, 
        h=.2, lty='dashed', 
        lwd=fustmp1$lwd
      )
      # different chr
      fustmp2 <- fustmp %>% filter(chr.5!=chr.3)
      bed1 <- fustmp2 %>% select(chr.3,bincenter3) 
      bed2 <- fustmp2 %>% select(chr.5,bincenter5) 
      circos.genomicLink(bed1, bed2, 
        col = colorspace::adjust_transparency('darkgoldenrod',0.8), 
        border = NA, h=.4, 
        lty='dashed', 
        lwd=fustmp2$lwd
      )
    }
  }else{plotfus <- NULL}

  if(!is.null(title)){
    title(title) 
  }
  circos.clear()
  dev.off()
  
  return(list(cnv=plotcnv,sv=plotsv,fusion=plotfus))
}
