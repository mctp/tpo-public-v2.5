.annotateFeat <- function(feat, ref, refdata) {
    ## prefer promoter overlap
    hit.prom <- as.data.table(findOverlaps(promoters(feat, 0, 1), ref))
    hit.feat <- as.data.table(findOverlaps(feat, ref))
    hit.feat <- hit.feat[!(queryHits %in% hit.prom$queryHits)]
    hit.feat <- hit.feat[,.SD[1],by=queryHits]
    hit <- rbind(hit.prom, hit.feat)
    ## use row of NAs for missing overlaps
    na.row <- head(refdata, 1)
    na.row[1,] <- NA
    refdata <- rbind(refdata, na.row)
    refidx <- rep(length(ref) + 1, length(feat))
    refidx[hit$queryHits] <- hit$subjectHits
    mcols(feat) <- cbind(mcols(feat), refdata[refidx,])
    return(feat)
}

#' @export
segOut <- function(digest) {
    seg <- digest$seg
    fit <- digest$fit
    mcols(seg) <- fit[,.(seg,C,K,sC=round(sC, 4),lr=round(lr, 4),nlr,naf,block)]
    return(seg)
}

#' Return a genes x C,K GRanges with most exotic C,K chosen for each gene
#'
#' @export
geneOut <- function(digest, gene) {
    fit <- digest$fit
    seg <- digest$seg
    seg.data <- fit[,.(seg,C,K,sC=round(sC,4),lr=round(lr,4),block)]
    hit.feat <- as.data.table(findOverlaps(gene, seg))
    hit.feat <- cbind(hit.feat, seg.data[hit.feat$subjectHits])
    hit <- hit.feat[order(C!=0, C!=1, K!=0, -C),.SD[1],.(queryHits)]
    na.row <- head(seg.data, 1)
    na.row[1,] <- NA
    seg.data <- rbind(seg.data, na.row)
    seg.idx <- rep(length(seg) + 1, length(gene))
    seg.idx[hit$queryHits] <- hit$subjectHits
    mcols(gene) <- cbind(mcols(gene), seg.data[seg.idx,])
    return(gene)
}

#' Return a table of C,K GRanges with all C,K (with duplicates if they exist)
#'
#' @export
allGeneOut <- function(digest, gene) {
    fit <- digest$fit
    seg <- digest$seg
    seg.data <- fit[,.(seg,C,K,sC=round(sC,4),lr=round(lr,4))]
    ol <- findOverlaps(gene,seg)
    out <- gene[queryHits(ol)]
    mcols(out) <- cbind(mcols(out),seg.data[subjectHits(ol),])
    return(out)
}

#' Return deduplicated Gene, C and K from all tiles (only covered genes)
#'
#' @export
geneTileOut <- function(digest, gene) {
    tile <- digest$tile
    #annotate tiles with gene names
    ol <- findOverlaps(tile,genes)
    tile$gene <- NA
    tile$gene[queryHits(ol)] <- genes[subjectHits(ol)]$gene_name
    #extract and dedup the gene,c,k table
    out <- unique(values(tile)[!is.na(tile$gene),c("C","K","gene")])
    #remove any uninteresting duplicates
    dups <- out$gene[duplicated(out$gene)]
    out <- out[-which(out$gene %in% dups & out$C==2 & out$K==1),]
    return(out)
}


.seg.data <- function(tile, var, seg, opts) {
    snp <- var[var$t.PASS]
    snp.segidx <- findOverlaps(snp, seg, select="first", maxgap=opts$tile.shoulder-1)
    snp.seg <- cbind(
        as.data.table(mcols(snp)[,c("SOMATIC", "mask.loose", "mask.strict", "t.GT", "t.AF", "t.DP", "t.PASS")]),
        seg=snp.segidx
    )
    snp.data <- snp.seg[,.(
        mzd=sum(ifelse(t.AF>0.5, t.AF * t.DP, (1-t.AF) * t.DP)) / sum(t.DP)
    ), seg]
    ## tiles
    tile.segidx <- findOverlaps(tile, seg, select="first")
    tile.seg <- cbind(
        as.data.table(mcols(tile)),
        seg=tile.segidx
    )
    tile.data <- tile.seg[,.(
        lr=mean(lr, na.rm=TRUE),
        gc=mean(gc, na.rm=TRUE),
        blacklist=mean(blacklist, na.rm=TRUE),
        unmasked=mean(unmasked, na.rm=TRUE),
        gap=mean(gap, na.rm=TRUE)
    ), seg]
    setkey(tile.data, seg)
    setkey(snp.data, seg)
    seg.data <- merge(tile.data, snp.data, all=TRUE)
    return(seg.data)
}

.annotateSegs <- function(seg, tile, var, opts) {
    tile <- .addBafTile(tile,var,opts)
    ## arms
    tmp <- findOverlaps(seg, tile, select="first")
    seg$arm <- tile[tmp]$arm

    ## tile num
    tmp <- data.table(
        segid = queryHits(findOverlaps(seg, tile, select="all"))
    )
    seg$tile.count <- tmp[, .(n = .N), keyby = .(segid)]$n

    ## seg lr-baf  and their variance
    tmp <- data.table(
        segid = queryHits(findOverlaps(seg, tile, select="all")),
        lr = tile$lr,
        baf = tile$baf,
        hq = tile$hq
    )
    tmp.mean <- tmp[, .(lr.mean = mean(lr,na.rm=TRUE),baf.mean = mean(baf,na.rm=TRUE)), keyby = .(segid)]
    tmp[hq==FALSE, lr:=NA_real_] #remove high variance regions
    tmp[hq==FALSE, baf:=NA_real_] #remove high variance regions
    tmp.var <- tmp[, .(lr.var = sd(lr, na.rm=TRUE),baf.var = sd(baf,na.rm=TRUE)), keyby = .(segid)]
    tmp.median <- tmp[, .(lr.median = median(lr,na.rm=TRUE),baf.median = median(baf,na.rm=TRUE)), keyby = .(segid)]
    seg$lr.median <- tmp.median$lr.median
    seg$baf.median <- tmp.median$baf.median
    seg$lr.var <- tmp.var$lr.var
    seg$baf.var <- tmp.var$baf.var

    ## hq
    tmp <- data.table(
        segid = queryHits(findOverlaps(seg, tile, select="all")),
        hq = tile$hq
    )
    seg$hq <- tmp[, .(hq = sum(hq, na.rm = TRUE)/.N), keyby = .(segid)]$hq

    return(seg)
}

.gisticSegFileInp <- function(cnv.fns, model.fns=NULL,segtype="logratio", opts) {
  #
  # segtype can be "logratio" or "absolute"
  #
  ## Extract each sample data
  if (segtype == "logratio") {
    rows <- foreach(file_num = 1:length(cnv.fns)) %dopar% {
      cnv <- readRDS(cnv.fns[file_num])
      sample <- rep(cnv.fns[file_num],length(cnv$seg))
      chromosome <- as.character(seqnames(cnv$seg))
      start.pos <- start(cnv$seg)
      end.pos <- end(cnv$seg)
      num.markers <- as.numeric(table(findOverlaps(cnv$tile,cnv$seg,select = "first")))
      seg.CN <- sapply((split(cnv$tile$lr,findOverlaps(cnv$tile,cnv$seg,select = "first"))),median,na.rm=TRUE)
      return(tibble(sample,chromosome,start.pos,end.pos,num.markers,seg.CN))
    }
  } else if (segtype == "absolute") {
    rows <- foreach(file_num = 1:length(cnv.fns)) %dopar% {
      cnv <- readRDS(cnv.fns[file_num])
      model <- readRDS(model.fns[file_num])
      cand <- model$fine[order(cand, -iter), .SD[1], cand][order(-L)]
      topcand <- cand[1,]$cand
      seg <- model$segs[[topcand]]
      sample <- rep(cnv.fns[file_num],length(cnv$seg))
      chromosome <- as.character(seqnames(cnv$seg))
      start.pos <- start(cnv$seg)
      end.pos <- end(cnv$seg)
      num.markers <- seg$nlr
      seg.CN <- seg$C
      return(tibble(sample,chromosome,start.pos,end.pos,num.markers,seg.CN))
    }
  } else {
    warning("Wrong segment value type!")
  }

  output <- do.call(rbind, rows)
  ## remove charachters from chr names
  output$chromosome <- .getChrName(output$chromosome)

  return(output)
}

.gisticMarkerFileInp <- function(cnv.fns, opts) { ## start pos of each tile, is the marker
  if (opts$output.same.marker) {
    cnv <- readRDS(cnv.fns[1])
    marker.name <- paste("tile",1:length(cnv$tile),sep=".")
    chromosome <- .getChrName(as.character(seqnames(cnv$tile)))
    marker.pos <- as.integer((end(cnv$tile)+start(cnv$tile))/2)
    return(tibble(marker.name,chromosome,marker.pos))
  } else {
    rows <- foreach(file_num = 1:length(cnv.fns)) %dopar% {
      cnv <- readRDS(cnv.fns[file_num])
      marker.name <- paste("tile",1:length(cnv$tile),sep=".")
      chromosome <- .getChrName(as.character(seqnames(cnv$tile)))
      marker.pos <- as.integer(start(cnv$tile))
      return(tibble(marker.name,chromosome,marker.pos))
    }
    return(do.call(rbind, rows))
  }
}

.gisticCNVFileInp <- function(germlineCNVS,opts) {
  ## after detecting germ cnvs, fill in this function
}

.getChrName <- function(chr.name) {
  chr.name <- substr(chr.name,4,nchar(chr.name))
  chr.name[chr.name == "X"] <- 23
  chr.name[chr.name == "Y"] <- 24
  chr.name <- as.integer(chr.name)
  return(chr.name)
}

#' @export
exportGistic <- function(cnv.fns, model.fns=NULL, segtype="logratio", dest, opts) {
    out.seg <- .gisticSegFileInp(cnv.fns, model.fns, segtype, opts)
    # out.marker <- .gisticMarkerFileInp(cnv.fns, opts)
    # out.cnv <- gisticCNVFileInp(germlineCNVS,opts)
    write.table(out.seg, paste0(dest,".segFile.txt"), sep = "\t", row.names = FALSE, col.names = opts$output.col.name, quote = FALSE,append = FALSE)
    # write.table(out.marker, paste0(dest,".markerFile.txt"), sep = "\t", row.names = FALSE, col.names = opts$output.col.name, quote = FALSE,append = FALSE)
    # write.table(out.cnv, paste0(dest,".cnvFile.txt"), sep = "\t", row.names = FALSE, col.names = opts$output.col.name, quote = FALSE)
}
