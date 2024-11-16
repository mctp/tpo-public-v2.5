
#### wgii
.cinGii <- function(digest) {
    ## calculating median copy number
    seg <- segOut(digest)
    seg <- seg[!(seqnames(seg) %in% c("chrX", "chrY"))]
    seg <- seg[!is.na(seg$C)]

    C.median <- weighted.median(seg$C, w = as.numeric(width(seg)))

    ## Aberrant copynumber segments
    abseg <- seg[seg$C != C.median,]

    ## calc gii
    gii <- sum(width(abseg), na.rm = TRUE)/sum(width(seg), na.rm = TRUE)

    return(gii)
}

.cinWGii <- function(segs) {
    ## calculating median copy number

    if("data.table" %in% class(segs)) {
      segs <- makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
    }
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]
    segs <- segs[!is.na(segs$C)]

    C.median <- weighted.median(segs$C, w = as.numeric(width(segs)))

    ## gii per chr
    seg.split <- (split(segs,seqnames(segs),drop = TRUE))
    chrgii <- sapply(seg.split, function(seg.chr) {

        abseg <- seg.chr[seg.chr$C != C.median,]
        localgii <- sum(width(abseg))/sum(width(seg.chr))
        return(localgii)
    })

    wgii <- mean(chrgii, na.rm=TRUE)

    return(wgii)
}

.LOHPortion <- function(digest) {
    ## calculating median copy number
    seg <- segOut(digest)
    seg <- seg[!(seqnames(seg) %in% c("chrX", "chrY"))]
    seg <- seg[!is.na(seg$K)]
    C.median <- weighted.median(seg$C, w = as.numeric(width(seg)))

    ## Aberrant copynumber segments
    abseg <- seg[seg$K == 0,]

    ## calc gii
    loh <- sum(width(abseg))/sum(width(seg))

    return(loh)
}




#### LST

.cinLST <- function(segs,gobj,opts) {
    if("data.table" %in% class(segs)) {
      segs <- makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
    }
    ## set up
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]
    segs <- segs[!is.na(segs$sC)]
    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    lst <- .cin.calcLST(segs,gobj,opts)

    return(lst)

}

.cin.calcLST <- function(segs, gobj, opts) {
    seg.split <- split((segs),segs$arm, drop = TRUE)

    LSTs <- lapply(seg.split, function(seg.chr) {
        ## filter for seg length
        seg.chr <- seg.chr[width(seg.chr) > opts$cin.lst.min.seglen]
        if (length(seg.chr) > 0) {
            ## find distance
            ## TODO: CHECK FOR SEG.CHR WITH ONE MEMBER
            seg.chr <- .cin.lst.distance(seg.chr)
            seg.chr <- .cin.lst.imbalance(seg.chr)
        }

        return(seg.chr)
    })

    segs <- unlist(as(LSTs, "GRangesList"))

    valid.bps <- .cin.findValidLSTbp(segs, opts)

    return(length(valid.bps))
}

.cin.lst.distance <- function(segs) { ## copied from gUtils
    # distance of row i from row i+1
    ix1 <- head(seq_along(segs), -1)
    ix2 <- tail(seq_along(segs), -1)
    dist <- GenomicRanges::distance(segs[ix1], segs[ix2])

    segs$dist <- NA
    segs$dist[ix1] <- dist

    return(segs)
}

.cin.lst.imbalance <- function(segs) { ## finding if two segments have different C or K (i.e AI)
    ix1 <- head(seq_along(segs), -1)
    ix2 <- tail(seq_along(segs), -1)

    AIs <- ((segs[ix1]$C != segs[ix2]$C) | (segs[ix1]$K != segs[ix2]$K))
    AIs[is.na(AIs)] <- FALSE

    segs$isLST <- FALSE
    segs$isLST[ix1] <- AIs

    return(segs)
}

.cin.findValidLSTbp <- function(segs, opts) {
    bps <- which((segs$dist < opts$cin.lst.min.segdist) &
                     (segs$isLST == TRUE))
    return(bps)
}




#### tAI

.cinNtAI <- function(segs, gobj, opts) {
    if("data.table" %in% class(segs)) {
      segs <- makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
    }
    ## set up
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]

    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    AIs <- .cin.calcAIs(segs, gobj, opts)

    return(sum(AIs[,tAI]))
}

.cin.calcAIs <- function(segs, gobj, opts) {
    C.median <- weighted.median(segs$C, w = as.numeric(width(segs)))
    seg.split <- split(as.data.table(segs),seqnames(segs), drop = TRUE)

    AIs <- lapply(seg.split, function(seg.chr) {
        if (C.median %% 2 == 0) {
            seg.chr[, AI := fifelse(K == C-K, FALSE, TRUE, na = FALSE)]
        } else {
            seg.chr[, AI := fifelse(K != 0 & C == C.median ,FALSE, TRUE, na = FALSE)]
        }
        seg.chr <- .cin.len.filter(seg.chr, opts$cin.NtAI.min.len)
        seg.chr <- .cin.isArmLevel(seg.chr, opts$cin.arm.level.portion)
        seg.chr <- .cin.NaAI.isOnTel(seg.chr, gobj, opts$cin.telomere.len)

        seg.chr[, tAI := len.filter.pass & !armLevel & onTel & AI]
        return(seg.chr)
    })
    return(rbindlist(AIs))
}

.cin.len.filter <- function(seg, min.len) {
    # filter the length of the segment
    if ("centromere" %in% colnames(seg)) {
        seg[,len.filter.pass := (width*(1-centromere) > min.len)]
    } else {
        seg[,len.filter.pass := (width > min.len)]
    }

    return(seg)
}

.cin.isArmLevel <- function(seg, cin.arm.level.portion) {
    seg[,armLevel := (width/sum(width)) > cin.arm.level.portion, by=arm]
    return(seg)
}

.cin.NaAI.isOnTel <- function(seg, gobj, cin.telomere.len) {
    seg[, onTel := (start < cin.telomere.len | seqlengths(gobj$seqi)[seqnames]-cin.telomere.len < end)]
    return(seg)
}




#### loh
.cinNLOH <- function(segs, gobj, opts) {
    if("data.table" %in% class(segs)) {
      segs <- makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
    }
    ## set up
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]

    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    LOHs <- .cin.calcLOH(segs, gobj, opts)

    return(sum(LOHs[,LOH]))
}

.cin.calcLOH <- function(segs, gobj, opts) {
    seg.split <- split(as.data.table(segs),seqnames(segs), drop = TRUE)
    LOHs <- lapply(seg.split, function(seg.chr) {
        seg.chr <- .cin.len.filter(seg.chr, opts$cin.LOH.min.len)
        seg.chr <- .cin.isArmLevel(seg.chr,opts$cin.arm.level.portion)
        seg.chr <- .cin.LOH.isLOH(seg.chr)

        seg.chr[, LOH := (isLOH & !armLevel & len.filter.pass)]

        return(seg.chr)
    })
    return(rbindlist(LOHs))
}

.cin.LOH.isLOH <- function(seg) {
    seg[,isLOH := K==0 & C!= 0]
    seg[is.na(isLOH), isLOH := FALSE]

    return(seg)
}




#### Shannon diversity score

.cinSI <- function(segs, gobj, opts) {
  ## set up
  if("data.table" %in% class(segs)) {
    segs <- makeGRangesFromDataFrame(segs, keep.extra.columns = TRUE)
  }
  segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]
  
  ## add refs
  segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)
  
  ## make buckets
  buckets <- .cin.SI.makeBucket(opts$cin.SI.bucket.min, opts$cin.SI.bucket.max, opts$cin.SI.bucket.len)
  
  ## calculating the measure
  SI <- .cin.SI.calcSI(segs, buckets, opts)
  
  return(SI)
}

.cin.SI.makeBucket <- function(start, end, width) {
    bps <- seq(start, end, width)
    buckets <- data.table(start = head(bps, -1), end = tail(bps, -1))
    return(buckets)
}

.cin.SI.calcSI <- function(segs, buckets, opts) {
    # sc %in% [start, end)
    segs <- segs[!is.na(segs$sC)]
    buckets[, p := sum(segs$sC >= start & segs$sC < end)/length(segs), by=1:nrow(buckets)]
    buckets[, sSI := (-1 * p * log(p))]

    SI <- sum(buckets$sSI, na.rm=TRUE)

    return(SI)
}



#### Focal gains

#' Computes the number of focal gains
#'
#' @param DIGEST digest/vault$tables$segment
#' @return integer of focal gains
#' 
#' @export
.cinFG <- function(DIGEST, PLOIDY){
  DIGEST[['tally']] <- ifelse( (DIGEST[['width']]<=7e6) & (DIGEST[['C']]>=(PLOIDY+0.75)) & (DIGEST[['C']]>=3), 1, 0 )
  fg = sum(DIGEST[['tally']])
  return(fg)
}



#### Misc

.cinAiCNA <- function(digest, gobj, opts) {
    ## set up
    segs <- segOut(digest)
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]

    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    s_aicna <- .cin.AiCNA.calcSAiCNA(segs, gobj, opts)

    return(s_aicna)
}

.cin.AiCNA.calcSAiCNA <- function(segs, gobj, opts) {
    AIs <- .cin.calcAIs(segs, gobj, opts)
    AIs <- AIs[!is.na(AIs$C)]

    AiCNA <- AIs[, AiCNA := (C > 0 & AI)] ## expanded version of paper's score

    AiCNA.p <- sum(AiCNA[AiCNA==TRUE,width])/sum(AiCNA[,width])
    AiCNA.n <- nrow(AiCNA[AiCNA==TRUE & width > opts$cin.AiCNA.min.len & !armLevel,])
    AiCNA.s <- AiCNA.p * AiCNA.n

    return(AiCNA.s)

}

.cinAbCNA <- function(digest, gobj, opts) {
    ## set up
    segs <- segOut(digest)
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]

    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    s_abcna <- .cin.AbCNA.calcSAbCNA(segs, gobj, opts)

    return(s_abcna)
}

.cin.AbCNA.calcSAbCNA <- function(segs, gobj, opts) {
    AIs <- .cin.calcAIs(segs, gobj, opts)
    AIs <- AIs[!is.na(AIs$C)]

    AbCNA <- AIs[, AbCNA := (!AI)] ## expanded version of paper's score

    AbCNA.n <- nrow(AbCNA[AbCNA==TRUE & width > opts$cin.AbCNA.min.len,])
    AbCNA.s <- AbCNA.n

    return(AbCNA.s)

}

.cinCnLOH <- function(digest, gobj, opts) {
    ## set up
    segs <- segOut(digest)
    segs <- segs[!(seqnames(segs) %in% c("chrX", "chrY"))]

    ## add refs
    segs <- .cin.addRefs(segs, gobj, addcentromere = TRUE, addarm = TRUE)

    s_cnloh <- .cin.CnLOH.calcSCnLOH(segs, gobj, opts)

    return(s_cnloh)
}

.cin.CnLOH.calcSCnLOH <- function(segs, gobj, opts) {
    AIs <- .cin.calcAIs(segs, gobj, opts)
    AIs <- AIs[!is.na(AIs$C)]

    median.C <- weighted.median(AIs$C, AIs$width,na.rm = TRUE)
    CnLOH <- AIs[, CnLOH := (C == median.C & K==0)] ## expanded version of paper's score
    CnLOH[is.na(CnLOH), CnLOH := FALSE]

    CnLOH.p <- sum(CnLOH[CnLOH==TRUE,width])/sum(CnLOH[,width])
    CnLOH.n <- nrow(CnLOH[CnLOH==TRUE & width > opts$cin.CnCNA.min.len,])
    CnLOH.s <- CnLOH.p * CnLOH.n

    return(CnLOH.s)

}

.cin.addRefs <- function(segs, gobj, addcentromere=TRUE, addarm=TRUE) {
    if (isTRUE(addcentromere)) {
        ## adding centromere to segs
        tmp <- findOverlaps(segs, gobj$refs$centromere)
        tmp <- data.table(
            seg=queryHits(tmp),
            cent=width(pintersect(segs[queryHits(tmp)], gobj$refs$centromere[subjectHits(tmp)]))
        )
        setkey(tmp, seg)
        tmp <- tmp[J(seq_along(segs))]
        tmp[is.na(cent), cent:=0]
        tmp <- tmp[,.(cent=sum(cent)),by=seg]
        segs$centromere <- tmp$cent / width(segs)
    }

    if (isTRUE(addarm)) {
        ## adding arms to the segs
        cyto <- gobj$refs$cytoband
        cyto$arm <- paste0(seqnames(cyto),str_sub(cyto$name, 1, 1))

        ovlap <- findOverlaps(segs, cyto, select = "first")
        segs$arm <- factor(cyto[ovlap]$arm, levels = unique(cyto[ovlap]$arm))

    }

    return(segs)

}
