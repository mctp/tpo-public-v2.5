sv2.tnscope <- function(var) {
    tmp.g <- geno(var)
    tmp.r <- rowRanges(var)
    tmp.i <- info(var)
    tmp.map <- data.table(
        bnd1=names(var),
        bnd2=info(var)$MATEID
    )
    tmp.map <- tmp.map[grepl("bnd_", bnd1) & grepl("bnd_", bnd2)]
    tmp.map[,
            ID:=ifelse(
                as.integer(str_replace(bnd1, "bnd_", "")) < as.integer(str_replace(bnd2, "bnd_", "")),
                paste(bnd1, bnd2, sep = ":"), paste(bnd2, bnd1, sep = ":")
            )]
    tmp.map$CHR1 <- as.character(seqnames(tmp.r[tmp.map$bnd1]))
    tmp.map$POS1 <- start(tmp.r[tmp.map$bnd1])
    tmp.map$ALT1 <- as.character(unlist(tmp.r[tmp.map$bnd1]$ALT))
    tmp.map$CHR2 <- as.character(seqnames(tmp.r[tmp.map$bnd2]))
    tmp.map$POS2 <- start(tmp.r[tmp.map$bnd2])
    tmp.map$ALT2 <- as.character(unlist(tmp.r[tmp.map$bnd2]$ALT))
    tmp.map$TYPE <- tmp.i[tmp.map$bnd1,]$SVTYPE
    tmp.map$QUAL <- tmp.r[tmp.map$bnd1,]$QUAL
    tmp.map$FILTER <- tmp.r[tmp.map$bnd1]$FILTER
    tmp.map$IMPRECISE <- tmp.i[tmp.map$bnd1,]$IMPRECISE
    SR <- do.call(rbind, tmp.g$AD[,1])
    tmp.map$TRS <- SR[tmp.map$bnd1,1] + SR[tmp.map$bnd2,1]
    tmp.map$TAS <- SR[tmp.map$bnd1,2] + SR[tmp.map$bnd2,2]
    tmp.map$NRS <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$NAS <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$TRP <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$TAP <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$NRP <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$NAP <- rep(NA_integer_, nrow(tmp.map))
    tmp.map$CONTIG <- rep(NA_character_, nrow(tmp.map))
    tbl.bnd <- tmp.map[,.SD[1],ID]
    tbl.bnd$bnd1 <- NULL
    tbl.bnd$bnd2 <- NULL
    sv2 <- tbl.bnd
    return(sv2)
}

sv2.manta <- function(var) {
    tmp.g <- geno(var)
    tmp.sr1 <- tmp.g$SR[,1]
    tmp.sr1[lengths(tmp.sr1)==0] <- list(c(0,0))
    tmp.sr2 <- tmp.g$SR[,2]
    tmp.sr2[lengths(tmp.sr2)==0] <- list(c(0,0))
    tmp.pr1 <- tmp.g$PR[,1]
    tmp.pr1[lengths(tmp.pr1)==0] <- list(c(0,0))
    tmp.pr2 <- tmp.g$PR[,2]
    tmp.pr2[lengths(tmp.pr2)==0] <- list(c(0,0))
    tmp.rr <- unname(cbind(
        do.call(rbind, tmp.sr1),
        do.call(rbind, tmp.sr2),
        do.call(rbind, tmp.pr1),
        do.call(rbind, tmp.pr2)
    ))
    suppressWarnings({
        info(var)$TRS <- tmp.rr[,3]
        info(var)$TAS <- tmp.rr[,4]
        info(var)$NRS <- tmp.rr[,1]
        info(var)$NAS <- tmp.rr[,2]
        info(var)$TRP <- tmp.rr[,7]
        info(var)$TAP <- tmp.rr[,8]
        info(var)$NRP <- tmp.rr[,5]
        info(var)$NAP <- tmp.rr[,6]
    })
    is.deldup <- info(var)$SVTYPE!="BND"
    var.deldup <- var[is.deldup]
    tmp.r <- rowRanges(var.deldup)
    tmp.i <- info(var.deldup)
    tbl.deldup <- data.table(
        ID=names(tmp.r),
        CHR1=as.character(seqnames(tmp.r)),
        ALT1="-",
        POS1=start(tmp.r),
        CHR2=as.character(seqnames(tmp.r)),
        POS2=tmp.i$END,
        ALT2="-",
        TYPE=tmp.i$SVTYPE,
        QUAL=tmp.i$SOMATICSCORE,
        FILTER=tmp.r$FILTER,
        IMPRECISE=tmp.i$IMPRECISE,
        TRS=tmp.i$TRS,
        TAS=tmp.i$TAS,
        NRS=tmp.i$NRS,
        NAS=tmp.i$NAS,
        TRP=tmp.i$TRP,
        TAP=tmp.i$TAP,
        NRP=tmp.i$NRP,
        NAP=tmp.i$NAP,
        CONTIG=tmp.i$CONTIG
    )
    is.bnd <- info(var)$SVTYPE=="BND"
    var.bnd <- var[is.bnd]
    tmp.r <- rowRanges(var.bnd)
    tmp.i <- info(var.bnd)
    tmp.g <- geno(var.bnd)
    tmp.map <- data.table(
        ID=str_sub(names(var.bnd), 1,-3),
        bnd1=names(var.bnd),
        bnd2=unlist(info(var.bnd)$MATEID)
    )
    tmp.map$CHR1 <- as.character(seqnames(tmp.r[tmp.map$bnd1]))
    tmp.map$POS1 <- start(tmp.r[tmp.map$bnd1])
    tmp.map$ALT1 <- unlist(tmp.r[tmp.map$bnd1]$ALT)
    tmp.map$CHR2 <- as.character(seqnames(tmp.r[tmp.map$bnd2]))
    tmp.map$POS2 <- start(tmp.r[tmp.map$bnd2])
    tmp.map$ALT2 <- unlist(tmp.r[tmp.map$bnd2]$ALT)
    tmp.map$TYPE <- tmp.i[tmp.map$bnd1,]$SVTYPE
    tmp.map$QUAL <- tmp.i[tmp.map$bnd1,]$SOMATICSCORE
    tmp.map$FILTER <- tmp.r[tmp.map$bnd1]$FILTER
    tmp.map$IMPRECISE <- tmp.i[tmp.map$bnd1,]$IMPRECISE
    tmp.map$TRS <- tmp.i[tmp.map$bnd1,]$TRS
    tmp.map$TAS <- tmp.i[tmp.map$bnd1,]$TAS
    tmp.map$NRS <- tmp.i[tmp.map$bnd1,]$NRS
    tmp.map$NAS <- tmp.i[tmp.map$bnd1,]$NAS
    tmp.map$TRP <- tmp.i[tmp.map$bnd1,]$TRP
    tmp.map$TAP <- tmp.i[tmp.map$bnd1,]$TAP
    tmp.map$NRP <- tmp.i[tmp.map$bnd1,]$NRP
    tmp.map$NAP <- tmp.i[tmp.map$bnd1,]$NAP
    tmp.map$CONTIG <- tmp.i[tmp.map$bnd1,]$CONTIG
    tbl.bnd <- tmp.map[,.SD[1],ID]
    tbl.bnd$bnd1 <- NULL
    tbl.bnd$bnd2 <- NULL
    sv2 <- rbind(tbl.deldup, tbl.bnd)
    return(sv2)
}
