.fixADTN <- function(var) {
    ## the expand function assumes that the AD field
    ## is populated the same for all samples
    if (ncol(geno(var)$AD)==2) {
        na1 <- is.na(geno(var)$AD[,1])
        len1 <- lengths(geno(var)$AD[,1])
        na2 <- is.na(geno(var)$AD[,2])
        len2 <- lengths(geno(var)$AD[,2])
        stopifnot(all(len1 == len2 | na1 | na2))
        geno(var)$AD[na2,2] <- lapply(len1[na2], function(n) rep(0, n))
        geno(var)$AD[na1,1] <- lapply(len2[na1], function(n) rep(0, n))
        len1 <- lengths(geno(var)$AD[,1])
        len2 <- lengths(geno(var)$AD[,2])
        stopifnot(all(len1 == len2))
    } else {
        var <- var[!is.na(geno(var)$AD)[,1]]
    }
    return(var)
}

.fixDPTN <- function(var) {
    ## sometimes DP value is NA seen only in DNAScope
    ## but adding to TNScope too.
    DPS <- intersect(names(geno(var)), c("DP", "AFDP"))
    if (ncol(geno(var)$AD)==2) {
        for (DP in DPS) {
            na1 <- is.na(geno(var)[[DP]][,1])
            na2 <- is.na(geno(var)[[DP]][,2])
            geno(var)[[DP]][na1,1] <- 0
            geno(var)[[DP]][na2,2] <- 0
        }
    } else {
        var <- var[!is.na(geno(var)[[DPS]])[,1]]
    }
    return(var)
}

.fixADTNE <- function(var) {
    if (ncol(geno(var)$AD)==2) {
        ADT <- geno(var)$AD[,1,2]
        ADT[is.na(ADT)] <- 0
        geno(var)$AD[,1,2] <- ADT
        ADN <- geno(var)$AD[,2,2]
        ADN[is.na(ADN)] <- 0
        geno(var)$AD[,2,2] <- ADN
    }
    return(var)
}

.fixAFTN <- function(var) {
    if (!("DP" %in% names(geno(var)))) {
        ## assumes AFDP is present
        dp.geno <- geno(header(var))["AFDP",]
        rownames(dp.geno) <- "DP"
        geno(header(var)) <- rbind(geno(header(var)), dp.geno)
        geno(var)$DP <- rowSums(geno(var)$AD, dims=2)
    }
    if (!("AF" %in% names(geno(var)))) {
        ## assumes AD and DP is present
        af.geno <- geno(header(var))["DP",]
        af.geno$Type <- "Float"
        af.geno$Description <- "Allele fraction of the event in the tumor"
        rownames(af.geno) <- "AF"
        geno(header(var)) <- rbind(geno(header(var)), af.geno)
        if (ncol(geno(var)$DP)==2) {
            ## Tumor-normal
            ADT <- geno(var)$AD[,1,2]
            ADN <- geno(var)$AD[,2,2]
            DPT <- geno(var)$DP[,1]
            DPN <- geno(var)$DP[,2]
            AF <- cbind(ADT / DPT, ADN / DPN)
            colnames(AF) <- colnames(geno(var)$DP)
            geno(var)$AF <- AF
        } else {
            ## Tumor-Only
            ADT <- geno(var)$AD[,1,2]
            DPT <- geno(var)$DP[,1]
            AF <- as.matrix(ADT / DPT)
            colnames(AF) <- colnames(geno(var)$DP)
            geno(var)$AF <- AF
        }
    }
    return(var)
}

.addContext <- function(var, gobj) {
    ## update header
    ctx.info <- DataFrame(Number=1, Type="String", Description="Trinucleotide Context", row.names="context")
    info(header(var)) <- rbind(info(header(var)), ctx.info)
    sbs96.info <- DataFrame(Number=1, Type="String", Description="SBS96 Class", row.names="SBS96")
    info(header(var)) <- rbind(info(header(var)), sbs96.info)
    if (length(var)==0) {
        info(var)$context <- character(0)
        info(var)$SBS96 <- character(0)
    } else {
        var.ref <- ref(var)
        var.ref[!(var.ref %in% c("A", "C", "G", "T"))] <- "N"
        var.ref <- DNAStringSet(var.ref) # DNAStringSet is required to work with some VCF files
        var.alt <- alt(var)
        var.alt[!(var.alt %in% c("A", "C", "G", "T"))] <- "N" # needed to work with structural
        var.alt <- DNAStringSet(var.alt)
        var.seq <- getSeq(gobj$seq, rowRanges(var)+1)
        flip <- var.ref %in% c("A","G")
        sbs <- width(var.ref)==1 & width(var.alt)==1 & var.ref!="N" & var.alt!="N"
        ## trinuc context
        context <- as.character(var.seq)
        ## Trinuc context not relevant for non-point mutations, so revert them to NA
        context[!sbs] <- NA
        info(var)$context <- context
        ## SBS96 context
        var.ref[(flip)] <- reverseComplement(var.ref[(flip)])
        var.alt[(flip)] <- reverseComplement(var.alt[(flip)])
        var.seq[(flip)] <- reverseComplement(var.seq[(flip)])
        var.seqL <- subseq(var.seq, 1, 1)
        var.seqR <- subseq(var.seq, 3, 3)
        var.ref <- BStringSet(var.ref)
        var.alt <- BStringSet(var.alt)
        var.seqL <- BStringSet(var.seqL)
        var.seqR <- BStringSet(var.seqR)
        var.sbs96 <- xscat(var.seqL, "[", var.ref,">", var.alt, "]", var.seqR)
        var.sbs96 <- as.character(var.sbs96)
        var.sbs96[!sbs] <- NA_character_
        info(var)$SBS96 <- var.sbs96
    }
    return(var)
}

.fixPV <- function(var) {
    if (!("PV" %in% names(info(var)))) {
        ## fisher exact test on read-counts
        pv.info <- DataFrame(Number=1, Type="Float", Description="Fisher's exact test p-value", row.names="PV")
        info(header(var)) <- rbind(info(header(var)), pv.info)
        if (ncol(geno(var)$DP)==2) {
            ADT <- geno(var)$AD[,1,2]
            ADN <- geno(var)$AD[,2,2]
            DPT <- geno(var)$DP[,1]
            DPN <- geno(var)$DP[,2]
            info(var)$PV <- apply(data.table(ADT,DPT-ADT,ADN,DPN-ADN), 1, row.fisher.test)[1,]
        } else {
            info(var)$PV <- 1
        }
    }
    return(var)
}

.csqTable <- function(var) {
    if (!all(unlist(is.na(info(var)$CSQ)))) {
        csq.desc <- info(header(var))["CSQ","Description"]
        csq.fields <- str_split(str_match(csq.desc, "Format: (.*)")[2], "\\|")[[1]]
        csq.vec <- str_split(unlist(info(var)$CSQ), "\\|")
        .csq.sel <- seq(1, length(csq.fields))
        .csq.fun <- function(x) x[.csq.sel]
        csq.tbl <- as.data.table(t(mapply(csq.vec, FUN=.csq.fun)))
        setnames(csq.tbl, csq.fields)
        csq.tbl[,ID:=rep(names(var), lengths(info(var)$CSQ))]
    } else {
        csq.tbl <- data.table(ID=names(var), ALLELE_NUM=1, CSQ=NA_character_)
    }
    return(csq.tbl)
}

#' @export
concatVars <- function(vars, cvcf) {
    vcfs <- sapply(vars, function(var) {
        if(!is.null(var)) {
          out <- tempfile()
          out.gz <- paste0(out, ".gz")
          writeVcf(sort(var), out)
          system2("bgzip", out)
          system2("tabix", out.gz)
          return(out.gz)
        }
    })
    vcfs <- vcfs[!sapply(vcfs,is.null)]
    system2("bcftools", c("concat", "-a -o", c(cvcf, vcfs)))
    unlink(vcfs)
    unlink(paste0(vcfs, ".tbi"))
}

#' @export
writeVcfGz <- function(var, fn.gz) {
    fn <- str_replace(fn.gz, ".gz$", "")
    writeVcf(sort(var), fn)
    system2("bgzip", fn)
    system2("tabix", fn.gz)
}

#' @export
concatVarsFromFile <- function(vcfs, cvcf) {
  system2("bcftools", c("concat", "-a -o", c(cvcf, vcfs)))
  unlink(vcfs)
  unlink(paste0(vcfs, ".tbi"))
}

#' @export
annotateVcfWithManta <- function(tn.var, mn.var, overlap.args=NULL) {
    ## convert TNScope and Manta VCF object to Breakpoint Ranges
    ## only attempt to annotate BND breakends
    tn.gr <- breakpointRanges(tn.var[info(tn.var)$SVTYPE=="BND"])
    mn.gr <- breakpointRanges(mn.var, info_columns=c("SOMATICSCORE","CONTIG"))
    ## find overlaps between TNScope and Manta ranges
    tnmn.hits <- do.call(findBreakpointOverlaps, c(list(tn.gr, mn.gr), overlap.args))
    mn.gr.match <- mn.gr[subjectHits(tnmn.hits),c("sourceId","FILTER","SOMATICSCORE","CONTIG")]
    tn.gr.match <- tn.gr[queryHits(tnmn.hits)]
    mn.gr.match$tnscopeId <- names(tn.gr.match)
    ## calculate Manta read-support
    mn.var.match <- mn.var[mn.gr.match$sourceId]
    sr1 <- geno(mn.var.match)$SR[,1]
    sr1[lengths(sr1)==0] <- list(c(0,0))
    nsr <- do.call(rbind, sr1)
    pr1 <- geno(mn.var.match)$PR[,1]
    pr1[lengths(pr1)==0] <- list(c(0,0))
    npr <- do.call(rbind, pr1)
    sr2 <- geno(mn.var.match)$SR[,2]
    sr2[lengths(sr2)==0] <- list(c(0,0))
    tsr <- do.call(rbind, sr2)
    pr2 <- geno(mn.var.match)$PR[,2]
    pr2[lengths(pr2)==0] <- list(c(0,0))
    tpr <- do.call(rbind, pr2)
    manta.reads <- paste(tsr[,2],tsr[,1],tpr[,2],tpr[,1],nsr[,2],nsr[,1],npr[,2],npr[,1],sep=",")
    ## update header
    manta.info.header <- DataFrame(Number=c(1,1,1,8),
                                   Type=c("String","String","String","String"),
                                   Description=c("Manta FILTER", "Manta SOMATICSCORE", "Manta CONTIG",
                                   "Manta READS [TSR_ALT,TSR_REF,TPR_ALT,TPR_REF,NSR_ALT,NSR_REF,NPR_ALT,NPR_REF]"),
                                   row.names=c("MANTA_FILTER", "MANTA_SOMATICSCORE", "MANTA_CONTIG", "MANTA_READS"))
    info(header(tn.var)) <- rbind(info(header(tn.var)), manta.info.header)
    if (length(tn.var)>0) {
        info(tn.var)$MANTA_FILTER <- NA_character_
        info(tn.var[mn.gr.match$tnscopeId])$MANTA_FILTER <- mn.gr.match$FILTER
        info(tn.var)$MANTA_SOMATICSCORE <- NA_character_
        info(tn.var[mn.gr.match$tnscopeId])$MANTA_SOMATICSCORE <- mn.gr.match$SOMATICSCORE
        info(tn.var)$MANTA_CONTIG <- NA_character_
        info(tn.var[mn.gr.match$tnscopeId])$MANTA_CONTIG <- mn.gr.match$CONTIG
        info(tn.var)$MANTA_READS <- NA_character_
        info(tn.var[mn.gr.match$tnscopeId])$MANTA_READS <- manta.reads
    } else {
        info(tn.var)$MANTA_FILTER <- character(0)
        info(tn.var)$MANTA_SOMATICSCORE <- character(0)
        info(tn.var)$MANTA_CONTIG <- character(0)
        info(tn.var)$MANTA_READS <- character(0)
    }
    return(tn.var)
}

#' @export
annotateVcfWithGridss <- function(tn.var, gripss.var, overlap.args) {
    return(tn.var)
}

#' @export
expand2 <- function(var, rename=TRUE) {
    ## This version of expand knows how to
    ## properly deal with AD fields from GATK
    ## and CSQ fields from VEP
    ## CSQs from VEP are actually Type=A, but encoded as Type=.
    csq.tbl <- .csqTable(var)

    ## pack CSQ by ALT
    csq.pack <- data.table(
        ID=csq.tbl$ID,
        ALLELE_NUM=as.integer(csq.tbl$ALLELE_NUM),
        CSQ=unlist(info(var)$CSQ)
    )
    ## for silent/synonymous mutations VEP sometimes provides only
    csq.pack <- csq.pack[is.na(ALLELE_NUM), ALLELE_NUM:=1]
    csq.pack <- csq.pack[,.(CSQA=list(CSQ)),by=.(ID, ALLELE_NUM)]

    ## All variant ALT combinations
    var.alt <- data.table(
        ID=rep(names(var), lengths(alt(var))),
        ALLELE_NUM=NA_integer_
    )
    var.alt[,ALLELE_NUM:=1:nrow(.SD),by=ID]
    var.alt[,idx:=1:nrow(var.alt)]
    setkey(csq.pack, ID, ALLELE_NUM)
    setkey(var.alt, ID, ALLELE_NUM)
    csq.alt <- csq.pack[var.alt]
    setkey(csq.alt, idx)

    ## expand and replace CSQ with per-alt CSQ
    vare <- expand(var)

    ## this function is really terribly broken it
    ## produces ExpandedVcf objects that does not preserve
    ## the data type of 'Type=A' info columns and the output
    ## differs (list or vector) depending on whether there
    ## are multiple ALTs in a VCF
    vinf <- info(header(vare))
    acols <- rownames(vinf[vinf$Number=="A",])
    for (acol in acols) {
        info(vare)[[acol]] <- as(info(vare)[[acol]], "List")
    }

    ## get the correct value for data type 'Type=R', RPA, MBQ, MFRL, MMQ
    ## [[1,2], [3,4,5,6], [7,8,9]] -> [[1,2], [3,4], [3,5], [3,6], [7,8], [7,9]]
    rcols <- rownames(vinf[vinf$Number=="R",])
    for (rcol in rcols) {
        info(vare)[[rcol]] <- unlist(lapply(as.list(info(var)[[rcol]]), function(e) {
            lapply(e[-1], function(i) c(e[1], i))
            }),recursive = FALSE)
    }

    ## make sure we have AF and DP
    info(vare)$CSQ <- CharacterList(csq.alt$CSQA)
    if (rename) {
        ## this step duplicates columsn in rowRanges(vare), really
        names(vare) <- sprintf("%s:%s_%s/%s",
                               as.character(seqnames(vare)),
                               as.character(start(vare)),
                               as.character(ref(vare)),
                               as.character(alt(vare))
                               )
    }
    return(vare)
}

#' @export
importVcf <- function(vcf, gobj, prune=TRUE) {
    var <- readVcf(vcf, genome = unique(genome(gobj$seqi)))
    ## readVcf is broken it does not preserve all of CSQ
    gzvcf <- gzfile(vcf)
    ln <- readLines(gzvcf)
    close(gzvcf)
    ln <- grep("^#", ln, value=TRUE, invert=TRUE)
    vcsq <- str_match(ln, "CSQ=([^\t]*)\t")[,2]
    if (!("CSQ" %in% rownames(info(header(var))))) {
        csq.row  <- new("DFrame", rownames = "CSQ", nrows = 1L,
                      listData = list(Number = ".", Type = "String", Description = "Consequence"))
        info(header(var)) <- rbind(info(header(var)), csq.row)
    }
    info(var)$CSQ <- CharacterList(str_split(vcsq, ","))
    if (prune) {
        shared.levels <- intersect(seqlevels(var), seqlevels(gobj$seqi))
        var <- keepSeqlevels(var, shared.levels, pruning.mode="coarse")
        seqlevels(var) <- seqlevels(gobj$seqi)
    }
    gc()
    return(var)
}

#' @export
importVcfReport <- function(vcf, gobj, context=TRUE, prune=TRUE, rename=TRUE) {
    var <- importVcf(vcf, gobj, prune)
    var <- .fixADTN(var)
    var <- .fixDPTN(var)
    if(length(var)>0) {
      var <- expand2(var, rename)
      var <- .fixADTNE(var)
      var <- .fixAFTN(var)
      var <- .fixPV(var)
    }
    if (context) {
        var <- .addContext(var, gobj)
    }
    return(var)
}
