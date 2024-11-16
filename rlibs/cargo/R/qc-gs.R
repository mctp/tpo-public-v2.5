.gsCordsAlignQC <- function(gs_align) {
    id <- basename(gs_align)
    alig.stat <- file.path(gs_align, paste0(id, "-aligstat.txt"))
    gc.summary <- file.path(gs_align, paste0(id, "-gcsummary.txt"))
    isize.txt <- file.path(gs_align, paste0(id, "-isize.txt"))
    qc <- list(align.id=id)
    stat <- system2("gsutil", c("cat", alig.stat), stdout=TRUE)
    tmp <- data.table::fread(paste(stat, collapse="\n"))[CATEGORY=="PAIR"]
    qc$total.reads <- tmp$TOTAL_READS
    qc$mapped.reads <- tmp$PF_READS_ALIGNED
    qc$mismatch.rate <- tmp$PF_MISMATCH_RATE
    qc$indel.rate <- tmp$PF_INDEL_RATE
    qc$mean.read.length <- tmp$MEAN_READ_LENGTH
    qc$pct.chimeras <- tmp$PCT_CHIMERAS
    qc$strand.balance <- tmp$STRAND_BALANCE
    gc <- system2("gsutil", c("cat", gc.summary), stdout=TRUE)
    tmp <- data.table::fread(paste(gc, collapse="\n"))
    qc$at.dropout <- tmp$AT_DROPOUT
    qc$gc.dropout <- tmp$GC_DROPOUT
    isize <- system2("gsutil", c("cat", isize.txt), stdout=TRUE)
    tmp <- data.table::fread(paste(isize[1:3], collapse="\n"))
    qc$median.insert.size <- tmp$MEDIAN_INSERT_SIZE
    qc$insert.size.mad <- tmp$MEDIAN_ABSOLUTE_DEVIATION
    return(as.data.table((qc)))
}

.gsCordsPostalignQC <- function(gs_postalign) {
    id <- basename(gs_postalign)
    dedup.metrics <- file.path(gs_postalign, paste0(id, "-dedup_metrics.txt"))
    qc <- list(id=id)
    dedup <- system2("gsutil", c("cat", dedup.metrics), stdout=TRUE)
    tmp <- data.table::fread(paste(dedup[1:3], collapse="\n"))
    qc$percent.duplication <- tmp$PERCENT_DUPLICATION
    qc$est.lib.size <- tmp$ESTIMATED_LIBRARY_SIZE
    return(as.data.table((qc)))
}

.gsCordsMiscQC <- function(gs_misc) {
    id <- basename(gs_misc)
    wgsmetrics.fn <- file.path(gs_misc, paste0(id, "-wgsmetrics.txt"))
    qc <- list(id=id)
    wgsmetrics <- system2("gsutil", c("cat", wgsmetrics.fn), stdout=TRUE)
    tmp <- data.table::fread(paste(wgsmetrics[1:3], collapse="\n"))
    qc$mean.coverage <- tmp$MEAN_COVERAGE
    qc$pct.1x  <- tmp$PCT_1X
    qc$pct.5x  <- tmp$PCT_5X
    qc$pct.15x  <- tmp$PCT_15X
    qc$pct.30x  <- tmp$PCT_30X
    qc$pct.60x  <- tmp$PCT_60X
    qc$pct.90x  <- tmp$PCT_90X
    qc$het.snp.sensitivity <- tmp$HET_SNP_SENSITIVITY
    return(as.data.table((qc)))
}

gsCordsQC <- function(id, alig.ids, gs_repo, align=TRUE, postalign=TRUE, misc=TRUE) {
    if (align) {
        align.qc <- rbindlist(lapply(alig.ids, function(alig.id) {
            cbind(id=id, .gsCordsAlignQC(file.path(gs_repo, "cords-align", alig.id)))
        }))
    } else {
        align.qc <- structure(list(align.id = character(0), total.reads = integer(0),
            mapped.reads = integer(0), mismatch.rate = numeric(0), indel.rate = numeric(0),
            mean.read.length = numeric(0), pct.chimeras = numeric(0),
            strand.balance = numeric(0), at.dropout = numeric(0), gc.dropout = numeric(0),
            median.insert.size = numeric(0), insert.size.mad = numeric(0)), row.names = c(NA,
            0L), class = c("data.table", "data.frame"))
    }
    if (postalign) {
        postalign.qc <- .gsCordsPostalignQC(file.path(gs_repo, "cords-postalign", id))
    } else {
        postalign.qc <- structure(list(align.id = character(0), percent.duplication = numeric(0),
            est.lib.size = integer(0)), row.names = c(NA, 0L), class = c("data.table",
            "data.frame"))
    }
    if (misc) {
        misc.qc <- .gsCordsMiscQC(file.path(gs_repo, "cords-misc", id))
    } else {
        misc.qc <- structure(list(align.id = character(0), mean.coverage = numeric(0),
            pct.1x = numeric(0), pct.5x = numeric(0), pct.15x = numeric(0),
            pct.30x = numeric(0), pct.60x = numeric(0), pct.90x = numeric(0),
            het.snp.sensitivity = numeric(0)), row.names = c(NA, 0L), class =
            c("data.table", "data.frame"))
    }
    setkey(align.qc, id)
    setkey(postalign.qc, id)
    setkey(misc.qc, id)
    cords.qc <- Reduce(function(...) merge(..., all = TRUE), list(align.qc, postalign.qc, misc.qc))
    return(cords.qc)
}

.gsCrispAlignQC <- function(gs_align) {
    id <- basename(gs_align)
    alig.log <- file.path(gs_align, paste0(id, "-alig.log"))
    cut.log <- file.path(gs_align, paste0(id, "-cut.log"))
    mrg.log <- file.path(gs_align, paste0(id, "-mrg.log"))
    cut.out <- file.path(gs_align, paste0(id, "-cut.out"))
    mrg.out <- file.path(gs_align, paste0(id, "-mrg.out"))
    qc <- list(align.id=id)
    alig.tmp <- system2("gsutil", c("cat", alig.log), stdout=TRUE)
    qc$total.reads <- as.integer(str_match(alig.tmp[6], "\\d+"))
    qc$mapped.reads <- as.integer(str_match(alig.tmp[9], "\\d+")[,1])
    qc$multi.reads <- as.integer(str_match(alig.tmp[24], "\\d+")[,1])
    qc$total.splices <- as.integer(str_match(alig.tmp[13], "\\d+")[,1])
    cut.tmp <- system2("gsutil", c("cat", cut.log), stdout=TRUE)
    qc$total.precut <- as.integer(str_match(cut.tmp[2], "\\d+")) / 2
    qc$rrna.reads <- as.integer(str_match(cut.tmp[3], "\\d+")) / 2
    cut.tmp2 <- system2("gsutil", c("cat", cut.out), stdout=TRUE)
    qc$trimmed.reads <- as.integer(str_match(cut.tmp2[22], "\\d+")) / 2
    mrg.tmp <- system2("gsutil", c("cat", mrg.log), stdout=TRUE)
    qc$median.insert <- as.integer(str_match(mrg.tmp[2], "\\d+"))
    mrg.tmp2 <- system2("gsutil", c("cat", mrg.out), stdout=TRUE)
    qc$unique.kmers <- as.integer(str_match(mrg.tmp2[30], "\\d+"))
    qc$joined.reads <- as.integer(str_match(mrg.tmp2[38], "\\d+"))
    return(as.data.table(qc))
}

.gsCrispCodacQC <- function(gs_codac) {
    id <- basename(gs_codac)
    qc <- list(id=id)
    ##
    stat.rds <- file.path(gs_codac, paste0(id, "-stat-rep.rds"))
    tmp <- tempfile()
    alig.tmp <- system2("gsutil", c("-q", "cp", stat.rds, tmp))
    stat <- readRDS(tmp)
    qc <- c(qc, stat$bpt.stat[c("art.rte", "dup.rte", "lig.rte")])
    return(as.data.table(qc))
}

.gsCrispQuasrQC <- function(gs_quasr) {
    id <- basename(gs_quasr)
    qc <- list(id=id)
    ##
    count.summary <- file.path(gs_quasr, paste0(id, "-count.summary"))
    count.tmp <- system2("gsutil", c("cat", count.summary), stdout=TRUE)
    qc$counted.reads <- as.integer(str_match(count.tmp[2], "\\d+"))
    qc$uncounted.multi.reads <- as.integer(str_match(count.tmp[8], "\\d+")) + as.integer(str_match(count.tmp[13], "\\d+"))
    qc$uncounted.nonexon.reads <- as.integer(str_match(count.tmp[11], "\\d+"))
    return(as.data.table(qc))
}

gsCrispQC <- function(id, alig.ids, gs_repo, align=TRUE, codac=TRUE, quasr=TRUE) {
    if (align) {
        align.qc <- rbindlist(lapply(alig.ids, function(alig.id) {
            cbind(id=id, .gsCrispAlignQC(file.path(gs_repo, "crisp-align", alig.id)))
        }))
    } else {
        align.qc <- structure(list(id = character(0), align.id=character(0), total.reads = integer(0),
            mapped.reads = integer(0), multi.reads = integer(0), total.splices = integer(0),
            total.precut = numeric(0), rrna.reads = numeric(0), trimmed.reads = numeric(0),median.insert = integer(0),
            unique.kmers = integer(0), joined.reads = integer(0)), row.names = c(NA,
            0L), class = c("data.table", "data.frame"))
    }
    if (codac) {
        codac.qc <- .gsCrispCodacQC(file.path(gs_repo, "crisp-codac", id))
    } else {
        codac.qc <- structure(list(id = character(0), art.rte = numeric(0), dup.rte = numeric(0),
            lig.rte = numeric(0)), row.names = c(NA, 0L), class = c("data.table", "data.frame"))
    }
    if (quasr) {
        quasr.qc <- .gsCrispQuasrQC(file.path(gs_repo, "crisp-quasr", id))
    } else {
        quasr.qc <- structure(list(id = character(0), counted.reads = integer(0),
            uncounted.multi.reads = integer(0), uncounted.nonexon.reads = integer(0)),
            row.names = c(NA, 0L), class = c("data.table", "data.frame"))
    }
    setkey(align.qc, id)
    setkey(codac.qc, id)
    setkey(quasr.qc, id)
    crisp.qc <- Reduce(function(...) merge(..., all = TRUE), list(align.qc, codac.qc, quasr.qc))
    return(crisp.qc)
}
