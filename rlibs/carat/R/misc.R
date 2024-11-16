.debug <- function() {
    if (debug) {
        ## Note: do not run in production
        geno.cols.omit <- c("AF", "AD", "DP", "MQRankSumPS", "ALT_F1R2", "ALT_F2R1")
        geno.cols <- setdiff(names(geno(var)), geno.cols.omit)
        geno.dts <- foreach(col=geno.cols) %do% {
            tmp <- as.data.table(geno(var)[[col]])
            setnames(tmp, paste0(col,c("T","N")))
        }
        geno.dt <- do.call(cbind, geno.dts)
        info.cols.omit <- c("SOMATIC", "GERMLINE", "PATHOGENIC", "CLINVAR_CLNSIG", "COSMIC_CNT", "STR", "SOR",
                            "GNOMAD_AF", "dbSNP_RS", "TLOD", "NLOD", "NLODF", "PV", "CSQ", "CLINVAR_ALLELEID",
                            "MUTATION_HOTSPOT", "CIEND", "CIPOS", "END", "dbSNP_SAO", "UNIPROT_FEATURE")
        info.cols <- setdiff(names(info(var)), info.cols.omit)
        info.dt <- as.data.table(info(var)[,info.cols])
        rep <- cbind(rep, info.dt, geno.dt)
    }
}

.gobj <- function(genome) {
    gname <- tolower(genome)
    if (gname %in% c("hg38", "mm10")) {
        gstyle <- "UCSC"
        gsexchr <- c(X="chrX", Y="chrY")
    } else if (gname %in% c("grch38", "grcm38")) {
        gstyle <- "NCBI"
        gsexchr <- c(X="X", Y="Y")
    } else {
        stop(sprintf("Genome: %s not supported", genome))
    }
    if (genome %in% c("hg38", "grch38")) {
        gseq <- BSgenome.Hsapiens.UCSC.hg38
        gseqi <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
    } else if (genome %in% c("mm10", "grcm38")) {
        gseq <- BSgenome.Mmusculus.UCSC.mm10
        gseqi <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
    }
    ## gseqi
    gseqi <- keepStandardChromosomes(gseqi)
    gseqi <- dropSeqlevels(gseqi, "chrM")
    suppressPackageStartupMessages({
        ## BioC bug
        seqlevelsStyle(gseqi) <- gstyle
        seqlevelsStyle(gseq) <- gstyle
    })
    ## gobj
    gobj <- list(
        seq=gseq,
        seqi=gseqi,
        sexchr=gsexchr
    )
    return(gobj)
}
