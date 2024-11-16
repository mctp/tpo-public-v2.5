.gobj <- function(genome, fasta, refs) {
    gname <- tolower(genome)
    if (gname %in% c("hg38", "mm10")) {
        gstyle <- "UCSC"
        gsexchr <- c(X="chrX", Y="chrY")
        gspecies  <- "human"
    } else if (grepl("grc", gname)) {
        gstyle <- "NCBI"
        gsexchr <- c(X="X", Y="Y")
        gspecies <- "mouse"
    } else {
        stop(sprintf("Genome: %s not supported", genome))
    }
    if (gspecies == "human") {
        gseq <- BSgenome.Hsapiens.UCSC.hg38
        gseqi <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
        gpar <- as(c("chrX:1-2781479", "chrX:155701383-156030895"), "GRanges")
    } else if (gspecies == "mouse") {
        gseq <- BSgenome.Mmusculus.UCSC.mm10
        gseqi <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
        gpar <- as(c("chrX:169969759-170931299"), "GRanges")
    }
    ## gseqi
    gseqi <- keepStandardChromosomes(gseqi)
    gseqi <- dropSeqlevels(gseqi, "chrM")
    ## TODO: this masks some bug in Bioconductor
    suppressPackageStartupMessages({
        seqlevelsStyle(gseqi) <- gstyle
    })
    seqlevelsStyle(gseq) <- gstyle
    seqlevelsStyle(gpar) <- gstyle
    ##
    genome(gpar) <- unique(genome(gseq))
    ## gobj
    gobj <- list(
        seq=gseq,
        seqi=gseqi,
        fasta=fasta,
        sexchr=gsexchr,
        par=gpar
    )
    if (refs) {
        gobj$refs <- .grefs(gname, gseqi)
    }
    return(gobj)
}

.grefs <- function(gname, gseqi) {
    if (gname %in% c("hg38", "grch38")) {
        grefns <- list(
            blacklist = system.file("extdata/hg38/blacklist.bed.gz", package="cnvex"),
            cytoband = system.file("extdata/hg38/cytoband.bed.gz", package="cnvex"),
            centromere = system.file("extdata/hg38/centromere.bed.gz", package="cnvex"),
            unmask.loose = system.file("extdata/hg38/1000G-loose-unmask.bed.gz", package="cnvex"),
            unmask.strict = system.file("extdata/hg38/1000G-strict-unmask.bed.gz", package="cnvex"),
            giab.isdifficult = system.file("extdata/hg38/GIAB_alldifficultregions.bed.gz", package="cnvex"),
            tumor.only.mask  = system.file("extdata/hg38/tumorOnly.mask.bed.gz", package="cnvex")
            ## meth = system.file("extdata/hg38/wgbs_clean4.bed.gz", package="cnvex"),
            ## repliseq = system.file("extdata/hg38/repliseq.bed.gz", package="cnvex")
        )
    } else if (gname %in% c("mm10", "grcm38")) {
        grefns <- list(
            blacklist = system.file("extdata/mm10/blacklist.bed.gz", package="cnvex"),
            cytoband = system.file("extdata/mm10/cytoband.bed.gz", package="cnvex"),
            centromere = system.file("extdata/mm10/centromere.bed.gz", package="cnvex"), ## TODO
            unmask.loose = system.file("extdata/mm10/mouse-loose-unmask.bed.gz", package="cnvex"),
            unmask.strict = system.file("extdata/mm10/mouse-strict-unmask.bed.gz", package="cnvex"),
            giab.isdifficult = system.file("extdata/mm10/GIAB_alldifficultregions.bed.gz", package="cnvex"),
            tumor.only.mask  = system.file("extdata/mm10/tumorOnly.mask.bed.gz", package="cnvex")
            ## meth = system.file("extdata/mm10/wgbs_clean4.bed.gz", package="cnvex"),
            ## repliseq = system.file("extdata/mm10/repliseq.bed.gz", package="cnvex")
        )
    }
    grefs <- lapply(grefns, function(refn) {
        .robust.import(refn, gseqi)
    })
    return(grefs)
}
