.runMosdepth <- function(bam, by, fasta, cores) {
    if (is.character(by) && file.exists(by)) {
        cov.col <- "V5"
    } else if (!is.na(suppressWarnings(as.integer(by)))) {
        cov.col <- "V4"
    } else {
        stop("by should be a file or window-size")
    }
    prefix <- tempfile("mosdepth_")
    ret <- system2("mosdepth", sprintf("-x -f %s -F 1796 -n -t%s -b %s %s %s", fasta, cores, by, prefix, bam))
    out.fn <- list.files(dirname(prefix), paste0(basename(prefix), ".regions.bed.gz$"),
                         full.names = TRUE)
    if (!ret && file.exists(out.fn)) {
        ## error-prone solution
        ## CRAM files can be sorted chr1 chr2 ... or 1 10 ...
        ## mosdepth preserves CRAM sort order which may be
        ## different than the order of the 'by' file
        ## for this to work correctly 'by' needs to be sorted
        ## in exactly the same way
        out.data <- with(fread(out.fn), GRanges(V1, IRanges(V2, V3), cov=get(cov.col)))
        cov <- sort(sortSeqlevels(out.data))$cov
    } else {
        stop(sprintf("mosdepth run failed: %s", ret))
    }
    out.fns <- list.files(dirname(prefix), paste0(basename(prefix)), full.names = TRUE)
    rets <- sapply(out.fns, unlink)
    if (any(rets)) {
        stop("could not remove temp. files")
    }
    return(cov)
}

.runMosdepthMulti <- function(bams, ...) {
    out <- rowSums(mapply(FUN=function(bam) {
        out1 <- .runMosdepth(bam, ...)
    }, bams))
    return(out)
}

.runMosdepthTile <- function(bams, gt, fasta, cores) {
    bed <- tempfile("mosdepth_", fileext=".bed")
    export(granges(gt), bed)
    out <- .runMosdepthMulti(bams, bed, fasta, cores)
    unlink(bed)
    return(out)
}

.runCBS <- function(y, y.weights=NULL, args) {
    n <- length(y)
    chrom <- rep(1, n)
    maploc <- 1:n
    genomdat <- y
    cna <- DNAcopy::CNA(genomdat, chrom, maploc)
    ## this is crazy
    if (!is.null(y.weights)) {
        y.weights[y.weights==0 | is.na(y.weights)] <- .Machine$double.eps
        inp <- c(list(cna, weights=y.weights), args)
    } else {
        inp <- c(list(cna), args)
    }
    capture.output(
        res <- do.call(DNAcopy::segment, inp)
    )
    bkp <- res$output$loc.end[-length(res$output$loc.end)]
    if (length(bkp)==0) {
        bkp <- integer()
    }
    return(bkp)
}

.runSmooth <- function(lr, gr) {
    obj <- CNA(lr, as.integer(seqnames(gr)), floor((start(gr)+end(gr))/2),
               data.type = "logratio",
               sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    return(adj)
}
