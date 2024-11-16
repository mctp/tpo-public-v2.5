.importSentieon <- function(vcf, tidx, nidx, gobj, opts) {
    header <- scanVcfHeader(vcf)
    ## if popaf is provided (CNVEX_POPAF) but missing from VCF file
    ## treat as if not provided (ignored)
    if (!(opts$popaf %in% vcfFields(header)[["info"]])) {
        opts$popaf <- NULL
    }
    param <-ScanVcfParam(
        info=c("AF", opts$popaf),
        geno=c("GT", "AD", "DP", "GQ")
    )
    if (!is.null(opts$which)) {
        vcfWhich(param) <- which
    }
    var <- .importVcf(vcf, param, gobj, opts)
    var.gr <- rowRanges(var)
    var.gr$mask.loose <- info(var)$mask.loose
    var.gr$mask.strict <- info(var)$mask.strict
    var.gr$mask.giab <- info(var)$mask.giab
    var.gr$tumor.only.mask <- info(var)$tumor.only.mask
    if (!is.null(opts$popaf)) {
        ## CNVEX_POPAF is provided and present
        ## population.af is the sum of all ALT's
        population.af <- sapply(info(var)[[opts$popaf]], sum)
        population.af[is.na(population.af)] <- 0
        var.gr$population.af <- population.af
    } else {
        ## CNVEX is not provided or missing filters will restrict
        ## variants population.af if too low or too high. Setting to 
        ## 0.5 means variants will not be filtered.
        var.gr$population.af <- 0.5
    }
    var.gr$SOMATIC <- FALSE
    var.gr$TYPE <- "SNV"
    if (!is.null(tidx)) {
        t.AD <- sapply(geno(var)$AD[,tidx], "[", 2)
        t.DP <- geno(var)$DP[,tidx]
        var.gr$t.GT <- geno(var)$GT[,tidx]
        var.gr$t.AF <- t.AD / t.DP
        var.gr$t.DP <- t.DP
    }
    if (!is.null(nidx)) {
        n.AD <- sapply(geno(var)$AD[,nidx], "[", 2)
        n.DP <- geno(var)$DP[,nidx]
        var.gr$n.GT <- geno(var)$GT[,nidx]
        var.gr$n.AF <- n.AD / n.DP
        var.gr$n.DP <- n.DP
    }
    return(var.gr)
}

getVcfImporter <- function(caller) {
    ## import
    if (is.null(caller)) {
        func <- NULL
    } else if (caller=="sentieon") {
        func <- .importSentieon
    } else {
        stop("Variant caller not supported.")
    }
    return(func)
}
