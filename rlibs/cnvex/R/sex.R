.detect.sex.coverage.n <- function(var, tile, gobj, opts) {
    ## based on the median coverage ratio between X and autosomes
    chrx.cov <- median(tile[tile$blacklist == 0 & tile$target & seqnames(tile) %in% gobj$sexchr["X"]]$n.cov.raw)
    chra.cov <- median(tile[tile$blacklist == 0 & tile$target & seqnames(tile) %ni% gobj$sexchr]$n.cov.raw)
    sex <- ifelse((chrx.cov / chra.cov > opts$sex.x.covratio.threshold), "female", "male")
    return(sex)
}

.detect.sex.coverage.t <- function(var, tile, gobj, opts) {
    ## based on the presence of Y chromosome
    chry.cov <- median(tile[tile$blacklist == 0 & tile$target & seqnames(tile) %in% gobj$sexchr["Y"]]$t.cov.raw)
    sex <- ifelse(chry.cov < opts$sex.y.covzero.threshold, "female", "male")
    return(sex)
}

.detect.sex.snp <- function(tn, var, tile, gobj, opts) {
    ## based on the density of heterozygous SNPs on X
    if (tn=="tumor") {
        sel <- mcols(var)[,"t.GT"] %in% c("0/0/0/1", "0/0/1/1", "0/1/1/1", "0/1") & var$t.PASS
    } else if (tn=="normal") {
        sel <- mcols(var)[,"n.GT"] %in% c("0/0/1/1", "0/1") & var$n.PASS
    }
    chrx.snp.n <- length(var[sel & seqnames(var) %in% gobj$sexchr["X"]])
    chra.snp.n <- length(var[sel & seqnames(var) %ni% gobj$sexchr])
    chrx.tile.n <- length(tile[tile$target & seqnames(tile) %in% gobj$sexchr["X"]])
    chra.tile.n <- length(tile[tile$target & seqnames(tile) %ni% gobj$sexchr])
    chrx.snp.density <- chrx.snp.n / chrx.tile.n
    chra.snp.density <- chra.snp.n / chra.tile.n
    sex <- ifelse((chrx.snp.density / chra.snp.density > opts$sex.x.snpratio.threshold), "female", "male")
    return(sex)
}

.detect.sex <- function(var, tile, gobj, opts) {
    min.x.snp <- 12
    have.n.gt  <- !all(is.na(var$n.GT))
    have.t.gt  <- !all(is.na(var$t.GT))
    if (have.n.gt && (sum(var$n.PASS & seqnames(var) %in% gobj$sexchr["X"])>min.x.snp)) {
        ## 1st preference, heterozygous SNPs in the normal
        sex <- .detect.sex.snp("normal", var, tile, gobj, opts)
    } else if (!all(is.na(tile$n.cov))) {
        ## 2nd preference, coverage in the normal
        sex <- .detect.sex.coverage.n(var, tile, gobj, opts)
    } else if (have.t.gt && (sum(var$t.PASS & seqnames(var) %in% gobj$sexchr["X"])>min.x.snp)) {
        ## 3rd preference, heterozygous SNPs in the tumor
        sex <- .detect.sex.snp("tumor", var, tile, gobj, opts)
    } else {
        ## 4th preference, coverage in the tumor
        sex <- .detect.sex.coverage.t(var, tile, gobj, opts)
    }
    return(sex)
}

.get.sexcopy <- function(sex, tile, gobj) {
    copy <- rep(2L, length(tile))
    if (sex=="male") {
        x.copy <- 1L
        y.copy <- 1L
    } else if (sex=="female") {
        x.copy <- 2L
        y.copy <- 0L
    } else {
        stop("wrong sex.")
    }
    ## keep tile in PAR on chrX diploid
    copy[as.logical(seqnames(tile) %in% gobj$sexchr["X"]) & !(tile %over% gobj$par)] <- x.copy
    copy[as.logical(seqnames(tile) %in% gobj$sexchr["Y"])] <- y.copy
    return(copy)
}
