.robust.import <- function(fn, seqi, ...) {
    tmp <- import(fn, ...)
    tmp <- keepStandardChromosomes(tmp, pruning.mode="coarse")
    if (length(tmp)>0) {
        tmp <- renameSeqlevels(tmp, mapSeqlevels(seqlevels(tmp), head(seqlevelsStyle(seqi), 1)))
    }
    shared.levels <- intersect(seqlevels(tmp), seqlevels(seqi))
    tmp <- keepSeqlevels(tmp, shared.levels, pruning.mode="coarse")
    seqlevels(tmp) <- seqlevels(seqi)
    seqinfo(tmp) <- seqi
    return(tmp)
}

.robust.import.vcf <- function(fn, param, seqi) {
    var <- readVcf(fn, param=param)
    shared.levels <- intersect(seqlevels(var), seqlevels(seqi))
    var <- keepSeqlevels(var, shared.levels, pruning.mode="coarse")
    seqlevels(var) <- seqlevels(seqi)
    seqinfo(var) <- seqi
    return(var)
}

.parse.overrides  <- function(overrides) {
    if (!is.null(overrides)) {
        tmp <- sapply(str_split(overrides, ",")[[1]], str_split, ":", USE.NAMES=FALSE)
        vals1 <- sapply(tmp, "[", 2)
        vals2 <- suppressWarnings(as.list(as.numeric(vals1)))
        vals2[is.na(vals2)] <- vals1[is.na(vals2)]
        names(vals2) <- sapply(tmp, "[", 1)
    } else {
        vals2 <- list()
    }
    return(vals2)
}

"%ni%" <- Negate("%in%")
