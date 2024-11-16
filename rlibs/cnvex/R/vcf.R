.preFilterVar <- function(var, opts) {
    var <- var[
        lengths(fixed(var)$ALT) == 1 &
        isSNV(var) &
        !is.na(qual(var))
    ]
    if (ncol(geno(var)$DP)>=1) {
        var <- var[!is.na(geno(var)$DP[,1])]
        var <- var[geno(var)$DP[,1] > opts$baf.min.dp]
    }
    if (ncol(geno(var)$DP)>=2) {
        var <- var[!is.na(geno(var)$DP[,2])]
        var <- var[geno(var)$DP[,2] > opts$baf.min.dp]
    }
    return(var)
}

.maskVar <- function(var, gobj) {
    ## masked regions 1000G
    loose           <- gobj$refs$unmask.loose
    strict          <- gobj$refs$unmask.strict
    giab            <- gobj$refs$giab.isdifficult
    tumor.only.mask <- gobj$refs$tumor.only.mask
    tmp <- DataFrame(Number=c("1","1","1","1"),
                     Type=c("Integer", "Integer","Integer", "Integer"),
                     Description=c("loose mask", "strict mask","giab mask", "tumor only mask"),
                     row.names=c("mask.loose", "mask.strict","mask.giab","tumor.only.mask"))
    info(header(var)) <- rbind(info(header(var)), tmp)
    info(var)$mask.loose <- as.integer(!(var %over% loose)) # 0 means it is good (not masked)
    info(var)$mask.strict <- as.integer(!(var %over% strict)) # 0 means it is good (not masked)
    info(var)$mask.giab <- as.integer((var %over% giab)) ## giab is difficult regions, 0 means it is good (not masked)
    info(var)$tumor.only.mask <- as.integer(!(var %over% tumor.only.mask)) # 0 means it is good (not masked)
    return(var)
}

.importVcf <- function(vcf, param, gobj, opts) {
    var <- .robust.import.vcf(vcf, param, gobj$seqi)
    var <- .preFilterVar(var, opts)
    var <- .maskVar(var, gobj)
    return(var)
}
