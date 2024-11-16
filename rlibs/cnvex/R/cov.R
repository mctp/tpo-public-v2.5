.normCoverage <- function(cov, tile, opts) {
    ## normalized coverage
    cov.n <- cov*width(tile)
    ## on-target
    cov.n.tar.total <- sum(cov.n[tile$target], na.rm = TRUE)
    cov.n.tar <- cov.n[tile$target]
    cov.n.tar <- cov.n.tar/(cov.n.tar.total/1e6)
    cov.n.tar <- 1000*cov.n.tar/width(tile[tile$target])
    cov.n[tile$target] <- cov.n.tar
    ## off-target
    cov.n.off.total <- sum(cov.n[!tile$target], na.rm = TRUE)
    cov.n.off <- cov.n[!tile$target]
    cov.n.off <- cov.n.off/(cov.n.off.total/1e6)
    cov.n.off <- 1000*cov.n.off/width(tile[!tile$target])
    cov.n[!tile$target] <- cov.n.off
    return(cov.n)
}

## .normCoverageSex <- function(cov, tile, sex, gobj, opts) {
##     ## normalized coverage
##     cov.n <- cov * width(tile)
##     ## female/male ratio
##     chrx <- as.logical(seqnames(tile) %in% gobj$sexchr["X"])
##     chry <- as.logical(seqnames(tile) %in% gobj$sexchr["Y"])
##     auto <- as.logical(seqnames(tile) %ni% gobj$sexchr)
##     cov.male <- sum(cov[auto])+sum(cov[chrx])+sum(cov[chry])
##     cov.female <- sum(cov[auto])+(2*sum(cov[chrx]))
##     cov.adjust <- ifelse(sex == "male", (cov.female/cov.male), 1)
##     ## on-target
##     cov.n.tar.total <- sum(cov.n[tile$target], na.rm = TRUE) * cov.adjust
##     cov.n.tar <- cov.n[tile$target]
##     cov.n.tar <- cov.n.tar/(cov.n.tar.total/1e6)
##     cov.n.tar <- 1000*cov.n.tar/width(tile[ tile$target])
##     cov.n[ tile$target] <- cov.n.tar
##     ## off-target
##     cov.n.off.total <- sum(cov.n[!tile$target], na.rm = TRUE) * cov.adjust
##     cov.n.off <- cov.n[!tile$target]
##     cov.n.off <- cov.n.off/(cov.n.off.total/1e6)
##     cov.n.off <- 1000*cov.n.off/width(tile[!tile$target])
##     cov.n[!tile$target] <- cov.n.off
##     return(cov.n)
## }
