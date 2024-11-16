#' @export
runCODAC <- function(dir, ann, par) {
    beg.time <- Sys.time()
    bam <- bamStat(dir, ann)
    spl <- readSplices(dir, ann)
    jnc <- readJunctions(dir, ann, par)
    bpt <- collapseBreakpoints(jnc, ann, par)
    bun <- makeBundle(bpt, jnc, "full")
    end.time <- Sys.time()
    log <- list(version=packageVersion("codac"), beg.time=beg.time, end.time=end.time, gc=gc())
    run <- list(dir=dir, spl=spl, bam=bam, bun=bun, log=log)
    return(run)
}
