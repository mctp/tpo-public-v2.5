bamStat <- function(dir, ann) {
    bam.stat <- list(
        chim.pe.reads = .aln.reads(dir$bam.pe.fn),
        chim.se.reads = .aln.reads(dir$bam.se.fn),
        alig.reads = .aln.reads(dir$bam.fn)
    )
    return(bam.stat)
}
