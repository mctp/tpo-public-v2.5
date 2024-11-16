.findAnnotationFiles <- function(genome) {
    if (genome=="hg38") {
        ann.files <- list(
            goi.fn = system.file("extdata", "hg38", "genes-of-interest.txt", package="codac"),
            loi.fn = system.file("extdata", "hg38", "loci-of-interest.txt", package="codac"),
            gtb.fn = system.file("extdata", "hg38", "genes-to-blacklist.txt", package="codac"),
            art.fn = system.file("extdata", "hg38", "artifacts-long.csv.gz", package="codac"),
            ars.fn = system.file("extdata", "hg38", "artifacts-short.csv.gz", package="codac"),
            rep.fn = system.file("extdata", "hg38", "rm.bed.gz", package="codac"),
            ret.fn = system.file("extdata", "hg38", "retro.bed.gz", package="codac"),
            cyt.fn = system.file("extdata", "hg38", "cytoband.txt", package="codac"),
            seg.fn = system.file("extdata", "hg38", "segdup.gz", package="codac"),
            alt.fn = system.file("extdata", "hg38", "scaffold-placement.txt", package="codac"),
            sv.fn  = system.file("extdata", "hg38", "1000G-sv-v2.vcf.gz", package="codac"),
            low.fn = system.file("extdata", "hg38", "dust50.bed.gz", package="codac"),
            igx.fn = system.file("extdata", "hg38", "igx.bed", package="codac"),
            sno.fn = system.file("extdata", "hg38", "sno.gff.gz", package="codac")
        )
    } else {
        write("Genome not supported.\n", stderr())
        quit("no", 1)
    }
    return(ann.files)
}

#' @export
makeDirectory <- function(dir.fn) {
    jnc.pe.fn <- list.files(dir.fn, pattern="*chim-pe.jnc.gz$", full.names=TRUE)
    jnc.se.fn <- list.files(dir.fn, pattern="*chim-se.jnc.gz$", full.names=TRUE)
    ## prefer BAM over CRAM
    bam.pe.fn <- head(list.files(dir.fn, pattern="*chim-pe.cram$|*chim-pe.bam$", full.names=TRUE), 1)
    bam.se.fn <- head(list.files(dir.fn, pattern="*chim-se.cram$|*chim-se.bam$", full.names=TRUE), 1)
    bam.fn <- head(list.files(dir.fn, pattern="*alig.cram$|*alig.bam$", full.names=TRUE), 1)
    sj.fn <- list.files(dir.fn, pattern="*-sj.tab.gz$", full.names=TRUE)
    if (
        length(jnc.se.fn)!=1 ||
        length(jnc.pe.fn)!=1 ||
        length(bam.pe.fn)!=1 ||
        length(bam.se.fn)!=1 ||
        length(sj.fn)!=1
    ) {
        write("Input directory does not conform.\n", stderr())
        quit("no", 1)
    }
    dir <- list(
      jnc.pe.fn = jnc.pe.fn,
      jnc.se.fn = jnc.se.fn,
      bam.pe.fn = bam.pe.fn,
      bam.se.fn = bam.se.fn,
      bam.fn = bam.fn,
      sj.fn = sj.fn,
      sid = basename(dir.fn)
    )
    return(dir)
}

#' @export
makeParams <- function(config.file) {
    ## presets
    if (!file.exists(config.file)) {
        config.file <- system.file("extdata", "conf", paste0(config.file, ".conf"), package="codac")
    }
    if (!file.exists(config.file)) {
        write(sprintf("Cannot find config file: %s.\n", config.file), stderr())
        quit("no", 1)
    }
    tmp <- fread(config.file)
    pres <- as.list(tmp$value)
    names(pres) <- tmp$param
    presn <- suppressWarnings(as.numeric(pres))
    presl <- suppressWarnings(as.logical(pres))
    pres[!is.na(presn)] <- presn[!is.na(presn)]
    pres[!is.na(presl)] <- presl[!is.na(presl)]
    return(pres)
}
