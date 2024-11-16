#!/usr/bin/env Rscript
#
# This script calls sex based on a chrX variant call matrix
#
suppressMessages(require(optparse))
suppressMessages(require(data.table))

bam.sex <- function(bam.fn) {
    tmp <- system2("samtools", c("idxstat", bam.fn), stdout = TRUE)
    tmp <- fread(paste(tmp, collapse="\n"))
    tmp <- tmp[V1 %in% c("chrX", "chrY")]
    tmp[,cov:=V3/V2]
    x.y <- tmp[V1=="chrX"]$cov/tmp[V1=="chrY"]$cov
    write(sprintf("X-Y ratio: %.2f", x.y), stderr())
    sex <- ifelse(x.y>3, "XX", ifelse(x.y<1.5, "XY", NA_character_))
    return(sex)
}

snp.sex <- function(gvcf.fn) {
    sex <- NA_character_
    return(sex)
}

option_list <- list(
  make_option(c("-s", "--snp"), type='character', default = NULL,
              help="Input genotype (GVCF file)"),
  make_option(c("-b", "--bam"), type='character', default = NULL,
              help="Input alignment (BAM file)"),
  make_option(c("-o", "--output"), type='character', default = "sex.txt",
              help="Output txt file name")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$bam) && is.null(opt$snp)) {
    write("Error: Neither BAM or GVCF file provided.", stderr())
    quit("no", 1)
}
if (!is.null(opt$bam) && !is.null(opt$snp)) {
    write("Error: Both BAM and GVCF file provided.", stderr())
    quit("no", 1)
}
if (!is.null(opt$bam)) {
    call <- bam.sex(opt$bam)
} else {
    call <- snp.sex(opt$snp)
}

call.str <- sprintf("SEX=%s", call)
writeLines(call.str, opt$output)
