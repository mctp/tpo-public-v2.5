suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(optparse))

transcriptsToGene <- function(gtf.fn) {
    gtf <- import(gtf.fn)
    gtf.tx <- gtf[gtf$type=="transcript"]
    dt.tx <- as.data.table(gtf.tx)
    t2gn <- dt.tx[,.(
        target_id=transcript_id,
        transcript_id=paste(transcript_id, transcript_version, sep="."),
        gene_id=paste(gene_id, gene_version, sep="."),
        gene_name
    )]
    return(t2gn)
}

expressionSummary <- function(kallisto.fn, gtf.fn) {
    t2gn <- transcriptsToGene(gtf.fn)
    setkey(t2gn, target_id)
    kallisto.out <- fread(kallisto.fn)
    setkey(kallisto.out, target_id)
    kallisto.out <- kallisto.out[t2gn]
    gene.out <- kallisto.out[,.(est_counts=sum(est_counts), tpm=sum(tpm)), .(gene_id, gene_name)]
    return(gene.out)
}

option_list = list(
    optparse::make_option(c("-g", "--gtf"), type="character",
                          default = NA_character_,
                          help="Location of annotation file [.gtf]"),
    optparse::make_option(c("-i", "--inp"), type="character",
                          default = NA_character_,
                          help="Input kallisto file [abundance.tsv]"),
    optparse::make_option(c("-o", "--out"), type="character",
                          default = NA_character_,
                          help="Output file [.csv]")
    )
parser = optparse::OptionParser(
                       "Rscript gene_summary.R -g [gtf_file] -i [kallisto] -o [csv_file]",
                       description=c("Summarize Expression at the gene level.\n"),
                       epilogue=c(
                           "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
                           "Michigan Center for Translational Pathology (c) 2019\n"),
                       option_list=option_list
                   )

args <- optparse::parse_args(parser, positional_arguments=TRUE)

if(!file.exists(args$options$gtf)) {
    optparse::print_help(parser)
    write("Annotation file not found.\n", stderr())
    quit("no", 1)
}

if(!file.exists(args$options$inp)) {
    optparse::print_help(parser)
    write("Input file not found.\n", stderr())
    quit("no", 1)
}

if(is.na(args$options$out)) {
    optparse::print_help(parser)
    write("Output file not provided.\n", stderr())
    quit("no", 1)
}

if(file.exists(args$options$out)) {
    optparse::print_help(parser)
    write("Output file exists.\n", stderr())
    quit("no", 1)
}

out <- expressionSummary(args$options$inp, args$options$gtf)
fwrite(out, args$options$out)
