suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))

nameMap <- function(bclout, bclsheet, format) {
    all.fq <- list.files(bclout, pattern="*.fastq.gz$")
    tmp1 <- as.data.table(do.call(rbind, str_split(all.fq, "_")))
    tmp1 <- tmp1[,.(Barcode=V1,
                    Lane=str_replace(V3, "L00", ""),
                    Read=str_replace(V4, "R", ""),
                    File=all.fq,
                    Size=file.size(file.path(bclout, all.fq)) / 1e6
                    )]
    tmp2 <- fread(bclsheet, colClasses="character", sep="\t", header=TRUE)
    setkey(tmp1, Barcode, Lane)
    setkey(tmp2, Barcode, Lane)
    tmp3 <- tmp2[tmp1]
    tmp3 <- tmp3[!is.na(Library),DestFile:=sprintf(format, Library, Flowcell, Lane, Read)]
    fin <- tmp3[,.(File, DestFile, Size)]
    return(fin)    
}

moveFiles <- function(bclout, bclsheet, format) {
    nm <- nameMap(bclout, bclsheet, format)
    nm.ok <- nm[!is.na(DestFile)]
    dir.create(file.path(bclout, "rename"), showWarnings=FALSE)
    fwrite(nm.ok, file.path(bclout, "file_rename.txt"))
    for (i in seq_len(nrow(nm.ok))) {
        system2("mv", c(file.path(bclout, nm.ok$File[i]), file.path(bclout, "rename", nm.ok$DestFile[i])))
    }
}

unmoveFiles <- function(bclout, bclsheet) {
    nm.ok <- fread(file.path(bclout, "file_rename.txt"))
    for (i in seq_len(nrow(nm.ok))) {
        system2("mv", c(file.path(bclout, "rename", nm.ok$DestFile[i]), file.path(bclout, nm.ok$File[i])))
    }
    file.remove(file.path(bclout, "file_rename.txt"))
}


option_list = list(
    optparse::make_option(c("-o", "--out"), type="character",
                          default = "",
                          help="bcl2fastq output folder"),
    optparse::make_option(c("-s", "--sheet"), type="character",
                          default = "",
                          help="library sheet"),
    optparse::make_option(c("-u", "--undo"), type="logical",
                          default = FALSE,
                          help="undo filemove"),
    optparse::make_option(c("-f", "--format"), type="character",
                          default = "%s-%s-%s-%s.fq.gz",
                          help="format string for output names, Lib, FC, Lane, Read")
    )
parser = optparse::OptionParser(
                       "Rscript fastq_rename.R [options]",
                       description=c("Rename FASTQ files.\n"),
                       epilogue=c(
                           "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
                           "Michigan Center for Translational Pathology (c) 2019\n"),
                       option_list=option_list
                   )

args <- optparse::parse_args(parser, positional_arguments=TRUE)

if(!file.exists(args$options$out)) {
    optparse::print_help(parser)
    write("bcl2fastq output folder not found.\n", stderr())
    quit("no", 1)
}

if(!file.exists(args$options$sheet)) {
    optparse::print_help(parser)
    write("library sheet not found.\n", stderr())
    quit("no", 1)
}

if(args$options$undo && !file.exists(file.path(args$options$out, "file_rename.txt"))) {
    optparse::print_help(parser)
    write("Nothing to undo.\n", stderr())
    quit("no", 1)
}

if(!args$options$undo && file.exists(file.path(args$options$out, "file_rename.txt"))) {
    optparse::print_help(parser)
    write("Files already renamed.\n", stderr())
    quit("no", 1)
}

if (args$options$undo) {
    unmoveFiles(args$options$out, args$options$sheet)
} else {
    moveFiles(args$options$out, args$options$sheet, args$options$format)
}
