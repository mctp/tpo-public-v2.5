# A script to read an illumina sample sheet and write the reverse complement of index2 for use with the illumina reverse complement workflow (i.e. chemistry 1.5)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Biostrings))

option_list = list(
    optparse::make_option(c("-i", "--input"), type="character",
                          default=NULL,
                          help="input Illumina sample sheet"),
    optparse::make_option(c("-o", "--output"), type="character",
                          default=NULL,
                          help="output Illumina sample sheet")
)
parser = optparse::OptionParser(
  "reverse_complement() [options]",
  description=c("convert a 1.0 sample sheet to 1.5"),
  epilogue=c(
      "Michigan Center for Translational Pathology (c) 2022\n"),
  option_list=option_list
  )


args <- optparse::parse_args(parser, positional_arguments=FALSE)

#args$input <- '/data/deployment/stage/tpo/pipe/bcl/novaseq_s4_di_merge.csv'
#args$output <- '/data/deployment/stage/mioncoseq2/foo'

str <- readLines(args$input)

#break up the file by the header
header <- grep('^Lane',str)

ds <- data.frame(raw=str[seq(header+1,length(str))]) %>%
        tidyr::separate(raw,into=strsplit(str[header],',')[[1]], sep=',') %>%
        mutate(
          index2=as.character(
            reverseComplement(
              DNAStringSet(index2)
            )
          )
        )

writeLines(str[1:header-1],con=args$output)
suppressWarnings( # a warning about writing headers to an existing file
  write.table(ds, file=args$output, sep=',', append=T,row.names=F, quote=F)
)

stopifnot(length(readLines(args$input))==length(readLines(args$output)))
