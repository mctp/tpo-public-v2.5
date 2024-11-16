# Takes as input a fasta file and constructs a set of genomic ranges
# for all poly-A sequences of length 6 or greater

# Code adapted from https://www.biostars.org/p/280234/

suppressPackageStartupMessages(require(Biostrings))
suppressPackageStartupMessages(require(GenomicRanges))
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(optparse))

#Args
option_list = list(
    optparse::make_option(c("-f", "--fasta"), type = "character",
                          default = NULL,
                          help = "Reference FASTA file"),
    optparse::make_option(c("-o", "--out"), type = "character",
                          default = NULL,
                          help = "Output filename for GenomicRanges RDS file")
)

parser = optparse::OptionParser(
    "build_msi_ref [options]",
    description=c("Construct a GenomicRanges object for all poly-A sequences of length 5 or greater\n"),
    epilogue=c(
        "Michigan Center for Translational Pathology (c) 2020\n"),
    option_list=option_list
)

args <- optparse::parse_args(parser, positional_arguments = FALSE)

## Avoid floating point numbers for larger genomic coordinates:
options(scipen=999)

## Function searches a given sequence (full.seq) for pattern matches (pattern.seq)
## Will catch all homopolymers of equal or greater length
seq.check <- function(pattern.seq, full.seq){
    ## Get all (also redundant) coordinates of the exact pattern match:
    tmp.match <- matchPattern(pattern.seq, full.seq[[1]])
    
    ## If no matches are found, output empty GRanges, else output GRanges with unique coordinates:
    if (length(tmp.match@ranges) == 0){
        return(
            GRanges()
        )
    } else{ 
        return( 
            suppressWarnings(GenomicRanges::reduce(GRanges(seqnames = full.seq@ranges@NAMES, ranges = ranges(tmp.match))))
        )
    }
}

####################################################################################################################
#We only really care about the primary autosome assemblies
ref_fasta <- import(con = FastaFile(args$fasta),
                    format = "fasta", type = "DNA")[1:22]

#Extract contig names in case FASTA names field has extra info
names(ref_fasta) <- sapply(names(ref_fasta), function(x) strsplit(x, split = " ")[[1]][1])

#Get all homopolymer sequences of length 5 or more
polya_ranges <- suppressWarnings(do.call("c",
                                         lapply(names(ref_fasta), function(x) seq.check(pattern.seq = "AAAAA",
                                                                                        full.seq = ref_fasta[x, ]))))

polyt_ranges <- suppressWarnings(do.call("c",
                                         lapply(names(ref_fasta), function(x) seq.check(pattern.seq = "TTTTT",
                                                                                        full.seq = ref_fasta[x, ]))))

polyc_ranges <- suppressWarnings(do.call("c",
                                         lapply(names(ref_fasta), function(x) seq.check(pattern.seq = "CCCCC",
                                                                                        full.seq = ref_fasta[x, ]))))

polyg_ranges <- suppressWarnings(do.call("c",
                                         lapply(names(ref_fasta), function(x) seq.check(pattern.seq = "GGGGG",
                                                                                        full.seq = ref_fasta[x, ]))))

homopolymer_ranges <- c(polya_ranges, polyt_ranges, polyc_ranges, polyg_ranges)

saveRDS(homopolymer_ranges, file = args$out)
