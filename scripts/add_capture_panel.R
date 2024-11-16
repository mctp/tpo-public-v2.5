#!/usr/bin/Rscript

option_list = list(
        optparse::make_option(c("--targets"), type="character",
                              default=NULL,
                              help="Targets bed file"),
        optparse::make_option(c("--genome"), type="character",
                              default=NULL,
                              help="Genome fasta (for assigning contigs)"),
        optparse::make_option(c("--exclude"), type="character",
                              default=NULL,
                              help="optional bed file of regions to exclude from annotation (genotyping or backbone SNPs etc.)"),
        optparse::make_option(c("--annotation"), type="character",
                              default=NULL,
                              help="existing TPO annotated_regions bed file"),
        optparse::make_option(c("--padding"), type="integer",
                              default=50,
                              help="padding added to both sides of target bed file for annotation (default:50)"),
        optparse::make_option(c("--picard"), type="character",
                              default=NULL,
                              help="path to picard.jar")
    )
parser = optparse::OptionParser("Rscript add_capture_panel.R [options]",
  description=c(
    "Generate necessary region files to add a new capture panel to TPO. Returns a MISC.interval_list for cords-misc, a CNV.bed for cords-cnvex and a ANNO.bed for cords-anno (annotation bed + padded [targets - excluded]). All output files are sorted and reduced. \n
\nInstructions\n
- Run this script
- Run 'moke refs_pull'
- copy MISC.interval_list -> refs/<build>/capture/
- copy ANNO.bed -> refs/<build>/custom/
- Bump the ROOT_VER in the mokefile and run 'moke refs_push' and 'moke gcp_root_build'
- CNV.bed -> tpo/rlibs/cnvex/inst/extdata/capture/
- Run 'moke code_build' (version bumping if needed)\n"
    ),
  epilogue=c("Michigan Center for Translational Pathology (c) 2021\n"),
  option_list=option_list
  )

args <- optparse::parse_args(parser, positional_arguments=FALSE)
###################################################################################
###################################################################################
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))

cat("\nMaking bed files...\n")
# Make the CNV bed file (just sort and reduce the targets)
target_regions <- import.bed(args$targets)
target_regions <- reduce(sort(target_regions))
export.bed(target_regions,'CNV.bed')
# print some sanity checks
sprintf("Target Contigs:%s", 
  paste0(
    stringr::str_sort(unique(
      as.character(
        seqnames(target_regions)
      )
    ), numeric=T),
    collapse=','
    )
)
sprintf('Target Size %.02fMb',sum(width(target_regions))/1e6)

# Next make the interval list for Hsmetrics (targets, all regions)
## Run Picard to make the interval_list
cat("\nWriting interval_list using picard...\n")
cmd=sprintf(
  "java -jar %s BedToIntervalList I=./CNV.bed O=MISC.interval_list SD=%s &>picard.log", 
  args$picard,args$genome
)
system(cmd)

# Make the annotated region bed file (old annotation file + [targets - excluded]
cat("\nMaking annotation bed file...\n")
old_anno <- import.bed(args$annotation)
strand(old_anno) <- "*"
if(!is.null(args$exclude)){
  excluded <- import.bed(args$exclude)
  to_add <- setdiff(target_regions,excluded) + args$padding
}else{
  to_add <- target_regions + args$padding
}
new_anno <- sort(reduce(union(old_anno,to_add)))

#some sanity checks
stopifnot(all(overlapsAny(to_add,new_anno))) # all regions to add must be in output
sprintf('(Padded) Regions to add: %.02fMb',sum(width(to_add))/1e6)
sprintf('Old Annotation Size %.02fMb',sum(width(old_anno))/1e6)
sprintf('New Annotation Size %.02fMb',sum(width(new_anno))/1e6)

export.bed(new_anno, 'ANNO.bed')

cat("\nDone.\n")

