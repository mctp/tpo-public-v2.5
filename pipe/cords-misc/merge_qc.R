# Combine the various QC metrics to one file and save it
# USAGE: Rscript merge_qc.R $BID $ID ${ALINID[@]}
# i.e. <name of the postalign> <name of the run> <array of input ALN ids>
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))

# parse args
args <- commandArgs(trailingOnly=TRUE)
bid <- args[1]
id <- args[2]

# get file paths
dedup <- sprintf("/input/%s/%s-dedup_metrics.txt",bid,bid)
hsmetrics <- sprintf("/output/%s-hsmetrics.txt",id)
aligstats <- sprintf("/output/%s-aligstat.txt",id)
geno <- sprintf("/output/%s-genotype.csv",id)
virus <- sprintf("/output/%s-virus.txt",id)

parse_geno <- function(gt) {
  gtds <- read.csv(gt)
  # convert the hets to IUPAC i.e. A/C -> M
  gt <- as.character(gtds$final_call)
  gt <- gsub("/","",gt)
  gt <- gsub("UNKNOWN","ACGT",gt)
  dict <- names(Biostrings::IUPAC_CODE_MAP)
  names(dict) <- Biostrings::IUPAC_CODE_MAP
  gtt <- as.character(dict[gt])
  gtt <- paste(gtt, collapse="")
  return(gtt)
}

if (all(c(
    file.exists(dedup),
    file.exists(aligstats),
    file.exists(geno),
    file.exists(virus)
    ))
) {

  # read the data
  #dedup metrics
  dedup <- fread(dedup, skip=1, nrows=1)
  dedup$ON <- bid
  dedup$FROM <- "dedup_metrics"

  out <- data.frame()
  #alignment statistics
  aligstats <- fread(aligstats, skip=1) %>%
      filter(CATEGORY=="PAIR")
  aligstats$ON <- id
  aligstats$FROM <- "BWA"

  #genotypes
  geno <- parse_geno(geno)
  geno <- data.frame(genotype=geno, ON=id, FROM="geno")

  #Virus Calls (may be empty)
  v <- try(read.delim(virus, skip=4, header=FALSE, stringsAsFactors=FALSE))
  if (class(v)=="try-error") {
    virus <- data.frame(virus="", ON=id, FROM="virus")
  } else {
    v %>%
      mutate(virus=gsub("SP:","",V6)) %>%
      mutate(call=sprintf("%s (%.0f reads, %s)", virus,V8,V9)) %>%
      select(call) %>% .$call %>%
      paste0(collapse=";") -> v
    virus <- data.frame(virus=v, ON=id, FROM="virus")
  }


  # transpose and merge
  qc <- rbind(
    melt(dedup,id.vars=c("ON","FROM")),
    melt(aligstats,id.vars=c("ON","FROM")),
    melt(geno,id.vars=c("ON","FROM")),
    melt(virus,id.vars=c("ON","FROM"))
  )

  # hsmetrics
  if(file.exists(hsmetrics)) {
    hsmetrics <- fread(hsmetrics, skip=1, nrows=1)
    hsmetrics$ON <- bid
    hsmetrics$FROM<- "hsmetrics"
    qc.hs <- melt(hsmetrics,id.vars=c("ON","FROM"))
    qc <- rbind(qc, qc.hs)
  }

  # save the output
  write.table(qc,file=sprintf("/output/%s-QC-summary.txt",id),sep="\t", row.names=FALSE, quote=FALSE)
  saveRDS(qc,sprintf("/output/%s-QC-summary.rds", id))

} else {
  warning("input files missing, QC not merged")
}
