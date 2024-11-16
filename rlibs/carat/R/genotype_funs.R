#' Convert tpo genotyping csv to iupac string
#'
#' @param tpocsv path to csv from tpo genotyping
#' @return A chr of the IUPAC string
convert_to_iupac <- function(tpocsv) {
  gtds <- read.csv(tpocsv)
  #convert the hets to IUPAC i.e. A/C -> M
  gt <- as.character(gtds$final_call)
  gt <- gsub("/","",gt)
  gt <- gsub("UNKNOWN","ACGT",gt)
  dict <- names(IUPAC_CODE_MAP)
  names(dict) <- IUPAC_CODE_MAP
  gtt <- as.character(dict[gt])
  gtt <- paste(gtt, collapse="")
  #return IUPAC string
  gtt
}

#' Compare one genotype to many
#'
#' @param goi a (named) genotype string in IUPAC format
#' @return out A dataframe with the names and scores of all comparisons
geno_comp <- function(goi, pop, filter=0) {
  sm <- nucleotideSubstitutionMatrix(match=1,mismatch=0,baseOnly=FALSE,type="DNA")
  diag(sm) <- 1
  pal <- pairwiseAlignment(pop, goi, type="global", substitutionMatrix=sm)
  out <- data.frame(score=score(pal), gt = names(pop)) %>%
    mutate(score=100*score/max(score)) %>%
    filter(score>filter) %>%
    arrange(desc(score))
 # return data frame of scores
 out
}


#examples
#library(Biostrings)
#r <- function(x){paste(sample(names(IUPAC_CODE_MAP), 166, replace=T), collapse='')}
#demo <- sapply(seq(10000), FUN=r)
#demo <- DNAStringSet(demo)
#names(demo) <- seq(10000)
#
#system.time(foo <- geno_comp(demo[1],demo[1:5000]))
#user  system elapsed
#1.682   0.030   1.715
