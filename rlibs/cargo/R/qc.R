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

make_geno_comp_table <- function(a, b) {
    a <- Biostrings::DNAString(a)
    b <- Biostrings::DNAString(b)
    sm <- Biostrings::nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=FALSE, type="DNA")
    diag(sm) <- 1

    #fix GOH as mismatch
    hets <- unlist(strsplit("RYSWKMBDHV",""))
    homs <- c("A","T","G","C")
    sm[homs,hets] <- 0

    al <- Biostrings::pairwiseAlignment(a, b, type="global", substitutionMatrix=sm)
    #compareStrings(al)
    df <- data.table(
        a=unlist(strsplit(as.character(a),"")),
        b=unlist(strsplit(as.character(b),"")),
        stringsAsFactors=F
        )
    df$score <- apply(df, MARGIN=1, function(x){sm[x[1], x[2]]})
    df$score <- df$score
    df$pos <- seq_len(nrow(df))
    #
    df
}

geno_comp_multi <- function(gt,ind=1) {
    dna <- DNAStringSet(gt$genotype)
    sm <- nucleotideSubstitutionMatrix(match=1, mismatch=0, baseOnly=FALSE, type="DNA")
    diag(sm) <- 1

    #fix GOH as mismatch
    hets <- unlist(strsplit("RYSWKMBDHV",""))
    homs <- c("A","T","G","C")
    sm[homs,hets] <- 0
    out <-  pairwiseAlignment(dna, dna[[ind]], type="global", substitutionMatrix=sm)
    sc <- score(out)
    gt$score <- sc/sc[ind]
    out <- as.data.table(arrange(gt,desc(score)))

    return(out)
}