#' carat
#'
#' @name carat
#' @docType package
#' @importFrom tpolib row.fisher.test
#' @importFrom parallel detectCores
#' @importFrom GenomeInfoDb keepStandardChromosomes Seqinfo dropSeqlevels seqlevels keepSeqlevels seqlevels<- seqinfo<- seqnames seqnames<- seqinfo seqlevelsStyle<- genome
#' @importFrom GenomicRanges tileGenome reduce granges findOverlaps width pintersect mcols mcols<- start start<- end end<- split gaps strand promoters GRangesList GRanges  width<- width makeGRangesFromDataFrame
#' @importFrom S4Vectors queryHits subjectHits DataFrame %in% endoapply elementNROWS split rbind List
#' @importFrom IRanges %over% CharacterList to IRanges start end ranges ranges<-
#' @importFrom data.table data.table setkey as.data.table fread setDT rbindlist dcast.data.table copy melt := setnames set SJ setcolorder
#' @importFrom stringr str_sub str_match str_replace str_split
#' @importFrom Biostrings getSeq IUPAC_CODE_MAP nucleotideSubstitutionMatrix pairwiseAlignment score
#' @importFrom VariantAnnotation readVcf writeVcf info geno header info<- geno<- header<- fixed fixed fixed<- alt ref
#' @importFrom StructuralVariantAnnotation breakpointRanges findBreakpointOverlaps
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom igraph clusters graph_from_edgelist
#' @importFrom Biostrings DNAStringSet DNAStringSetList BStringSet xscat subseq reverseComplement
#' @importFrom BiocGenerics unlist
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel stopImplicitCluster
NULL
