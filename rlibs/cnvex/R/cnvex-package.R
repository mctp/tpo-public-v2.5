#' cnvex
#'
#' @name cnvex
#' @docType package
#' @importFrom tpolib absMedDiff get.opts plog_sum_exp max_na_rm sigmoid minSum
#' @importFrom rtracklayer import export
#' @importFrom matrixStats rowMedians rowMads rowMins rowMaxs colMedians
#' @importFrom Rsamtools BamFile
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom GenomeInfoDb keepStandardChromosomes Seqinfo dropSeqlevels seqlevels keepSeqlevels seqlevels<- seqinfo<- seqnames seqnames<- seqinfo mapSeqlevels renameSeqlevels seqlengths seqlevelsStyle<- seqlevelsStyle sortSeqlevels genome
#' @importFrom GenomicRanges tileGenome reduce granges findOverlaps width pintersect mcols mcols<- start end start<- end<- split gaps strand promoters GRangesList GRanges nearest
#' @importFrom S4Vectors queryHits subjectHits DataFrame %in% endoapply elementNROWS
#' @importFrom IRanges %over% IRanges disjoin
#' @importFrom data.table data.table setkey as.data.table fread fwrite setDT rbindlist dcast.data.table copy melt := fcase setnames
#' @importFrom stringr str_sub str_match str_replace str_split str_detect str_remove
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom VariantAnnotation readVcf ScanVcfParam info geno header info<- header<- info<- geno<- vcfWhich<- vcfWhich qual qual<- fixed isSNV vcfFields scanVcfHeader
#' @importFrom limma loessFit weighted.median
#' @importFrom DNAcopy CNA smooth.CNA
#' @importFrom DelayedArray rowRanges
#' @importFrom jointseg jointSeg estimateSd
#' @importFrom raster raster focal Which
#' @importFrom robfilter hybrid.filter
#' @importFrom fastICA fastICA
#' @importFrom MASS ginv
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom parallel detectCores parLapply
#' @importFrom logger log_debug log_info log_warn log_error
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom karyoploteR getDefaultPlotParams getCytobandColors plotKaryotype kpAddCytobands getChromosomeNamesBoundingBox kpAddChromosomeNames kpRect kpAxis kpPoints kpAddMainTitle getMainTitleBoundingBox autotrack kpPlotRegions kpAddLabels kpArea
#' @importFrom regioneR toGRanges
#' @importFrom cowplot plot_grid
#' @importFrom ggpubr theme_pubr ggparagraph
#' @importFrom egg ggarrange
#' @import ggplot2
NULL
