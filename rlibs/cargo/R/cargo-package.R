#' cargo
#'
#' @name cargo
#' @docType package
#' @importFrom tpolib list.files.null get.opts pickExisting
#' @importFrom parallel detectCores
#' @importFrom GenomeInfoDb keepStandardChromosomes Seqinfo dropSeqlevels seqlevels keepSeqlevels seqlevels<- seqinfo<- seqnames seqnames<- seqinfo seqlevelsStyle<- genome
#' @importFrom GenomicRanges tileGenome reduce granges findOverlaps width pintersect mcols mcols<- start start<- end end<- split gaps strand promoters GRangesList GRanges  width<- width makeGRangesFromDataFrame elementMetadata elementMetadata<- disjoin
#' @importFrom S4Vectors queryHits subjectHits DataFrame %in% endoapply elementNROWS split rbind List
#' @importFrom rtracklayer import
#' @importFrom cnvex cinNLOH cinWgii cinNtAI cinLST
#' @importFrom carat fusionCSQ structuralCSQ
#' @importFrom codac svFormat
#' @importFrom IRanges findOverlapPairs to IRanges ranges<- subsetByOverlaps nearest
#' @importFrom limma weighted.median
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom tibble column_to_rownames
#' @importFrom purrr accumulate
#' @importFrom edgeR cpm rpkm
#' @importFrom magrittr set_colnames '%$%' %<>%
#' @importFrom dplyr '%>%' select filter mutate arrange rename bind_rows distinct relocate tbl collect inner_join pull rows_update group_by n_distinct left_join tibble transmute semi_join ungroup summarize tally
#' @importFrom data.table data.table as.data.table fread copy := setcolorder fifelse rbindlist SJ setcolorder setnames melt setkey key setkeyv foverlaps
#' @importFrom VariantAnnotation readVcf
#' @importFrom stringr str_sub str_match str_replace str_split str_length
#' @importFrom modeest mfv
#' @importFrom tidyr spread unite unnest
#' @importFrom grid gpar unit grid.text
#' @importFrom circlize colorRamp2
#' @importFrom trackViewer lolliplot
#' @importFrom maftools gisticChromPlot readGistic
#' @importFrom duckdb dbConnect dbDisconnect dbSendQuery duckdb_unregister duckdb_register_arrow duckdb_shutdown dbListTables
#' @importFrom arrow write_parquet
#' @importFrom dtplyr lazy_dt
#' @importFrom modeest mlv
#' @importFrom diptest dip.test dip
#'
NULL
