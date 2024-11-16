
#' Returns a matrix usable for oncoprints
#'
#' @param FEATURES A vector of gene names
#' @param SAMPLES A vector of sample IDs
#' @param HITS A 3 col dataframe of hits
#' @return A df usable to make oncoprint
#'
#' @export
createOncoMatrix <- function(FEATURES, ANNOTATED, FIELD, SUBFIELD){
  
  ## prepare
  hits = .prepareOncoHits(ANNOTATED, FIELD, SUBFIELD, FEATURES)
  non_hits = .prepareOncoNonHits(ANNOTATED, FIELD, FEATURES)
  missing = .prepareOncoMissing(ANNOTATED, FIELD, FEATURES)
  ## add non-hits/missing
  hits = rbind(hits,non_hits)
  hits = rbind(hits,missing,fill=T)
  hits = hits[!duplicated(paste(hits$id, hits$gene_name, sep = '..')), ]
  hits = hits[!is.na(hits$id),]
  ## reformat
  hits = tidyr::spread(hits, key = 'id', value = 'annotation')
  hits = tibble::column_to_rownames(hits, 'gene_name')
  hits = hits[order(rownames(hits)) , order(colnames(hits))]
  hits = as.matrix(hits)
  return(hits)
}

#' Returns a df usable for making oncoMatrix
#'
#' @param ANNOTATED An annotated metavault
#' @param FIELD Part of annotated metavault
#' @param SUBFIELD Type of FIELD (if necessary)
#' @param FEATURES the features of interest
#' @return A df of hits of specified type
#'
.prepareOncoHits <- function(ANNOTATED, FIELD, SUBFIELD, FEATURES){
  hits = ANNOTATED[[FIELD]][ (annotation %in% SUBFIELD) & gene_name %in% FEATURES, c('id', 'gene_name', 'annotation')]
  if(length(hits$id)>0){
    hits[['annotation']] = 1
    hits = distinct(hits)
    return(hits)
  }
  return(NULL)
  
}

#' Returns a df usable for making oncoMatrix
#'
#' @param ANNOTATED An annotated metavault
#' @param FIELD Part of annotated metavault
#' @param FEATURES the features of interest
#' @return A df of hits of specified type
#'
.prepareOncoNonHits <- function(ANNOTATED, FIELD, FEATURES){

  if(FIELD %in% c('somatic','structural','gene.copy', 'germline')){
    ids = ANNOTATED$meta$id[!is.na(ANNOTATED$meta$tlib)]
    non_hits = expand.grid(id = ids, gene_name = FEATURES)
    non_hits$annotation = 0
  }
  if(FIELD %in% c('fusions')){
    ids = ANNOTATED$meta$id[!is.na(ANNOTATED$meta$rlib)]
    non_hits = expand.grid(id = ids, gene_name = FEATURES)
    non_hits$annotation = 0
  }
  return(non_hits)
}

#' Returns a df usable for making oncoMatrix
#'
#' @param ANNOTATED An annotated metavault
#' @param FIELD Part of annotated metavault
#' @param FEATURES the features of interest
#' @return A df of hits of specified type
#'
.prepareOncoMissing <- function(ANNOTATED, FIELD, FEATURES){
  missing=NULL
  if(FIELD %in% c('somatic','structural','gene.copy','germline')){
    if( anyNA(ANNOTATED$meta$tlib) ){
      ids = ANNOTATED$meta$id[is.na(ANNOTATED$meta$tlib)]
      missing = expand.grid(id = ids, gene_name = FEATURES)
      missing$annotation = 0
    }
  }
  if(FIELD %in% c('fusions')){
    if( anyNA(ANNOTATED$meta$rlib) ){
      ids = ANNOTATED$meta$id[is.na(ANNOTATED$meta$rlib)]
      missing = expand.grid(id = ids, gene_name = FEATURES)
      missing$annotation = 0 
    }
  }
  return(missing)
}

#' Returns a df usable for making oncoMatrix
#'
#' @param ANNOTATED An annotated metavault
#' @param FIELD Part of annotated metavault
#' @return A matrix of CIN values
#'
#' @export
createOncoCin <- function(ANNOTATED, FIELD){
  ## get cin values
  cin = ANNOTATED$cin
  ## add missing
  tmp = data.table('id' = ANNOTATED$meta$id[!ANNOTATED$meta$id %in% cin$id])
  cin = rbind(cin, tmp, fill=T)
  ## reformat
  cin = cin[order(cin$id), ]
  cin = cin[,c('id',FIELD), with=F]
  cin = tibble::column_to_rownames(cin, 'id')
  cin = t(cin)
  return(cin)
}

#' Returns a df usable for making oncoMatrix
#'
#' @param ANNOTATED An annotated metavault
#' @param FIELD Part of annotated metavault
#' @return A matrix of chr arm gain/loss
#'
#' @export
createOncoArms <- function(ANNOTATED, ARMS){
  ## get cin values
  arms = ANNOTATED$arms
  ## filter
  arms = arms[arm %in% ARMS,]
  ## add missing
  tmp = data.table('id' = ANNOTATED$meta$id[!ANNOTATED$meta$id %in% arms$id])
  arms = rbind(arms, tmp, fill=T)
  ## reformat
  arms = arms[,c('id','arm','annotation'), with=F]
  arms = tidyr::spread(arms, key=arm,value=annotation,fill=NA)
  arms = arms[order(arms$id)]
  arms = arms[,'<NA>':=NULL]
  arms = tibble::column_to_rownames(arms, 'id')
  arms = t(arms)
  return(arms)
}




#### Misc

#' Returns a df usable for making oncoMatrix
#'
#' @param MATRIX A prepared OncoMatrix
#' @param ORDER Vector of IDs in correct order
#' @return A matrix of chr arm gain/loss
#'
reorderOncoMatrix <- function(MATRIX, ORDER){
  ncol=ncol(MATRIX)
  nms = rownames(MATRIX)
  matrix = matrix(MATRIX[,ORDER], ncol = ncol, dimnames = list(nms,ORDER))
  return(matrix)
}

#' Simplifies multiple oncomatrix features
#'
#' @param OM A prepared OncoMatrix
#' @param NF Character value for new feature
#' @param OFS vector of old feature names from OM
#' @param VALUE (optional) value to set feature as
#' @return A simplified matrix with old features removed and new added
#'
#' @export
simplifyOncoMatrixFeature <- function(OM, NF, OFS, VALUE=NULL){
  ## get
  old = OM[rownames(OM)%in% OFS, ]
  om = OM[!rownames(OM)%in% OFS, ]
  new = colSums(old)
  ## set PRN
  if(!is.null(VALUE)){
    new[] = VALUE
  }
  ## combine
  nom = rbind(new, om)
  rownames(nom)[1] = NF
  return(nom)
}