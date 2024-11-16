
#### Create 

#' Given the metavault, annotates. 
#'
#' @param METAVAULT The metavault
#' @param GOBJ The gobj
#' @param OPTS The opts
#' @return list of annotated data.table
#' @export
annotateMetaVault <- function(METAVAULT, GOBJ, OPTS){
  ## remove failed variants
  METAVAULT = vaultFilter(METAVAULT)
  ## extract and annotate vaults
  ids = METAVAULT$meta |> pull(id)
  vaults = lapply(ids, dxVaultSubsetById, dxv=METAVAULT)
  annot_vaults = lapply(vaults, annotateVault, GOBJ, OPTS)
  names(annot_vaults) <- ids
  
  ## format
  annot_metaVault = reformatAnnotVaultList(annot_vaults)
  return(annot_metaVault)
}




#' Given the annotated list, squishes list of lists to list of dfs
#'
#' @param ALIST The annotated list of vaults
#' @return list of annotated data.table
#' @export
reformatAnnotVaultList <- function(ALIST) {
  ## squish
  som = getAnnotListTable(ALIST, 'somatic')
  grm = getAnnotListTable(ALIST, 'germline')
  str = getAnnotListTable(ALIST, 'structural')
  gene.copy = getAnnotListTable(ALIST, 'gene.copy')
  gene.expression = getAnnotListTable(ALIST, 'gene.expression')
  fus = getAnnotListTable(ALIST, 'fusions')
  cin <- getAnnotListTable(ALIST, 'cin')
  arms <- getAnnotListTable(ALIST, 'arms')
  chromosomes <- getAnnotListTable(ALIST, 'chromosomes')
  meta <- getAnnotListTable(ALIST, 'meta')
  ## organize
  annotation <- list(
    somatic = som,
    germline = grm,
    structural = str,
    gene.copy = gene.copy,
    gene.expression = gene.expression,
    fusions = fus,
    cin = cin,
    arms = arms,
    chromosomes = chromosomes,
    meta = meta
  )
  return(annotation)
}

#' Gets tables from annotated list
#'
#' @param ALIST The annotated list
#' @param FIELD the table to get
#' @return list of annotated data.table
getAnnotListTable <- function(ALIST, FIELD) {
  nms = names(ALIST)
  int = lapply(nms, .getAnnotListTable, ALIST, FIELD) %>%
    rbindlist(.) %>%
    return()
}

#' Gets tables from annotated list for a id
#'
#' @param id THe id to grab
#' @param ALIST The annotated list
#' @param FIELD the table to get
#' @return annotated data.table
.getAnnotListTable <- function(ID, ALIST, FIELD) {
  tmp = ALIST[[ID]][[FIELD]]
  if(!is.null(tmp)){
    if(dim(tmp)[1]>0){
      return(tmp)
    }
  }
  return(NULL)
}




#### Misc

#' Given a vector of ids, subsets vault
#'
#' @param MVA An annotated metavault
#' @param IDS a vector of IDs
#' @return metavault
#' @export 
subsetAnnotatedMetaVault <- function(MVA, IDS) {
  
  MVA = lapply(MVA, function(X){
    X = X[id %in% IDS,]
  })
  return(MVA)
}
