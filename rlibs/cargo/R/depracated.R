#### Keys

LIKELY_BENIGN <- c("Benign", "Benign/Likely_benign", "Likely_benign")
UNLIKELY_BENIGN <- c("Pathogenic/Likely_pathogenic", "Pathogenic", "Likely_pathogenic", "drug_response")

#### Variant Triage/filters

#' Adds pass/fail value to vault tables based on annotation, evidence, consequence,
#' and artifact filters
#'
#' @param vault vault/metavault
#' @param anno anno object (list); contains gene info
#' @param fset filter set (name); optional name of the custom filters
#' @param annotation boolean; whether filter on annotation
#' @param evidence boolean; whether filter on evidence
#' @param consequence boolean; whether filter on consequence
#' @param artifact boolean; whether filter on artifact
#' @param pid boolean; whether to filter on PID
#' @param GOI optional vector of genes of interest
#'
#' @return vault w/ pass/fail info for features
#'
triageVault <- function(vault, anno=NULL, fset=NULL, annotation=TRUE, evidence=TRUE, consequence=TRUE, artifact=TRUE, pid=TRUE, GOI=NULL) {
  ## import default filter functions
  source(system.file("extdata/filters/default.R", package="cargo"))
  ## if provided import filtering functions which will override defaults
  if (!is.null(fset)) {
    if (file.exists(fset)) {
      source(fset)
    } else {
      source(system.file(sprintf("extdata/filters/%s.R", fset), package="cargo"))
    }
  }
  if (!is.null(anno)) {
    vault$anno <- anno
  }
  if (annotation) {
    vault <- .triageAnnotationVault(vault)
  }
  if (artifact) {
    vault <- .triageArtifactVault(vault)
  }
  if (consequence) {
    vault <- .triageConsequenceVault(vault)
  }
  if (evidence) {
    vault <- .triageEvidenceVault(vault, GOI)
  }
  if (pid) {
    vault <- .triagePhaseVault(vault)
  }
  return(vault)
}

#' Filter variants based on "pass" filter
#'
#' @param V - triaged vault
#' @return vault with failed variants (somatic, germline, structural, fusion) removed
filterVault <- function(V) {
  V <- .removeFailedVariantTable(V, "somatic")
  V <- .removeFailedVariantTable(V, "germline")
  V <- .removeFailedVariantTable(V, "structural")
  V <- .removeFailedVariantTable(V, "fusion")
  return(V)
}

#' Changes pass/fail value of variants based on annotation filter
#'
#' @param vault - vault
#' @param anno anno object (list); contains gene info
#'
#' @return filtered vault
.triageAnnotationVault <- function(vault) {
  gids <- vault$anno$genes$gene_id
  ## somatic
  tbl <- getVaultTable(vault, 'somatic')
  if (!is.null(tbl)) {
    if( !'pass' %in% colnames(tbl) ){ tbl=cbind(data.table(pass='pass'), tbl) }
    new_tbl <- tbl %>%
      dplyr::mutate(pass = ifelse(pass=='pass' & !GENE %in% gids, 'F-annotation', pass)) %>%
      dplyr::relocate(pass)
    vault <- setVaultTable(vault, 'somatic', new_tbl)
  }
  ## germline
  tbl <- getVaultTable(vault, 'germline')
  if (!is.null(tbl)) {
    if( !'pass' %in% colnames(tbl) ){ tbl=cbind(data.table(pass='pass'), tbl) }
    new_tbl <- tbl %>%
      dplyr::mutate(pass = ifelse(pass=='pass' & !GENE %in% gids, 'F-annotation', pass)) %>%
      dplyr::relocate(pass)
    vault <- setVaultTable(vault, 'germline', new_tbl)
  }
  ## structural
  tbl <- getVaultTable(vault, 'structural')
  if (!is.null(tbl)) {
    if( !'pass' %in% colnames(tbl) ){ tbl=cbind(data.table(pass='pass'), tbl) }
    new_tbl <- tbl %>%
      dplyr::mutate(pass = ifelse(pass=='pass' & !(GENE1 %in% gids | GENE2 %in% gids), 'F-annotation', pass)) %>%
      dplyr::relocate(pass)
    vault <- setVaultTable(vault, 'structural', new_tbl)
  }
  ## fusion
  tbl <- getVaultTable(vault, 'fusion')
  if (!is.null(tbl)) {
    if( !'pass' %in% colnames(tbl) ){ tbl=cbind(data.table(pass='pass'), tbl) }
    in_gids <- function(x) any(x %in% gids)
    new_tbl <- tbl %>%
      dplyr::mutate(pass = ifelse(pass=='pass' &
                                    !(sapply(str_split(gene_ids.5.1, ":"), in_gids) |
                                      sapply(str_split(gene_ids.3.1, ":"), in_gids) |
                                      sapply(str_split(gene_ids.5.2, ":"), in_gids) |
                                      sapply(str_split(gene_ids.3.2, ":"), in_gids)), 'F-annotation', pass)) %>%
      dplyr::relocate(pass)
    vault <- setVaultTable(vault, 'fusion', new_tbl)
  }
  return(vault)
}

#' Changes pass/fail value of variants based on evidence filter
#'
#' @param vault - vault
#' @param GOI - optional vector of genes of interst
#'
#' @return filtered vault
.triageEvidenceVault <- function(vault, GOI=NULL) {
  ## somatic
  if (!is.null(vault$tables$somatic)) {
    if( !'pass' %in% colnames(vault$tables$somatic) ){ vault[['tables']][['somatic']] = cbind(data.table(pass='pass'),vault[['tables']][['somatic']]) }
    vault$tables$somatic <- somatic.filter.fun(vault$tables$somatic, GOI)
  }
  ## germline
  if (!is.null(vault$tables$germline)) {
    if( !'pass' %in% colnames(vault$tables$germline) ){ vault[['tables']][['germline']] = cbind(data.table(pass='pass'),vault[['tables']][['germline']]) }
    vault$tables$germline <- germline.filter.fun(vault$tables$germline, GOI)
  }
  ## structural
  if (!is.null(vault$tables$structural)) {
    if( !'pass' %in% colnames(vault$tables$structural) ){ vault[['tables']][['structural']] = cbind(data.table(pass='pass'),vault[['tables']][['structural']]) }
    vault$tables$structural <-structural.filter.fun(vault$tables$structural)
  }
  ## fusion
  if (!is.null(vault$tables$fusion)) {
    if( !'pass' %in% colnames(vault$tables$fusion) ){ vault[['tables']][['fusion']] = cbind(data.table(pass='pass'),vault[['tables']][['fusion']]) }
    vault$tables$fusion <- fusion.filter.fun(vault$tables$fusion)
  }
  ## segment
  if (!is.null(vault$tables$segment)) {
    vault$tables$segment <- segment.filter.fun(vault$tables$segment)
  }
  return(vault)
}

#' Changes pass/fail value of variants based on consequence filter
#'
#' @param vault - vault
#'
#' @return filtered vault
.triageConsequenceVault <- function(vault) {
  ## somatic
  if (!is.null(vault$tables$somatic)) {
    if( !'pass' %in% colnames(vault$tables$somatic) ){ vault[['tables']][['somatic']] = cbind(data.table(pass='pass'),vault[['tables']][['somatic']]) }
    vault$tables$somatic <- consequence.filter.fun(vault$tables$somatic)
  }
  ## germline
  if (!is.null(vault$tables$germline)) {
    if( !'pass' %in% colnames(vault$tables$germline) ){ vault[['tables']][['germline']] = cbind(data.table(pass='pass'),vault[['tables']][['germline']]) }
    vault$tables$germline <- consequence.filter.fun(vault$tables$germline)
  }
  return(vault)
}

#' Changes pass/fail value of variants based on artifact filter
#'
#' @param vault - vault
#'
#' @return filtered vault
.triageArtifactVault <- function(vault) {
  ## somatic
  if (!is.null(vault$tables$somatic)) {
    if( !'pass' %in% colnames(vault$tables$somatic) ){ vault[['tables']][['somatic']] = cbind(data.table(pass='pass'),vault[['tables']][['somatic']]) }
    vault$tables$somatic <- artifact.filter.fun(vault$tables$somatic)
  }
  return(vault)
}

#' Changes pass/fail value of variants based on phase-id (pid) filter
#'
#' @param vault - vault
#'
#' @return filtered vault
.triagePhaseVault <- function(vault) {
  ## somatic
  if (!is.null(vault$tables$somatic)) {
    if( !'pass' %in% colnames(vault$tables$somatic) ){ vault[['tables']][['somatic']] = cbind(data.table(pass='pass'),vault[['tables']][['somatic']]) }
    vault$tables$somatic <- phase.filter.fun(vault$tables$somatic)
  }
  return(vault)
}

#' Removes variants that do not pass filters from a specific table
#'
#' @param V - filtered vault
#' @param TABLE - the table to remove failed variants
#' @return updated vault
.removeFailedVariantTable <- function(V, TABLE){
  V[['tables']][[TABLE]] <- V[['tables']][[TABLE]][triage=='PASS',]
  V[['maps']][[TABLE]] <- V[['maps']][[TABLE]][V[['tables']][[TABLE]],on=.(id, var_id), nomatch=0]
  return(V)
}




#### Rescue variant functions

#' Changes pass/fail value of variants based on provided table
#'
#' @param V - vault
#' @param TBL a table containing id,var_id,table,&pass value
#' @return updated vault
rescueVariants <- function(V, TBL) {
  ## Error handling
  if( is.null(TBL) ){ stop("No variant table provided; consider using .changeVariantPass for single variant") }
  if( !all(c('pass','id','var_id','table') %in% colnames(TBL)) ){ stop("The following cols needed in TBL: pass,id,var_id,table") }
  tmp = c(paste0(V$tables$somatic$id, '--', V$tables$somatic$var_id), paste0(V$tables$germline$id, '--', V$tables$germline$var_id),
          paste0(V$tables$structural$id, '--', V$tables$structural$var_id), paste0(V$tables$fusion$id, '--', V$tables$fusion$var_id))
  keys = paste0(TBL$id, '--', TBL$var_id)
  notInVault = keys[!keys %in% tmp]
  if( length(notInVault)>0 ){
    stop( paste0("The following id--var_id pairs are not in the vault:  ", paste0( notInVault, collapse = ', ')) )
    }
  if( !all(TBL$table %in% c('somatic', 'germline', 'structural', 'fusion')) ){ stop("A provided table is not present in vault") }
  ## change/return
  for(row in seq_len(nrow(TBL))) {
    V = rescueVariant(V, TBL$table[row], TBL$id[row], TBL$var_id[row], TBL$pass[row])
  }
  return(V)
}

#' Changes pass/fail value of a single variant based on provided values
#'
#' @param V - vault
#' @param TABLE The table the variant is in
#' @param IDD - The sample ID
#' @param VAR_IDD - the variant id
#' @param VALUE The P/F value to change to
#' @return updated vault
rescueVariant <- function(V, TABLE, IDD, VAR_IDD, VALUE){
  V[['tables']][[TABLE]][id==IDD & var_id==VAR_IDD, pass:=VALUE]
  return(V)
}




#### Inject variants

#' Inject new variants into the vault
#'
#' @param V a vault
#' @param TBL table of variants to add
#' @param FIELD which table to add variants to
#' @return updated vault
#'
injectVariants <- function(V, TBL, FIELD){
  ## make pass
  TBL[['pass']] <- "pass"
  ## get
  prev_table = getVaultTable(V, FIELD=FIELD)
  prev_map = getVaultMap(V, FIELD=FIELD)
  ## remove excess columns
  TBL <- TBL %>% subset(select = colnames(TBL) %in% colnames(prev_table))
  ## key
  prev_table$id.varid <- paste0(prev_table$id, '--', prev_table$var_id)
  TBL$id.varid <- paste0(TBL$id, '--', TBL$var_id)
  ## check
  # if( sum(colnames(prev_table) != colnames(TBL))>0 ){ print('TBL colnames do not match specified vault table FIELD'); break }
  if( sum(TBL$id.varid %in% prev_table$id.varid)>0 ){
    print(paste0("The following id-var_id pairs already present in table: ",
                 paste0(TBL$id.varid[TBL$id.varid %in% prev_table$id.varid], collapse = "; ")))
    prev_table$pass[prev_table$id.varid %in% TBL$id.varid] <- "pass"
    TBL <- TBL[!(TBL$id.varid %in% prev_table$id.varid),]
  }
  ## combine table
  new_table = plyr::rbind.fill(TBL, prev_table)
  new_table[['id.varid']] <- NULL
  ## combine map
  tmp = data.frame(id=TBL$id, var_id=TBL$var_id, GENE=TBL$GENE)
  new_map <- plyr::rbind.fill(tmp, prev_map)
  ## make into datatable
  new_table <- data.table(new_table)
  new_map <- data.table(new_map)
  ## set
  V <- setVaultTable(V, FIELD=FIELD, new_table)
  V <- setVaultMap(V, FIELD=FIELD, new_map)
  return(V)
}

#' Given a table, makes a map
#'
#' @param TBL vault table
#' @return vault map
#'
.makeMapFromTable <- function(TBL) {
  map_cols = c('id','var_id', 'GENE')
  if( sum(!map_cols %in% colnames(TBL))>0 ){ print(paste0('provided table missing 1+ required cols: ', paste0(map_cols, collapse = ', '))) }
  tmp = TBL[, ..map_cols]
  ## temporary
  colnames(tmp)[3] <- 'gene_id'
  return(tmp)
}

#' Create metaVault. Main function.
#'
#' @param VAULTS List of vaults
#' @param ANNO Object returned by vaultAnno
#' @return metaVault
metaVault <- function(VAULTS, ANNO=NULL) {

  # load vaults
  ids <- sapply(VAULTS, function(v){ return( getVaultMeta(v)[['id']] )} )
  names(VAULTS) <- ids
  print('Vaults loaded')

  # combine maps
  somatic.maps <- combineMapsTables(VAULTS, 'maps', 'somatic')
  germline.maps <- combineMapsTables(VAULTS, 'maps', 'germline')
  structural.maps <- combineMapsTables(VAULTS, 'maps', 'structural')
  fusion.maps <- combineMapsTables(VAULTS, 'maps', 'fusion')
  segment.maps <- combineMapsTables(VAULTS, 'maps', 'segment')
  if(nrow(somatic.maps)>0){setkey(somatic.maps, id, var_id)}
  if(nrow(germline.maps)>0){setkey(germline.maps, id, var_id)}
  if(nrow(structural.maps)>0){setkey(structural.maps, id, var_id, end)}
  if(nrow(fusion.maps)>0){setkey(fusion.maps, id, var_id, end)}
  if(nrow(segment.maps)>0){setkey(segment.maps, id, var_id)}
  print('Maps combined')

  # combine tables
  somatic.tables <- combineMapsTables(VAULTS, 'tables', 'somatic')
  germline.tables <- combineMapsTables(VAULTS, 'tables', 'germline')
  structural.tables <- combineMapsTables(VAULTS, 'tables', 'structural')
  fusion.tables <- combineMapsTables(VAULTS, 'tables', 'fusion')
  segment.tables <- combineMapsTables(VAULTS, 'tables', 'segment')
  gene.copy <- combineMapsTables(VAULTS, 'tables', 'gene.copy')
  gene.expression <- combineMapsTables(VAULTS, 'tables', 'gene.expression')
  if(nrow(somatic.tables)>0){setkey(somatic.tables, id, var_id)}
  if(nrow(germline.tables)>0){setkey(germline.tables, id, var_id)}
  if(nrow(structural.tables)>0){setkey(structural.tables, id, var_id)}
  if(nrow(fusion.tables)>0){setkey(fusion.tables, id, var_id)}
  if(nrow(segment.tables)>0){setkey(segment.tables, id, var_id)}
  if(nrow(gene.copy)>0){setkey(gene.copy, id, GENE)}
  if(nrow(gene.expression)>0){setkey(gene.expression, id, GENE)}
  print('Tables combined')

  # combine meta
  meta <- combineMeta(VAULTS)
  print('Meta combined')

  # combine QC
  qc <- combineQC(VAULTS)

  # set null PRN
  if( nrow(somatic.tables)==0 ){ somatic.tables = somatic.maps = NULL }
  if( nrow(germline.tables)==0 ){ germline.tables = germline.maps = NULL }
  if( nrow(structural.tables)==0 ){ structural.tables = structural.maps = NULL }
  if( nrow(fusion.tables)==0 ){ fusion.tables = fusion.maps = NULL }
  if( nrow(segment.tables)==0 ){ segment.tables = segment.maps = NULL }
  if( nrow(gene.copy)==0 ){ gene.copy = NULL }
  if( nrow(gene.expression)==0 ){ gene.expression = NULL }

  # create metaVault
  mv <- list(
    maps = list(somatic = somatic.maps,
                germline = germline.maps,
                structural = structural.maps,
                fusion = fusion.maps,
                segment = segment.maps),
    tables = list(somatic = somatic.tables,
                  germline = germline.tables,
                  structural = structural.tables,
                  fusion = fusion.tables,
                  segment = segment.tables,
                  gene.copy = gene.copy,
                  gene.expression = gene.expression),
    meta = meta,
    qc = qc,
    anno = ANNO
  )

  # return
  return(mv)
}

#' Combines tables from lists of vaults. Wrapper fxn.
#'
#' @param VAULTS A list of filtered vault object
#' @param FIELD1 Table/map of vault
#' @param FIELD2 Which table/map
#' @return data.table of combined TABLE with id ID
combineMapsTables <- function(VAULTS, FIELD1, FIELD2){
  if( FIELD1=='maps' ) { tmp <- lapply(VAULTS, getVaultMap, FIELD=FIELD2) }
  if( FIELD1=='tables' ){ tmp <- lapply(VAULTS, getVaultTable, FIELD=FIELD2) }
  tmp <- tmp[!sapply(tmp, is.null)] # this should not be necessary bug segfault
  tmp <- rbindlist(tmp)
  return(tmp)
}

#' Creates data.table from metadata of a list of vaults
#'
#' @param VAULTS List of vaults
#' @return data.table of meta data from all vaults
combineMeta <- function(VAULTS) {
  lapply(VAULTS, getVaultMeta) %>%
    rbindlist(.) %>%
    return()
}

#' Creates data.table from qc of a list of vaults
#'
#' @param VAULTS List of vaults
#' @return data.table of qc data from all vaults
combineQC <- function(VAULTS) {
  qc <- lapply(VAULTS, getVaultQC) ## TODO: turn into a data.table
  return(qc)
}

#' Subsets vault(s) from a metavault
#'
#' @param METAVAULT A metaVault
#' @param IDENTIFIER A ID within the metaVault
#' @return vault for that id
extractVault <- function(IDENTIFIERS, METAVAULT) {

  # extract meta
  meta <- METAVAULT$meta[id %in% IDENTIFIERS]
  
  # set to null at baseline
  somatic.maps <- germline.maps <- structural.maps <- fusion.maps <- segment.maps <- NULL
  somatic.tables <- germline.tables <- structural.tables <- fusion.tables <- segment.tables <- gene.copy <- gene.expression <- NULL
  
  # extract maps if not null
  if( !is.null(METAVAULT$maps$somatic) ){ somatic.maps <- METAVAULT$maps$somatic[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$maps$germline) ){ germline.maps <- METAVAULT$maps$germline[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$maps$structural) ){ structural.maps <- METAVAULT$maps$structural[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$maps$fusion) ){ fusion.maps <- METAVAULT$maps$fusion[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$maps$segment) ){ segment.maps <- METAVAULT$maps$segment[id %in% IDENTIFIERS] }

  # extract tables
  if( !is.null(METAVAULT$tables$somatic) ){ somatic.tables <- METAVAULT$tables$somatic[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$germline) ){ germline.tables <- METAVAULT$tables$germline[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$structural) ){ structural.tables <- METAVAULT$tables$structural[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$fusion) ){ fusion.tables <- METAVAULT$tables$fusion[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$segment) ){ segment.tables <- METAVAULT$tables$segment[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$gene.copy) ){ gene.copy <- METAVAULT$tables$gene.copy[id %in% IDENTIFIERS] }
  if( !is.null(METAVAULT$tables$gene.expression) ){ gene.expression <- METAVAULT$tables$gene.expression[id %in% IDENTIFIERS] }
    
  # handle empty maps/tables
  if(length(somatic.tables$id) == 0){ somatic.maps = somatic.tables = NULL }
  if(length(germline.tables$id) == 0){ germline.maps = germline.tables = NULL }
  if(length(structural.tables$id) == 0){ structural.maps = structural.tables = NULL }
  if(length(fusion.tables$id) == 0){ fusion.maps = fusion.tables = NULL }
  if(length(segment.tables$id) == 0){ segment.maps = segment.tables = NULL }

  # export
  v <- list(
    maps = list(somatic = somatic.maps,
                germline = germline.maps,
                structural = structural.maps,
                fusion = fusion.maps,
                segment = segment.maps),
    tables = list(somatic = somatic.tables,
                  germline = germline.tables,
                  structural = structural.tables,
                  fusion = fusion.tables,
                  segment = segment.tables,
                  gene.copy = gene.copy,
                  gene.expression = gene.expression),
    meta = meta
  )
  return(v)
}

#' Subsets a metavault to only contain select genes by id
#'
#' @param GENEIDS vector of gene ids
#' @param METAVAULT metaVault
#' @param SYMBOLS vector of gene symbols (optional)
#' @param ANNO vault annotation object
#' @return vault for select ids
extractGene <- function(METAVAULT, GENEIDS=NULL, SYMBOLS=NULL, ANNO=NULL) {

  if (!is.null(SYMBOLS) && !is.null(ANNO)) {
    converted <- ANNO$genes[ANNO$genes$gene_name %in% SYMBOLS]$gene_id
    GENEIDS <- unique(c(GENEIDS, converted))
  }

  msom <- tsom <- NULL
  if (!is.null(METAVAULT$maps$somatic) && !is.null(METAVAULT$tables$somatic)) {
    msom <- METAVAULT$maps$somatic[gene_id %in% GENEIDS]
    tsom <- METAVAULT$tables$somatic[msom]
  }

  mgrm <- tgrm <- NULL
  if (!is.null(METAVAULT$maps$germline) && !is.null(METAVAULT$tables$germline)) {
    mgrm <- METAVAULT$maps$germline[gene_id %in% GENEIDS]
    tgrm <- METAVAULT$tables$germline[mgrm]
  }

  mstr <- tstr <- NULL
  if (!is.null(METAVAULT$maps$structural) && !is.null(METAVAULT$tables$structural)) {
    mstr <- METAVAULT$maps$structural[gene_id %in% GENEIDS]
    tstr <- METAVAULT$tables$structural[mstr]
  }

  mfus <- tfus <- NULL
  if (!is.null(METAVAULT$maps$fusion) && !is.null(METAVAULT$tables$fusion)) {
    mfus <- METAVAULT$maps$fusion[gene_id %in% GENEIDS]
    tfus <- METAVAULT$tables$fusion[mfus]
  }

  mseg <- tseg <- NULL
  if (!is.null(METAVAULT$maps$segment) && !is.null(METAVAULT$tables$segment)) {
    mseg <- METAVAULT$maps$segment[gene_id %in% GENEIDS]
    tseg <- METAVAULT$tables$segment[mfus]
  }

  if (!is.null(METAVAULT$tables$gene.copy)) {
    tgc <- METAVAULT$tables$gene.copy[GENE %in% GENEIDS]
  }

  if (!is.null(METAVAULT$tables$gene.expression)) {
    tge <- METAVAULT$tables$gene.expression[GENE %in% GENEIDS]
  }

  # export
  mv <- list(
    maps = list(somatic = msom,
                germline = mgrm,
                structural = mstr,
                fusion = mfus,
                segment = mseg),
    tables = list(somatic = tsom,
                  germline = tgrm,
                  structural = tstr,
                  fusion = tfus,
                  segment = tseg,
                  gene.copy = tgc,
                  gene.expression = tge),
    meta = METAVAULT$meta
  )
  return(mv)
}

#' Creates a list of vaults from a metavault
#'
#' @param METAVAULT A metaVault
#' @return list of vaults
splitMetaVault <- function(METAVAULT) {
    ids <- as.list(METAVAULT$meta$id)
    vs <- lapply(ids, extractVault, METAVAULT=METAVAULT)
    return(vs)
}
