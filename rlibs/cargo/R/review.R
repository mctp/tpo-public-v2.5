
som_cols <- c('id', 'triage', 'gene_name', 'var_id', 'HGVSp', 'Consequence', 'IMPACT',
              'AFT', 'AFN', 'ADT', 'ADN', 'DPN', 'ADT_FWD', 'ADT_REV', 'str', 'tlod', 'nlod', 
              'dbsnp', 'clinid', 'clinvar', 'cosmic_cnt', 'gnomad_af', 'kg_af', 'pid', 'ccf')
grm_cols <- c('id', 'triage', 'gene_name', 'var_id', 'HGVSp', 'Consequence', 'IMPACT',
              'AFT', 'AFN', 'ADT', 'ADN', 'DPN', 'str',
              'dbsnp', 'clinid', 'clinvar', 'cosmic_cnt', 'gnomad_af', 'kg_af', 'pid')
str_cols <- c('id', 'triage', 'gene_name.1', 'gene_id.1', 'gene_name.2', 'gene_id.2', 
              'CSQ1', 'CSQ2', 'var_id', 'bnd.id', 'topo', 'insert', 
              'QUAL', 'AFT', "ADT", 'DPT', 'MANTA.FILTER', 'MANTA.SCORE','MANTA.CONTIG')
fus_cols <- c('id', 'triage', 'gene.5', 'gene.3', 'var_id', 'junction.5','junction.3',
              "spanning.reads", "encompassing.reads", "total.reads", "breakpoints",
              "spliced","HQ","TS","inframe","distance","topology","mm2.valid","gmap.valid",        
              "read fraction 5'", "read fraction 3'", "recurrent(5';3')", "repetitive(5';3')",
              "chain", "contig")




# variant
###########

#' generates table for review; subsets cols from somatic table & adds from meta
#'
#' @param MV A metavault 
#' @param GENE_NAME single gene name
#' @return A table for manual review of variants
#' 
#' @export
getReviewTableSomatic <- function(MV, GENE_NAME){
  tbl = MV$tables$somatic %>%
    ## filter
    filter(gene_name==GENE_NAME) %>%
    filter(!grepl('artifact|consequence|annotation', triage)) %>%
    ## merge
    merge(., MV$meta, by='id') %>%
    ## reformat
    select(any_of(som_cols), purity, ploidy) %>%
    mutate(flag = ifelse( (grepl('path', clinvar, ignore.case=T) | cosmic_cnt>5) & AFT>0.02, 'FLAG', '')) %>%
    mutate(decision='', notes='') %>%
    arrange(triage, desc(flag), Consequence)
  return(tbl)
}

#' generates table for review; subsets cols from germline table & adds from meta
#'
#' @param MV A metavault 
#' @param GENE_NAME single gene name
#' @return A table for manual review of variants
#' 
#' @export
getReviewTableGermline <- function(MV, GENE_NAME){
  tbl = MV$tables$germline %>%
    ## filter
    filter(gene_name==GENE_NAME) %>%
    filter(!grepl('artifact|consequence|annotation', triage)) %>%
    ## merge
    merge(., MV$meta, by='id') %>%
    ## reformat
    select(any_of(grm_cols), purity, ploidy) %>%
    mutate(flag = ifelse( (grepl('path', clinvar, ignore.case=T) | cosmic_cnt>5) & AFT>0.2 & AFN>0.2, 'FLAG', '')) %>%
    mutate(decision='', notes='') %>%
    arrange(triage, desc(flag), Consequence)
  return(tbl)
}

#' generates table for review; subsets cols from fusion table & adds from meta
#'
#' @param MV A metavault 
#' @param GENE_NAME single gene name
#' @return A table for manual review of variants
#' 
#' @export
getReviewTableFusion <- function(MV, GENE_NAME){
  ## filter string
  str <- paste0( paste0('^',GENE_NAME, '$'), '|', paste0('^', GENE_NAME, ':'), '|', paste0(':', GENE_NAME, '$'), '|', paste0(':', GENE_NAME, ':'))
  ## 
  tbl = MV$tables$fusion %>%
    ## filter
    filter( grepl(str, gene_names.5.1) |  grepl(str, gene_names.5.2) | grepl(str, gene_names.3.1) | grepl(str, gene_names.3.2)) %>%
    filter(!grepl('annotation', triage)) %>%
    ## reformat
    .fusionFormat(., REORDER=T) %>%
    merge(., MV$meta, by='id') %>%
    select(any_of(fus_cols), purity, ploidy) %>%
    mutate(flag = ifelse( (HQ=='yes'&mm2.valid=='yes'&gmap.valid=='yes'&total.reads>5), 'FLAG', '')) %>%
    mutate(decision='', notes='') %>%
    arrange(id, triage, desc(flag), desc(inframe), desc(spliced))
  return(tbl)
}

#' generates table for review; subsets cols from structural table & adds from meta
#'
#' @param MV A metavault 
#' @param GENE_NAME single gene name
#' @return A table for manual review of variants
#' 
#' @export
getReviewTableStructural <- function(MV, GENE_NAME){
  ## 
  tbl = MV$tables$structural %>%
    ## filter
    filter( grepl(paste0('^',GENE_NAME, '$'), gene_name.1) |  grepl(paste0('^',GENE_NAME, '$'), gene_name.2) ) %>%
    filter(!grepl('annotation', triage)) %>%
    ## reformat
    merge(., MV$meta, by='id') %>%
    select(any_of(str_cols), purity, ploidy) %>%
    mutate(flag = ifelse( (MANTA.FILTER=='PASS' & AFT>0.1 & ADT>5), 'FLAG', '')) %>%
    mutate(decision='', notes='') %>%
    arrange(triage, desc(flag))
  return(tbl)
}

#' reformats review table for storage
#'
#' @param TBL output of getReviewTable fxns; post manual review
#' @param FIELD which table the variant originates from 
#' @return A table that can be used with 
#' 
#' @export
formatReviewTable <- function(TBL, FIELD){
  tbl <- TBL %>%
    ## reformat
    mutate(table=FIELD) %>%
    mutate(triage = ifelse(grepl('pass|keep|override', decision, ignore.case=T), 'keep', '')) %>%
    mutate(triage = ifelse(grepl('fail|manual|remove', decision, ignore.case=T), 'manual', triage)) %>%
    select(triage, id, var_id, table, notes) 
  return(tbl)
}




# Feature
###########

#' compiles data types for manual review of sample-level features
#' i.e. is a sample SPOP-mutant or WNT-altered
#'
#' @param MVF filtered metavault
#' @param FEATURES vector of features; gene_names
#' @return A table that can be reviewed by hand
#' 
#' @export
reviewFeatures <- function(MVF, FEATURES){
  
  ## filter vault
  mvf <- vaultFilter(MVF)
  
  ## gather
  som <- .reviewFeaturesSomatic(mvf, FEATURES)
  grm <- .reviewFeaturesGermline(mvf, FEATURES)
  gc <- .reviewFeaturesGeneCopy(mvf, FEATURES)
  exp <- .reviewFeaturesExpression(mvf, FEATURES)
  meta <- MVF$meta %>% select(id, purity, ploidy)
  
  ## hits
  tbl <- mvf %>%  
    tableHitsMetavault(., 
                       HIT_FIELDS = c('somatic', 'germline', 'structural', 'fusion', 'gene.copy'), 
                       GENES = FEATURES) %>%
    ## add
    mutate(key = paste0(id, '--', gene_id)) %>%
    merge(., som, by='key', all.x=T) %>%
    merge(., grm, by='key', all.x=T) %>%
    merge(., meta, by='id', all.x=T) %>%
    merge(., gc, by='key', all.x=T) %>%
    merge(., exp, by='key', all.x=T) %>%
    ## filter
    filter( hits>0 | nseg>1 ) %>%
    filter( !(gene.copy==1 & hits==1) | is.na(gene.copy) ) %>%
    ## reformat
    mutate(decision='', notes='') %>%
    arrange(gene_name, gene.copy, somatic, germline, structural, fusion) 
  
  return(tbl)
  
}

#' Compiles variants
#'
#' @param MVF filtered metavault
#' @param FEATURES vector of features; gene_names
#' @return A table of pertinent variant cols
#' 
.reviewFeaturesSomatic <- function(MVF, FEATURES){
  tbl <- MVF$tables$somatic %>%
    ## filter
    filter(gene_name %in% FEATURES) %>%
    ## reformat
    mutate(key = paste0(id, '--', gene_id)) %>%
    distinct(key, Consequence, HGVSp, AFT, ccf, clinvar, cosmic_cnt, gnomad_af) %>%
    group_by(key) %>%
    mutate(som_consequence = paste0(Consequence, collapse='; ')) %>%
    mutate(som_hgvsp = paste0(HGVSp, collapse='; ')) %>%
    mutate(som_aft = paste0(AFT, collapse='; ')) %>%
    mutate(som_ccf = paste0(ccf, collapse='; ')) %>%
    mutate(som_clinvar = paste0(clinvar, collapse='; ')) %>%
    mutate(som_cosmic_cnt = paste0(cosmic_cnt, collapse='; ')) %>%
    mutate(som_gnomad_af = paste0(gnomad_af, collapse='; ')) %>%
    ungroup() %>%
    distinct(key, som_hgvsp, som_consequence, som_aft, som_ccf, som_clinvar, som_cosmic_cnt, som_gnomad_af)
  return(tbl)
}

#' Compiles variants
#'
#' @param MVF filtered metavault
#' @param FEATURES vector of features; gene_names
#' @return A table of pertinent variant cols
#' 
.reviewFeaturesGermline <- function(MVF, FEATURES){
  tbl <- MVF$tables$germline %>%
    ## filter
    filter(gene_name %in% FEATURES) %>%
    ## reformat
    mutate(key = paste0(id, '--', gene_id)) %>%
    distinct(key, Consequence, HGVSp, clinvar, cosmic_cnt, gnomad_af) %>%
    group_by(key) %>%
    mutate(grm_consequence = paste0(Consequence, collapse='; ')) %>%
    mutate(grm_hgvsp = paste0(HGVSp, collapse='; ')) %>%
    mutate(grm_clinvar = paste0(clinvar, collapse='; ')) %>%
    mutate(grm_cosmic_cnt = paste0(cosmic_cnt, collapse='; ')) %>%
    mutate(grm_gnomad_af = paste0(gnomad_af, collapse='; ')) %>%
    ungroup() %>%
    distinct(key, grm_hgvsp, grm_consequence, grm_clinvar, grm_cosmic_cnt, grm_gnomad_af)
}

#' Compiles variants
#'
#' @param MVF filtered metavault
#' @param FEATURES vector of features; gene_names
#' @return A table of pertinent variant cols
#' 
.reviewFeaturesGeneCopy <- function(MVF, FEATURES){
  tbl <- MVF$tables$gene.copy %>%
    ## filter
    filter(gene_name %in% FEATURES) %>%
    ## reformat
    mutate(key = paste0(id, '--', gene_id)) %>%
    distinct(key, nseg)
}

#' Compiles expression data
#'
#' @param MVF filtered metavault
#' @param FEATURES vector of features; gene_names
#' @return trpkm for features
#' 
.reviewFeaturesExpression <- function(MVF, FEATURES){
  if( is.null(MVF$tables$gene.expression) ){ return(NULL) }
  tbl <- MVF$tables$gene.expression %>%
    filter(gene_name %in% FEATURES) %>%
    mutate(key=paste0(id, '--', gene_id)) %>%
    distinct(key, trpkm)
  return(tbl)
}



# misc
###########

#' reformats fusion table for easier viewing; 
#' slightly modified version of svFormat
#'
#' @param TBL fusion table
#' @param REORDER should the table be reordered based on 'chain'
#' @return A user-friendly fusion table
#' 
#' @export
.fusionFormat <- function(TBL, REORDER=TRUE) {
  
  rep <- TBL[,.(
    id,
    triage,
    var_id,
    gene.5=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
                  ifelse(gene_names.3.2!="", sprintf("(%s)", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
    junction.5=sprintf("%s:%s:%s", chr.5, pos.5, str.5),
    gene.3=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
                  ifelse(gene_names.5.2!="", sprintf("(%s)", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
    junction.3=sprintf("%s:%s:%s", chr.3, pos.3, str.3),
    spanning.reads=sum.jnc * (l1=="spn"),
    encompassing.reads=tot.enc,
    total.reads=tot.jnc+tot.enc,
    breakpoints=sum.bpt,
    spliced=ifelse(d2a, "yes", ""),
    HQ=ifelse(hq.bpt, "yes", ""),
    TS=ifelse(ts.warn, "yes", ""),
    inframe=ifelse(orf, "yes", ""),
    distance=dst,
    topology=topo,
    mm2.valid=paste(ifelse(mm2.valid, "yes", "")),
    gmap.valid=paste(ifelse(gmap.valid, "yes", "")),
    "read fraction 5'" = format(round(sum.jnc / (tot.jnc.5 + tot.sp.jnc.5), 2), nsmall=2, scientific=999),
    "read fraction 3'" = format(round(sum.jnc / (tot.jnc.3 + tot.sp.jnc.3), 2), nsmall=2, scientific=999),
    "recurrent(5';3')" = paste(unq.rec.5, unq.rec.3, sep=";"),
    "repetitive(5';3')" = paste(art.5, art.3, sep=";"),
    chain=sv.chain,
    contig=sapply(lapply(ctg.seq, as.character), "[", 1)
  )]
  if (REORDER) {
    rep <- rep[order(chain)]
  }
  return(rep)
}

