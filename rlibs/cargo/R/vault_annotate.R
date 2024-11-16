
#### combine

#' Annotates vault
#'
#' @param VAULT The vault object
#' @param GOBJ the cnvex Gobj via getGobj()
#' @param OPTS cnvex opts via getOpts()
#' @return List of annotated data
#' @export
#'
annotateVault <- function(VAULT, GOBJ, OPTS) {
  id = VAULT$meta |> pull(id);print(id)
  # collect
  #########
  ## somatic mutations
  if( nrow(VAULT[['tables']][['somatic']])!=0 ) {
    som = annotateMutations(VAULT, 'somatic')
  } else{som = NULL}
  ## germline mutations
  if( nrow(VAULT[['tables']][['germline']])!=0 ) {
    grm = annotateMutations(VAULT, 'germline')
  } else{grm = NULL}
  ## structural mutations
  if( nrow(VAULT[['tables']][['structural']])!=0 ) {
    str = annotateMutationsStructural(VAULT)
  } else{str = NULL}
  ## gene copy number
  if( nrow(VAULT[['tables']][['gene.copy']])!=0 ) {
    gene.copy = annotateGeneCopyNumber(VAULT)
  } else{gene.copy = NULL}
  ## gene expression
  if( nrow(VAULT[['tables']][['gene.expression']])!=0 ) {
    exp = annotateGeneExpression(VAULT)
  } else{exp = NULL}
  ## fusion
  if( nrow(VAULT[['tables']][['fusion']])!=0 ) {
    fus = annotateFusions(VAULT)
  } else{fus = NULL}
  ## Chromosome gain/loss
  if( nrow(VAULT[['tables']][['segment']])!=0 ) {
    arms = annotateArms(VAULT)
    chromosomes = annotateChromosomes(arms)
  } else{
    arms = NULL
    chromosomes = NULL
  }
  ## genomic instability
  if( nrow(VAULT[['tables']][['segment']])!=0 ) {
    cin = annotateCin(VAULT, GOBJ, OPTS)
    ## arm level totals
    cin[['chr_whole']] = .tallyChrAlterations(arms, chromosomes, 'chr_whole')
    cin[['chr_whole_gain']] = .tallyChrAlterations(arms, chromosomes, 'chr_whole_gain')
    cin[['chr_whole_loss']] = .tallyChrAlterations(arms, chromosomes, 'chr_whole_loss')
    cin[['arm_only']] = .tallyChrAlterations(arms, chromosomes, 'arm_only')
    cin[['arm_only_gain']] = .tallyChrAlterations(arms, chromosomes, 'arm_only_gain')
    cin[['arm_only_loss']] = .tallyChrAlterations(arms, chromosomes, 'arm_only_loss')
    cin[['arm_total']] = .tallyChrAlterations(arms, chromosomes, 'arm_total')
    cin[['arm_total_gain']] = .tallyChrAlterations(arms, chromosomes, 'arm_total_gain')
    cin[['arm_total_loss']] = .tallyChrAlterations(arms, chromosomes, 'arm_total_loss')
  } else{cin = NULL}
  #########

  # combine
  #########
  annotation <- list(
    somatic = som,
    germline = grm,
    structural = str,
    gene.copy = gene.copy,
    gene.expression = exp,
    fusions = fus,
    cin = cin,
    arms = arms,
    chromosomes = chromosomes,
    meta = VAULT[['meta']]
  )
  return(annotation)
  #########
}




#### gene copy number

#' Given the vault, gets the average basepair-weighted C value for each chromosome arm.
#'
#' @param VAULT The vault object
#' @return Data table with gene, CN annotation, and required information
#' @export
annotateGeneCopyNumber <- function(VAULT){
  if(!is.null(VAULT$tables$gene.copy)){
    if( !all(is.na(VAULT$tables$segment$C)) ) {

      ## get C/K
      tmp = VAULT$tables$gene.copy
      gene = GRanges(tmp[!is.na(tmp$chr), c('chr', 'start', 'end', 'id', 'gene_id', 'gene_name', 'Cmin', 'Kmin')])

      ## get arm average
      aa <- getArmAverage(VAULT)

      ## ensure gene overlaps with single chr arm (per ENSG00000276128)
      hits <- GenomicRanges::findOverlaps(gene, aa, ignore.strand = T, type = 'any')
      hits_df <- as.data.table(hits)
      hits_df[['width']] <- width(pintersect(gene[queryHits(hits)], aa[subjectHits(hits)]))
      data.table::setorder(hits_df, cols = "queryHits", -"width")
      overlaps <- hits_df[!duplicated(hits_df[['queryHits']])]

      ## add arm average to genes
      tmp <- as.data.table( aa[ overlaps[['subjectHits']] , ] )[,c('chr_arm', 'armAvg')]
      gaa <- cbind(as.data.table(gene), tmp)

      ## weighted median
      gaa[['weighted_median']] <- getWeightedMedian(VAULT)

      ## segment length
      tmp <- getGeneSegLength(VAULT)
      gaasl <- merge(gaa,tmp, by.x='gene_id', by.y='gene_id', all.x = T)

      ## annotate
      tmp <- labelCn(gaasl, VAULT$meta$sex)
      gaasl[['annotation']] = tmp

      ## reformat
      out = gaasl[,c('id','seqnames','start','end','width','strand','gene_id','gene_name','annotation','Cmin','Kmin','chr_arm','armAvg','weighted_median','seg_width')]
      colnames(out)[c(2,10,11,13)] = c('chr','c_min','k_min','arm_avg')
      return(out)
    }
  }
  return(NULL)
}

#' Given the vault, gets the average basepair-weighted C value for each chromosome arm.
#'
#' @param VAULT The vault object
#' @return data.table with gene and arm average
#' @export
getArmAverage <- function(VAULT){

  # get arm positional info
  ###########################
  assembly <- VAULT$meta$assembly
  chr_arms <- getArmInfo(assembly)
  ###########################

  # get weighted arm average
  ###########################
  segments <- GenomicRanges::GRanges(VAULT$tables$segment) # segment data from vault
  chr_arms$armAvg <- NA
  # iterate through chromosome arms
  for(arm in 1:length(chr_arms)){
    ## overlaps and their C values
    overlaps <- to(GenomicRanges::findOverlaps(chr_arms[arm], segments, ignore.strand = T, type = 'any'))
    c_vals <- as.data.table(segments[overlaps, 'C'])$C
    ## get weighted C value
    if(length(overlaps) > 0 & length(c_vals) > 0){
      ## include only overlapping bp
      overlaps <- as.data.table(suppressWarnings(GenomicRanges::pintersect(IRanges::findOverlapPairs(chr_arms[arm], segments))))
      ## weighted average C
      overlaps$len = overlaps$end - overlaps$start
      overlaps$cvals <- c_vals
      overlaps$weightedC <- overlaps$cvals*overlaps$len
      aa <- sum(overlaps$weightedC, na.rm = T)/sum(overlaps$len, na.rm = T)
    } else{aa <- NA}
    chr_arms[arm]$armAvg <- aa
  }
  return(chr_arms)
  ###########################
}

#' Pulls the coordinates for chromosome arms given an assembly
#'
#' @param ASSEMBLY The vault's assembly information
#' @return A GRanges object with the chromosome arm information
#' TODO: Handle other assemblies
getArmInfo <- function(ASSEMBLY){
  # hg38
  if(grepl('hg38', ASSEMBLY)){
    chr_arms <- read.table(system.file("extdata/hg38/arm.bed.gz", package="cnvex")) %>%
      select(-c(V5,V6)) %>%
      set_colnames(c("chrom","chromStart","chromEnd", 'chr_arm')) %>%
      GenomicRanges::GRanges()
    return(chr_arms)
  }
  # mm10
  if(grepl('mm10', ASSEMBLY)){
    print('mm10 not yet supported - see Ryan/Noshad about including a chr arms file')
    return(NA)
  }
}

#' Gets weighted median of sample.
#'
#' @param VAULT The vault object
#' @return the weighted median of the library
getWeightedMedian <- function(VAULT) {
  if(!is.null( VAULT[['tables']][['segment']] )) {
    if(!all(is.na( VAULT[['tables']][['segment']][['C']] ))) {
      wm = weighted.median(x =VAULT$tables$segment$C, w=as.numeric(VAULT$tables$segment$width), na.rm = T)
      return(wm)
    }
  }
  return(NULL)
}

#' Gets weighted median of sample.
#'
#' @param VAULT The vault object
#' @return the weighted median of the library
#' TODO: This might be leaking
getGeneSegLength <- function(VAULT) {
  VAULT$tables$segment %>%
    ## combine genes w/ segment lengths
    select(var_id, width) %>%
    rename(seg_width = width) %>%
    merge(., VAULT$maps$segment, by = 'var_id', all.y = T) %>%
    ## take shortest length per gene
    arrange(desc(seg_width)) %>%
    filter(!duplicated(gene_id)) %>%
    select(gene_id, seg_width) %>%
    return()
}

#' Wrapper script to annotate cnvex output.
#'
#' @param INTERMEDIATE An intermediate object from annotateGeneCopyNumber()
#' @param SEX Sex of sample from vault object
#' @return The copy number labels
labelCn <- function(INTERMEDIATE, SEX) {
  apply(INTERMEDIATE, MARGIN = 1, FUN = .labelCN, SEX) %>%
    return()
}

#' Annotates copy number of a row based on
#' C,K,sex, and chromosome arm average C value
#'
#' @param ROW The unified output corresponding to a single gene
#' @param SEX Sex of sample from vault object
#' @return The annotation for the gene
.labelCN <- function(ROW, SEX){
  S = SEX
  C = as.numeric(ROW[['Cmin']])
  K = as.numeric(ROW[['Kmin']])
  A = ROW[['chr_arm']]
  AA = as.numeric(ROW[['armAvg']])
  WM = as.numeric(ROW[['weighted_median']])
  SW = as.numeric(ROW[['seg_width']])
  if(is.na(C)){return(NA)}
  if(is.na(SW)){SW = -1}
  if(C==0 & !grepl('X|Y', A)) {
    return('Homo-del')
  } else if(C==0 & ((S=='female' & grepl('X', A)) | S=='male') ) {
    return('Homo-del')
  } else if(C>=8 & C >= (AA*2)) {
    return('Amplification')
  } else if(C>3 & C>WM & SW<=5e6 & SW >-1){
    return('Gain')
  } else if ( (C==1 | (!is.na(K)&K==0)) & (S=='female' | (S=='male' & !grepl('X', A)) ) ) {
    return('LOH')
  } else {
    return('None')
  }
}

#' Tallies chromosomes gained/lost in several ways
#'
#' @param ARMS The output of annotateArms()
#' @param VALUE which tally to return
#' @return integer of total number of whole chromosomes gained/lost in defined way
.tallyChrAlterations <- function(ARMS, CHROMOSOMES, VALUE){
  ## calculate
  tmp = ARMS %>%
    dplyr::filter(annotation != 'none') %>%
    dplyr::group_by(chr, annotation) %>%
    dplyr::mutate(n = dplyr::n_distinct(arm))
  ## return
  val = 0
  ## whole
  if( VALUE=='chr_whole' ){ try( (val=table(tmp$n==2)[[2]]/2), silent = T) }
  if( VALUE=='chr_whole_gain' ){ try( (val=table(tmp$annotation=='gain' & tmp$n==2)[[2]]/2), silent=T) }
  if( VALUE=='chr_whole_loss' ){ try( (val=table(tmp$annotation=='loss' & tmp$n==2)[[2]]/2), silent=T) }
  ## arm only
  if( VALUE=='arm_only' ){ try( (val=table(tmp$n==1)[[2]] ), silent = T) }
  if( VALUE=='arm_only_gain' ){ try( (val=table(tmp$annotation=='gain' & tmp$n==1)[[2]] ), silent = T) }
  if( VALUE=='arm_only_loss' ){ try( (val=table(tmp$annotation=='loss' & tmp$n==1)[[2]] ), silent = T) }
  ## arm all
  if( VALUE=='arm_total' ){ try( (val=length(tmp$arm)), silent = T) }
  if( VALUE=='arm_total_gain' ){ try( (val=table(tmp$annotation=='gain')[[2]] ), silent = T) }
  if( VALUE=='arm_total_loss' ){ try( (val=table(tmp$annotation=='loss')[[2]] ), silent = T) }
  return(val)
}




#### fusions

#' Annotates fusions
#'
#' @param VAULT The vault object
#' @return The annotated fusions with ENSG/gene names
#' @export
annotateFusions <- function(VAULT){
  if(!is.null(VAULT[['tables']][['fusion']])) {

    ## replace blank ids with chr arms
    int = VAULT[['tables']][['fusion']]
    int[['gene_ids.5.1']] = ifelse(int[['gene_ids.5.1']]=='' & int[['gene_ids.3.1']]!='', int[['cyt.5.1']], int[['gene_ids.5.1']])
    int[['gene_ids.3.1']] = ifelse(int[['gene_ids.3.1']]=='' & int[['gene_ids.5.1']]!='', int[['cyt.3.1']], int[['gene_ids.3.1']])
    int[['gene_ids.5.2']] = ifelse(int[['gene_ids.5.2']]=='' & int[['gene_ids.3.2']]!='', int[['cyt.5.2']], int[['gene_ids.5.2']])
    int[['gene_ids.3.2']] = ifelse(int[['gene_ids.3.2']]=='' & int[['gene_ids.5.2']]!='', int[['cyt.3.2']], int[['gene_ids.3.2']])
    ## replace blank symbols with chr arms
    int[['gene_names.5.1']] = ifelse(int[['gene_names.5.1']]=='' & int[['gene_names.3.1']]!='', int[['cyt.5.1']], int[['gene_names.5.1']])
    int[['gene_names.3.1']] = ifelse(int[['gene_names.3.1']]=='' & int[['gene_names.5.1']]!='', int[['cyt.3.1']], int[['gene_names.3.1']])
    int[['gene_names.5.2']] = ifelse(int[['gene_names.5.2']]=='' & int[['gene_names.3.2']]!='', int[['cyt.5.2']], int[['gene_names.5.2']])
    int[['gene_names.3.2']] = ifelse(int[['gene_names.3.2']]=='' & int[['gene_names.5.2']]!='', int[['cyt.3.2']], int[['gene_names.3.2']])

    # reformat
    int_51 = .annotateFusions(int, '5.1')
    int_31 = .annotateFusions(int, '3.1')
    int_52 = .annotateFusions(int, '5.2')
    int_32 = .annotateFusions(int, '3.2')
    int = rbind(int_51, int_31, int_52, int_32)
    if(!is.null(int)){ distinct(int) }
    return(int)
  }
  return(NULL)
}

#' Annotates fusions for specified field
#'
#' @param VAULT The vault object
#' @param FIELD The portion of the fusion to reformat
#' @return The annotated fusions for specified field
.annotateFusions <- function(INT, FIELD){
  ## set params
  if(FIELD=='5.1'){
    tmp = c('id', 'gene_ids.5.1', 'gene_names.5.1', 'gene_ids.3.1', 'gene_names.3.1', 'CSQF5')
    csq = 'CSQF5'
    other = '3.1'
  } else if(FIELD=='3.1'){
    tmp = c('id', 'gene_ids.3.1', 'gene_names.3.1', 'gene_ids.5.1', 'gene_names.5.1', 'CSQF3')
    csq = 'CSQF3'
    other = '5.1'
  } else if(FIELD=='5.2') {
    tmp = c('id', 'gene_ids.5.2', 'gene_names.5.2', 'gene_ids.3.2', 'gene_names.3.2', 'CSQR5')
    csq = 'CSQR5'
    other = '3.2'
  } else if(FIELD=='3.2'){
    tmp = c('id', 'gene_ids.3.2', 'gene_names.3.2', 'gene_ids.5.2', 'gene_names.5.2', 'CSQR3')
    csq = 'CSQR3'
    other = '5.2'
  } else { print('incorrect field chosen'); break }

  ## select and modify
  int = INT[ , ..tmp]
  int = int[!is.na(int[[csq]]), ]
  if(length(int$id)==0){ return(NULL) }

  int$tmp <- ifelse(!int[[csq]] %in% c('fusion3', 'fusion5'), 'structural', int[[csq]])
  int[['TYPE_PARTNER_GENE']] = paste0(int[['tmp']], '-', int[[paste0('gene_ids.', other)]])
  int[['TYPE_PARTNER_SYMBOL']] = paste0(int[['tmp']], '-', int[[paste0('gene_names.', other)]])

  ## reformat and return
  colnames(int)[c(2:3,6)] = c('gene_id', 'gene_name', 'annotation')
  int = int[, c('id', 'gene_id', 'gene_name', 'annotation', 'TYPE_PARTNER_GENE', 'TYPE_PARTNER_SYMBOL') ]
  int$gene_id = strsplit(int$gene_id, split = ':')
  int$gene_name = strsplit(int$gene_name, split = ':')
  int <- tidyr::unnest(int, cols = c('gene_id', 'gene_name'))
  return(int)
}

#### Genomic instability

#' Annotates measures of CIN
#'
#' @param VAULT The vault object
#' @param GOBJ the cnvex Gobj via getGobj()
#' @param OPTS cnvex opts via getOpts()
#' @return The segments
annotateCin <- function(VAULT, GOBJ, OPTS){
  if(!is.null( VAULT[['tables']][['segment']] )) {
    if(!all(is.na( VAULT[['tables']][['segment']][['C']] ))) {
      ## prepare
      tmp = data.table(id = VAULT[['meta']][['id']], wgii = NA, lst = NA, nloh = NA, ntai = NA)
      seg = VAULT[['tables']][['segment']]
      ## calculate
      tmp[['wgii']] = cinWgii(seg)
      tmp[['lst']] = cinLST(seg, GOBJ, OPTS)
      tmp[['nloh']] = cinNLOH(seg, GOBJ, OPTS)
      tmp[['ntai']] = cinNtAI(seg, GOBJ, OPTS)
      tmp[['fg']] = cinFG(seg, VAULT$meta$ploidy)
      tmp[['si']] = cinSI(seg, GOBJ, OPTS)
      ##
      tmp[['si']] = cinSI(seg, GOBJ, OPTS)
      tmp[['si']] = cinSI(seg, GOBJ, OPTS)
      return(tmp)
    }
  }
  return(NULL)
}




#### Mutations

#' Annotates the mutations for provided type
#'
#' @param VAULT The vault object
#' @param TYPE somatic/germline
#' @return The mutations
annotateMutations <- function(VAULT, TYPE){
  if(!is.null(VAULT[['tables']][[TYPE]])) {
    ## annotate
    tmp = VAULT[['tables']][[TYPE]]
    tmp[['annotation']] = ifelse(grepl('frameshift|stop_gained|stop_lost|start_lost|splice', tmp[['Consequence']]), 'truncating', tmp[['Consequence']])
    tmp[['annotation']] = ifelse(grepl('missense|inframe|protein_altering_variant', tmp[['annotation']]), 'missense', tmp[['annotation']])
    tmp[['annotation']] = ifelse(grepl('synonymous|stop_retained_variant', tmp[['annotation']]), 'synonymous', tmp[['annotation']])
    tmp[['annotation']] = ifelse(grepl('UTR|stream|intergenic', tmp[['annotation']]), 'noncoding', tmp[['annotation']])
    ## AA pos
    tmp[['AA_pos']] = as.numeric(gsub('^p\\.[A-Z]{1,6}([0-9]{1,4}).*', '\\1', tmp[['HGVSp']], ignore.case = TRUE))
    ## reformat
    tmp = distinct(tmp[, c('id', 'gene_id', 'gene_name', 'annotation', 'AA_pos')])
    return(tmp)
  }
  return(NULL)
}

#' Annotates the structural mutations for provided type
#'
#' @param VAULT The vault object
#' @param TYPE structural
#' @return The mutations
#'
#' TODO: this currently oversimplifies and misses variants in different alleles
#'
annotateMutationsStructural <- function(VAULT){
  tmp = VAULT[['tables']][['structural']]
  ## gene 1
  g1 = distinct(tmp[!is.na(tmp$gene_id.1), c('id', 'gene_id.1', 'gene_name.1', 'CSQ1')])
  colnames(g1) <- c('id', 'gene_id', 'gene_name', 'CSQ')
  ## gene 2
  g2 = distinct(tmp[!is.na(tmp$GENE2), c('id', 'gene_id.2', 'gene_name.2', 'CSQ2')])
  colnames(g2) <- c('id', 'gene_id', 'gene_name', 'CSQ')
  ## annotate
  g <- rbind(g1, g2)
  g$tmp <- ifelse(!g$CSQ %in% c('fusion3', 'fusion5'), 'structural', g$CSQ)
  g$TYPE_PARTNER_GENE = paste(g$tmp, g$gene_id, sep = '-')
  g$TYPE_PARTNER_SYMBOL = paste(g$tmp, g$gene_name, sep = '-')
  g$annotation = g$CSQ
  g=g[,c('id','gene_id','gene_name','annotation','TYPE_PARTNER_GENE','TYPE_PARTNER_SYMBOL')]
  g <- distinct(g)
  return(g)
}




#### Gene expression

#' Annotates measures of CIN
#'
#' @param VAULT The vault object
#' @param GOBJ the cnvex Gobj via getGobj()
#' @param OPTS cnvex opts via getOpts()
#' @return The segments
annotateGeneExpression <- function(VAULT, GOBJ, OPTS){
  exp = VAULT[['tables']][['gene.expression']]
  colnames(exp)[2] = c('gene_id')
  return(exp)
}




#### Chr arm gains/losses

#' Annotates measures of CIN
#'
#' @param ARNS The output of annotateArms()
#' @return Annotated chr gains/losses
#' @export
#'
annotateChromosomes <- function(ARMS){
  ## keep
  tmp <- ARMS %>%
    group_by(id, chr) %>%
    mutate(start = min(start, na.rm = T)) %>%
    mutate(end = max(end, na.rm = T)) %>%
    distinct(id, chr, start, end, wm, ploidy)
  ## whole
  chr <- ARMS %>%
    mutate(int = ifelse(annotation %in% c('gain'), 1, 0)) %>%
    mutate(int = ifelse(annotation %in% c('loss'), -1, int)) %>%
    group_by(id, chr) %>%
    mutate(sum = sum(int, na.rm = T)) %>%
    mutate(annotation = ifelse(sum==2, 'gain', 'none')) %>%
    mutate(annotation = ifelse(sum== -2, 'loss', annotation)) %>%
    distinct(id, chr, annotation) %>%
    merge(., tmp, by=c('id','chr')) 
  return(chr)
}


#' Annotates measures of CIN
#'
#' @param VAULT The vault object
#' @return Annotated chr arm gains/losses
#' @export
#'
annotateArms <- function(VAULT){

  if(!is.null( VAULT[['tables']][['segment']] )) {
    if(!all(is.na( VAULT[['tables']][['segment']][['C']] ))) {

      # get
      pos = getArms(VAULT)
      pos = pos[!grepl('X|Y', pos$arm),]
      wm <- getWeightedMedian(VAULT)
      ploidy <- VAULT$meta$ploidy
      id = VAULT$meta$id

      # annotate
      pos$id <- id
      arms <- as.list(pos$arm)
      labels <- lapply(arms, .annotateArms, VAULT=VAULT, POS=pos, WM=wm)
      labels = as.data.table(do.call(rbind, labels))
      colnames(labels) = c('annotation','proportion')

      ## reformat
      labels$proportion <- as.numeric(labels$proportion)
      labels$wm <- wm
      labels$ploidy <- ploidy
      labpos = cbind(pos,labels)
      labpos <- labpos[,c('id','arm', 'annotation', 'chr', 'start', 'end', 'wm', 'ploidy', 'proportion')]
      return(labpos)
    }
  }
  return(NULL)
}

#' Given the arm returns, compresses region, overlaps with segments and returns
#' the annotation.
#'
#' @param ARM the arm you're annotating
#' @param DIGEST The digest file
#' @param ARMS The coordinates for the arms
#' @param WM The sample's weighted median copy number
#' @return The annotated unified object.
.annotateArms <- function(ARM, VAULT, POS, WM) {

  # current info
  cur_arm <- GRanges(POS[POS$arm==ARM,])
  seg_meta <- GRanges(VAULT$tables$segment)
  ploidy = VAULT$meta$ploidy

  # overlap
  ol <- to(findOverlaps(cur_arm, seg_meta, ignore.strand = T, type = 'any'))
  cvals <- as.data.table(seg_meta[ol, 'C'])$C

  # if there are non-NA values, assign, else set to NA
  if(length(ol) > 0 & sum(!is.na(cvals)) > 0){

    ## mode of segments
    c_mode = getMode(VAULT)

    ## get segments of chr arm
    ol <- as.data.table(suppressWarnings(pintersect(findOverlapPairs(seg_meta, cur_arm))))
    ol <- ol[,c('seqnames','start','end','C')]
    ol <- ol[order(start),]
    ol$C <- ifelse(is.na(ol$C), -1, ol$C)

    ## annotate
    ol$annotation = ifelse(ol$C>WM & ol$C>c_mode & ol$C>(ploidy+0.5) & ol$C>2, 'gain', 'none')
    ol$annotation = ifelse(ol$C<WM & ol$C<c_mode & ol$C<(ploidy-0.5) & ol$C<4, 'loss', ol$annotation)

    ## compress
    ol = compressOverlapSegments(OVERLAP=ol)
    label = labelOverlap(ol)
    return(label)
  }
  return(NA)
}

#' Given overlap returns the compressed overlap
#'
#' @param OVERLAP dataframe of overlap of segments with arm
#' @return The compressed overlap df
#'
compressOverlapSegments <- function(OVERLAP){

  ## initialize
  i <- 1
  new_ol <- NULL

  # while position is w/in overlap
  while (i <= dim(OVERLAP)[1]) {
    pos = 0
    end <- OVERLAP$end[i]
    cur_c <- OVERLAP$C[i]
    cur_annotation <- OVERLAP$annotation[i]
    ## combine segments w/ same annotation if any
    same_annotation <- (OVERLAP$annotation == cur_annotation)[-(1:i)]
    if(length(same_annotation) != 0){
      same_annotation <- append(same_annotation, F)
      ## before first False position
      pos <- Position(function(x) x == F, same_annotation)-1
      end <- OVERLAP$end[(pos+i)]
      cur_c = paste(OVERLAP$C[i:(pos+i)], collapse = ';')
    }

    # add to new overlap list
    new_ol <- rbind(new_ol, data.table(seqnames = OVERLAP$seqnames[i], start = OVERLAP$start[i], end = end, C = cur_c, annotation = cur_annotation))
    # reindex
    i <- i+pos+1
  }
  return(new_ol)
}

#' Given compressed overlap returns the C, label, and proportion of arm
#'
#' @param OVERLAP dataframe of overlap of segments with arm
#' @return df containing the C values, label, and proportino of arm
#'
labelOverlap <- function(OVERLAP){
  OVERLAP$proplen <- ( (OVERLAP$end-OVERLAP$start) / sum(OVERLAP$end-OVERLAP$start) )
  OVERLAP = OVERLAP[,sum:=sum(proplen),by="annotation"]
  OVERLAP = distinct(OVERLAP[,c('annotation','sum')])
  OVERLAP = OVERLAP[order(-sum)]
  annotation = unlist(OVERLAP[1,1],use.names = F)
  proportion = round(unlist(OVERLAP[1,2],use.names = F),2)
  if(proportion<0.5){ annotation = 'none'}
  return( c(annotation, proportion) )
}




#### Table hits

#' Counts hits/gene for a single vault
#'
#' @param MVF A vault/metavault
#' @param HIT_FIELDS the tables to count towards hits
#' @param RULES variant filtering rules
#' @return DF of hits/gene
#'
#' @export
tableHitsMetavault <- function(MVF, HIT_FIELDS=c('somatic','germline','structural','gene.copy'), RULES=NULL, GENES=NULL){

  ## remove failed variants
  mvf = vaultFilter(MVF, RULES)

  ## get tables
  som = getVaultTable(mvf, 'somatic') |> as.data.table()
  grm = getVaultTable(mvf, 'germline') |> as.data.table()
  str = getVaultTable(mvf, 'structural') |> as.data.table()
  fus = getVaultTable(mvf, 'fusion') |> as.data.table()
  gc = getVaultTable(mvf, 'gene.copy') |> as.data.table()

  ## get genes/ids; break if none
  genes <- getGenes(SOM = som,
                    GRM = grm,
                    STR = str,
                    FUS = fus,
                    GENE.COPY = gc)
  if(!is.null(GENES)) { genes <- genes[V2 %in% GENES, ] }
  ids = getVaultId(mvf)
  if(nrow(genes)==0){ return(NULL) }
  if(is.null(ids)){ print('No IDs in meta');break }

  # make hits table
  hits = makeHitsTable(ids, genes)

  ## tally
  hits[['gene.copy']] = .tallyCN(hits, gc)
  hits[['somatic']] = .tallyMutations(hits, som)
  hits[['germline']] = .tallyMutations(hits, grm)
  hits[['structural']] = .tallyMutations(hits, str)

  ## fusions
  tmp = .getMutations(fus)
  hits[['fusion']] = NA
  if(!is.null(fus)){ hits[['fusion']] = ifelse(hits[['id.gene']] %in% tmp, 1, 0) }

  ## sum hits
  hits[['hits']] = rowSums(hits[,c(..HIT_FIELDS)],na.rm = T)
  hits[['id.gene']] = NULL
  return(hits)

}


#' Makes a table of hits
#'
#' @param IDS vector of ids
#' @param GENES DF of GENES and SYMBOLS
#' @return DF of hits
#'
makeHitsTable <- function(IDS,GENES){
  hits = data.table('id' = list(IDS), gene_id = GENES$V1, gene_name = GENES$V2)
  hits = tidyr::unnest(hits, cols = 'id')
  hits[['id.gene']] = paste0(hits[['id']], '--', hits[['gene_id']])
  hits = as.data.table(hits)
  return(hits)
}

#' tallies copy number hits
#'
#' @param HITS the hits table
#' @param GC gene copy numebr table
#' @return vector of cn losses
#'
.tallyCN <- function(HITS, GC){
  HITS[['tmp']] = NA
  if(!is.null(GC)){
    GC[['id.gene']] = paste0(GC[['id']], '--', GC[['gene_id']])
    loss1 = GC[Kmin==0,]
    loss2 = GC[Cmin==0,]
    HITS[['tmp']] = ifelse(HITS[['id']] %in% GC[['id']], 0, NA)
    HITS[['tmp']] = ifelse(HITS[['id.gene']] %in% loss1[['id.gene']], 1, HITS[['tmp']])
    HITS[['tmp']] = ifelse(HITS[['id.gene']] %in% loss2[['id.gene']], 2, HITS[['tmp']])
  }
  return(HITS[['tmp']])
}

#' tallies mutations
#'
#' @param HITS the hits table
#' @param TBL the mutation type
#' @return table
#'
.tallyMutations <- function(HITS, TBL){
  if(!is.null(TBL)){
    muts = .getMutations(TBL)
    cnt = as.data.table(table( muts ))
    cnt=cnt[N>1]
    tally = ifelse(HITS[['id']] %in% TBL[['id']], 0, NA)
    tally = ifelse(HITS[['id.gene']] %in% muts, 1, tally)
    tally = ifelse(HITS[['id.gene']] %in% cnt[['muts']], 2, tally)
    return(tally)
  }
  return(NA)
}

#' gets id/GENE key
#'
#' @param TABLE which table to get
#' @return vector of id--gene values
#'
.getMutations <- function(TBL){
  ## som/germ
  if( 'gene_id' %in% colnames(TBL) ){
    return( paste0(TBL[['id']], '--', TBL[['gene_id']]) )
  }
  ## str
  if( 'gene_id.1' %in% colnames(TBL) ){
    g1 = paste0(TBL[['id']], '--', TBL[['gene_id.1']])
    g2 = paste0(TBL[['id']], '--', TBL[['gene_id.2']])
    g = c(g1,g2)[grepl('--ENSG', c(g1,g2))]
    return(g)
  }
  ## fus
  if( all(c('gene_names.5.1', 'gene_names.3.1', 'gene_names.5.2', 'gene_names.3.2') %in% colnames(TBL) ) ){
    ## get
    g51 = TBL[,c('id', 'gene_ids.5.1')] %>% set_colnames(c('id','gene_id'))
    g31 = TBL[,c('id', 'gene_ids.3.1')] %>% set_colnames(c('id','gene_id'))
    g52 = TBL[,c('id', 'gene_ids.5.2')] %>% set_colnames(c('id','gene_id'))
    g32 = TBL[,c('id', 'gene_ids.3.2')] %>% set_colnames(c('id','gene_id'))
    g = rbindlist(list(g51,g52,g31,g32))
    ## unnest
    g$gene_id = strsplit(g$gene_id, split = ':')
    g <- tidyr::unnest(g, cols = c('gene_id'))
    ## return
    return( paste0(g[['id']], '--', g[['gene_id']]) )
  }
}




#### Misc

#' Sets metadata field
#'
#' @param VAULT The vault object
#' @param FIELD The metadata field to change
#' @param VALUE The value to set the metadata to
#'
#'  @export
setMetaField <- function(VAULT, FIELD, VALUE){
  VAULT[['meta']][[FIELD]] = VALUE
  return(VAULT)
}

#' get all available genes for a sample
#'
#' @param SOM annotated somatic
#' @param GRM annotated germline
#' @param STR annotated structural
#' @param FUS annotated fusions
#' @param GENE.COPY annotated gene copy number
#' @return DF ensembl/symbol
#'
getGenes <- function(SOM, GRM, STR, FUS, GENE.COPY){
  genes <- list(data.table(SOM$gene_id, SOM$gene_name),
                data.table(GRM$gene_id, GRM$gene_name),
                data.table(STR$gene_id, STR$gene_name),
                data.table(FUS$gene_id, FUS$gene_name),
                data.table(GENE.COPY$gene_id, GENE.COPY$gene_name))
  genes <- genes[ unlist(lapply(genes, function(X){ dim(X)[1]>0 })) ]
  genes <- rbindlist(genes)
  genes <- distinct(genes)
  return(genes)
}

#' gets id(s) from vault/metavault/annotated metavault
#'
#' @param V The vault-like object
#' @return vector of ids
#'
#'  @export
getVaultId <- function(V){
  ids <- V$meta |> pull(id)
  return(ids)
}

#' gets the chr arms for the vault's assembly
#'
#' @param VAULT The vault object
#' @return df of chr arms
#'
#'  @export
getArms <- function(VAULT){
  assembly = unique(VAULT$meta$assembly)
  if(assembly=='hg38'){
    arms <- read.table(system.file("extdata/hg38/arm.bed.gz", package="cnvex"))
    colnames(arms) <- c('chr','start','end','arm','score','strand')
    return(arms)
  }
  print(paste0("Invalid assembly: ", assembly))
  return(NULL)
}

#' gets the mode
#'
#' @param VAULT The vault object
#' @return df of chr arms
#'
getMode <- function(VAULT){
  c_mode <- modeest::mfv(VAULT$tables$segment$C)
  c_mode <- mean(c_mode, na.rm = T)
  return(c_mode)
}
