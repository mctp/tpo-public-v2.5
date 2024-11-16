
#### CN heatmap

#' Plots CN heatmaps for samples w/ segment data
#'
#' @param MV The metavault
#' @param GOBJ the Gobj
#' @param BIN width of bin to use on genome; default 1MB
#' @param TYPE the type of CN heatmap; c('Annotated', 'RelativeC', 'ActualC')
#' @param META optional DF of meta to order/split. 2 columns ('order'/'split') with ids as rownames.
#' @return plot
#' @export
#' 
plotCnHeatmap <- function(MV, GOBJ, BIN=1e6, TYPE, META=NULL){
  
  ## prepare
  annot = prepareCnHeatmap(MV,GOBJ,BIN)
  
  ## plot
  if(TYPE=='Annotated'){
    p = .plotCnHeatmapAnnotatedC(ANNOTATED = annot, META)
    return(p)
  }
  if(TYPE=='RelativeC'){
    p = .plotCnHeatmapRelativeC(ANNOTATED = annot, META)
    return(p)
  }
  if(TYPE=='ActualC'){
    p = .plotCnHeatmapActualC(ANNOTATED = annot, META)
    return(p)
  }
  if(TYPE=='BinaryK'){
    p = .plotCnHeatmapBinaryK(ANNOTATED = annot, META)
    return(p)
  }
  
  ## error catch
  print('Invalid type. Please select from: "Annotated", "RelativeC", "ActualC".')
  return(NULL)
  
}

#' prepares data for plotting CN heatmaps
#'
#' @param MV The metavault
#' @param GOBJ the Gobj
#' @param BIN width of bin to use on genome; default 1MB
#' @return Df of prepare data
#' 
prepareCnHeatmap <- function(MV, GOBJ, BIN=1e6){
  
  # get starting data
  tiles <- getTiles(GOBJ, BIN)
  segments = getSegments(MV)
  arms = GRanges( getArms(MV) )
  
  ## overlap tiles, segments, and arms
  overlap_tiles = getOverlapWidth(QUERY=tiles, SUBJECT=segments)
  overlap_arms = getOverlapWidth(QUERY=tiles, SUBJECT=arms)[,1:2]
  
  ## merge segments/tiles together
  segments_dt = as.data.table(segments)[,c('id','C','K')]
  segments_dt$subjectHits <- 1:nrow(segments_dt)
  tiles_ol = merge(tiles, overlap_tiles, by='queryHits', all=T)
  tiles_ol = merge(tiles_ol, segments_dt, by='subjectHits', all.x = T)
  
  ## merge arms
  arms_dt = as.data.table(arms)[,c('arm')]
  arms_dt$armHits <- 1:nrow(arms_dt)
  colnames(overlap_arms) <- c('queryHits', 'armHits')
  tiles_ola = merge(tiles_ol, overlap_arms, by = 'queryHits', all = T)
  tiles_ola = merge(tiles_ola, arms_dt, by='armHits',all.x = T)
  tiles_ola = data.table(tiles_ola[,!(names(tiles_ola) %in% 'armHits')])
  
  ## handle uncovered regions 
  uncov = tiles_ola[is.na(id),]
  tmp = list(unique(segments_dt$id))
  uncov$id <- rep(tmp, length(uncov$id))
  uncov <- tidyr::unnest(uncov, cols = id)
  uncov$id <- as.character(uncov$id)
  tiles_olu = rbind( tiles_ola[!is.na(id)], uncov )
  
  
  ## handle multiple matches
  tiles_olud <- tiles_olu %>%
    mutate(weightedC:= C*overlap_width) %>%
    group_by(id, queryHits) %>%
    mutate(sumWidth = sum(overlap_width, na.rm = T)) %>%
    mutate(sumWeightedC = sum(weightedC, na.rm = T)) %>%
    mutate(tileC = sumWeightedC/sumWidth) %>%
    mutate(maxC = ifelse( all(is.na(C)), NA, max(C, na.rm = T) )) %>%
    mutate(minC = ifelse( all(is.na(C)), NA, min(C, na.rm = T) )) %>%
    mutate(K = ifelse( all(is.na(K)), NA, min(K, na.rm = T) )) %>%
    ungroup() %>%
    select(queryHits, seqnames, start, end, tileC, maxC, minC, id, arm, K) %>%
    distinct() %>%
    ## remove dups d/t arms
    arrange(arm) %>%
    filter(!duplicated(paste0(id, '..', queryHits))) %>%
    ## order
    arrange(queryHits, id) %>%
    as.data.table()
  
  ## add weighted medians, ploidies, and relative C values
  ids = MV[['meta']][['id']]
  vaults = lapply(ids, extractVault, MV)
  wms = lapply(vaults, getWeightedMedian)
  wms = data.table('id'=ids,'wm'=as.numeric(as.character(wms)))
  pwm = merge(tiles_olud, wms, by='id', all.x = T)
  cmodes = lapply(vaults, getMode)
  cmodes = data.table('id'=ids,'c_mode'=as.numeric(as.character(cmodes)))
  pwmc = merge(pwm, cmodes, by='id')
  pwmcm = merge(pwmc, MV$meta[,c('id','ploidy')], by='id', all.x = T)
  pwmcm[, relC:= tileC/ploidy]
  sex = getSex(MV)
  pwmcms = merge(pwmcm, sex, by='id', all.x = T)
  
  ## annotate autosomes
  annot = pwmcms
  annot$annot = ifelse( (annot$maxC>annot$ploidy+3) & (annot$maxC>annot$wm+3) & (annot$maxC>8), 'amplification', 'none')
  annot$annot = ifelse( (annot$annot=='none') & (annot$maxC>annot$ploidy+0.5) & (annot$maxC>annot$wm) & (annot$maxC>2) & annot$maxC>annot$c_mode, 'gain', annot$annot)
  annot$annot = ifelse( (annot$annot=='none') & (annot$minC<annot$ploidy-0.5) & (annot$minC<annot$wm) & (annot$minC<4) & annot$minC<annot$c_mode, 'loss', annot$annot)
  annot$annot = ifelse( (annot$minC == 0), 'homo-del', annot$annot)
  annot$annot = ifelse(is.na(annot$tileC), 'missing', annot$annot)
  ## annotate X chromosome for men
  annot$annot <- ifelse(annot$annot %in% c('gain', 'loss') & annot$seqnames=='chrX' & annot$sex=='male', 'none', annot$annot)
  annot$annot = ifelse(
    annot$seqnames=='chrX' & annot$sex=='male' & (!annot$annot %in% c('amplification','homo-del','missing')) & 
      (annot$tileC>(annot$ploidy/2)+0.5) & (annot$tileC>annot$wm/2) & (annot$tileC>1) & annot$tileC>annot$c_mode/2,
    'gain',annot$annot)
  annot$annot = ifelse(
    annot$seqnames=='chrX' & annot$sex=='male' & (!annot$annot %in% c('amplification','homo-del','missing')) & 
      (annot$tileC<(annot$ploidy/2)-0.5) & (annot$tileC<annot$wm/2) & (annot$tileC<3) & annot$tileC<annot$c_mode/2,
    'loss',annot$annot)
  
  return(annot)
  
}

#' Plots Annotated heatmap (amplification, gain, ect.)
#'
#' @param ANNOTATED Output of prepareCnHeatmap
#' @param META optional DF of meta to order/split. 2 columns ('order'/'split') with ids as rownames.
#' @return Annotated plot
#' 
.plotCnHeatmapAnnotatedC <- function(ANNOTATED, META){
  
  ## format for plot
  tp = ANNOTATED[,c('id','queryHits','annot')]
  tp = tidyr::spread(tp, key = queryHits, value = annot)
  tp = as.matrix(tp)
  rownames(tp) <- tp[,1]
  tp = tp[,-1]
  
  ## Row split/order
  ROW_SPLIT <- rep('CNA plot: Annotated',nrow(tp))
  ROW_ORDER <- NULL
  if(!is.null(META)){
    rn = rownames(META)[rownames(META) %in% ANNOTATED$id]
    ## split
    if( 'split' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'split'] )
      colnames(tmp) <- 'split'
      tmp$id <- rn
      ROW_SPLIT <- tmp$split
    }
    ## order
    if( 'order' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'order'] )
      colnames(tmp) <- 'order'
      tmp$id <- rn
      ROW_ORDER <- tmp$order
    }
    ## match
    rows = match(tmp$id, rownames(tp))
    tp=tp[rows,]
  }
  
  ## Col order
  ORDER_COL = dplyr::distinct(ANNOTATED[,c('queryHits','arm')])
  ORDER_COL = ORDER_COL[order(queryHits),][['arm']]
  tmp = c(paste0('chr',rep(1:22,each=2), rep(c('p','q'), 22)), 'chrXp', 'chrXq')
  ORDER_COL = factor(ORDER_COL, levels = tmp)
  
  ## plot
  COLORS = structure(c('#000000','#0072B5FF','grey89','#BC3C29FF','#E18727FF','grey60'), names = c('homo-del','loss','none','gain','amplification','missing'))
  COLUMN_TITLES = c( paste0('chr',rep(1:22,each=2)), 'chrX', 'chrX')
  COLUMN_TITLES[seq_along(COLUMN_TITLES)%%2==0] = ''
  p <- Heatmap(tp,
               ## rows
               show_row_dend = F,
               show_row_names = F,
               row_split = ROW_SPLIT,
               row_order =  ROW_ORDER,
               cluster_rows = T,
               cluster_row_slices = F,
               row_title_rot = 90,
               row_gap = grid::unit(3, "mm"), 
               row_title_gp = grid::gpar(fontsize = 10), 
               ## cols
               show_column_dend = F, 
               show_column_names = F,
               cluster_columns = F,
               column_title = COLUMN_TITLES,
               column_split = ORDER_COL,
               column_title_gp = grid::gpar(fontsize = 10), 
               column_title_side = 'top',
               column_title_rot = 90,
               column_gap = unit(c(rep(c(0.5, 3),22),0.5), "mm"),
               ## other
               col = COLORS,
               na_col = '#000000',
               name = 'Annotation'
  );p
  return(p)
}

#' Plots relative C (C:ploidy)
#'
#' @param ANNOTATED Output of prepareCnHeatmap
#' @param META optional DF of meta to order/split. 2 columns ('order'/'split') with ids as rownames.
#' @return Relative C plot
#' 
.plotCnHeatmapRelativeC <- function(ANNOTATED, META){
  
  ## format for plot
  tp = ANNOTATED[,c('id','queryHits','relC')]
  MAX = ceiling(max(ANNOTATED$relC, na.rm = T))
  tp = tidyr::spread(tp, key = queryHits, value = relC)
  tp = as.matrix(tp)
  nms <- tp[,1]
  tp = tp[,-1]
  tp = matrix(as.numeric(tp), ncol = ncol(tp))
  rownames(tp) <- nms

  
  ## Row split/order
  ROW_SPLIT <- rep('CNA plot: Annotated',nrow(tp))
  ROW_ORDER <- NULL
  if(!is.null(META)){
    rn = rownames(META)[rownames(META) %in% ANNOTATED$id]
    ## split
    if( 'split' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'split'] )
      colnames(tmp) <- 'split'
      tmp$id <- rn
      ROW_SPLIT <- tmp$split
    }
    ## order
    if( 'order' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'order'] )
      colnames(tmp) <- 'order'
      tmp$id <- rn
      ROW_ORDER <- tmp$order
    }
    ## match
    rows = match(tmp$id, rownames(tp))
    tp=tp[rows,]
  }
  
  ## Col order
  ORDER_COL = dplyr::distinct(ANNOTATED[,c('queryHits','arm')])
  ORDER_COL = ORDER_COL[order(queryHits),][['arm']]
  tmp = c(paste0('chr',rep(1:22,each=2), rep(c('p','q'), 22)), 'chrXp', 'chrXq')
  ORDER_COL = factor(ORDER_COL, levels = tmp)
  
  ## plot
  MAX= ifelse(MAX <3, 3, MAX)
  COLORS = circlize::colorRamp2(c(0, 0.5, 1, 2, MAX), c('#022d46',"#0571B0", "white", "#CA0020", '#650010'))
  COLUMN_TITLES = c( paste0('chr',rep(1:22,each=2)), 'chrX', 'chrX')
  COLUMN_TITLES[seq_along(COLUMN_TITLES)%%2==0] = ''
  p <- Heatmap(tp,
               ## rows
               show_row_dend = F,
               show_row_names = F,
               row_split = ROW_SPLIT,
               row_order = ROW_ORDER,
               cluster_rows = T,
               cluster_row_slices = F,
               row_title_rot = 90,
               row_gap = grid::unit(3, "mm"),
               row_title_gp = grid::gpar(fontsize = 10), 
               ## cols
               show_column_dend = F, 
               show_column_names = F,
               cluster_columns = F,
               column_title = COLUMN_TITLES,
               column_split = ORDER_COL,
               column_title_gp = grid::gpar(fontsize = 10), 
               column_title_side = 'top',
               column_title_rot = 90,
               column_gap = unit(c(rep(c(0.5, 3),22),0.5), "mm"),
               ## other
               col = COLORS,
               name = 'Relative C'
  ); return(p)
}


#' Plots actual C (unadjusted C)
#'
#' @param ANNOTATED Output of prepareCnHeatmap
#' @param META optional DF of meta to order/split. 2 columns ('order'/'split') with ids as rownames.
#' @return Actual C plot
#' 
.plotCnHeatmapActualC <- function(ANNOTATED, META){
  
  ## format for plot
  tp = ANNOTATED[,c('id','queryHits','tileC')]
  MAX = ceiling(max(ANNOTATED$tileC, na.rm = T))
  tp = tidyr::spread(tp, key = queryHits, value = tileC)
  tp = as.matrix(tp)
  nms <- tp[,1]
  tp = tp[,-1]
  tp = matrix(as.numeric(tp), ncol = ncol(tp))
  rownames(tp) <- nms
  
  
  ## Row split/order
  ROW_SPLIT <- rep('CNA plot: Actual',nrow(tp))
  ROW_ORDER <- NULL
  if(!is.null(META)){
    rn = rownames(META)[rownames(META) %in% ANNOTATED$id]
    ## split
    if( 'split' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'split'] )
      colnames(tmp) <- 'split'
      tmp$id <- rn
      ROW_SPLIT <- tmp$split
    }
    ## order
    if( 'order' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'order'] )
      colnames(tmp) <- 'order'
      tmp$id <- rn
      ROW_ORDER <- tmp$order
    }
    ## match
    rows = match(tmp$id, rownames(tp))
    tp=tp[rows,]
  }
  
  ## Col order
  ORDER_COL = dplyr::distinct(ANNOTATED[,c('queryHits','arm')])
  ORDER_COL = ORDER_COL[order(queryHits),][['arm']]
  tmp = c(paste0('chr',rep(1:22,each=2), rep(c('p','q'), 22)), 'chrXp', 'chrXq')
  ORDER_COL = factor(ORDER_COL, levels = tmp)
  
  ## plot
  COLORS = circlize::colorRamp2(c(0, 1, MAX), c("#0571B0", "white", "#CA0020"))
  COLUMN_TITLES = c( paste0('chr',rep(1:22,each=2)), 'chrX', 'chrX')
  COLUMN_TITLES[seq_along(COLUMN_TITLES)%%2==0] = ''
  p <- Heatmap(tp,
               ## rows
               show_row_dend = F,
               show_row_names = F,
               row_split = ROW_SPLIT,
               row_order = ROW_ORDER,
               cluster_rows = T,
               cluster_row_slices = F,
               row_title_rot = 90,
               row_gap = grid::unit(3, "mm"),
               row_title_gp = grid::gpar(fontsize = 10), 
               ## cols
               show_column_dend = F, 
               show_column_names = F,
               cluster_columns = F,
               column_title = COLUMN_TITLES,
               column_split = ORDER_COL,
               column_title_gp = grid::gpar(fontsize = 10), 
               column_title_side = 'top',
               column_title_rot = 90,
               column_gap = unit(c(rep(c(0.5, 3),22),0.5), "mm"),
               ## other
               col = COLORS,
               name = 'Actual C'
  ); return(p)
}


#' Plots binary K (LOH)
#'
#' @param ANNOTATED Output of prepareCnHeatmap
#' @param META optional DF of meta to order/split. 2 columns ('order'/'split') with ids as rownames.
#' @return Binary K plot
#' 
.plotCnHeatmapBinaryK <- function(ANNOTATED, META){
  
  ## format for plot
  tp <- ANNOTATED %>%
    mutate(LOH = ifelse(K>0, 'Heterozygous', 'LOH')) %>%
    distinct(id,queryHits,LOH) %>%
    tidyr::spread(., key = queryHits, value = LOH) %>%
    tibble::column_to_rownames('id') %>%
    as.matrix()
  #MAX = ceiling(max(ANNOTATED$tileC, na.rm = T))
  #nms <- tp[,1]
  #tp = tp[,-1]
  #tp = matrix(as.numeric(tp), ncol = ncol(tp))
  #rownames(tp) <- nms
  
  
  ## Row split/order
  ROW_SPLIT <- rep('LOH plot: ',nrow(tp))
  ROW_ORDER <- NULL
  if(!is.null(META)){
    rn = rownames(META)[rownames(META) %in% ANNOTATED$id]
    ## split
    if( 'split' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'split'] )
      colnames(tmp) <- 'split'
      tmp$id <- rn
      ROW_SPLIT <- tmp$split
    }
    ## order
    if( 'order' %in% colnames(META) ){
      tmp = as.data.frame( META[rownames(META) %in% ANNOTATED$id, 'order'] )
      colnames(tmp) <- 'order'
      tmp$id <- rn
      ROW_ORDER <- tmp$order
    }
    ## match
    rows = match(tmp$id, rownames(tp))
    tp=tp[rows,]
  }
  
  ## Col order
  ORDER_COL = dplyr::distinct(ANNOTATED[,c('queryHits','arm')])
  ORDER_COL = ORDER_COL[order(queryHits),][['arm']]
  tmp = c(paste0('chr',rep(1:22,each=2), rep(c('p','q'), 22)), 'chrXp', 'chrXq')
  ORDER_COL = factor(ORDER_COL, levels = tmp)
  
  ## plot
  COLORS = c('Heterozygous'="grey89", 'LOH'="black")
  COLUMN_TITLES = c( paste0('chr',rep(1:22,each=2)), 'chrX', 'chrX')
  COLUMN_TITLES[seq_along(COLUMN_TITLES)%%2==0] = ''
  p <- Heatmap(tp,
               ## rows
               show_row_dend = F,
               show_row_names = F,
               row_split = ROW_SPLIT,
               row_order = ROW_ORDER,
               cluster_rows = T,
               cluster_row_slices = F,
               row_title_rot = 90,
               row_gap = grid::unit(3, "mm"),
               row_title_gp = grid::gpar(fontsize = 10), 
               ## cols
               show_column_dend = F, 
               show_column_names = F,
               cluster_columns = F,
               column_title = COLUMN_TITLES,
               column_split = ORDER_COL,
               column_title_gp = grid::gpar(fontsize = 10), 
               column_title_side = 'top',
               column_title_rot = 90,
               column_gap = unit(c(rep(c(0.5, 3),22),0.5), "mm"),
               ## other
               col = COLORS, 
               na_col = 'grey60',
               name = 'Binary K'
  ); return(p)
}




#### Lollipop
#' Plots Lolliplots from metavaults
#'
#'
#'
#' @param MV The metavault or .maf file. If using .maf file then change FORMAT to "MAF"
#' @param SYMBOL Gene symbol
#' @param TRANSCRIPT Which trascript to use from the GTF. If not specified, it will use the most common transcript found for the gene.
#' @param PROTEIN indicating whether to plot protein or gene info. Choices = c("protein", "gene")
#' @param GTF GTF table containing exon information
#' @param STACK_SAME How to deal with mutations at same location?; c("stack", "unique", "none")
#' @param REMOVE_EXCESS Remove exons at either end that contain no information?
#' @param REMOVE_INTRON Remove intronic regions so exons are easier to see? This will cause issues if there are expected mutations in intron regions!
#' @param INTRON_GAP if shrinking introns, how much gap to leave between exons?
#' @param ZOOM If you want to zoom into a specific region to see; c("chrom", "chromstart", "chromend"), else = NULL
#' @param COLOR If you want to add colors to the exons; c("None", "Stripes", "Rainbow")
#' @param LABEL_ON_FEATURE puts exon numbers on the bars instead. But will remove x axis scale.
#' @param BREAKPOINT T/F indicating whether to plot breakpoints based on fusion information. Only works when REMOVE_INTRON is FALSE!
#' @param PROTEIN_DOMAINS Add file to plot protein domains on track. Default = F, need to supply protein file as gff3 format.
#' @param SEPARATE allows separation of lollipops by up or down if provided a column name. Must be only 2 categories.
#' @param FORMAT pick input format as metavault "MVF" or MAF format "MAF"
#' @param PATIENTS a list of patients to filter by
#' @return plot
#' 
#' @export
plotLollipop <- function(MV, SYMBOL, GTF, PROTEIN_DOMAINS=F, FORMAT="MVF",
                         PATIENTS = NULL,
                         STACK_SAME="stack", REMOVE_INTRON = T, FIELD='somatic', PROTEIN="protein", 
                         REMOVE_EXCESS=F, BREAKPOINT=F, COLOR="rainbow", SHOW_EXON_NAME=F,
                         INTRON_GAP = NULL, LABEL_ON_FEATURE=F, SEPARATE=NULL, TRANSCRIPT="common"){
  
  # Initialize
  ###############
  ## params
  proteinlib <- c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
                  'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N', 
                  'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W', 
                  'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M')
  ## subset
  if(FORMAT == "MVF"){
    tbl = getVaultTable(MV, FIELD) %>% dplyr::filter(SYMBOL==!!SYMBOL & pass == "pass")
  } else if (FORMAT == "MAF"){
    #assume all are filtered and pass
    tbl = MV %>% dplyr::filter(Hugo_Symbol==!!SYMBOL)
    names(tbl)[names(tbl) == 'Transcript_ID'] <- 'TRANSCRIPT'
    names(tbl)[names(tbl) == 'Patient_ID'] <- 'id'
    names(tbl)[names(tbl) == 'Chromosome'] <- 'chr'
    names(tbl)[names(tbl) == 'Start_Position'] <- 'pos'
    
  }
  
  if(!is.null(PATIENTS)){
    tbl <- tbl %>% dplyr::filter(id %in% PATIENTS)
  }
  ## get most common transcript
  if(TRANSCRIPT=="common"){
    tt <- data.frame(table(tbl[['TRANSCRIPT']]))
    tt <- tt[['Var1']][tt[['Freq']]==max(tt[['Freq']])]
  } else { tt <- TRANSCRIPT }
  ## subset GTF file for blocks to be plotted
  ft <- GTF %>% dplyr::filter(geneName==SYMBOL & type=='CDS' & transcript==tt)
  ftt <- GRanges( paste0("chr",ft[['contig']]), IRanges(ft[['start']], ft[['end']], names=paste0( ft[['exonNumber']])))
  ftt$featureLayerID <- "gene"
  tbl[['score']] = 1
  ###############
  
  ## Change Protein AA code
  if(PROTEIN!="protein"){
    gp <- "var_id"
  } else {
    gp <- "HGVSp"
    tbl <- tbl[!is.na(tbl$HGVSp),]
    tbl[['HGVSp']] <- gsub("p[.]", "", tbl[['HGVSp']])
    for(i in 1:nrow(tbl)){
      for(p in names(proteinlib)){
        if( grepl(p, tbl[['HGVSp']][i], ignore.case = T) ){
          tbl[['HGVSp']][i] <- gsub(p, proteinlib[p], tbl[['HGVSp']][i])
        }
      }
    }
  }
  
  ## Check if there is a separation parameter
  if(!is.null(SEPARATE)){
    tbl <- tbl[!is.na(tbl[[SEPARATE]])]
    if(length(unique(tbl[[SEPARATE]])) != 2){
      warning("More than 2 categories in SEPARATE parameter")
    }
  }
  
  ## determine height (stack) of each location
  if(STACK_SAME == "stack"){
    tbl2 <- tbl %>% dplyr::count(score, chr, pos, !!!syms(gp), !!!syms(SEPARATE))
  } else if (STACK_SAME == "unique") {
    tbl2 <- tbl %>% dplyr::transmute(chr, pos, !!!syms(gp), !!!syms(SEPARATE), n=1) %>% dplyr::distinct()
  } else {
    tbl2 <- tbl %>% dplyr::transmute(chr, pos, !!!syms(gp), !!!syms(SEPARATE), n=1)
  }
  ## feature
  snpside <- ifelse(rep(!is.null(SEPARATE), nrow(tbl2)),
                    c("top", "bottom")[1+as.numeric(tbl2[[SEPARATE]] == unique(tbl2[[SEPARATE]])[1])],
                    "top")
  feature <- GRanges(tbl2[['chr']], IRanges(tbl2[['pos']], width = 1, names = tbl2[[gp]]), score=tbl2[['n']], 
                     SNPsideID = snpside)
  
  ###############
  
  ## load protein domains that overlap w/ feature track
  ###############
  if(!is.null(PROTEIN_DOMAINS)){
    #make sure each column is there
    if(F %in% (c("contig", "start", "end") %in% PROTEIN_DOMAINS)){
      PROTEIN_DOMAINS.gr <- GRanges(paste0("chr", PROTEIN_DOMAINS[['contig']]), 
                                    IRanges(PROTEIN_DOMAINS[['start']], PROTEIN_DOMAINS[['end']], names=PROTEIN_DOMAINS[['proteinDomainName']]))
      PROTEIN_DOMAINS.gr <- subsetByOverlaps(PROTEIN_DOMAINS.gr, ftt)
      if(length(PROTEIN_DOMAINS.gr) > 0){
        PROTEIN_DOMAINS.gr$'featureLayerID' <- "protein"
        ftt <- c(ftt, PROTEIN_DOMAINS.gr)
      }
    } else{
      warning("Protein file format incorrect")
      break
    }
  }
  ###############
  
  ## Remove Introns from feature track, features, and protein track
  if(REMOVE_INTRON){
    ftt.df <- data.frame(ftt[ftt$featureLayerID == "gene"])
    ftt.df$id <- names(ftt[ftt$featureLayerID == "gene"])
    feature.df <- data.frame(feature)
    feature.df$id <- names(feature)
    if(!is.null(PROTEIN_DOMAINS)){
      protein_track.df <- data.frame(ftt[ftt$featureLayerID == "protein"])
      protein_track.df$id <- names(ftt[ftt$featureLayerID == "protein"])
    } 
    ftt.df <- ftt.df[order(ftt.df$start, decreasing=F),]
    protein_track.df <- protein_track.df[order(protein_track.df$start, decreasing=F),]
    
    for(i in 1:nrow(ftt.df)){
      #get gap between previous end and current start
      if(i == 1){
        to_subtract <- ftt.df$start[i]
      } else {
        to_subtract <- ftt.df$start[i] - ftt.df$end[i-1]
      }
      ftt.df$start[i:nrow(ftt.df)] <- ftt.df$start[i:nrow(ftt.df)] - to_subtract
      ftt.df$end[i:nrow(ftt.df)] <- ftt.df$end[i:nrow(ftt.df)] - to_subtract
      feature.df$start[feature.df$start > ftt.df$end[i]] <- feature.df$start[feature.df$start>ftt.df$end[i]] - to_subtract
      feature.df$end[feature.df$end > ftt.df$end[i]] <- feature.df$end[feature.df$end>ftt.df$end[i]] - to_subtract
      if(!is.null(PROTEIN_DOMAINS)){
        protein_track.df$start[protein_track.df$start > ftt.df$end[i]] <- protein_track.df$start[protein_track.df$start>ftt.df$end[i]] - to_subtract
        protein_track.df$end[protein_track.df$end > ftt.df$end[i]] <- protein_track.df$end[protein_track.df$end>ftt.df$end[i]] - to_subtract
      }
    }
  }
  
  #get accurate amino acid location
  aapos <- as.numeric( sub("\\D*(\\d+).*", "\\1", names(feature)) )
  
  #rebuild features
  feature <- GRanges(feature.df$seqnames, IRanges(if(PROTEIN=="protein"){aapos}else{feature.df$start}, 
                                                  if(PROTEIN=="protein"){aapos}else{feature.df$end}, 
                                                  names=feature.df$id), 
                     SNPsideID = feature.df$SNPsideID, 
                     score=feature.df$score)
  
  #one more merging in case of overlapping
  f <- data.frame(feature)
  f$names <- names(feature)
  f <- f %>% dplyr::group_by(names, SNPsideID) %>% dplyr::summarise(score = sum(score), start=start, end=end, seqnames=seqnames) %>% distinct()
  feature <- GRanges(f$seqnames, IRanges(f$start, f$end, names=f$names), SNPsideID = f$SNPsideID, score=f$score)
  
  
  #rebuild feature track and protein track
  ftt <- GRanges(ftt.df$seqnames, IRanges(if(PROTEIN=="protein"){ftt.df$start/3}else{ftt.df$start}, 
                                          if(PROTEIN=="protein"){ftt.df$end/3}else{ftt.df$end}, 
                                          names=ftt.df$id),featureLayerID = "gene")
  
  
  protein_track <- GRanges(protein_track.df$seqnames, IRanges(if(PROTEIN=="protein"){protein_track.df$start/3}else{protein_track.df$start},
                                                              if(PROTEIN=="protein"){protein_track.df$end/3}else{protein_track.df$end},
                                                              names=protein_track.df$id),featureLayerID = "protein")
  
  ftt <- c(protein_track, ftt)
  ## remove ends of feature track not within regions with somatic calls
  ###############
  if(REMOVE_EXCESS){
    ftt <- subsetByOverlaps(ftt, GRanges(unique(seqnames(feature)), IRanges(min(start(feature)), max(end(feature)))))
  }
  ftt$fill <- c("#FF8833")
  ###############
  
  
  #get breakpoints from fusions if available, future stuff...
  ###############
  if(BREAKPOINT){
    breakpoints <- MV$tables$fusion %>% dplyr::transmute(
      pass=pass,
      gene1=ifelse(gene_names.5.1!="",                 str_match(gene_names.5.1, "[^:]*")[,1],
                   ifelse(gene_names.3.2!="", sprintf("%s", str_match(gene_names.3.2, "[^:]*")[,1]), cyt.5.1)),
      gene2=ifelse(gene_names.3.1!="",                 str_match(gene_names.3.1, "[^:]*")[,1],
                   ifelse(gene_names.5.2!="", sprintf("%s", str_match(gene_names.5.2, "[^:]*")[,1]), cyt.3.1)),
      chr.5 = chr.5,pos.5 = pos.5,chr.3=chr.3,pos.3=pos.3) %>% dplyr::filter(
        pass=="pass" & (gene1 == SYMBOL) | gene2 == SYMBOL)
    breakpoints <- breakpoints %>% dplyr::select(gene1, chr.5, pos.5) %>% dplyr::filter(gene1==SYMBOL)
  }
  ###############
  
  
  # Add colors
  ###############
  if(grepl("rainbow", COLOR, ignore.case = T)){
    ggcolors <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    ftt$fill[ftt$featureLayerID=="gene"] <- ggcolors(length(ftt[ftt$featureLayerID=="gene"]))
    if(!is.null(PROTEIN_DOMAINS)){
      unique_domain <- unique(names(ftt)[ftt$featureLayerID=="protein"])
      protein_col <- ggcolors(length(unique_domain))[match(names(ftt), unique_domain)]
      protein_col <- protein_col[!is.na(protein_col)]
      ftt$fill[ftt$featureLayerID=="protein"] <- protein_col
    }
  }
  ###############
  
  if(start(ftt[names(ftt) == min(as.numeric(names(ftt)[ftt$featureLayerID=="gene"]))]) > 1){
    tot.len <- max(end(ftt))
    ftt <- GRanges(seqnames(ftt), IRanges(start=tot.len-end(ftt), 
                                          end=tot.len-start(ftt),
                                          names=names(ftt)),
                   featureLayerID=ftt$featureLayerID,
                   fill=ftt$fill)
  }
  if(!SHOW_EXON_NAME){
    names(ftt)[ftt$featureLayerID=="gene"] <- ""
  }
  
  
  
  # Miscs and labels
  ###############
  chrom = seqnames(ftt)[1]
  ftt$height=0.05
  
  ## adjusting chromosome label height
  if(LABEL_ON_FEATURE){
    ylheight=0.03
  } else {
    ylheight=0.15
  }
  
  lolliplot(feature, ftt, label_on_feature=LABEL_ON_FEATURE)
  # grid.text("label of x-axis here", x=.5, y=.01, just="bottom")
  grid.text(SYMBOL, x=.5, y=.98, just="top", 
            gp=gpar(cex=1.5, fontface="bold"))
  grid.text(chrom, x=0.05, y=ylheight, just="bottom")
  ###############
  
}



#### Misc

#' Gets tiles given GOBJ and bin width
#'
#' @param GOBJ the Gobj
#' @param BIN width of bin to use on genome; default 1MB
#' @return granges of tiles
#' 
getTiles <- function(GOBJ, BIN){
  tiles = tileGenome(GOBJ$seqi, tilewidth = BIN, cut.last.tile.in.chrom = T)
  tiles = tiles[!grepl('Y', seqnames(tiles))]
  tiles$queryHits <- 1:length(tiles)
  return(tiles)
}

#' Returns segments for all samples that have some non-NA C values
#'
#' @param MV the metavault
#' @return granges of segments
#' 
getSegments <- function(MV){
  ## get
  segments <- MV$tables$segment
  ## remove samples w/ all NAs
  segments[, allNA:= all(is.na((C))), by=list(id)]
  segments = segments[allNA==F,]
  segments[,allNA := NULL]  
  segments = GRanges(segments)
  return(segments)
}

#' Given 2 granges objects, returns the overlap w/ widths of overlap
#'
#' @param QUERY typically the tiles
#' @param SUBJECT typically the segment data
#' @return overlap df with widths of overlap
#' 
getOverlapWidth <- function(QUERY, SUBJECT){
  overlap = findOverlaps(query = QUERY, subject = SUBJECT, select='all')
  overlap_width = width(pintersect(QUERY[queryHits(overlap)],SUBJECT[subjectHits(overlap)]))
  overlap = cbind(as.data.table(overlap), overlap_width)
}

#' Given meta vault get sex for each id
#'
#' @param MV metavault
#' @return df of id/sex
#' 
getSex <- function(MV){
  tmp=MV$meta[,c('id','sex')]
  return(tmp)
}

## 
# #Diagnostic plot
# homopolymer_freq_plot_subset <- ggplot(data = homopolymer_freqs %>%
#                                     filter(length %in% 7:20),
#                                 aes(x = length, y = freq)) +
#     ylab("INDEL frequency") +
#     xlab("Homopolymer length") +
#     geom_bar(stat = "identity") +
#     annotate(geom = "text", label = strsplit(basename(args$out), '-')[[1]][1], size = 2.5,
#             x = Inf, y = Inf, hjust = 1, vjust = 1) +
#     annotate(geom = "text", label = "% sites mutated, len 7-20", size = 2.5,
#             x = Inf, y = Inf, hjust = 1, vjust = 6.5) +
#     annotate(geom = "text", label = sprintf("%0.2f", (homopolymer_freqs[homopolymer_freqs$length %in% 7:20, ]$indels %>% sum(na.rm = TRUE))/(homopolymer_freqs[homopolymer_freqs$length %in% 7:20, ]$sites %>% sum())*100), size = 2.5,
#             x = Inf, y = Inf, hjust = 1, vjust = 8) +
#     coord_cartesian(xlim = c(6.5, 20.5), ylim = c(0, .25))
# homopolymer_freq_plot_subset
# ggsave(filename = paste(args$out, "-plot.png", sep = ""),
#     plot = homopolymer_freq_plot_subset,
#     height = 3, width = 3, units = "in", device = "png")

