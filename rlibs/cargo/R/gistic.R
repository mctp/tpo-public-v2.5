
#### Prepare

#' Prep gistic table
#'
#' @param MV The metavault
#' @param SEG_CN_METHOD How the seg.CN parameter is calculated
#' @return Table for gistic
#' 
#' @export
prepGisticTable <- function(MV, SEG_CN_METHOD='centered_c'){
  
  ## get
  ids = getVaultId(MV)
  tbl = getVaultTable(V=MV, FIELD='segment')
  
  ## run
  lapply(ids, .prepGisticTable, tbl, SEG_CN_METHOD) %>%
    rbindlist() %>%
    setkey(., sample, chromosome, start.pos) %>%
    return()
}

#' Prep gistic table for single sample
#'
#' @param ID The sample ID
#' @param TBL The segment table
#' @param SEG_CN_METHOD How the seg.CN parameter is calculated
#' @return Table for gistic for a single sample
#' 
.prepGisticTable <- function(ID, TBL, SEG_CN_METHOD){
  
  ## get
  seg <- TBL[id == ID][order(chr)]
  
  ## error handling
  if( nrow(seg) == 0 ){ return(NULL) }
  if( all(is.na(seg[['C']])) ){ return(NULL) }
  
  ## reformat
  data.table(sample = seg[['id']], 
             chromosome = gsub("^chr", "", seg[['chr']]), 
             start.pos = seg[['start']], 
             end.pos = seg[['end']],
             num.markers = seg[['nlr']],
             seg.CN = .segCN(seg, SEG_CN_METHOD)
  ) %>%
    return()
}

#' Compute seg.CN
#'
#' @param SEGMENTS The vault segment table
#' @param SEG_CN_METHOD How the seg.CN parameter is calculated
#' @return vector of seg.CN values

.segCN <- function(SEGMENTS, SEG_CN_METHOD='centered_c'){
  ## default
  if(SEG_CN_METHOD=='centered_c'){
    return( (SEGMENTS[['C']] - weighted.median( SEGMENTS[['C']], as.numeric(SEGMENTS[['width']]), na.rm = TRUE)) )  
  }
  ## Add others PRN
  if(SEG_CN_METHOD=='lr'){
    return(SEGMENTS[['lr']])  
  }
  
  ## invalid method
  print(paste0('Invalid SEG_CN_METHOD: `', SEG_CN_METHOD, '` selected'))
  return(NA)
}




#### Plot

#' Plot MAF's Chrom plot from GISTIC data
#'
#' @param ALL_LESTIONS path to all_lesions.conf_99.txt
#' @param AMP.GENES path to amp_genes.conf_99.txt
#' @param DEL.GENES path to del_genes.conf_99.txt
#' @param SCORES.GISTIC path to scores.gistic
#' @return plot
#' 
#' @export
plotMafGisticChromPlot <- function(ALL_LESIONS, AMP.GENES, DEL.GENES, SCORES.GISTIC,
                                   FDR=0.1, 
                                   REF_BUILD='hg38', 
                                   CYTOBAND_TEXT_SIZE=1.3, 
                                   CYTOBAND_OFFSET=0.05,
                                   GENES_TEXT_SIZE=10){
  
  maftools::readGistic(
    gisticAllLesionsFile = ALL_LESIONS, 
    gisticAmpGenesFile = AMP.GENES, 
    gisticDelGenesFile = DEL.GENES, 
    gisticScoresFile = SCORES.GISTIC
  ) %>%
    maftools::gisticChromPlot(
      ., 
      fdrCutOff = FDR, 
      ref.build = REF_BUILD, 
      cytobandTxtSize = CYTOBAND_TEXT_SIZE, 
      cytobandOffset = CYTOBAND_OFFSET, 
      mutGenesTxtSize = GENES_TEXT_SIZE
    ) %>%
    return()
}

#### Run Gisitc

#' Run GISTIC using the docker image at https://github.com/ShixiangWang/install_GISTIC
#'
#' @param MV Metavault containing all the samples
#' @param samples samples to include in analysis (if NULL it is all the samples)
#' @param out.fn output directory to write the results
#' @param SEG_CN_METHOD CNV value to be used for segment table of Gistic
#' @return 0
#' 
#' @export

runGistic <- function(MV, samples=NULL, out.fn, SEG_CN_METHOD = "centered_c") {
  gistictbl <- prepGisticTable(MV, SEG_CN_METHOD)
  
  ## pick samples
  if(!is.null(samples)) {
    gistictbl <- gistictbl[sample %in% samples]
  }
  
  ## make the input file
  inp <- tempfile("gisitc_")
  fwrite(gistictbl, inp,sep = "\t", col.names = FALSE)
  
  ## prep the docker image
  image.exist <- system("docker manifest inspect docker.io/shixiangwang/gistic:latest > /dev/null ; echo $?", intern = TRUE)
  if(image.exist != "0") {
    docker.install <- system("docker pull shixiangwang/gistic")  
  }
  
  ## cp input to docker directory
  system(sprintf("cp %s %s/USER_INPUT.txt", inp, out.fn))
  
  ## prep the docker script
  GISTIC_LOC="/opt/GISTIC"
  refgenefile=sprintf("%s/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat", GISTIC_LOC)
  #sprintf("SegFile at: %s", segfile)
  sprintf("Output at: %s", out.fn)
  sprintf("Reference at: %s", refgenefile)
  DOCKER_INPDIR <- sprintf("%s/run_result/USER_INPUT.txt", GISTIC_LOC)
    
  ## run GISITC
  print("--- running GISTIC ---")
  print("[1] Starting gistic...")
  DOCKER_OUTDIR=sprintf("%s/run_result", GISTIC_LOC)
  docker_scripts = sprintf("docker run --rm -v %s:%s shixiangwang/gistic -b %s -seg %s -refgene %s -rx 0 -js 4 -broad 0.98 -cap 1.5 -twoside 1 -genegistic 1 -smallmem 0 -conf 0.99 -armpeel 1 -savegene 1 -saveseg 1 -qvt 0.1",
                           out.fn,
                           DOCKER_OUTDIR,
                           DOCKER_OUTDIR,
                           DOCKER_INPDIR,
                           refgenefile)
  
  system(docker_scripts)
  system(sprintf("rm %s", inp))
  
  # ## make the plot
  # ALL_LESIONS    <- inp
  # AMP.GENES
  # DEL.GENES
  # SCORES.GISTIC
  # 
  # plotMafGisticChromPlot(ALL_LESIONS, AMP.GENES, DEL.GENES, SCORES.GISTIC,
  #                                    FDR=0.1, 
  #                                    REF_BUILD='hg38', 
  #                                    CYTOBAND_TEXT_SIZE=1.3, 
  #                                    CYTOBAND_OFFSET=0.05,
  #                                    GENES_TEXT_SIZE=10)
  
  return(0)
}
