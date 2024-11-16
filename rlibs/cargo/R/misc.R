#' @export
cargoConfig <- function(settings, args=list()) {
    cfg <- get.opts(settings, args)
    return(cfg)
}

.configToGroup <- function(vdir) {
    cfg.fn <- list.files(vdir, pattern="config.txt$", full.names=TRUE)
    cfg <- fread(cfg.fn, sep="=", quote='"', header=FALSE)
    setnames(cfg, c("key", "value"))
    ## patient-tumor.normal.rna
    tmp1 <- str_split(cfg[key=="ID"]$value, "\\-")
    tmp2 <- str_split(tmp1[[1]][2], "\\.")
    group <- list(
        case = tmp1[[1]][1],
        cohort = NULL,
        tumor.dna.sample=tmp2[[1]][1],
        normal.dna.sample=tmp2[[1]][2],
        tumor.rna.sample=tmp2[[1]][3],
        normal.rna.sample=NULL
    )
    return(group)
}

.genesClipped <- function(tile, genes) {

    ## target regions
    tile <- tile[tile$target == TRUE]

    hits <- as.data.table(findOverlaps(genes, tile, minoverlap = 1))

    ## clip boundry
    hits <- hits[, .(start.tile = min(subjectHits), end.tile = max(subjectHits)), by = queryHits]

    ## pick targeted genes
    genes.clip <- genes[hits$queryHits]

    ## clip targeted genes
    ranges(genes.clip) <- IRanges(start = pmax(start(tile[hits$start.tile]), start(genes.clip)),
                                  end = pmin(end(tile[hits$end.tile]), end(genes.clip)))

    ## add off targeted genes
    offtarget.genes <- genes[!(genes$gene_id %in% genes.clip$gene_id)]
    genes.out <- c(offtarget.genes, genes.clip)
    genes.out <- sort(genes.out)

    return(genes.out)
}

#' Loads GTF file from location and adds useful columns
#'
#'
#' @param GTF path to gtf file, default currently is my file
#' @return GTF object
#' 
#' @export
load_GTF <- function(GTF){
  
  ## load
  message("Loading annotation")
  gtf <- scan(GTF, what=list(contig="",src="",type="",start=0,end=0,score="",strand="",frame="",attributes=""), sep="\t", comment.char="#", quote='"', multi.line=F)
  attr(gtf, "row.names") <- .set_row_names(length(gtf[[1]]))
  class(gtf) <- "data.frame"
  gtf[['contig']] <- gsub('^chr', '', gtf[['contig']])
  
  ## parse
  gtf[['geneID']] <- .parseGtfAttribute("gene_id", gtf)
  gtf[['geneName']] <- .parseGtfAttribute("gene_name", gtf)
  gtf[['geneName']] <- ifelse(gtf[['geneName']] == "", gtf[['geneID']], gtf[['geneName']] )
  gtf[['transcript']] <- .parseGtfAttribute("transcript_id", gtf)
  gtf[['exonNumber']] <- .parseGtfAttribute("exon_number", gtf)
  gtf[['gene_version']] <- .parseGtfAttribute("gene_version", gtf)
  gtf[['transcript_version']] <- .parseGtfAttribute("transcript_version", gtf)
  
  ## return
  return(gtf)
  
}

#' Parses GTF attribute
#'
#' @param ATTRIBUTE - the attribute to parse
#' @param GTF path to gtf file, default currently is my file
#' @return GTF object
#' 
.parseGtfAttribute <- function(ATTRIBUTE, GTF) {
  parsed <- sub(paste0(".*", ATTRIBUTE, "[ =]([^;]+).*"), "\\1", GTF[['attributes']], perl=T)
  parsed <- ifelse(parsed == GTF[['attributes']], "", parsed)
  return(parsed)
}
