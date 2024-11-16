#### directory listing

.quasrDir <- function(pth, sample) {
    if (is.null(pth)) {
        dir <- NULL
    } else {
        dir <- list(
            config=list.files.null(pth, "-config.txt$", full.names = TRUE),
            count=list.files.null(pth, "-count$", full.names = TRUE),
            count.log=list.files.null(pth, "-count.log$", full.names = TRUE),
            count.summary=list.files.null(pth, "-count.summary$", full.names = TRUE)
        )
    }
    return(dir)
}

.cnvexDir <- function(pth) {
    if (is.null(pth)) {
        dir <- NULL
    } else {
        dir <- list(
            config=list.files.null(pth, "-config.txt$", full.names = TRUE),
            digest=list.files.null(pth, "-somatic-digest.*rds$", full.names = TRUE, n=1)
        )
    }
    return(dir)
}

.codacDir <- function(pth) {
    if (is.null(pth)) {
        dir <- NULL
    } else {
        dir <- list(
            config=list.files.null(pth, "-config.txt$", full.names = TRUE),
            fusion=list.files.null(pth, "-rep-sv.rds$", full.names = TRUE)
        )
    }
    return(dir)
}

.caratDir <- function(pth) {
    if (is.null(pth)) {
        dir <- NULL
    } else {
        dir <- list(
            config=list.files.null(pth, "-config.txt$", full.names = TRUE),
            somatic=list.files.null(pth, "-somatic-tnscope-report.rds$|-somatic-mutect2-report.rds$", full.names = TRUE),
            somatic.vcf=list.files.null(pth, "-somatic-tnscope-annotated.vcf.gz$|-somatic-mutect2-annotated.vcf.gz$", full.names = TRUE),
            germline=list.files.null(pth, "-somatic-dnascope-report.rds$|-somatic-haplotypecaller-report.rds$", full.names = TRUE),
            germline.vcf=list.files.null(pth, "-somatic-dnascope-annotated.vcf.gz$|-somatic-haplotypecaller-annotated.vcf.gz$", full.names = TRUE),
            structural=list.files.null(pth, "-structural-tnscope-report.rds$", full.names = TRUE),
            msi=list.files.null(pth, "-msi-score.csv$", full.names = TRUE)
        )
    }
    return(dir)
}

.miscDir <- function(pth) {
    if (is.null(pth)) {
        dir <- NULL
    } else {
        dir <- list(
            config=list.files.null(pth, "-config.txt$", full.names = TRUE),
            summary=list.files.null(pth, "-QC-summary.rds$", full.names = TRUE)
        )
    }
    return(dir)
}




#### data imports

.import.config <- function(vdir, input) {
    if (is.null(vdir[[input]]$config)) {
        config <- NULL
    } else {
        config <- fread(vdir[[input]]$config[1], sep="=", header=FALSE)
        setnames(config, c("variable", "value"))
    }
    return(config)
}

importQuasr <- function(vdir) {
    ## tumor
    if (is.null(vdir$quasr_t$count)) {
        summary <- NULL
        gene.count <- NULL
    } else {
        ## summary
        summary <- melt(fread(vdir$quasr_t$count.summary), id="Status")[,.(variable=Status, read.pairs=value)]
        ## counts
        tmp <- fread(vdir$quasr_t[["count"]])
        gene.count <- tmp[,c(1,6,7)]
        setnames(gene.count, c("gene_id", "gene_length", "count"))
    }
    tlist <- list(
        tconfig=.import.config(vdir, "quasr_t"),
        tsummary=summary,
        tgene.count=gene.count
    )
    ## normal
    if (is.null(vdir$quasr_n$count)) {
      summary <- NULL
      gene.count <- NULL
    } else {
      ## summary
      summary <- melt(fread(vdir$quasr_n$count.summary), id="Status")[,.(variable=Status, read.pairs=value)]
      ## counts
      tmp <- fread(vdir$quasr_n[["count"]])
      gene.count <- tmp[,c(1,6,7)]
      setnames(gene.count, c("gene_id", "gene_length", "count"))
    }
    nlist <- list(
      nconfig=.import.config(vdir, "quasr_n"),
      nsummary=summary,
      ngene.count=gene.count
    )

    return(c(tlist, nlist))
}

importCnvex <- function(vdir) {
    if (is.null(vdir$cnvex$digest)) {
        digest <- NULL
    } else {
        digest <- readRDS(vdir$cnvex$digest)
    }
    list(
        config=.import.config(vdir, "cnvex"),
        digest=digest
    )
}

importCodac <- function(vdir) {
    if (is.null(vdir$codac$fusion)) {
        fusion <- NULL
    } else {
        fusion <- readRDS(vdir$codac$fusion)
        if (dim(fusion)[1] < 1) {
          fusion <- NULL
        }
    }
    codac.reps <- filterCodac(list(
        config=.import.config(vdir, "codac"),
        fusion=fusion
    ))
    return(codac.reps)
}

importCarat <- function(vdir) {
    if (is.null(vdir$carat$somatic)) {
        somatic <- NULL
    } else {
        somatic <- readRDS(vdir$carat$somatic)
        if (dim(somatic)[1] < 1) {
          somatic <- NULL
        }
    }
    if (is.null(vdir$carat$germline)) {
        germline <- NULL
    } else {
        germline <- readRDS(vdir$carat$germline)
        if (dim(germline)[1] < 1) {
          germline <- NULL
        }
    }
    if (is.null(vdir$carat$structural)) {
        structural <- NULL
    } else {
        structural <- readRDS(vdir$carat$structural)
        if (dim(structural)[1] < 1) {
          structural <- NULL
        }
    }
    if (is.null(vdir$carat$msi)) {
        msi <- NULL
    } else {
        msi <- fread(vdir$carat$msi)
        if (dim(msi)[1] < 1) {
          msi <- NULL
        }
    }

    cfg <- .import.config(vdir, "carat")
    libs <- c(NA, NA)
    if (!is.null(vdir$carat$somatic.vcf)) {
      tmp <- colnames(readVcf(vdir$carat$somatic.vcf))
      if (length(tmp)==1) {
        libs[1] <- tmp
      }
      if(length(tmp)==2) {
        libs <- tmp
      }
    }
    cfg <- rbind(cfg, data.table(variable=c("TUMOR.SAMPLE", "NORMAL.SAMPLE"), value=libs))
    carat.reps <- filterCarat(tumor_only=length(tmp)==1,
      list(
        config=cfg,
        somatic=somatic,
        germline=germline,
        structural=structural,
        msi=msi
    ))
    return(carat.reps)
}

importMisc <- function(vdir, type) {
    if (type=="tumor") {
      tmpdir <- vdir$misc_t$summary
    } else {
      tmpdir <- vdir$misc_n$summary
    }
    if (is.null(tmpdir)) {
        qc <- NULL
    } else {
        qc <- readRDS(tmpdir)
        if (dim(qc)[1] < 1) {
          qc <- NULL
        }
    }
    return(qc)
}




#### data pre-filters

filterCodac <- function(codac.reps) {
    if (!is.null(codac.reps$fusion)) {
        codac.reps$fusion <- fusion.prefilter(codac.reps$fusion)
    }
    return(codac.reps)
}

filterCarat <- function(carat.reps,tumor_only=FALSE) {
    ## somatic
    if (!is.null(carat.reps$somatic)) {
        if(tumor_only) {
          carat.reps$somatic <- somatic.to.prefilter(carat.reps$somatic)
          print("using tumor only prefilters for somatic")
        } else {
          carat.reps$somatic <- somatic.prefilter(carat.reps$somatic)
        }
    }
    ## germline
    if (!is.null(carat.reps$germline)) {
        if(tumor_only) {
          carat.reps$germline <- germline.to.prefilter(carat.reps$germline)
          print("using tumor only prefilters for germline")
        } else {
          carat.reps$germline <- germline.prefilter(carat.reps$germline)
        }
    }
    ## structural
    if (!is.null(carat.reps$structural)) {
        carat.reps$structural <- structural.prefilter(carat.reps$structural)
    }
    return(carat.reps)
}




#### tables / maps

createMaps <- function(data, anno, cfg) {
    maps <- list(
        somatic=somaticMap(data$carat, cfg),
        germline=germlineMap(data$carat, cfg),
        structural=structuralMap(data$carat, cfg),
        fusion=fusionMap(data$codac, cfg),
        segment=segmentMap(data$cnvex, anno, cfg)
    )
    return(maps)
}

createTables <- function(data, anno, cfg) {
    tbls <- list(
        somatic=somaticTable(data$carat, anno, cfg),
        germline=germlineTable(data$carat, anno, cfg),
        structural=structuralTable(data$carat, anno, cfg),
        fusion=fusionTable(data$codac, anno, cfg),
        segment=segmentTable(data$cnvex, anno, cfg),
        gene.copy=geneCopyTable(data$cnvex, anno, cfg),
        gene.expression=geneExpressionTable(data$quasr, anno, cfg)
    )
    return(tbls)
}

## structural map/table

structuralMap <- function(carat.reps, cfg) {
    structural.map <- NULL
    if (!is.null(carat.reps$structural)) {
        tmp <- carat.reps$structural
        tmp.1 <- tmp[,.(gene_id=GENE1, id=cfg$meta$group$id, var_id=bnd.id, end="e1")]
        tmp.2 <- tmp[,.(gene_id=GENE2, id=cfg$meta$group$id, var_id=bnd.id, end="e2")]
        structural.map <- unique(rbind(tmp.1, tmp.2))
        setkey(structural.map, id, var_id, end)
        setcolorder(structural.map)
    }
    return(structural.map)
}

structuralTable <- function(carat.reps, anno, cfg) {
    structural.tbl <- NULL
    if (!is.null(carat.reps$structural)) {
        structural.tbl <- structuralCSQ(carat.reps$structural)
        structural.tbl[, id:=cfg$meta$group$id]
        structural.tbl[, var_id:=bnd.id]
        setkey(structural.tbl, id, var_id)
        setcolorder(structural.tbl)
    }
    return(structural.tbl)
}

## somatic map/table

somaticMap <- function(carat.reps, cfg) {
    somatic.map <- NULL
    if (!is.null(carat.reps$somatic)) {
        tmp <- carat.reps$somatic
        somatic.map <- tmp[,.(gene_id=GENE, id=cfg$meta$group$id, var_id=ID)]
        setkey(somatic.map, id, var_id)
        setcolorder(somatic.map)
    }
    return(somatic.map)
}

somaticTable <- function(carat.reps, anno, cfg) {
    somatic.tbl <- NULL
    if (!is.null(carat.reps$somatic)) {
        somatic.tbl <- carat.reps$somatic
        somatic.tbl[,id:=cfg$meta$group$id]
        somatic.tbl[,var_id:=ID]
        setkey(somatic.tbl, id, var_id)
        setcolorder(somatic.tbl)
    }
    return(somatic.tbl)
}

## germline map/table

germlineMap <- function(carat.reps, cfg) {
    germline.map <- NULL
    if (!is.null(carat.reps$germline)) {
        tmp <- carat.reps$germline
        germline.map <- tmp[,.(gene_id=GENE, id=cfg$meta$group$id, var_id=ID)]
        setkey(germline.map, id, var_id)
        setcolorder(germline.map)
    }
    return(germline.map)
}

germlineTable <- function(carat.reps, anno, cfg) {
    germline.tbl <- NULL
    if (!is.null(carat.reps$germline)) {
        germline.tbl <- carat.reps$germline
        germline.tbl[, id:=cfg$meta$group$id]
        germline.tbl[, var_id:=ID]
        setkey(germline.tbl, id, var_id)
        setcolorder(germline.tbl)
    }
    return(germline.tbl)
}

## fusion

fusionMap <- function(codac.reps, cfg) {
    ## gene to junction
    fusion.map <- NULL
    if (!is.null(codac.reps$fusion)) {
        tmp <- codac.reps$fusion
        tmp[,var_id:=sprintf("%s:%s:%s::%s:%s:%s", chr.5, pos.5, str.5, chr.3, pos.3, str.3)]
        tmp.5.1 <- tmp[,.(
            gene_id=unlist(str_split(gene_ids.5.1, ":")), var_id),
            .(gene_ids.5.1,var_id)][,.(var_id, gene_id, end="f5")]
        tmp.3.1 <- tmp[,.(
            gene_id=unlist(str_split(gene_ids.3.1, ":")), var_id),
            .(gene_ids.3.1,var_id)][,.(var_id, gene_id, end="f3")]
        tmp.5.2 <- tmp[,.(
            gene_id=unlist(str_split(gene_ids.5.2, ":")), var_id),
            .(gene_ids.5.2,var_id)][,.(var_id, gene_id, end="r5")]
        tmp.3.2 <- tmp[,.(
            gene_id=unlist(str_split(gene_ids.3.2, ":")), var_id),
            .(gene_ids.3.2,var_id)][,.(var_id, gene_id, end="r3")]
        tmp.53 <- unique(rbind(tmp.5.1, tmp.3.1, tmp.5.2, tmp.3.2)[gene_id!=""])
        tmp.53[,id:=cfg$meta$group$id]
        fusion.map <- tmp.53[,.(gene_id, id, var_id, end)]
        setkey(fusion.map, id, var_id, end)
        setcolorder(fusion.map)
    }
    return(fusion.map)
}

fusionTable <- function(codac.reps, anno, cfg) {
    fusion.tbl <- NULL
    if (!is.null(codac.reps$fusion)) {
        fusion.tbl <- fusionCSQ(codac.reps$fusion)
        fusion.tbl[,id:=cfg$meta$group$id]
        fusion.tbl[,var_id:=sprintf("%s:%s:%s::%s:%s:%s", chr.5, pos.5, str.5, chr.3, pos.3, str.3)]
        setkey(fusion.tbl, id, var_id)
        setcolorder(fusion.tbl)
    }
    return(fusion.tbl)
}

## CNV

segmentTable <- function(cnvex.reps, anno, cfg) {
    segment.tbl <- NULL
    if (!is.null(cnvex.reps$digest)) {
        seg <- cnvex.reps$digest$seg
        fit <- cnvex.reps$digest$fit
        tmp <- cbind(as.data.table(seg), fit)
        setnames(tmp, "seqnames", "chr")
        segment.tbl <- cbind(id=cfg$meta$group$id, var_id=sprintf("%s:%s-%s", tmp$chr, tmp$start, tmp$end), tmp)
        setkey(segment.tbl, id, var_id)
        setcolorder(segment.tbl)
    }
    return(segment.tbl)
}

segmentMap <- function(cnvex.reps, anno, cfg) {
    segment.map <- NULL
    if (!is.null(cnvex.reps$digest)) {
        seg <- cnvex.reps$digest$seg
        mcols(seg) <- cnvex.reps$digest$fit
        tmp <- as.data.table(findOverlaps(seg, anno$genes))
        setkey(tmp, queryHits)
        tmp <- cbind(tmp, as.data.table(mcols(anno$genes[tmp$subjectHits])))
        tmp <- cbind(tmp, as.data.table(seg[tmp$queryHits]))
        segment.map <- tmp[,.(gene_id, id=cfg$meta$group$id, var_id=sprintf("%s:%s-%s", seqnames, start, end))]
        setkey(segment.map, id, var_id)
        setcolorder(segment.map)
    }
    return(segment.map)
}

geneCopyTable <- function(cnvex.reps, anno, cfg) {
    ## merge segments, genes, and fit information
    genecopy.tbl <- NULL
    if (!is.null(cnvex.reps$digest)) {
        genes.clipped <- .genesClipped(cnvex.reps$digest$tile, anno$genes)
        seg <- cnvex.reps$digest$seg
        mcols(seg) <- cnvex.reps$digest$fit
        tmp <- as.data.table(findOverlaps(seg, genes.clipped))
        setkey(tmp, queryHits)
        tmp <- cbind(tmp, as.data.table(genes.clipped[tmp$subjectHits]))
        tmp <- cbind(tmp, as.data.table(mcols(seg[tmp$queryHits])))
        tmp <- tmp[,.(seqnames,start,end,gene_id, gene_name, C, sC, K)]
        tmp[C==0,K:=0]
        tmp[C==1,K:=0]
        genecopy.tbl <- suppressWarnings(tmp[,.(
            nseg=.N,
            Kmin=min(K, na.rm = TRUE),
            Kmax=max(K, na.rm = TRUE),
            Cmin=min(C, na.rm = TRUE),
            Cmax=max(C, na.rm = TRUE),
            Csub=mean(sC, na.rm = TRUE)
        ),.(seqnames,start,end,gene_id,gene_name)]) ## genes with only NA's raise warning
        genecopy.tbl[is.infinite(Kmin) | is.infinite(Kmax), ":="(Kmin=NA, Kmax=NA)]
        setkey(genecopy.tbl, gene_id)
        genecopy.tbl <- genecopy.tbl[SJ(genes.clipped$gene_id)] ## add genes on chrM
        genecopy.tbl[,id:=cfg$meta$group$id]
        setnames(genecopy.tbl, c("gene_id", "gene_name"), c("GENE","SYMBOL"))
        setkey(genecopy.tbl, id, GENE)
    }
    return(genecopy.tbl)
}

## RNA

geneExpressionTable <- function(quasr.reps, anno, cfg) {
    expression.tbl <- NULL
    ## tumor
    if (!is.null(quasr.reps$tgene.count)) {
        expression.tbl <- copy(quasr.reps$tgene.count)
        setnames(expression.tbl, "count", "tcount")
        tmp <- as.data.table(anno[["genes"]])[, c("gene_id", "gene_name")]
        expression.tbl <- merge(tmp, expression.tbl, by="gene_id", all.y=TRUE)
        expression.tbl$tcpm <- cpm(expression.tbl$tcount)
        expression.tbl$trpkm <- rpkm(expression.tbl$tcount, expression.tbl$gene_length)
        expression.tbl <- expression.tbl[,.(gene_id, gene_name, gene_length, tcount, tcpm, trpkm)]
        setkey(expression.tbl, gene_id)
        expression.tbl <- expression.tbl[SJ(anno$genes$gene_id)] ## add genes on chrM
        expression.tbl[,id:=cfg$meta$group$id]
    }
    texpression.tbl <- copy(expression.tbl)
    ## normal
    expression.tbl <- NULL
    if (!is.null(quasr.reps$ngene.count)) {
      expression.tbl <- copy(quasr.reps$ngene.count)
      setnames(expression.tbl, "count", "ncount")
      tmp <- as.data.table(anno[["genes"]])[, c("gene_id", "gene_name")]
      expression.tbl <- merge(tmp, expression.tbl, by="gene_id", all.y=TRUE)
      expression.tbl$ncpm <- cpm(expression.tbl$ncount)
      expression.tbl$nrpkm <- rpkm(expression.tbl$ncount, expression.tbl$gene_length)
      expression.tbl <- expression.tbl[,.(gene_id, gene_name, gene_length, ncount, ncpm, nrpkm)]
      setkey(expression.tbl, gene_id)
      expression.tbl <- expression.tbl[SJ(anno$genes$gene_id)] ## add genes on chrM
      expression.tbl[,id:=cfg$meta$group$id]
    }
    nexpression.tbl <- copy(expression.tbl)

    if(!is.null(quasr.reps$ngene.count) && !is.null(quasr.reps$ngene.count)) {
      expression.tbl <- base::merge(texpression.tbl, nexpression.tbl[, c("gene_id", "ncount", "ncpm", "nrpkm")], by = "gene_id", all.y = TRUE)
      setnames(expression.tbl, c("gene_id", "gene_name"), c("GENE","SYMBOL"))
    } else if(!is.null(quasr.reps$tgene.count)) {
      expression.tbl <- texpression.tbl
      setnames(expression.tbl, c("gene_id", "gene_name"), c("GENE","SYMBOL"))
      expression.tbl$ncount <- NA_real_
      expression.tbl$ncpm <- NA_real_
      expression.tbl$nrpkm <- NA_real_
    } else if (!is.null(quasr.reps$ngene.count)) {
      expression.tbl <- nexpression.tbl
      setnames(expression.tbl, c("gene_id", "gene_name"), c("GENE","SYMBOL"))
      expression.tbl$tcount <- NA_real_
      expression.tbl$tcpm <- NA_real_
      expression.tbl$trpkm <- NA_real_
    }
    if (!is.null(expression.tbl)) {
        setcolorder(expression.tbl, c("id", "GENE", "SYMBOL", "gene_length",
                                      "tcount", "tcpm", "trpkm",
                                      "ncount", "ncpm", "nrpkm"))
        setkey(expression.tbl, id, GENE)
    }
    return(expression.tbl)
}




#### vault

#' @export
vaultDir <- function(cfg) {
    vdir <- list(
        "quasr_t"=.quasrDir(cfg$data$tquasr),
        "quasr_n"=.quasrDir(cfg$data$nquasr),
        "codac"=.codacDir(cfg$data$codac),
        "cnvex"=.cnvexDir(cfg$data$cnvex),
        "carat"=.caratDir(cfg$data$carat),
        "misc_t"=.miscDir(cfg$data$misc_t),
        "misc_n"=.miscDir(cfg$data$misc_n)
    )
    return(vdir)
}

#' @export
vaultDump <- function(dir) {
    impt <- list(
        cnvex = importCnvex(dir),
        quasr = importQuasr(dir),
        codac = importCodac(dir),
        carat = importCarat(dir),
        misc_t = importMisc(dir, "tumor"),
        misc_n = importMisc(dir, "normal")
    )
    return(impt)
}

#' @export
vaultAnno <- function(cfg) {
    if (!is.null(cfg$annotation$gene.model)) {
        genes <- import(cfg$annotation$gene.model, feature.type="gene")
        if (cfg$drop.chrM) {
            genes <- dropSeqlevels(genes, "chrM", pruning.mode = "coarse")
        }
    } else {
        genes <- NULL
    }
    if (!is.null(cfg$annotation$homopolymers)) {
        homopolymers <- readRDS(cfg$annotation$homopolymers)
    } else {
        homopolymers <- NULL
    }
    anno <- list(
        genes=genes,
        homopolymers=homopolymers
    )
    return(anno)
}

#' @export
vaultData <- function(dump, anno, cfg) {
    ## output
    maps <- createMaps(dump, anno, cfg)
    tables <- createTables(dump, anno, cfg)
    data <- list(
        maps=maps,
        tables=tables
    )
    return(data)
}

#' @export
vaultMeta <- function(dump, cfg) {

    ## initialize
    meta <- data.table(id=NA, case=NA, cohort=NA)

    ## case/cohort/id
    tmp <- cfg[["meta"]][["group"]][["id"]]
    meta[["id"]] <- ifelse(length(tmp)>0, tmp, NA)
    tmp <- cfg[["meta"]][["group"]][["case"]]
    meta[["case"]] <- ifelse(length(tmp)>0, tmp, NA)
    tmp <- cfg[["meta"]][["group"]][["cohort"]]
    meta[["cohort"]] <- ifelse(length(tmp)>0, tmp, NA)

    ## tlib/nlib
    tmp <- ifelse(!is.null(dump$carat$config), dump$carat$config[variable=="TUMOR.SAMPLE",value], NA)
    meta[["tumor.dna.sample"]] <- ifelse(length(tmp)>0, tmp, NA)
    tmp <- ifelse(!is.null(dump$carat$config), dump$carat$config[variable=="NORMAL.SAMPLE",value], NA)
    meta[["normal.dna.sample"]] <- ifelse(length(tmp)>0, tmp, NA)

    ## rlib
    tmp <- ifelse(!is.null(dump$codac$config), dump$codac$config[variable=="ID",value], NA)
    tmp2 <- ifelse(!is.null(dump$quasr$tconfig), dump$quasr$tconfig[variable=="ID",value], NA)
    meta[["tumor.rna.sample"]] <-  pickExisting(tmp, tmp2)
    ### TODO: EMPTY FOR CODAC NORMAL CODE ###
    tmp3 <- ifelse(!is.null(dump$quasr$nconfig), dump$quasr$nconfig[variable=="ID",value], NA)
    meta[["normal.rna.sample"]] <- tmp3

    ## derived - cnvex
    tmp <- dump$cnvex$digest$sex
    meta[["sex"]] <- ifelse(length(tmp)>0, tmp, NA)
    tmp <- dump$cnvex$digest$purity
    meta[["purity"]] <- ifelse(length(tmp)>0, tmp, NA)
    tmp <- dump$cnvex$digest$ploidy
    meta[["ploidy"]] <- ifelse(length(tmp)>0, tmp, NA)

    setkey(meta, id)

    return(meta)
}

#' @export
vaultQC <- function(dump) {
    ## QC
    tmp <- dump$misc_t
    if (!is.null(tmp)) {
        qc <- tmp
    } else {
        qc <- NULL
    }
    tmp <- dump$misc_n
    if (!is.null(tmp)) {
        qc <- rbind(tmp,qc)
    }
    return(qc)
}

#' @export
vault <- function(data, anno, meta, qc) {
    ## output
    vault <- data
    vault$meta <- meta
    vault$qc <- qc
    return(vault)
}

#### Get/set

#' Gets specific vault table
#'
#' @param V vault
#' @param FIELD which table
#' @return vault table
#'
#' @export
getVaultTable <- function(V, FIELD) {
  return(V[["tables"]][[FIELD]])
}

#' Gets specific vault map
#'
#' @param V vault
#' @param FIELD which map
#' @return vault map
#'
#' @export
getVaultMap <- function(V, FIELD) {
  return(V[["maps"]][[FIELD]])
}

#' Gets specific vault meta
#'
#' @param V vault
#' @return vault meta
#'
#' @export
getVaultMeta <- function(V) {
  return(V[["meta"]])
}

#' Gets specific vault qc
#'
#' @param V vault
#' @return vault qc
#'
#' @export
getVaultQC <- function(V) {
  return(V[["qc"]])
}

#' Sets specific vault table
#'
#' @param V vault
#' @param FIELD which table
#' @param TBL new table
#' @return vault
#'
#' @export
setVaultTable <- function(V, FIELD, TBL) {
  V[["tables"]][[FIELD]] <- TBL
  return(V)
}

#' Sets specific vault map
#'
#' @param V vault
#' @param FIELD which map
#' @param MAP new map
#' @return vault
#'
#' @export
setVaultMap <- function(V, FIELD, MAP) {
  V[["maps"]][[FIELD]] <- MAP
  return(V)
}

#' Sets specific vault meta
#'
#' @param V vault
#' @param META new meta
#' @return vault
#'
#' @export
setVaultMeta <- function(V, META) {
  V[["meta"]] <- META
  return(V)
}

#' Sets specific vault QC
#'
#' @param V vault
#' @param QC new qc
#' @return vault
#'
#' @export
setVaultQC <- function(V, QC) {
  V[["qc"]] <- QC
  return(V)
}
