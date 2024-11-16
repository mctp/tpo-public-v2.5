#' Create dtVault object from legacy vault or metavault
#'
#' @param vault legacy vault
#' @param anno dtVaultAnno object
#' @param gene_model GTF file
#' @return dtVault
#'
#' @export
dtVault <- function(vault, anno=NULL, gene_model=NULL, ...) {
    if (is.null(anno)) {
        anno <- dtVaultAnno(gene_model, ...)
    }
    dtv <- dtVaultUpdate(vault, anno)
    return(dtv)
}

#' Create dtVaultAnno object
#'
#' @param gene_model GTF file
#' @param somatic_goi custom list of somatic genes of interest [default: goi-somatic-tier1]
#' @param germline_goi custom list of germline genes of interest [default: goi-germline]
#' @param somatic_art custom list of somatic artifact genes [default: art-somatic]
#' @param drop_chrM remove genes on chrM
#' @return dtVaultAnno
#'
#' @export
dtVaultAnno <- function(gene_model, somatic_goi=NULL, germline_goi=NULL, somatic_art=NULL, drop_chrM=TRUE) {
    ## somatic GOI
    if (is.null(somatic_goi)) {
        somatic_goi <- readLines(system.file("extdata/goi/goi-somatic-tier1.txt", package="cargo"))
    }
    ## germline GOI
    if (is.null(germline_goi)) {
        germline_goi <- readLines(system.file("extdata/goi/goi-germline.txt", package="cargo"))
    }
    ## artifact genes
    if (is.null(somatic_art)) {
        somatic_art <- readLines(system.file("extdata/goi/art-somatic.txt", package="cargo"))
    }
    genes <- import(gene_model, feature.type="gene")
    if (drop_chrM) {
        genes <- dropSeqlevels(genes, "chrM", pruning.mode = "coarse")
    }
    tmp <- as.data.table(genes)
    vanno <- tmp[,.(
        chr=as.character(seqnames), start, end, strand=as.character(strand), gene_id, gene_name, anno_goi=TRUE,
        somatic_goi=gene_name %in% somatic_goi, germline_goi=gene_name %in% germline_goi,
        somatic_art=gene_name %in% somatic_art
        )]
    setkey(vanno, gene_id)
    return(vanno)
}

#' Update vault to newer dtVault format
#'
#' Note: input vault is modified
#'
#' @param vault legacy vault
#' @param anno dtVaultAnno object
#' @param format dtVault format [v2,v3]
#' @return dtVault
#'
#' @export
dtVaultUpdate <- function(vault, anno=NULL, format="v2") {
    STANDARD.CHR <- paste0("chr", c(1:22, "X", "Y", "M"))
    ## replace anno if provided
    if (!is.null(anno)) {
        vault$anno <- anno
    }
    ## fail if no vault$anno
    if (is.null(vault$anno)) {
        stop("No anno provided or present in vault.")
    }
    ## check if already converted
    if (!is.null(vault$format) && vault$format==format) {
        NULL
    } else if (!is.null(vault$format) && vault$format=="v2" && format=="v3") {
        #### Vault v2 to v3 update
        ##  - add 'block' in tables$segment
        ##  - add 'ccf' and 'm' in tables$somatic
        ## just fill in NAs so it works
        vault$tables$somatic <- rbind(vault0$tables$somatic, vault$tables$somatic, fill=TRUE)
        vault$tables$segment <- rbind(vault0$tables$segment, vault$tables$segment, fill=TRUE)
    } else {
        fn <- system.file(sprintf("extdata/schema/vault-%s.schema", format), package="cargo")
        vault0 <- eval(parse(fn))
        #### legacy vault
        ## remove ID columns
        if ("ID" %in% names(vault$tables$somatic)) {
            vault$tables$somatic[,ID:=NULL]
        }
        if ("ID" %in% names(vault$tables$germline)) {
            vault$tables$germline[,ID:=NULL]
        }
        ## remove qc (not supported yet)
        vault$qc <- NULL
        ## fix column types
        if (class(vault$tables$fusion$ctg.seq)=="list") {
            tmp <- sapply(vault$tables$fusion$ctg.seq, as.character)
            tmp[lengths(tmp)==0] <- NA_character_
            tmp <- unlist(tmp)
            vault$tables$fusion$ctg.seq <- tmp
        }
        for (tbl_name in c("germline", "somatic")) {
            for (col_name in c("rmsk_hit", "pon")) {
                col <- vault$tables[[tbl_name]][[col_name]]
                if (class(col) == "list") {
                    tmp <- sapply(col, paste, collapse=";")
                    tmp[tmp==""]  <- NA_character_
                    vault$tables[[tbl_name]][[col_name]] <- tmp
                }
            }
            ## add useful filtering variables
            tbl <- vault$tables[[tbl_name]]
            if (!is.null(tbl)) {
                tbl[, ":="(
                    is.indel = str_length(ref)!=str_length(alt),
                    sblr = abs(log2((ADT_FWD+1)/(ADT_REV+1))),
                    all.pon = ((gnomad_af==0) + (is.na(pon)) + (kg_af==0))
                )]
                vault$tables[[tbl_name]] <- tbl
            }
        }
        ## add filter column
        for (tbl_name in c("germline", "somatic", "fusion", "structural", "segment")) {
            tbl <- vault$tables[[tbl_name]]
            if (!is.null(tbl)) {
                if ("filter" %in% names(tbl)) {
                    tbl[,vcf_filter:=filter]
                    tbl[,filter:=NULL]
                } else if (tbl_name %in% c("germline", "somatic")) {
                    tbl[,vcf_filter:="."]
                }
                tbl[,triage:=NA_character_]
                # move filter to front
                tbl <- setcolorder(tbl, c(tail(names(tbl),n=1), head(names(tbl), -1)))
                vault$tables[[tbl_name]] <- tbl
            }
        }
        ## make gene id columns cosistent
        if (!is.null(vault$tables$germline)) {
            if (!("pid" %in% colnames(vault$tables$germline))) {
                vault$tables$germline$pid <- NA_character_
            }
            setnames(vault$tables$germline, "GENE", "gene_id")
            setnames(vault$tables$germline, "SYMBOL", "gene_name")
            suppressWarnings(vault$tables$germline[,":="(
                Tn=NULL,
                Nn=NULL,
                Taf90=NULL,
                Naf90=NULL
            )])
            setcolorder(vault$tables$germline, colnames(vault0$tables$germline))
        } else {
            vault$tables$germline <- vault0$tables$germline
        }
        if (!is.null(vault$tables$somatic)) {
            if (!("pid" %in% colnames(vault$tables$somatic))) {
                vault$tables$somatic$pid <- NA_character_
            }
            setnames(vault$tables$somatic, "GENE", "gene_id")
            setnames(vault$tables$somatic, "SYMBOL", "gene_name")
            suppressWarnings(vault$tables$somatic[,":="(
                Tn=NULL,
                Nn=NULL,
                Taf90=NULL,
                Naf90=NULL,
                ccf=NULL,
                m=NULL
            )])
            setcolorder(vault$tables$somatic, colnames(vault0$tables$somatic))
        } else {
            vault$tables$somatic <- vault0$tables$somatic
        }
        if (!is.null(vault$tables$gene.expression)) {
            setnames(vault$tables$gene.expression, "GENE", "gene_id")
            setnames(vault$tables$gene.expression, "SYMBOL", "gene_name")
            if (!("tcount" %in% colnames(vault$tables$gene.expression))) {
                vault$tables$gene.expression$tcount <- vault$tables$gene.expression$count
                vault$tables$gene.expression$tcpm <- vault$tables$gene.expression$cpm
                vault$tables$gene.expression$trpkm <- vault$tables$gene.expression$rpkm
                vault$tables$gene.expression$ncount <- NA_integer_
                vault$tables$gene.expression$ncpm <- NA_real_
                vault$tables$gene.expression$nrpkm <- NA_real_
                vault$tables$gene.expression$count <- NULL
                vault$tables$gene.expression$cpm <- NULL
                vault$tables$gene.expression$rpkm <- NULL
            }
            setcolorder(vault$tables$gene.expression, colnames(vault0$tables$gene.expression))
            setkey(vault$tables$gene.expression, id, gene_id)
        } else {
            vault$tables$gene.expression <- vault0$tables$gene.expression
        }
        if (!is.null(vault$tables$gene.copy)) {
            vault$tables$gene.copy$chr <- as.character(vault$tables$gene.copy$seqnames)
            vault$tables$gene.copy$seqnames <- NULL
            setnames(vault$tables$gene.copy, "GENE", "gene_id")
            setnames(vault$tables$gene.copy, "SYMBOL", "gene_name")
            setcolorder(vault$tables$gene.copy, colnames(vault0$tables$gene.copy))
            setkey(vault$tables$gene.copy, id, gene_id)
        } else {
            vault$tables$gene.copy <- vault0$tables$gene.copy
        }
        if (!is.null(vault$tables$structural)) {
            setnames(vault$tables$structural, "GENE1", "gene_id.1")
            setnames(vault$tables$structural, "GENE2", "gene_id.2")
            setnames(vault$tables$structural, "SYMBOL1", "gene_name.1")
            setnames(vault$tables$structural, "SYMBOL2", "gene_name.2")
            vault$tables$structural$chr1 <- as.character(vault$tables$structural$chr1)
            vault$tables$structural$chr2 <- as.character(vault$tables$structural$chr2)
            vault$tables$structural$STRAND1 <- as.character(vault$tables$structural$STRAND1)
            vault$tables$structural$STRAND2 <- as.character(vault$tables$structural$STRAND2)
        } else {
            vault$tables$structural <- vault0$tables$structural
        }
        if (!is.null(vault$tables$fusion)) {
            vault$tables$fusion$sid <- as.character(vault$tables$fusion$sid)
            vault$tables$fusion$chr.5 <- as.character(vault$tables$fusion$chr.5)
            vault$tables$fusion$chr.3 <- as.character(vault$tables$fusion$chr.3)
            vault$tables$fusion$str.5 <- as.character(vault$tables$fusion$str.5)
            vault$tables$fusion$str.3 <- as.character(vault$tables$fusion$str.3)
            vault$tables$fusion$art.5 <- as.character(vault$tables$fusion$art.5)
            vault$tables$fusion$art.3 <- as.character(vault$tables$fusion$art.3)
            vault$tables$fusion$dst <- as.character(vault$tables$fusion$dst)
            vault$tables$fusion$topo <- as.character(vault$tables$fusion$topo)
            vault$tables$fusion$l3 <- as.character(vault$tables$fusion$l3)
        } else {
            vault$tables$fusion <- vault0$tables$fusion
        }
        if (!is.null(vault$tables$segment)) {
            vault$tables$segment$chr <- as.character(vault$tables$segment$chr)
            vault$tables$segment$strand <- as.character(vault$tables$segment$strand)
            setkey(vault$tables$segment, id, var_id)
        } else {
            vault$tables$segment <- vault0$tables$segment
        }
        ## meta boolean NA
        vault$meta$tumor.dna.sample <- as.character(vault$meta$tumor.dna.sample)
        vault$meta$tumor.rna.sample <- as.character(vault$meta$tumor.rna.sample)
        vault$meta$normal.dna.sample <- as.character(vault$meta$normal.dna.sample)
        vault$meta$normal.rna.sample <- as.character(vault$meta$normal.rna.sample)
        setkey(vault$meta, id)
        ## add anno
        vault$anno <- anno
        ## add storage/format
        vault$storage <- "data.table"
        vault$format <- "v2"
    }
    return(vault)
}

#' Stack (rbind) list of dtVaults
#'
#' @param vaults list of dtVaults
#' @param anno dtAnno object [optional]
#' @return dtVault
#'
#' @export
dtVaultStack <- function(vaults, anno=NULL) {
    ## use versioned vault schema
    format <- unique(sapply(vaults, "[[", "format"))
    fn <- system.file(sprintf("extdata/schema/vault-%s.schema", format), package="cargo")
    vault0 <- eval(parse(fn))
    for (group in c("maps", "tables")) {
        for (table_name in names(vault0[[group]])) {
            table <- rbindlist(lapply(lapply(vaults, "[[", group), "[[", table_name), fill = T)
            vault0[[group]][[table_name]] <- rbind(vault0[[group]][[table_name]], table, fill=T)
        }
    }
    vault0$meta <- rbind(vault0$meta, rbindlist(lapply(vaults, "[[", "meta")), fill=T)
    if (!is.null(anno)) {
        vault0$anno <- anno
    } else {
        vault0$anno <- vaults[[1]]$anno
    }
    vault0$format <- format
    vault0$storage <- "data.table"
    vault0 <- dtVaultIndex(vault0)
    return(vault0)
}

#' Restore dtVault default indexes
#'
#' Note: this opperation occurs in-place
#'
#' @param dbv dtVault
#' @return dtVault
#'
#' @export
dtVaultIndex <- function(dtv) {
    ## retrieve schema
    schema.fn <- system.file(sprintf("extdata/schema/vault-%s.schema", dtv$format), package="cargo")
    dts <- eval(parse(schema.fn))
    ## collect all vault tables
    for (dt_name in names(dtv$tables)) {
        setkeyv(dtv$tables[[dt_name]], key(dts$tables[[dt_name]]))
    }
    for (dt_name in names(dtv$maps)) {
        setkeyv(dtv$map[[dt_name]], key(dts$map[[dt_name]]))
    }
    setkeyv(dtv$meta, key(dts$meta))
    setkeyv(dtv$anno, key(dts$anno))
    return(dtv)
}

.dtVaultSchema <- function(dtv) {
    #### get empty table
    dtv0 <- list()
    ## maps
    dtv0$maps  <- lapply(dtv$maps, function(dbt) {
        dbt0 <- dbt[0]
        attr(dbt0, ".internal.selfref") <- NULL
        dbt0
    })
    ## tables
    dtv0$tables  <- lapply(dtv$tables, function(dbt) {
        dbt0 <- dbt[0]
        attr(dbt0, ".internal.selfref") <- NULL
        dbt0
    })
    ## meta
    dtv0$meta <- dtv$meta[0]
    attr(dtv0$meta, ".internal.selfref") <- NULL
    ## anno
    dtv0$anno <- dtv$anno[0]
    attr(dtv0$anno, ".internal.selfref") <- NULL
    dtv0$storage <- character(0)
    dtv0$format <- character(0)
    dtv.schema <- capture.output(dput(dtv0))
    return(dtv.schema)
}
