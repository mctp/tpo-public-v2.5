#' @export
process <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CNVEX file"),
        optparse::make_option(c("-t", "--tumor"), type="character",
                              default=NULL,
                              help="Tumor BAM/CRAM file(s) [required]"),
        optparse::make_option(c("-n", "--normal"), type="character",
                              default=NULL,
                              help="Normal BAM/CRAM file(s)"),
        optparse::make_option(c("-e", "--tvcf"), type="character",
                              default=NULL,
                              help="Target Somatic VCF file"),
        optparse::make_option(c("-w", "--gvcf"), type="character",
                              default=NULL,
                              help="Genome Somatic VCF file"),
        optparse::make_option(c("-l", "--pool"), type="character",
                              default=NULL,
                              help="Pool of normals"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default=NULL,
                              help="Settings Preset file [{exome/panel}-{pair/pool},<file>]"),
        optparse::make_option(c("-c", "--capture"), type="character",
                              default=NULL,
                              help="Capture target file [bed]"),
        optparse::make_option(c("-f", "--fasta"), type="character",
                              default=NULL,
                              help="Genome FASTA file for CRAM [fa]"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]"),
        optparse::make_option(c("-v", "--caller"), type="character",
                              default="sentieon",
                              help="Variant caller for VCF [sentieon]"),
        optparse::make_option(c("-p", "--popaf"), type="character",
                              default=NULL,
                              help="Population allele frequency  [KG_AF,GNOMAD_AF]"),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=4L,
                              help="Number of parallel threads"),
        optparse::make_option(c("-x", "--sex"), type="character",
                              default=NULL,
                              help="Patients sex [male,female]")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'cnvex::process()' [options]",
      description=c("Load BAM/CRAM/VCF and optional pool to create CNVEX object.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Required output file name not provided.\n", stderr())
        quit("no", 1)
    }
    if (!grepl("\\.rds$", args$out)) {
        optparse::print_help(parser)
        write("Output file does not end with '.rds'.\n", stderr())
        quit("no", 1)
    }

    ## Tumor
    if (!is.null(args$tumor)) {
        args$tumor <- str_split(str_remove(args$tumor, "[,;]$"), ",|;")[[1]]
        if (!all(file.exists(args$tumor))) {
            optparse::print_help(parser)
            write("Input tumor file(s) provided but do not exist.\n", stderr())
            quit("no", 1)
        }
    }

    ## Normal
    if (!is.null(args$normal)) {
        args$normal <- str_split(str_remove(args$normal, "[,;]$"), ",|;")[[1]]
        if (!all(file.exists(args$normal))) {
            optparse::print_help(parser)
            write("Input normal file(s) provided but do not exist.\n", stderr())
            quit("no", 1)
        }
    }

    ## VCF
    if (!is.null(c(args$gvcf, args$tvcf))) {
        if (!all(file.exists(c(args$gvcf, args$tvcf)))) {
            optparse::print_help(parser)
            write("Input VCF file(s) provided but do not exist.\n", stderr())
            quit("no", 1)
        }
    }

    ## Pool
    if (!is.null(args$pool) && !file.exists(args$pool)) {
        optparse::print_help(parser)
        write("Pool file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Setttings
    if (!is.null(args$settings) && !file.exists(args$settings)) {
        args$settings <- system.file(sprintf("extdata/settings/%s.R", args$settings), package="cnvex")
        if (!file.exists(args$settings)) {
            optparse::print_help(parser)
            write("Settings preset or file not provided.\n", stderr())
            quit("no", 1)
        }
    }

    ## Capture
    if (!is.null(args$capture) && !file.exists(args$capture) && !is.null(args$vcf["target"])) {
        optparse::print_help(parser)
        write("Target VCF provided without capture file.\n", stderr())
        quit("no", 1)
    }

    ## FASTA
    if (!is.null(args$fasta) && !file.exists(args$fasta)) {
        optparse::print_help(parser)
        write("Input FASTA file provided but does not exist.\n", stderr())
        quit("no", 1)
    }

    ## Input Combinations
    if (!is.null(args$inp) && (!is.null(args$tumor) || !is.null(args$normal) || !is.null(args$tvcf) || !is.null(args$gvcf))) {
        optparse::print_help(parser)
        write("Unsupported combination of input files (--inp and --tumor/--normal/--tvcf/--gvcf).\n", stderr())
        quit("no", 1)
    }
    if (is.null(args$inp) && is.null(args$tumor) && is.null(args$normal)) {
        optparse::print_help(parser)
        write("Required input files not provided (--inp or --tumor or --normal).\n", stderr())
        quit("no", 1)
    }

    ## valid sex
    if (!is.null(args$sex)) {
        if (!(args$sex %in% c("male","female"))) {
            optparse::print_help(parser)
            write("Provided sex should be either male or female.\n", stderr())
            quit("no", 1)
      }
    }

    ## get options object
    settings <- args$settings
    args$settings <- NULL
    capture <- args$capture
    args$capture <- NULL
    args$vcf <- c(genome=args$gvcf, target=args$tvcf)
    opts <- getOpts(settings, args)

    ## Pool
    if (is.null(args$pool) && (opts$lr.tumor == "pool" || opts$lr.normal == "pool" || opts$prune.nvar)) {
        optparse::print_help(parser)
        write("Selected 'lr.tumor', 'lr.normal' or 'prune.nvar' options require a pool.\n", stderr())
        quit("no", 1)
    }

    ## import pool
    if (!is.null(args$pool)) {
        pool <- readRDS(args$pool)
        if ((pool$method=="pca" && opts$pool.method %in% c("ica", "ica-nosex")) ||
            (pool$method=="ica" && opts$pool.method %in% c("pca", "pca-nosex"))
            ) {
            write("Provided 'pool' not compatible with options 'pool.method'.\n", stderr())
            quit("no", 1)
        }
    } else {
        pool <- NULL
    }

    ## create Genome object
    gobj <- getGobj(args$genome, args$fasta, TRUE)

    ## import data
    if (is.null(args$inp)) {
        cnv <- createCnv(args$vcf, capture, args$tumor, args$normal, gobj, opts)
    } else {
        cnv <- readRDS(args$inp)
        cnv <- upgradeCnv(cnv, opts)
    }

    ## normals for pair normalization
    if ((all(is.na(cnv$tile$n.cov.raw)) && opts$lr.tumor %in% c("pair"))) {
        optparse::print_help(parser)
        write("Required normal sample not provided for 'lr.tumor' option pair.\n", stderr())
        quit("no", 1)
    }

    ## options properties
    opts.out <- str_replace(args$out, "\\.rds$", "-opts.rds")
    opts.save <- copy(opts)
    opts.save$inp <- NULL
    opts.save$out <- NULL
    saveRDS(opts.save, opts.out)

    ## update CNV object
    cnv <- addHqTile(cnv, opts)
    cnv <- addPassVariant(cnv, opts)
    cnv <- addSex(cnv, gobj, opts)
    cnv <- addNormCoverage(cnv, opts)
    cnv <- addLogRatio(cnv, pool, opts)

    ## save
    saveRDS(cnv, args$out)
}

#' @export
segment <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output file"),
        optparse::make_option(c("-l", "--pool"), type="character",
                              default=NULL,
                              help="Pool of normals"),
        optparse::make_option(c("-a", "--sample"), type="character",
                              default="tumor",
                              help="Sample to analyze [tumor,normal]"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default=NULL,
                              help="Settings Preset file [{exome/panel}-{pair/pool},<file>]"),
        optparse::make_option(c("-j", "--cores"), type="integer",
                            default=detectCores(),
                            help="Number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'cnvex::segment()' [options]",
      description=c("Segment copy number data.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (is.null(args$inp)) {
        optparse::print_help(parser)
        write("Input file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Output file name not provided.\n", stderr())
        quit("no", 1)
    }

    ## Pool
    if (!is.null(args$pool) && !file.exists(args$pool)) {
        optparse::print_help(parser)
        write("Pool file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Setttings
    if (is.null(args$settings)) {
        opts.rds <- str_replace(args$inp, "\\.rds$", "-opts.rds")
        opts <- readRDS(opts.rds)
    } else {
        if (!file.exists(args$settings)) {
            args$settings <- system.file(sprintf("extdata/settings/%s.R", args$settings), package="cnvex")
        }
        if (!file.exists(args$settings)) {
            optparse::print_help(parser)
            write("Settings not found or not provided.\n", stderr())
            quit("no", 1)
        }
        opts <- getOpts(args$settings, list())
    }

    ## Pool
    if (is.null(args$pool) && (opts$prune.nvar)) {
        optparse::print_help(parser)
        write("Selected 'prune.nvar' option requires a pool.\n", stderr())
        quit("no", 1)
    }

    ## import pool
    if (!is.null(args$pool)) {
        pool <- readRDS(args$pool)
        if ((pool$method=="pca" && opts$pool.method %in% c("ica", "ica-nosex")) ||
            (pool$method=="ica" && opts$pool.method %in% c("pca", "pca-nosex"))
            ) {
            write("Provided 'pool' not compatible with options 'pool.method'.\n", stderr())
            quit("no", 1)
        }
    } else {
        pool <- NULL
    }

    ## read CNVEX object
    cnv <- readRDS(args$inp)

    ## create Genome object
    gobj <- getGobj(unique(genome(cnv$tile)), NULL, FALSE)

    ## segmentation
    registerDoParallel(cores = args$cores)
    mcnv <- modelCnv(args$sample, cnv, pool, gobj, opts)
    seg <- modelSeg(mcnv, pool, opts)
    stopImplicitCluster()

    ## save CNVEX OPT object
    saveRDS(seg, args$out)
}

#' @export
modelsearch <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output file"),
        optparse::make_option(c("-l", "--pool"), type="character",
                              default=NULL,
                              help="Pool of normals"),
        optparse::make_option(c("-a", "--sample"), type="character",
                              default="tumor",
                              help="Sample to analyze [tumor,normal]"),
        optparse::make_option(c("-e", "--segment"), type="character",
                              default=NULL,
                              help="Segmentation file"),
        optparse::make_option(c("-r", "--nogrid"), action="store_true",
                              default=FALSE,
                              help="Disable grid-search"),
        optparse::make_option(c("-f", "--nofine"), action="store_true",
                              default=FALSE,
                              help="Disable fine-search"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default=NULL,
                              help="Settings Preset file [{exome/panel}-{pair/pool},<file>]"),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'cnvex::modelsearch()' [options]",
      description=c("Search absolute copy number models.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (is.null(args$inp)) {
        optparse::print_help(parser)
        write("Input file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Output file name not provided.\n", stderr())
        quit("no", 1)
    }

    ## Pool
    if (!is.null(args$pool) && !file.exists(args$pool)) {
        optparse::print_help(parser)
        write("Pool file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Setttings
    if (is.null(args$settings)) {
        opts.rds <- str_replace(args$inp, "\\.rds$", "-opts.rds")
        opts <- readRDS(opts.rds)
    } else {
        if (!file.exists(args$settings)) {
            args$settings <- system.file(sprintf("extdata/settings/%s.R", args$settings), package="cnvex")
        }
        if (!file.exists(args$settings)) {
            optparse::print_help(parser)
            write("Settings not found or not provided.\n", stderr())
            quit("no", 1)
        }
        opts <- getOpts(args$settings, list())
    }

    ## Segment
    if (is.null(args$segment) || !file.exists(args$segment)) {
        optparse::print_help(parser)
        write("Segment file not provided or not found.\n", stderr())
        quit("no", 1)
    }

    ## get options object
    opts <- getOpts(args$settings, list())

    ## Pool
    if (is.null(args$pool) && (opts$prune.nvar)) {
        optparse::print_help(parser)
        write("Selected 'prune.nvar' option requires a pool.\n", stderr())
        quit("no", 1)
    }

    ## import pool
    if (!is.null(args$pool)) {
        pool <- readRDS(args$pool)
        if ((pool$method=="pca" && opts$pool.method %in% c("ica", "ica-nosex")) ||
            (pool$method=="ica" && opts$pool.method %in% c("pca", "pca-nosex"))
            ) {
            write("Provided 'pool' not compatible with options 'pool.method'.\n", stderr())
            quit("no", 1)
        }
    } else {
        pool <- NULL
    }

    ## read CNVEX object
    cnv <- readRDS(args$inp)
    seg <- readRDS(args$segment)

    ## create Genome object
    gobj <- getGobj(unique(genome(cnv$tile)), NULL, FALSE)

    ## segmentation and model-fitting
    registerDoParallel(cores = args$cores)
    mcnv <- modelCnv(args$sample, cnv, pool, gobj, opts)
    opt <- modelOpt(mcnv, seg, args$nogrid, args$nofine, opts)
    stopImplicitCluster()

    ## save CNVEX OPT object
    saveRDS(opt, args$out)
}

#' @export
digest <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output file"),
        optparse::make_option(c("-l", "--pool"), type="character",
                              default=NULL,
                              help="Pool of normals"),
        optparse::make_option(c("-a", "--sample"), type="character",
                              default="tumor",
                              help="Sample to analyze [tumor,normal]"),
        optparse::make_option(c("-m", "--model"), type="character",
                              default=NULL,
                              help="Sample model"),
        optparse::make_option(c("-e", "--segment"), type="character",
                              default=NULL,
                              help="Sample segment"),
        optparse::make_option(c("-p", "--pick"), type="character",
                              default=NULL,
                              help="Model pick"),
        optparse::make_option(c("-u", "--user"), type="character",
                              default=NULL,
                              help="Sign-off user"),
        optparse::make_option(c("-n", "--notes"), type="character",
                              default=NULL,
                              help="Sign-off notes"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default=NULL,
                              help="Settings Preset file [{exome/panel}-{pair/pool},<file>]"),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores")
    )
    parser <- optparse::OptionParser(
        "Rscript -e 'cnvex::digest()' [options]",
        description=c("Create CNVEX digest.\n"),
        epilogue=c(
            "Michigan Center for Translational Pathology (c) 2022\n"),
        option_list=option_list
        )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (is.null(args$inp)) {
        optparse::print_help(parser)
        write("Input file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Required output file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!dir.exists(dirname(args$out))) {
        optparse::print_help(parser)
        write("File directory does not exist.\n", stderr())
        quit("no", 1)
    }

    ## Pool
    if (!is.null(args$pool) && !file.exists(args$pool)) {
        optparse::print_help(parser)
        write("Pool file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Setttings
    if (is.null(args$settings)) {
        args$settings <- str_replace(args$inp, "\\.rds$", "-opts.rds")
    } else {
        if (!file.exists(args$settings)) {
            args$settings <- system.file(sprintf("extdata/settings/%s.R", args$settings), package="cnvex")
        }
    }
    if (!file.exists(args$settings)) {
        optparse::print_help(parser)
        write("Settings not found or not provided.\n", stderr())
        quit("no", 1)
    }

    ## get options object
    opts <- getOpts(args$settings, list())

    ## Pool
    if (is.null(args$pool) && (opts$prune.nvar)) {
        optparse::print_help(parser)
        write("Selected 'prune.nvar' option requires a pool.\n", stderr())
        quit("no", 1)
    }

    ## import pool
    if (!is.null(args$pool)) {
        pool <- readRDS(args$pool)
        if ((pool$method=="pca" && opts$pool.method %in% c("ica", "ica-nosex")) ||
            (pool$method=="ica" && opts$pool.method %in% c("pca", "pca-nosex"))
            ) {
            write("Provided 'pool' not compatible with options 'pool.method'.\n", stderr())
            quit("no", 1)
        }
    } else {
        pool <- NULL
    }

    ## Segment
    if (is.null(args$segment) || !file.exists(args$segment)) {
        optparse::print_help(parser)
        write("Segment file not provided or not found.\n", stderr())
        quit("no", 1)
    }

    ## pick
    if (str_detect(args$pick, ":")) {
        tmp <- str_split(args$pick, ":")
        purity <- as.numeric(tmp[[1]][1])
        ploidy <- as.numeric(tmp[[1]][2])
        eval <- CNVEX_EMPTY_EVAL
        eval$p <- purity
        eval$P <- ploidy
        topn <- 1L
        top0 <- TRUE
    } else if (startsWith(tolower(args$pick), "top") && !is.null(args$model)) {
        opt <- readRDS(args$model)
        eval <- opt$eval
        topn <- as.integer(str_replace(tolower(args$pick), "top", ""))
        top0 <- FALSE
    } else {
        optparse::print_help(parser)
        write("Unsupported pick setting.\n", stderr())
        quit("no", 1)
    }

    ## read CNVEX objects
    cnv <- readRDS(args$inp)
    seg <- readRDS(args$segment)

    ## create Genome object
    gobj <- getGobj(unique(genome(cnv$tile)), NULL, FALSE)
    log <- list(user=args$user, notes=args$notes)

    ## segmentation and model-fitting
    registerDoParallel(cores = args$cores)
    mcnv <- modelCnv(args$sample, cnv, pool, gobj, opts)
    mods <- modelGetEval(mcnv, seg, eval, topn, opts)
    foreach(i=seq_len(length(mods))) %dopar% {
        i_mod <- mods[[i]]
        i_eval <- as.list(eval[i])
        i_digest <- modelDigest(mcnv, seg, i_mod, i_eval$p, i_eval$P, i_eval, opts, log)
        i_pick <- ifelse(top0, i - 1, i)
        i_out <- paste0(args$out, sprintf("-digest-%02d.rds", i_pick))
        saveRDS(i_digest, i_out)
    }
    stopImplicitCluster()
}

#' @export
ggplots <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX digest file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output prefix"),
        optparse::make_option(c("-t", "--type"), type="character",
                              default="lr-baf",
                              help="Plot type"),
        optparse::make_option(c("-s", "--ggopts"), type="character",
                              default="gg-default",
                              help="GGplot preset file [gg-default]"),
        optparse::make_option(c("-r", "--override"), type="character",
                              default=NULL,
                              help="Preset overrides."),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores")
    )
    parser <- optparse::OptionParser(
        "Rscript -e 'cnvex::ggplots()' [options]",
        description=c("Plot CNVEX digest file.\n"),
        epilogue=c(
            "Michigan Center for Translational Pathology (c) 2022\n"),
        option_list=option_list
        )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (is.null(args$inp)) {
        optparse::print_help(parser)
        write("Input file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        args$inp <- list.files(dirname(args$inp), paste0(basename(args$inp), ".*rds"), full.names = TRUE)
        if (length(args$inp)==0) {
            optparse::print_help(parser)
            write("Input file prefix provided but files not found.\n", stderr())
            quit("no", 1)
        }
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Required output prefix not provided.\n", stderr())
        quit("no", 1)
    }
    if (!dir.exists(dirname(args$out))) {
        optparse::print_help(parser)
        write("Prefix directory does not exist.\n", stderr())
        quit("no", 1)
    }

    ## ggplot settings
    if (!is.null(args$ggopts) && !file.exists(args$ggopts)) {
        args$ggopts <- system.file(sprintf("extdata/plot/%s.R", args$ggopts), package="cnvex")
        if (!file.exists(args$ggopts)) {
            optparse::print_help(parser)
            write("GGplot preset or file not found.\n", stderr())
            quit("no", 1)
        }
    }

    ## get overrides
    over <- .parse.overrides(args$override)
    ## get gg.opts
    gg.opts <- getGgOpts(args$ggopts, over)
    ## get genes
    genes <- GRanges()

    registerDoParallel(cores = args$options$cores)
    foreach(inp=args$inp) %dopar% {
        ## read CNVEX digest object
        library(ggplot2)
        digest <- readRDS(inp)
        if (args$type=="arranged-plot") {
            plt <- ggArrangedPlot(digest, genes, gg.opts)
            out <- str_replace(str_replace(inp, "-digest-", "-arranged-"), ".rds", ".pdf")
            ggsave(out, plt)
        }
    }
}

#' @export
kpplots <- function() {

    option_list <- list(
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input ONE CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output prefix"),
        optparse::make_option(c("-t", "--type"), type="character",
                              default="lr-baf",
                              help="Plot type"),
        optparse::make_option(c("-k", "--kpopts"), type="character",
                              default="kp-default",
                              help="KaryoploteR preset file [kp-default]"),
        optparse::make_option(c("-r", "--override"), type="character",
                              default=NULL,
                              help="Preset overrides."),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores")
    )
    parser <- optparse::OptionParser(
        "Rscript -e 'cnvex::kpplots()' [options]",
        description=c("Plot CNVEX digest file.\n"),
        epilogue=c(
            "Michigan Center for Translational Pathology (c) 2022\n"),
        option_list=option_list
        )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)

    ## Input
    if (is.null(args$inp)) {
        optparse::print_help(parser)
        write("Input file not provided.\n", stderr())
        quit("no", 1)
    }
    if (!is.null(args$inp) && !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file provided but not found.\n", stderr())
        quit("no", 1)
    }

    ## Output
    if (is.null(args$out)) {
        optparse::print_help(parser)
        write("Required output prefix not provided.\n", stderr())
        quit("no", 1)
    }
    if (!dir.exists(dirname(args$out))) {
        optparse::print_help(parser)
        write("Prefix directory does not exist.\n", stderr())
        quit("no", 1)
    }

    ## KaryoploteR settings
    if (!is.null(args$kpopts) && !file.exists(args$kpopts)) {
        args$kpopts <- system.file(sprintf("extdata/plot/%s.R", args$kpopts), package="cnvex")
        if (!file.exists(args$kpopts)) {
            optparse::print_help(parser)
            write("KaryoploteR preset or file not found.\n", stderr())
            quit("no", 1)
        }
    }

    ## get overrides
    over <- .parse.overrides(args$override)

    ## get kp.opts
    kp.opts <- getKpOpts(args$kpopts, over)

    ## read CNVEX digest object
    digest <- readRDS(args$inp)
    if ("lr-baf" %in% args$type) {
        pdf(args$out, width=10, height=3)
        kp <- kpInit(genome=unique(genome(digest$tile)), kp.opts=kp.opts)
        kpLrBafPlot(digest, NULL, kp, kp.opts)
        invisible(dev.off())
    }
}


#' @export
pool <- function() {

    # arguments
    option_list <- list(
        optparse::make_option(c("-o", "--out"), type="character",
                              default="./",
                              help="Output folder address"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default=NULL,
                              help="Settings Preset file [{exome/panel}-{pair/pool},<file>]"),
        optparse::make_option(c("-j", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores"),
        optparse::make_option(c("-f", "--file"), type="character",
                              default=NULL,
                              help="a file containing all the sample addresses"),
        optparse::make_option(c("-d", "--pd"), type="character",
                              default=NULL,
                              help="Pool Data (pd) file"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'cnvex::pool()' [options] files",
      description=c("Create pool of normals CNVEX object.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=TRUE)


    ## Setttings
    if (!is.null(args$options$settings) && !file.exists(args$options$settings)) {
      args$options$settings <- system.file(sprintf("extdata/settings/%s.R", args$options$settings), package="cnvex")
      if (!file.exists(args$options$settings)) {
        optparse::print_help(parser)
        write("Settings preset or file not provided.\n", stderr())
        quit("no", 1)
      }
    }

    ## Multiple file sources
    if ((length(args$args) > 1) && !is.null(args$options$file)) {
      optparse::print_help(parser)
      write("Can not have both --file and files address list, only one should be provided.\n", stderr())
      quit("no", 1)
    }

    ## Not enough files
    if ((length(args$args) < 2) && is.null(args$options$file) && is.null(args$options$pd)) {
      optparse::print_help(parser)
      write("At least 2 CNVEX input files required.\n", stderr())
      quit("no", 1)
    }

    ## File not exist
    if (!is.null(args$options$file) && !file.exists(args$options$file)) {
      optparse::print_help(parser)
      write("Provided file address does not exist.\n", stderr())
      quit("no", 1)
    }

    ## Read the addresses
    if (is.null(args$options$file)) {
      cnv.fns <- as.character(args$args)
      args$args <- NULL
    } else {
      cnv.fns <- as.character(read.csv(args$options$file,header = FALSE)[,1])
      args$options$file <- NULL
      ## File have less than 2 addresses
      if (length(cnv.fns) < 2) {
        optparse::print_help(parser)
        write("At least 2 CNVEX input files required.\n", stderr())
        quit("no", 1)
      }
    }

    genome <- args$options$genome
    settings <- args$options$settings
    args$options$settings <- NULL
    opts <- getOpts(settings, args)

    gobj <- getGobj(genome,NULL,refs=TRUE)
    registerDoParallel(cores = args$options$cores)
    if (is.null(args$options$pd)) {
      pd <- createPoolData(cnv.fns, gobj, opts)
    } else {
      pd <- readRDS(args$options$pd)
    }
    pool <- createPool(pd, opts)
    stopImplicitCluster()

    saveRDS(pd, paste(args$options$out,"pd.rds",sep = ""))
    saveRDS(pool, paste(args$options$out,"pool.rds",sep = ""))
}
