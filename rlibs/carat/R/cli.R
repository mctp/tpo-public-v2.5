#' @export
reportVariants <- function() {

    option_list <- list(
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores"),
        optparse::make_option(c("-v", "--variant"), type="character",
                              default=NULL,
                              help="VCF variant call file"),
        optparse::make_option(c("-t", "--transcripts"), type="character",
                              default="generic",
                              help="transcript preference: [onco1500-v6a,generic]"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output RDS file"),
        optparse::make_option(c("-s", "--structural"), type="logical",
                              default=FALSE, action="store_true",
                              help="Toggle structural variant reporting")
    )
    parser <- optparse::OptionParser(
      "reportVariants() [options]",
      description=c("Report variants in a VCF file.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    if (
        is.null(args$variant) ||
        !file.exists(args$variant)
    ) {
        optparse::print_help(parser)
        write("Input file(s) not found.\n", stderr())
        quit("no", 1)
    }
    if (
        is.null(args$out) ||
        !dir.exists(dirname(args$out))
    ) {
        optparse::print_help(parser)
        write("Output file not provided or output directory not found.\n", stderr())
        quit("no", 1)
    }
    if (
        !(args$transcripts %in% c("onco1500-v6a", "generic"))
    ) {
        optparse::print_help(parser)
        write("Invalid transcripts preferences.\n", stderr())
        quit("no", 1)
    }

    gobj <- .gobj(args$genome)

    if (args$structural) {
        var1 <- importVcfReport(args$variant, gobj, context=FALSE, prune=FALSE, rename=FALSE)
        rep1 <- makeStructuralReport(var1, args$transcripts)
    } else {
        var1 <- importVcfReport(args$variant, gobj, context=TRUE, prune=TRUE, rename=TRUE)
        rep1 <- makeReport(var1, args$transcripts)
    }
    saveRDS(rep1, args$out)

}

#' @export
mnv <- function() {

    option_list <- list(
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=4L,
                              help="Number of cores"),
        optparse::make_option(c("-v", "--variant"), type="character",
                              default=NULL,
                              help="TNScope variant calls"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output VCF file name")

    )
    parser <- optparse::OptionParser(
      "mnv() [options]",
      description=c("Transform blocks of phased variants.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    if (
        is.null(args$variant) ||
        !file.exists(args$variant)
    ) {
        optparse::print_help(parser)
        write("Input file(s) not found.\n", stderr())
        quit("no", 1)
    }
    if (
        is.null(args$out) ||
        !dir.exists(dirname(args$out))
    ) {
        optparse::print_help(parser)
        write("Output file or directory not found.\n", stderr())
        quit("no", 1)
    }

    ## create Genome object
    gobj <- .gobj(args$genome)

    var <- importVcf(args$variant, gobj)
    registerDoParallel(cores = args$cores)
    vars <- phaseMNV(var, maxgap=1, gobj)
    stopImplicitCluster()
    concatVarsFromFile(vars, args$out)
}

#' @export
addManta <- function() {

    option_list <- list(
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=4L,
                              help="Number of cores"),
        optparse::make_option(c("-t", "--tnscope"), type="character",
                              default=NULL,
                              help="TNScope structural variant calls"),
        optparse::make_option(c("-m", "--manta"), type="character",
                              default=NULL,
                              help="Manta 'somaticsv' variant calls"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output VCF file name")

    )
    parser <- optparse::OptionParser(
      "addManta() [options]",
      description=c("Add Manta information.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    if (
        is.null(args$tnscope) ||
        !file.exists(args$tnscope) ||
        is.null(args$manta) ||
        !file.exists(args$manta)

    ) {
        optparse::print_help(parser)
        write("Input file(s) not found.\n", stderr())
        quit("no", 1)
    }
    if (
        is.null(args$out) ||
        !dir.exists(dirname(args$out))
    ) {
        optparse::print_help(parser)
        write("Output file or directory not found.\n", stderr())
        quit("no", 1)
    }
    oargs <- list(maxgap=100, sizemargin=0.5, restrictMarginToSizeMultiple=0.5)

    tn.var <- readVcf(args$tnscope, genome=args$genome)
    mn.var <- readVcf(args$manta, genome=args$genome)
    tn.var <- annotateVcfWithManta(tn.var, mn.var, overlap.args=oargs)
    writeVcf(tn.var, args$out)
}

#' @export
addMantaAndGridss <- function() {

    option_list <- list(
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=4L,
                              help="Number of cores"),
        optparse::make_option(c("-t", "--tnscope"), type="character",
                              default=NULL,
                              help="TNScope structural variant calls"),
        optparse::make_option(c("-m", "--manta"), type="character",
                              default=NULL,
                              help="Manta 'somaticsv' variant calls"),
        optparse::make_option(c("-d", "--gridss"), type="character",
                              default=NULL,
                              help="GRIDSS 'somatic' variant calls"),
        optparse::make_option(c("-g", "--genome"), type="character",
                              default="hg38",
                              help="Reference Genome [hg38,grch38,mm10,grcm38]"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output VCF file name")

    )
    parser <- optparse::OptionParser(
      "addMantaAndGridss() [options]",
      description=c("Add Manta information.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    if (
        is.null(args$tnscope) ||
        !file.exists(args$tnscope)
    ) {
        optparse::print_help(parser)
        write("Input TNSCOPE file not found.\n", stderr())
        quit("no", 1)
    }
    if (
        is.null(args$out) ||
        !dir.exists(dirname(args$out))
    ) {
        optparse::print_help(parser)
        write("Output file or directory not found.\n", stderr())
        quit("no", 1)
    }
    oargs <- list(maxgap=100, sizemargin=0.5, restrictMarginToSizeMultiple=0.5)

    tn.var <- readVcf(args$tnscope, genome=args$genome)
    if (!is.null(args$manta) && file.exists(args$manta)) {
        mn.var <- readVcf(args$manta, genome=args$genome)
        tn.var <- annotateVcfWithManta(tn.var, mn.var, overlap.args=oargs)
    }
    if (!is.null(args$gridss) && file.exists(args$gridss)) {
        gr.var <- readVcf(args$gridss, genome=args$genome)
        tn.var <- annotateVcfWithGridss(tn.var, gr.var, overlap.args=oargs)
    }

    writeVcf(tn.var, args$out)
}
