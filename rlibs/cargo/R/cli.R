#' @export
makeVault <- function() {
    option_list <- list(
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CARGO vault file"),
        optparse::make_option(c("-d", "--dump"), type="character",
                              default=NULL,
                              help="Output full CARGO import file"),
        optparse::make_option(c("-i", "--id"), type="character",
                              default=NULL,
                              help="unique id"),
        optparse::make_option(c("-a", "--case"), type="character",
                              default=NULL,
                              help="case id"),
        optparse::make_option(c("-t", "--cohort"), type="character",
                              default=NULL,
                              help="cohort name"),
        optparse::make_option(c("-r", "--carat"), type="character",
                              default=NULL,
                              help="carat-anno output (variants)"),
        optparse::make_option("--misc_t", type="character",
                              default=NULL,
                              help="cords-misc output tumor (QC, etc.)"),
        optparse::make_option("--misc_n", type="character",
                              default=NULL,
                              help="cords-misc output normal (QC, etc.)"),
        optparse::make_option(c("-c", "--cnvex"), type="character",
                              default=NULL,
                              help="cords-cnvex output (copy-number)"),
        optparse::make_option("--tquasr", type="character",
                              default=NULL,
                              help="crisp-quasr-tumor output (expression)"),
        optparse::make_option("--nquasr", type="character",
                              default=NULL,
                              help="crisp-quasr-normal output (expression)"),
        optparse::make_option(c("-f", "--codac"), type="character",
                              default=NULL,
                              help="crisp-codac output (fusion)"),
        optparse::make_option(c("-g", "--genes"), type="character",
                              default=NULL,
                              help="annotation GTF file"),
        optparse::make_option(c("--homo"), type="character",
                              default=NULL,
                              help="Homopolymers RDS file"),
        optparse::make_option(c("--targets"), type="character",
                              default=NULL,
                              help="Capture targets BED file"),
        optparse::make_option(c("-s", "--settings"), type="character",
                              default="default",
                              help="Settings Preset file [{preset_name},<file>]"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=detectCores(),
                              help="Number of cores")
    )
    parser <- optparse::OptionParser(
      "makeVault() [options]",
      description=c("Create CARGO vault.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2022\n"),
      option_list=option_list
      )
    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    ## Input GTF
    if (
        !file.exists(args$genes)
    ) {
        optparse::print_help(parser)
        write("Gene GTF input file(s) not found.\n", stderr())
        quit("no", 1)
    }
    ## Output
    if (
        is.null(args$out) ||
        !dir.exists(dirname(args$out))
    ) {
        optparse::print_help(parser)
        write("Output file not provided or output directory not found.\n", stderr())
        quit("no", 1)
    }

    ## Setttings
    if (!is.null(args$settings) && !file.exists(args$settings)) {
        args$settings <- system.file(sprintf("extdata/settings/%s.R", args$settings), package="cargo")
        if (!file.exists(args$settings)) {
            optparse::print_help(parser)
            write("Settings preset or file not provided.\n", stderr())
            quit("no", 1)
        }
    }

    ## Config
    .args <- list(
        meta=list(group=list(id=args$id, case=args$case, cohort=args$cohort, args$targets)),
        data=list(carat=args$carat, misc=args$misc, cnvex=args$cnvex, tquasr=args$tquasr,
            nquasr=args$nquasr, codac=args$codac,misc_t=args$misc_t,misc_n=args$misc_n),
        annotation=list(
            gene.model=args$genes,
            homopolymers=args$homo
        )
    )
    .vcfg <- cargoConfig(args$settings, .args)

    ## Data
    .vdir <- vaultDir(.vcfg)
    .vdump <- vaultDump(.vdir)
    if (!is.null(args$dump)) {
        saveRDS(.vdump, args$dump)
    }
    .vanno <- vaultAnno(.vcfg)
    .vdata <- vaultData(.vdump, .vanno, .vcfg)
    .vmeta <- vaultMeta(.vdump, .vcfg)
    .qc <- vaultQC(.vdump)
    .v <- vault(.vdata, .vanno, .vmeta, .qc)
    saveRDS(.v, args$out)
}
