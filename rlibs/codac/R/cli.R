#' @export
annotation <- function() {
    option_list <- list(
      optparse::make_option(c("-g", "--genome"), type="character",
                            default="hg38",
                            help="set genome version: hg38"),
      optparse::make_option(c("-j", "--cores"), type="integer",
                            default=detectCores(),
                            help="set the number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'library(methods);codac::annotation()' [options] gtf_file output_file",
      description=c("Create annotation file.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2019\n"),
      option_list=option_list
    )
    opt <- optparse::parse_args(parser, positional_arguments=TRUE)

    options(mc.cores=opt$options$cores)

    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("input_directory or output_file are missing.\n", stderr())
        quit("no", 1)
    }

    gtf.fn <- opt$args[1]
    out.fn <- opt$args[2]

    ann <- makeAnnotations(gtf.fn, opt$options$genome)

    saveRDS(ann, out.fn)
}

#' @export
run <- function() {
    option_list <- list(
      optparse::make_option(c("-a", "--annotation"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-c", "--config"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-j", "--cores"), type="integer",
                            default=detectCores(),
                            help="set the number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'library(methods);codac::run()' [options] input_directory output_prefix",
      description=c("Run CODAC to detect all types of chimeric RNAs.\n"),
      epilogue=c(
          "Michigan Center for Translational Pathology (c) 2019\n"),
      option_list=option_list
    )
    opt <- optparse::parse_args(parser, positional_arguments=TRUE)

    options(mc.cores=opt$options$cores)

    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("input_directory or output_prefix are missing.\n", stderr())
        quit("no", 1)
    }

    inp.pth <- opt$args[1]
    out.pfx <- opt$args[2]

    dir <- makeDirectory(inp.pth)
    ann <- readRDS(opt$options$annotation)
    par <- makeParams(opt$options$config)

    run <- runCODAC(dir, ann, par)

    saveRDS(run, paste0(out.pfx, "-run.rds"))
}

#' @export
callsv <- function() {
    option_list <- list(
      optparse::make_option(c("-a", "--annotation"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-c", "--config"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-g", "--gmap"), type="character",
                            default=NULL,
                            help="Location of gmap index"),
      optparse::make_option(c("-m", "--minimap"), type="character",
                            default=NULL,
                            help="Location of minimap2 index"),
      optparse::make_option(c("-j", "--cores"), type="integer",
                            default=detectCores(),
                            help="set the number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'library(methods);codac::callsv()' [options] codac_file output_prefix",
      description=c("Call structural variants.\n"),
      epilogue=c(
        "Michigan Center for Translational Pathology (c) 2019\n"),
      option_list=option_list
    )
    opt <- optparse::parse_args(parser, positional_arguments=TRUE)

    options(mc.cores=opt$options$cores)

    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("input_directory or output_prefix are missing.\n", stderr())
        quit("no", 1)
    }

    inp.pth <- opt$args[1]
    out.pfx <- opt$args[2]

    ann <- readRDS(opt$options$annotation)
    par <- makeParams(opt$options$config)

    run <- readRDS(inp.pth)
    sv.bun <- svBundle(run$bun, ann)
    sv.bun <- assembleBreakpoints(sv.bun, ann, par)
    sv.bun <- alignBreakpoints(sv.bun, ann, par, opt$options$minimap, opt$options$gmap)
    sv.bun <- validateBreakpoints(sv.bun, ann, par)
    sv.rep <- svReport(sv.bun, run$spl, ann, par)
    sv.fmt <- svFormat(sv.rep)

    ##
    saveRDS(sv.bun, paste0(out.pfx, "-bun-sv.rds"))
    saveRDS(sv.rep, paste0(out.pfx, "-rep-sv.rds"))
    fwrite(sv.fmt, paste0(out.pfx, "-rep-sv.csv"))

}

#' @export
stat <- function() {
    option_list <- list(
      optparse::make_option(c("-a", "--annotation"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-c", "--config"), type="character",
                            default="built-in",
                            help="Algorithm configuration file"),
      optparse::make_option(c("-j", "--cores"), type="integer",
                            default=detectCores(),
                            help="set the number of cores")
    )
    parser <- optparse::OptionParser(
      "Rscript -e 'library(methods);codac::stat()' [options] codac_file output_prefix",
      description=c("Calculate Chimeric Statistics.\n"),
      epilogue=c(
        "Michigan Center for Translational Pathology (c) 2019\n"),
      option_list=option_list
    )
    opt <- optparse::parse_args(parser, positional_arguments=TRUE)

    options(mc.cores=opt$options$cores)

    if (length(opt$args) < 2) {
        optparse::print_help(parser)
        write("input_directory or output_prefix are missing.\n", stderr())
        quit("no", 1)
    }

    inp.pth <- opt$args[1]
    out.pfx <- opt$args[2]

    ann <- readRDS(opt$options$annotation)
    par <- makeParams(opt$options$config)

    run <- readRDS(inp.pth)
    stat.rep <- statReport(run$bun, run$spl, run$bam, ann, par)

    ##
    saveRDS(stat.rep, paste0(out.pfx, "-stat-rep.rds"))

}
