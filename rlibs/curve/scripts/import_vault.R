#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
getCurrentFile <- function() {
    this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
    return(this_file)
}
TPO_RLIBS <- dirname(dirname(dirname(getCurrentFile())))
devtools::load_all(file.path(TPO_RLIBS, "tpolib"), quiet=TRUE)
devtools::load_all(file.path(TPO_RLIBS, "carat"), quiet=TRUE)
devtools::load_all(file.path(TPO_RLIBS, "cnvex"), quiet=TRUE)
devtools::load_all(file.path(TPO_RLIBS, "codac"), quiet=TRUE)
devtools::load_all(file.path(TPO_RLIBS, "cargo"), quiet=TRUE)
devtools::load_all(file.path(TPO_RLIBS, "curve"), quiet=TRUE)
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(RPostgreSQL))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

#parse arguments
option_list = list(
    optparse::make_option(c("-v", "--vault"), type="character",
                          help="file path of TPO vault"),
    optparse::make_option(c("-c", "--config"), type="character",
                          default='./rlibs/cargo/inst/extdata/settings/default.R',
                          help="cargo config file"),
    optparse::make_option(c("-g", "--gtf"), type="character",
                          default="grch38.108.clean.gtf",
                          help="File path of GTF Gene Model"),
    optparse::make_option(c("-f", "--fset"), type="character",
                          default=NULL,
                          help="fset passed to vaultTriage")
)
parser = optparse::OptionParser(
  "Rscript import_vault.R [options]",
  description=c("import a dtVault to a curve database"),
  epilogue=c(
      "Michigan Center for Translational Pathology (c) 2022\n"),
  option_list=option_list
  )

args <- optparse::parse_args(parser, positional_arguments=FALSE)

#populate the gxp_genes_table if needed
init_gxp_genes(args$gtf)

# load
stopifnot(file.exists(args$vault))
vault <- readRDS(args$vault)

#check
dxVaultCheck(vault)

#triage
vault.triage <- vaultTriage(vault, args$fset)

#loop through cases and import
ids <- vault.triage$meta$id
for(a in ids){
  print(sprintf("importing %s",a))
  vs <- dxVaultSubsetById(vault.triage, a)
  #add the grouping and metadata
  gid <- import_meta(vs)

  #import the data
  if(nrow(vs$tables$somatic)>0){
    import_variants(vs,'somatic',gid)
  }
  if(nrow(vs$tables$germline)>0){
    import_variants(vs,'germline',gid)
  }
  if(nrow(vs$tables$structural)>0){
    import_structural(vs,gid)
  }
  if(nrow(vs$tables$segment)>0){
    import_cnv(vs,gid,type='som')
  }
  if(nrow(vs$tables$gene.expression)>0){
    import_gxp(vs,gid,type='tumor')
  }
  if(nrow(vs$tables$fusion)>0){
    import_fus(vs,gid,type='tumor')
  }
}
warnings()
