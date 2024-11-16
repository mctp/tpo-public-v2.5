library(data.table)
library(stringr)
library(rtracklayer)
library(readxl)
setwd("/mnt/share/code/tpo/rlibs/carat/inst/extdata/known_transcripts")

####
## library(biomaRt)
## library(httr)
## httr::set_config(config(ssl_verifypeer = FALSE))
## ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version=97)
## res <- getBM(attributes=c("refseq_mrna", "ensembl_transcript_id"), mart= ensembl)
## res <- as.data.table(res)
## res <- unique(res[refseq_mrna!="" & ensembl_transcript_id!=""])
## fwrite(res, "/mnt/share/code/tpo/rlibs/carat/inst/extdata/known_transcripts/170620_refseq_to_ensembl.txt")
## ### this is final design document 
## tmp <- as.data.table(read_excel("/mnt/share/code/tpo-devel/carat/meta/onco1500_final.xlsx", sheet=1))
## mo.raw = unique(tmp[grep("NM_|NR_|ENST",`LOCUS or TYPE`)][,.(Feature_ID=`LOCUS or TYPE`, Gene=GENE)])[order(Gene)]
## fwrite(mo.raw, "210221_refseq_tx_orderings_v6.txt", sep="\t")


ccds.raw <- fread("./170620_CCDS2Sequence.current.txt")
ccds.short <- unique(ccds.raw[,.(id.short=str_match(nucleotide_ID, "[^\\.]*")[,1])])

civic.raw <- fread("./170620_Ensembl-v75_build37-hg19_UcscGenePred_CIViC-Genes.ensGene")
civic.short <- unique(civic.raw[,.(id.short=str_match(hg19.ensGene.name, "[^\\.]*")[,1])])

lrg.raw <- fread("./170620_list_LRGs_transcripts_xrefs.txt", skip=1)
lrg.short <- rbind(
    unique(lrg.raw[,.(id.short=str_match(ENSEMBL_TRANSCRIPT, "[^\\.]*")[,1])])[id.short!="-"],
    unique(lrg.raw[,.(id.short=str_match(REFSEQ_TRANSCRIPT, "[^\\.]*")[,1])])[id.short!="-"]
)

mo.raw.v4 <- fread("./170620_refseq_tx_orderings_v4.txt")
mo.short.v4 <- unique(mo.raw.v4[,.(id.short=str_match(Feature_ID, "[^\\.]*")[,1])])

mo.raw.v6 <- fread("./210221_refseq_tx_orderings_v6.txt")
mo.short.v6 <- unique(mo.raw.v6[,.(id.short=str_match(Feature_ID, "[^\\.]*")[,1])])

refseq2enst <- fread("./170620_refseq_to_ensembl.txt")
setkey(mo.short.v4, id.short)
setkey(refseq2enst, refseq_mrna)
mo.short.ensembl.v4 <- unique(refseq2enst[mo.short.v4,.(id.short=ensembl_transcript_id),nomatch=0])
setkey(mo.short.v6, id.short)
setkey(refseq2enst, refseq_mrna)
mo.short.ensembl.v6 <- unique(refseq2enst[mo.short.v6,.(id.short=ensembl_transcript_id),nomatch=0])

gencode.gtf <- import.gff("./170620_gencode.v34.annotation.tx.tags.gtf.gz")
gencode.raw <- as.data.table(mcols(gencode.gtf)[,c("transcript_id", "level", "transcript_support_level", "tags")])
gencode.raw <- unique(gencode.raw[grepl("basic", tags),
                   .(
                       id.short=str_match(transcript_id, "[^\\.]*")[,1], xbasic=TRUE, xvalid=level%in%1,
                       xtsl1=transcript_support_level%in%1, xappris1=grepl("appris_principal_1", tags)
                   )])
gencode.short <- gencode.raw[,.(id.short)]

refseq.gtf <- import.gff("./170620_GCF_000001405.39_GRCh38.p13_genomic.tx.gff.gz")
refseq.best.gtf <- refseq.gtf[refseq.gtf$source=="BestRefSeq"]
refseq.best.raw <- as.data.table(mcols(refseq.best.gtf)[,c("Name", "tag")])[,.(id.short=str_match(Name, "[^\\.]*")[,1], xselect=!is.na(tag))]
refseq.short <- refseq.best.raw[,.(id.short)]

canon <- unique(rbind(ccds.short, civic.short, lrg.short, mo.short.v4, mo.short.v6, gencode.short, refseq.short))
canon$known <- TRUE
canon$ccds <- FALSE
canon$civic <- FALSE
canon$lrg <- FALSE
canon$mo.v4 <- FALSE
canon$mo.v6 <- FALSE
canon$basic <- FALSE
canon$valid <- FALSE
canon$tsl1 <- FALSE
canon$appris1 <- FALSE
canon$select <- FALSE

setkey(canon, id.short)
setkey(gencode.raw, id.short)
setkey(mo.short.v4, id.short)
setkey(mo.short.ensembl.v4, id.short)
setkey(mo.short.v6, id.short)
setkey(mo.short.ensembl.v6, id.short)
setkey(lrg.short, id.short)
setkey(civic.short, id.short)
setkey(ccds.short, id.short)
setkey(refseq.best.raw, id.short)
canon[gencode.raw, ":="(
                       basic=xbasic,
                       valid=xvalid,
                       tsl1=xtsl1,
                       appris1=xappris1
)]
canon[refseq.best.raw, ":="(
                           select=xselect
)]
canon[mo.short.v4, mo.v4:=TRUE]
canon[mo.short.ensembl.v4, mo.v4:=TRUE]
canon[mo.short.v6, mo.v6:=TRUE]
canon[mo.short.ensembl.v6, mo.v6:=TRUE]
canon[ccds.short, ccds:=TRUE]
canon[civic.short, civic:=TRUE]
canon[lrg.short, lrg:=TRUE]

fwrite(canon, "canonical_transcripts2.txt", sep=",")
