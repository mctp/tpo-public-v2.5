suppressMessages({
    library(rtracklayer)
    library(data.table)
    library(stringr)
    library(BSgenome)
    library(GenomicFeatures)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(BSgenome.Mmusculus.UCSC.mm10)
})

include_immune_genes <- FALSE
only_preferred <- FALSE
tpopath <- '/data/deployment/dev/tpo'

## artifacts/input/refs/fasta (mostly for VEP)
## 1. download fa sequence from Ensembl
## 2. gunzip bgzip samtools faidx
## 3. add ribosomal RNA

## artifacts/input/refs/genomes (mostly for alignment)
## 1. download fa sequence from some source probably with chr prefix
## 2. gunzip samtools faidx
## 3. generate dict with picard
## 4. get alt??

## download full GTF file
## Human:
## Mouse: ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

## fix tags
#### zcat <file> | ./make_gtf_tags.py fix_tags > <file-tags>

##

inp.dir <- "/work/tmp/artifacts/input/refs/gtf/ensembl.GRCh38.97"
inp.dir <- "/mnt/share/code/tpo-devel/ensembl-108.1"
out.dir <- "/work/tmp/artifacts/input/refs/gtf/ensembl.GRCh38.97"
out.dir <- "/mnt/share/code/tpo-devel/ensembl-108.1"
gtf.inp.fn <- file.path(inp.dir, "Homo_sapiens.GRCh38.108-tags.gtf")
if(include_immune_genes){
  gtf.out.fn <- file.path(out.dir, "Homo_sapiens.GRCh38.108.clean.imm.gtf")
}else{
  gtf.out.fn <- file.path(out.dir, "Homo_sapiens.GRCh38.108.clean.gtf")
}
fa.out.fn <- file.path(out.dir, "Homo_sapiens.GRCh38.108.cdna.clean.fa")
genome <- BSgenome.Hsapiens.UCSC.hg38

inp.dir <- "/tpo/refs/grcm38/ensembl"
out.dir <- "/mnt/share/code/tpo-dev/mouse"
gtf.inp.fn <- file.path(inp.dir, "grcm38.102.all.gtf")
gtf.out.fn <- file.path(out.dir, "grcm38.102.clean.gtf")
fa.out.fn <- file.path(out.dir, "grcm38.102.clean.cdna.fa")
genome <- BSgenome.Mmusculus.UCSC.mm10


{

gen.gr <- import(gtf.inp.fn)
gen.dt <- as.data.table(gen.gr)
gen.dt.tx <- gen.dt[type=="transcript"]

## protein-coding trancripts of protein coding genes
prot.tx <- gen.dt.tx[
        gen.dt.tx$gene_biotype=="protein_coding" &
        gen.dt.tx$transcript_biotype=="protein_coding" &
        !grepl("\\.", gen.dt.tx$gene_name) &
        (
            ## remove annotated readthrough genes
            !grepl(".*[A-Z].*-.*[A-Z].*", gen.dt.tx$gene_name) |
            grepl("RP11-|MT-|HLA-", gen.dt.tx$gene_name)
        )
]$transcript_id

## non-coding transcripts and 'novel' protein-coding genes
nonc.tx <- gen.dt.tx[
    grepl("pseudogene|lincRNA", gen.dt.tx$gene_biotype) |
    (gen.dt.tx$gene_biotype=="protein_coding" & grepl("\\.", gen.dt.tx$gene_name))
]$transcript_id

imm.tx <- gen.dt.tx[
  grepl("IG_[CVDJ]_gene|TR_[CVDJ]_gene",gen.dt.tx$transcript_biotype)
]$transcript_id

## remove non-coding and immune transcripts which overlap a protein-conding transcript
gen.gr.prot <- gen.gr[(gen.gr$transcript_id %in% prot.tx)]
gen.gr.nonc <- gen.gr[(gen.gr$transcript_id %in% nonc.tx)]
gen.gr.nonc.tx <- gen.gr.nonc[gen.gr.nonc$type=="transcript"]
gen.gr.prot.tx <- gen.gr.prot[gen.gr.prot$type=="transcript"]
hits <- findOverlaps(gen.gr.nonc.tx, gen.gr.prot.tx)
gen.gr.nonc <- gen.gr.nonc[!(gen.gr.nonc$transcript_id %in% gen.gr.nonc.tx[queryHits(hits)]$transcript_id)]
#
gen.gr.imm <- gen.gr[(gen.gr$transcript_id %in% imm.tx)]
gen.gr.imm.tx <- gen.gr.imm[gen.gr.imm$type=="transcript"]
hits <- findOverlaps(gen.gr.imm.tx, gen.gr.prot.tx)
gen.gr.imm <- gen.gr.imm[!(gen.gr.imm$transcript_id %in% gen.gr.imm.tx[queryHits(hits)]$transcript_id)]


if(include_immune_genes){
  gen.gr.full <- c(gen.gr.prot, gen.gr.nonc, gen.gr.imm)
}else{
  gen.gr.full <- c(gen.gr.prot, gen.gr.nonc)
}

#
if(only_preferred){
  #load the preferred transcripts
  pref <- file.path(tpopath,'rlibs/carat/inst/extdata/onco1500-v6a/transcript_preference.csv')
  pref <- read.csv(pref,stringsAsFactors=F) %>% filter(rank==1)
  #map them to ENST
  map <- file.path(tpopath,'rlibs/carat/inst/extdata/ensembl97/refseq2ensembl.csv')
  map <- read.csv(map,stringsAsFactors=F)
  pref_enst <- merge(pref,map,
                  by.x=c('ensg','transcript_id'),
                  by.y=c('ensembl_gene_id','refseq'),
                  all.x=T
                  )
  #bring along any preferred ENST
  w <- which(is.na(pref_enst$ensembl_transcript_id))
  pref_enst$ensembl_transcript_id[w] <- pref_enst$transcript_id[w]
  #subset to preferred where appropriate
  out <-  gen.gr.full[
              !(gen.gr.full$gene_id %in% pref_enst$ensg) |
              (gen.gr.full$transcript_id %in% pref_enst$ensembl_transcript_id) | 
              (gen.gr.full$gene_id=='ENSG00000163098') # this won't map, including all 
            ]
 #make sure a transcript is chosen for all genes 
 stopifnot(all(gen.gr.full$gene_id %in% out$gene_id))
 gen.gr.full <- sort(out)
}

## create full GTF
gen.gr.full.gn <- unlist(range(split(gen.gr.full, gen.gr.full$gene_id)))
gen.gr.gn <- gen.gr[gen.gr$type=="gene"]
names(gen.gr.gn) <- gen.gr.gn$gene_id
mcols(gen.gr.full.gn) <- mcols(gen.gr.gn[names(gen.gr.full.gn)])
fin.gr <- sort(c(unname(gen.gr.full.gn), gen.gr.full))

## rename to UCSC main chromosomes
fin.gr.drop <- keepStandardChromosomes(fin.gr, pruning.mode="coarse")
fin <- sort(fin.gr.drop)
seqlevelsStyle(fin) <- 'UCSC'

## update from 97->106 has genes without gene_names
fin[is.na(fin$gene_name)]$gene_name <- fin[is.na(fin$gene_name)]$gene_id

#save the gtf
export(fin, gtf.out.fn)

## decide if needed
txdb <- makeTxDbFromGFF(gtf.out.fn)
seqlevelsStyle(txdb) <- "UCSC"
txseq <- extractTranscriptSeqs(genome, txdb)
names(txseq) <- transcripts(txdb)$tx_name
writeXStringSet(txseq, fa.out.fn)

}
