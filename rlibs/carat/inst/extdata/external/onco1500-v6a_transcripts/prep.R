library(stringr)
library(data.table)
library(biomaRt)
library("readxl")

## https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/GRCh38_interim_annotation/interim_GRCh38.p12_top_level_2019-01-25.gff3.gz

tmp <- data.table(read_excel("raw/OncoSeqv6a_Preferred_Annotations.xlsx"))
tmp[, transcript_id:=`LOCUS or TYPE`]
tmp <- tmp[grepl("N[MR]_|ENST", transcript_id)]
tmp[`EXON or SNP`=='NM_007313_0', transcript_id:="NM_007313"]
tmp[transcript_id=="NM_007313,NM_005157", transcript_id:="NM_005157"]
tmp <- tmp[,transcript_id:=str_match(transcript_id, "[^.]*")[,1]]

## fix refseq to refseq
tmp[transcript_id=="NM_001365420", transcript_id:="NM_001365417"]
tmp[transcript_id=="NM_020732", transcript_id:="NM_001374820"]
tmp[transcript_id=="NM_001304745", transcript_id:="NM_018283"]
tmp[transcript_id=="NM_001267568", transcript_id:="NM_182718"]
tmp[transcript_id=="NM_001267566", transcript_id:="NM_182769"]
tmp[transcript_id=="NM_1290039", transcript_id:="NM_001290039"]
tmp[transcript_id=="NM_001214903", transcript_id:="NM_001437"]
tmp[transcript_id=="NM_001146186", transcript_id:="NM_006210"]
tmp[transcript_id=="NM_01164603", transcript_id:="NM_001164603"]

## fix refseq to ensembl
tmp[transcript_id=="NM_001369461", transcript_id:="ENST00000397575"]
tmp[transcript_id=="NM_001369472", transcript_id:="ENST00000606230"]
tmp[transcript_id=="NM_001369476", transcript_id:="ENST00000636057"]
tmp[transcript_id=="NM_001374820", transcript_id:="ENST00000346085"]
tmp[transcript_id=="NM_001385557", transcript_id:="ENST00000433722"]

## 
tmp[,n.exons:=.N,transcript_id]
v6a <- unique(tmp[,.(transcript_id, symbol=GENE, n.exons)])

MAIN_CHROMS <- c(1:22, "X", "Y")
refseq2ensembl <- fread("../ensembl97/refseq2ensembl.csv")
refseq2ensembl <- refseq2ensembl[chromosome_name %in% MAIN_CHROMS]

setkey(v6a, transcript_id)
setkey(refseq2ensembl, refseq)
v6a.refseq <- v6a[refseq2ensembl, ensg:=ensembl_gene_id]

enst2ensg <- fread("enst2ensg.csv")
setkey(enst2ensg, ensembl_transcript_id)
v6a[enst2ensg, ensg:=ensembl_gene_id]

v6a[order(-n.exons),rank:=(1:nrow(.SD)),ensg]
v6a[is.na(ensg), rank:=NA]
v6a <- v6a[order(ensg, rank),.(ensg, transcript_id, rank)]

fwrite(v6a, "transcript_preference.csv")
