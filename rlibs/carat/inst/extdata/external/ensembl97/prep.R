library(biomaRt)
library(stringr)
library(data.table)

ensembl97 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', version=97)

attrs <- c("ensembl_gene_id", "refseq_mrna", "ensembl_transcript_id", 'chromosome_name')
tmp <- getBM(attributes=attrs, mart=ensembl97)
tmp.mrna <- data.table(tmp)
attrs <- c("ensembl_gene_id", "refseq_ncrna", "ensembl_transcript_id", 'chromosome_name')
tmp <- getBM(attributes=attrs, mart=ensembl97)
tmp.ncrna <- data.table(tmp)
tmp <- rbind(tmp.mrna, tmp.ncrna, use.names=FALSE)

refseq2ensembl <- tmp[refseq_mrna!="", .(ensembl_gene_id, ensembl_transcript_id, refseq=refseq_mrna, chromosome_name)]
fwrite(refseq2ensembl, "refseq2ensembl.csv")

attrs <- c("ensembl_gene_id", "ensembl_transcript_id")
tmp <- getBM(attributes=attrs, mart=ensembl97)
enst2ensg <- data.table(tmp)
fwrite(enst2ensg, "enst2ensg.csv")
