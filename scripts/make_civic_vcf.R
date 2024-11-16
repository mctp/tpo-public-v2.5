library(data.table)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)


version <-'01-Feb-2019'

vrt <- fread('01-Feb-2019-VariantSummaries.tsv', sep = '\t', 
          header = T)
evi <- fread('01-Feb-2019-ClinicalEvidenceSummaries.tsv', sep = '\t', 
          header = T)


# remove variants with no variant bases 
vrt <- vrt[ vrt$variant_bases!='',]
# remove long insertion or deletion
vrt <- vrt[nchar(vrt$reference_bases)<10 & nchar(vrt$variant_bases)<10,]

# keep columns we want
vrt <- vrt[,c('variant_id','chromosome','start','stop','reference_bases','variant_bases','civic_actionability_score'), with=FALSE]


# add evidence info
vids <- vrt$variant_id[vrt$variant_id %in% evi$variant_id] # some variant ids not found in evidence file because there are no accepted evidence

evi.dis <- sapply(vids, function(x) { # summarize diseases
  d=evi[evi$variant_id==x,]
  dis=table(d$disease)
  dis=paste(names(dis), dis, sep = ':')
  dis=paste(dis, collapse = '|')
})
names(evi.dis) <- vids

# summarize evidence type
type <- rep(FALSE, 4) 
names(type) <- c('Diagnostic','Predictive','Predisposing','Prognostic')
evi.type <- lapply(vids, function(x) {
  d=evi[evi$variant_id==x,]
  evi.type=unique(d$evidence_type[d$evidence_direction=='Supports'])
  out=type
  out[evi.type]=TRUE
  return(out)
})
names(evi.type) <- vids
evi.type <- do.call('rbind', evi.type)

# add evidence info to data frame
vrt$Disease <- evi.dis[as.character(vrt$variant_id)]
tmp <- matrix(FALSE, ncol = 4, nrow = nrow(vrt))
colnames(tmp) <- colnames(evi.type)
rownames(tmp) <- vrt$variant_id
tmp[rownames(evi.type),] <- evi.type
vrt <- cbind.data.frame(as.data.frame(vrt), tmp)

# add required columns for vcf
vrt$chromosome <- paste0('chr', vrt$chromosome)
colnames(vrt)[c(1,5:6,7)] <- c('ID','REF','ALT','Score')
vrt$QUAL <- '.'
vrt$FILTER <- 'PASS'


# deal with multiple alternative alleles
exp <- lapply(vrt$ALT, function(x) {
  alleles=strsplit(x, '[/]')[[1]]
})
vrt$ALT <- exp

# lift to hg38
gr <- GRanges(vrt)
tmp  <-  liftOver(gr, chain=import.chain("/mctp/users/yupingz/refs/hg19ToHg38.over.chain"))
gr.hg38  <-  unlist(tmp[lengths(tmp) == 1])
names(gr.hg38) <- NULL

# check and fill in missing reference bases
ref.bases <- getSeq(Hsapiens, gr.hg38)
disagree.i <- which(as.character(ref.bases)!=gr.hg38$REF & gr.hg38$REF!='')
gr.hg38 <- gr.hg38[-disagree.i]
gr.hg38$REF <- ref.bases[-disagree.i]
gr.hg38$ALT <- DNAStringSetList(gr.hg38$ALT)

# make vcf
vcf <- VCF(granges(gr.hg38))
ref(vcf) <- gr.hg38$REF
alt(vcf) <- gr.hg38$ALT
info(vcf) <- elementMetadata(gr.hg38)[,c('ID','Score','Disease','Diagnostic','Predictive','Predisposing','Prognostic')]
info(header(vcf)) <- DataFrame(Number=c(1,1,'.',0,0,0,0),
                            Type=c('Integer','Float','String',rep('Flag',4)),
                            Description=c('CIVIC variant id','CIVIC actionability score',
                                          'Diseases with evidence',
                                          'CIVIC evidence type-Diagnostic',
                                          'CIVIC evidence type-Predictive',
                                          'CIVIC evidence type-Predisposing',
                                          'CIVIC evidence type-Prognostic'),
                            row.names = colnames(info(vcf)))

# add header 
meta <- DataFrameList(DataFrame(Value='hg38', row.names = 'assembly'),
                                DataFrame(Value=version, row.names = 'CIVIC_version'))
names(meta) <- c('assembly','CIVIC_version')
meta(header(vcf)) <- meta

# add geno to save vcf (otherwise it gives error)
vcf0 <- readVcf('/mctp/users/yupingz/others/Civic/civic.vcf')
geno(vcf) <- geno(vcf0) 
writeVcf(vcf, paste0(version, "_civic.vcf"))
