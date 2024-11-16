library(VariantAnnotation)

#load the old cosmic vcf
old <- readVcf('./cosmic.vcf.gz')

#get the ranges and which are delins
rr <- rowRanges(old)
f <- fixed(old)
#these are DNAStringSet with all(length==1)
rr$ALT <- unlist(rr$ALT)
f$ALT <- unlist(f$ALT)
w <- which((nchar(rr$REF)==nchar(rr$ALT)) & nchar(rr$REF)>1)

length(w)
#[1] 22338
rr
as.data.frame(rr[w][1])
#seqnames  start    end width strand paramRangeID REF ALT QUAL FILTER
#COSM1492151        1 961920 961922     3      *         <NA> TGG TTT   NA      .

#remove the anchor base
rr[w]$REF <- subseq(rr[w]$REF,start=2) 
rr[w]$ALT <- subseq(rr[w]$ALT,start=2) 
f$REF[w] <- subseq(f$REF[w],start=2) 
f$ALT[w] <- subseq(f$ALT[w],start=2) 

#bump the start by one
start(rr[w]) <- start(rr[w])+1

#sanity check
all(rr$REF==f$REF)
all(rr$ALT==f$ALT)
mcols(rr) <- NULL

#make a new VCF object with the altered rowRanges
#for some reason rowRanges() works as an accessor but not for assignment?

new <- VCF(
  rowRanges=rr,
  colData=colData(old),
  fixed=f,
  info=info(old),
  geno=geno(old),
  collapsed=F
)

as.data.frame(rowRanges(new)[w][1])
#seqnames  start    end width strand paramRangeID REF ALT QUAL FILTER
#COSM1492151        1 961921 961922     2      *         <NA>  GG  TT   NA      .


#some sanity checks
length(old)==length(new)
#[1] TRUE
all(as.integer(info(old)$CNT)==as.integer(info(new)$CNT))
#[1] TRUE
all(end(rowRanges(old))==end(rowRanges(new)))
#[1] TRUE
all(start(rowRanges(old))==start(rowRanges(new)))
#[1] FALSE

writeVcf(new,'./cosmic.v87.norm.vcf')

