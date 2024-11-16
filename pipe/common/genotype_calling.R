#!/usr/bin/env Rscript
#
# This script is to generate SNP count matrix from a GVCF file
#
suppressMessages(require(optparse))
suppressMessages(require(VariantAnnotation))
suppressMessages(require(reshape2))

option_list <- list(
  make_option(c("-i", "--input"), type='character', default = NULL,
              help="Input gvcf file name"),
  make_option(c("-o", "--output"), type='character', default = NULL,
              help="Output csv file name")
)
opt <- parse_args(OptionParser(option_list=option_list))

##############
# Functions

# Follow Bob's criteria:
# to make a heterozygous call, require 10 reads of each allele, variant fractions >= 30% for each, 
# and >= 90% for the two bases combined.
# to make a homozygous call, require 15 reads of the allele and variant fraction >= 90%.

is_homozygous <- function(ct, tot) {
  return(ct >= 3 & ct/tot >= .8)
}
  
is_heterozygous <- function(ct1,ct2,tot) {
  return(tot >= 13 & (ct1/tot) >= .2 & (ct2/tot) >= .2 &
         (ct1+ct2)/tot >= .9)
}

gt_call=function(vcf_file) {
    
  vcf <- readVcf(vcf_file)
  evcf <- expand(vcf) # get expanded vcf

  dp <-geno(evcf)$DP # read depth
  gt <-geno(evcf)$GT # genotype
  ad <-geno(evcf)$AD # allele depths. it's an array with dimensions of 253 * 1 * 2. the last dimension, which includes 2 elements, are corresponding to ref and alt alleles.
  gr <- rowRanges(evcf) # snps in granges
  
  ref <- ref(evcf) # reference alleles
  alt <- alt(evcf) # alternative alleles
  
  snps <- paste(seqnames(gr), start(gr), sep = '_') # name snps by coordinates
  
  # compile data to use
  df <- data.frame(name=snps, ref=ref, alt=alt, ref.ct=ad[,1,1], alt.ct=ad[, 1,2], 
                   depth=dp[,1], genotype=gt[,1], stringsAsFactors = F)
  
  df[is.na(df)] <- 0
  df$ref.ct[df$genotype=='0/0'] <- df$depth[df$genotype=='0/0']
  
  # generate allele count matrix
  tmp1 <- df[,c(1,2,4)]
  tmp2 <- df[,c(1,3,5)]
  colnames(tmp1)[2:3] = colnames(tmp2)[2:3] = c('allele','count')
  tmp <- rbind(tmp1, tmp2)
  tmp <- tmp[!duplicated(tmp),] # remove duplicated reference allele rows
  snp.sum <- dcast(tmp, name ~ allele, value.var = 'count')
  snp.sum[is.na(snp.sum)] <- 0
  snp.sum <- snp.sum[,c("name", "<NON_REF>", "A", "C", "G", "T")] # drop possible indels
  snp.sum$tot <- rowSums(snp.sum[,c('A','C','G','T')])
  snp.sum$genotype <- df[match(snp.sum$name, df$name),]$genotype
  snp.sum$genotype[snp.sum$tot<3 & snp.sum$genotype!='0/1']='./.'
  snp.sum$genotype[snp.sum$tot<13 & snp.sum$genotype=='0/1']='./.'
  
  # make calls based on Bob's python script
  calls <- sapply(1:nrow(snp.sum), function(i) {
    ctA <- snp.sum[i, 'A']
    ctC <- snp.sum[i, 'C']
    ctG <- snp.sum[i, 'G']
    ctT <- snp.sum[i, 'T']
    ctTot <- snp.sum[i, 'tot']
    
    final_call <-  'UNKNOWN'
    
    if (is_homozygous(ctA, ctTot)) {
      final_call = 'A'
    }
    
    if (is_homozygous(ctC,ctTot)) {
      final_call = 'C'
    }
    
    if (is_homozygous(ctG,ctTot)) {
      final_call = 'G'
    }
    
    if (is_homozygous(ctT,ctTot)) {
      final_call = 'T'
    }
    
    if (is_heterozygous(ctA,ctC,ctTot)) {
      final_call = 'A/C'
    }
    
    if (is_heterozygous(ctA,ctG,ctTot)) {
      final_call = 'A/G'
    }
    
    if (is_heterozygous(ctA,ctT,ctTot)) {
      final_call = 'A/T'
    }
    
    if (is_heterozygous(ctC,ctG,ctTot)) {
      final_call = 'C/G'
    }
    
    if (is_heterozygous(ctC,ctT,ctTot)) {
      final_call = 'C/T'
    }
    
    if (is_heterozygous(ctG,ctT,ctTot)) {
      final_call = 'G/T'
    }
    return(final_call)
  })
  
  snp.sum$final_call=calls
  
  return(snp.sum)
}


################
gt_res=gt_call(opt$input)
gt_res=gt_res[,colnames(gt_res)!="<NON_REF>"]
write.csv(gt_res,  opt$output, row.names = FALSE)
