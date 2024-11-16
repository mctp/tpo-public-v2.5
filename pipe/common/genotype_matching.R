#!/usr/bin/env Rscript
#
# This script is to determine if a pair of samples (usually tumor and matched normal) 
#       have matching SNP profiles
#

suppressMessages(require(optparse))

option_list <- list( 
  make_option(c("-t", "--tumor"), type='character', default = NULL,
              help="Genotyping csv file of tumor or sample#1"),
  make_option(c("-n", "--normal"), type='character', default = NULL, 
              help="Genotyping csv file of normal or sample#2"),
  make_option(c("-o", "--output"), type='character', default = NULL, 
              help="Prefix for output files")
)
opt <- parse_args(OptionParser(option_list=option_list))

############################################
#1, read genotyping results of normal and tumor
gt_n=read.csv(opt$normal, stringsAsFactors = F)
gt_t=read.csv(opt$tumor, stringsAsFactors = F)
colnames(gt_n)[-1]=paste0("Normal_",colnames(gt_n)[-1])
colnames(gt_t)[-1]=paste0("Tumor_",colnames(gt_t)[-1])
comb=merge(gt_n, gt_t, by="name")
write.csv(comb,  paste0(opt$output,"-genotype.csv"), row.names = F)


############################################
#2, calculate T and N matching scores

# weighted pearson correlation
pearson_cor=function(X, Y) {
  rhos=unlist(lapply(1:nrow(X), function(i) cor(X[i,], Y[i,])))
  s=sum(c(rowSums(Nx),rowSums(Ny)))
  wt=rowSums(Nx)+rowSums(Ny)
  wt=wt/s
  res.wt=sum(rhos*wt) # weighted (by total sequencing depth) average correlation
  return(res.wt)
}


# X and Y are probability matrix of sample A and B, rows are positions, columns are A, C, G, T
# Nx and Ny are count matrix of sample A and B, rows are positions, columns are A, C, G, T
Nx=as.matrix(comb[, c("Normal_A","Normal_C","Normal_G","Normal_T")])
Ny=as.matrix(comb[, c("Tumor_A","Tumor_C","Tumor_G","Tumor_T")])
X=(Nx+1)/rowSums(Nx+4)
Y=(Ny+1)/rowSums(Ny+4)
k=which(rowSums(Nx)==0|rowSums(Ny)==0) # positions with no coverage in either N or T sample
Nx=Nx[-k,]
Ny=Ny[-k,]
X=X[-k,]
Y=Y[-k,]

Pearson_cor=pearson_cor(X, Y)
match_call= Pearson_cor>0.95 # arbitary cutoff

out = file(paste0(opt$output, "-genotype.results.txt"),"w")
writeLines(paste("#" ,nrow(X), "out of 166 SNPs with coverage in both Tumor and Normal"), con=out, sep="\n")
writeLines(paste("GT_SCORE_PC =", round(Pearson_cor,4)), con=out, sep="\n")
writeLines(paste("GT_MATCH =", match_call), con=out, sep="\n")
close(out)
