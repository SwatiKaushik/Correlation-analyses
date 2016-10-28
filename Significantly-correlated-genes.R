#Calculate correlation and pvalue of correlation of the given matrix

library(Hmisc)
library(MASS)

#Read data matrix
args <- commandArgs(TRUE)
data <- as.matrix(read.table(args[1], header =TRUE))
#bass = scale(t(data))
#attr(bass, "scaled:scale") =NULL
#attr(bass, "scaled:center") = NULL

#Calculate gene-gene correlation of the matrix
correlation.matrix = cor(t(data), method ="spearman")

#Round the matrix to 2 digits
rounded.matrix = round(correlation.matrix, digits=2)

#Calculate the pvalue of the matrix using package Hmisc
pvalue.matrix= rcorr(t(data), type="spearman")

#write output data files in the matrix
write.table(rounded.matrix, paste (args[1],".spearman",sep=""), sep ="\t")
write.table(pvalue.matrix$P, paste(args[1],".spearman.pvalue",sep=""), sep="\t")

