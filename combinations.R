#Make combinations of all the genes given as input
#Oct 12 2016

rm(list=ls())

args <- commandArgs(trailingOnly=TRUE)
nargs = length(args)

print.usage <- function(){
	 cat("Rscript combinations.R [input file] [output file]\n")
	 cat ("read input gene list\n")
}

if (nargs < 2) {  # no. of arguments
  print.usage()
  q(save="no",status=1)
}

gene.list <- args[[1]] # genelist
if (! file.exists(gene.list)) {
  cat("gene file ", gene.list,"does not exist\n")
  q(save="no",status=1)
}

#define output file
outfile <- args[2]

load('gene.list.RData')

gene.vector <- as.vector(gene.list$V1)
#Make combinations of all the genes
comb.gene <- t(combn(gene.vector, 2))
write.table(comb.gene, file =outfile, sep="\t", row.names = FALSE, quote = FALSE)

