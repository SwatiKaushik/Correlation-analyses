# Program to extract specific gene list from the correlation matrix and plot a heatmap of correlation of interactors
# Oct 5 2016

rm(list=ls())
library("gplots")

args <- commandArgs(trailingOnly=TRUE)
nargs = length(args)

# Print usage
print.usage <- function(){
	cat("Rscript plot.icorrelation.R [corrData] [interactorname][outputFile]\n")
	cat("Plots correlation matrix\n")
}

if (nargs < 2) {
  print.usage()
  q(save="no",status=1)
}

luad.corr <- args[[1]] # correlation matrix
if (! file.exists(luad.corr)) {
  cat("Correlation Data file ", luad.corr,"does not exist\n")
  q(save="no",status=1)
}

interact.list <- args[[2]] #interactor file
if (! file.exists(interact.list)) {
  cat("Interactor Data file ", interact.list,"does not exist\n")
  q(save="no",status=1)
}

source(luad.corr)
source(interact.list)

int.list <- as.vector(interact.list[,1])  #make vector of interactor file
write.table(int.list, file="int.txt", sep="\t")

# replace all ..* with nothing in column
#colnames(luad.corr) =gsub("\\..*","",colnames(luad.corr))  
#rownames(luad.corr) =gsub("\\..*","",rownames(luad.corr))

data.interactors <- luad.corr[,int.list]   #extract data of interactors from CorrData
data.interactors.log <- -log10(data.interactors)
is.na(data.interactors.log) <- do.call(cbind,lapply(data.interactors.log, is.infinite)) #convert Inf to NA
data.interactors.log[is.na(data.interactors.log)] <- 0  # assign NA to 0

#write the matrix
write.table(data.interactors.log, file="data.interactors.log.txt", sep="\t")

#plot heatmap
pdf(file="heatmap-correlation.pdf")
col.corr =colorRampPalette(c("white","red"),(10))
heatmap.2 (data.interactors.log,
           Rowv = TRUE,
           Colv= TRUE,
           distfun = dist,
           hclustfun = hclust,
           dendrogram = "none",
           scale = "none",
           na.rm = TRUE,
           col= col.corr,
           trace = "none"
           cexRow = 0.2, 
		   cexCol = 0.2)
		   
dev.off()		   
