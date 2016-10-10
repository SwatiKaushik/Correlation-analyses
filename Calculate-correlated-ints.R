#program to calculate percentage of coexpressed interactions across different tumor types
#These calculations were done using interactors at score of 1.5
#Total no. of interactions tested were 4664 (at Z 1.5) 
#These are -log10 (pvalue) of correlation and NOT correlation values

rm(list=ls())

#function to calculate number of interactors with pvalue more than 1.2
correlation.count <-function(file, filename){
	
	#take log10 of pvalue
	file$V4 <- -log10(file$V3)
	#assign NA, Inf to 0
	file$V4[is.infinite(file$V4)] <- 0 
	file$V4[is.na(file$V4)] <- 0 
	#count if pvalue is more than 1.2
	corr.sum <- sum(file$V4>1.2)
	write.table(file, file=paste(filename,".new", sep=""), sep="\t")
	return (corr.sum)
	
}	

args <- commandArgs(trailingOnly=TRUE)
nargs = length(args)

# Print usage
print.usage <- function(){
	 cat("Rscript Calculate-correlated-ints.R [output file]\n")
	 cat ("read files with .Corr extension\n")

}

if (nargs < 1) {  # no. of arguments
  print.usage()
  q(save="no",status=1)
}

#define output file
outfile <- args[1]

tumor.name <- numeric()
corr.value <- numeric()

#Read the files with extension "corr" in a working directory
correlation.files <- Sys.glob("*.Corr")
correlation.file.count = length(correlation.files)
cat ("Number of correlation files =", correlation.file.count,"\n" )

#For loop to read all the correlation files
for (i in 1:length(correlation.files)) {
	
	#print (correlation.files[i])
	if (! file.exists(correlation.files[i])) {
 		 cat("correlation.files[i] does not exist\n")	 
	}
	
	correlation.tumor <- read.table(correlation.files[i], header=F, sep="\t", fill=TRUE)
	correlation.out <- correlation.files[i]
	
	#Parse through the function to count interactors with significant pvalues of correlations
	int.corr <- correlation.count(correlation.tumor, correlation.out)
	#output on the terminal
	cat (correlation.files[i],"\t", int.corr,"\n")
	tumor.name <- append (tumor.name, correlation.files[i])
	corr.value <- append (corr.value, int.corr)
	
}

#calculate percent significantly correlated interactors

sign.corr.int <- sapply(corr.value,function(x) (x*100)/4664)
percent.sign <- cbind(tumor.name, corr.value, sign.corr.int)
write.table(percent.sign, file =outfile, sep="\t")
percentsign <- as.data.frame(percent.sign)

#plot barplot

pdf(file="test.pdf")
par(las=2)
percent.sign <- percentsign[order(percentsign$sign.corr.int),]
reduced.name <- sapply(percent.sign$tumor.name, function(x) substr(x,1,4)[[1]][1])
write.table(percent.sign, file="percentsign", sep="\t")
barplot(percent.sign$sign.corr.int, col = rainbow(10), names.arg = reduced.name)
dev.off()
