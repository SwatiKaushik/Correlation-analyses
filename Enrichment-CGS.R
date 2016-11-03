# compares two gene list to calculate the significance of the enrichment
# statistics: Fischer Exact test (one sided)
# Swati Kaushik, Nov 3 2016

rm(list=ls())
args <- commandArgs(TRUE)

#Calculate the significance of the enrichment
significance.calculation <- function (genename){
	
	#get data for the inpit gene
	kinase <- interaction.file[grep(genename, interaction.file$Bait),]
	ints <- subset(kinase, kinase$zscore> interaction.score, drop=FALSE)
	preys <- levels(droplevels(ints$Prey))
	
	#calculate overlap betwee two datsets
	overlap <- intersect(preys, cancer.census.genes)
	overlap.length <- length(overlap)
	preys.length <- length(preys)
	ccg.length <- length(cancer.census.genes)
	
	#calculate pvalue of overlap
	significance.test <- rbind(c(overlap.length, preys.length), c(ccg.length, 20000))
	significance.pvalue <- fisher.test(significance.test, alternative = "greater")$p.value
	
	#return data
	gene.names <- paste(overlap, sep=" ", collapse =",")
	data <- paste(significance.pvalue, overlap.length, preys.length, ccg.length, sep="\t")
	data2 <- paste(data, gene.names, sep="\t")
	return(data2)
	#return(significance.pvalue)
	
}

# print usage
print.usage <- function(){
	cat("Rscript Gene-enrichment.R [first gene file] [interesed-gene-list] [cutoff][outputFile]\n")
	cat("Performs enrichment analysis of the given genes\n")
}

# verify that at least required parameters are passed 
if(length(args)<4) {
	cat(print.usage)
	stop("Incorrect or missing required input!")
}

# verify if input file name is present in directory
cancer.census <- args[[1]] 
if (! file.exists(cancer.census)) {
	cat("gene file ", cancer.census,"does not exist\n")
}

# verify if input file name is present in directory
interaction.file <- args[[2]] 
if (! file.exists(interaction.file )) {
	cat("gene file ", interaction.file ,"does not exist\n")
}

# verify if input file name is present in directory
interaction.score <- args[[3]] 
if (interaction.score <0) {
	cat("warning: may not be true interactions\n")
}

#assign output file
output <- args[4]

#cancer gene census
cancer.census <- read.table(args[1], header =T, sep="\t", fill=TRUE)
cancer.census.genes <- levels(droplevels(cancer.census$Gene.Symbol))

interaction.file <- read.table(args[2], header=T)

#get unique gene names
unq.kinase <- unique(interaction.file$Bait)
unq <- levels(droplevels(unq.kinase))

#perform enrichment calculations
data <- lapply(unq, FUN = significance.calculation )
output <- (cbind(unq, data))

#write output
write.table(output, file=args[4], sep="\t")

