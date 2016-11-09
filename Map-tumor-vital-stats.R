# Separates patients id according to the vital status and generate two different groups to calculate the correlation of genes
# Swati Kaushik Nov 8 2016

#! /usr/bin/env Rscript 

library(Hmisc)

rm(list=ls())
args <- commandArgs(TRUE)

#extract tumorcode for different vital_status
extract.tumors <- function(type){

	type.file <- clinical.file[grepl(type, clinical.file$vital_status),]
	type.tumorcode <- substr(type.file$tumorcode,1,15)
	type.tumorcode.for <- gsub("-",'.',type.tumorcode, fixed=T)
	type.tumorcode <- toupper(type.tumorcode.for)
	return(type.tumorcode)
	
}	

#calculate overlap of identified tumorcode with the expression datasets
calculate.overlap <- function(expfile, type, vital_status) {

	deceased.datafile <- expression.file[,intersect(colnames(expfile), type)]
	write.table(deceased.datafile, file = paste(vital_status ,".patients", sep="") , sep ="\t")	
	return(deceased.datafile)
}

#calculate correlation of genesets
calculate.correlation <- function(expfile,vital_status) {

	type.corr <- cor(t(expfile), method = "spearman")
	type.corr.round = round(type.corr, digits=2)
	type.hmis <- rcorr(t(expfile), type="spearman")
	write.table(type.corr.round, file= paste(vital_status,".corr", sep=""), sep="\t")
	write.table(type.hmis$P, file= paste(vital_status,".patients.corr.pvalue", sep=""), sep="\t")
	
}

# print usage
print.usage <- function(){
	cat("Rscript Map-tumor-vital-stats.R [Expression file] [clinical data]\n")
}

#verify if two data files are passed
if(length(args)<2) {
	cat(print.usage)
	stop("Incorrect or missing required input!")
}

# verify if expression file is located in the working directory
expression.file <- args[[1]] 
if (! file.exists(expression.file)) {
	cat("Expression file ", expression.file,"does not exist\n")
}

# verify if clinical file is located in the working directory
clinical.file <- args[[2]] 
if (! file.exists(clinical.file)) {
	cat("Clinical file ", clinical.file,"does not exist\n")
}

#open expression file
expression.file <- read.table(args[1], header =T, fill = TRUE, sep="\t")
exp.colnames <- colnames(expression.file)
exp.tumorcode <- substr(exp.colnames,1,15)
exp.tumorcode <- toupper(exp.tumorcode)
colnames(expression.file) <- exp.tumorcode

#open clinical data file
clinical.file <- read.table(args[2], header=T, fill=TRUE, sep="\t")

#possible vital status
vital_status <- c("DECEASED", "LIVING")

for (i in 1:length(vital_status)) {
	
	deceased <- extract.tumors(vital_status[i])
	cat ("Computing overlap",vital_status[i],"\n")
	overlapped.genes <- calculate.overlap(expression.file, deceased, vital_status[i])
	cat ("Computing correlation for",vital_status[i],"\n")
	gene.corr <- calculate.correlation(overlapped.genes, vital_status[i])

}
