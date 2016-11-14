# Separates patients id according to the vital status and generate two different groups to calculate the correlation of genes
# Swati Kaushik Nov 14 2016 
# modified from Map-tumor-vital-stats.R

#! /usr/bin/env Rscript 

library(Hmisc)

rm(list=ls())
args <- commandArgs(TRUE)

#extract tumorcode for different vital_status
extract.tumors <- function(type, tumortype){

	type.file1 <- clinical.file[grepl(type, clinical.file$vital_status),]
	type.file <- type.file1[grepl(tumortype, type.file1$person_neoplasm_cancer_status),]
	type.tumorcode <- substr(type.file$tumorcode,1,15)
	type.tumorcode.for <- gsub("-",'.',type.tumorcode, fixed=T)
	type.tumorcode <- toupper(type.tumorcode.for)
	return(type.tumorcode)
	
}	

#calculate overlap of identified tumorcode with the expression datasets
calculate.overlap <- function(expfile, type, vital_status, tumorstatus) {

	deceased.datafile <- expression.file[,intersect(colnames(expfile), type)]
	write.table(deceased.datafile, file = paste(vital_status,tumorstatus ,".patients", sep="") , sep ="\t")	
	return(deceased.datafile)
}

#calculate correlation of genesets
calculate.correlation <- function(expfile,vital_status,tumorstatus) {

	type.corr <- cor(t(expfile), method = "spearman")
	type.corr.round = round(type.corr, digits=2)
	type.hmis <- rcorr(t(expfile), type="spearman")
	write.table(type.corr.round, file= paste(vital_status,tumorstatus,".corr", sep=""), sep="\t")
	write.table(type.hmis$P, file= paste(vital_status,tumorstatus,".patients.corr.pvalue", sep=""), sep="\t")
	
}

# print usage
print.usage <- function(){
	cat("Rscript Map-tumor-vital-stats.R [Expression file] [clinical data]\n")
}

#verify if two data files are passed
if(length(args)<3) {
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

# verify if selected gene file is located in the working directory
subset.genes <- args[[3]] 
if (! file.exists(subset.genes)) {
	cat("subset.genes ", subset.genes,"does not exist\n")
}

#open expression file
expression.file1 <- read.table(args[1], header =T, fill = TRUE, sep="\t")

#get the name of selected gene sets. This is to reduce the computation time
subset.genes <- read.table(args[3], header= T)
unq <- levels(droplevels(subset.genes$genes))
expression.file <- expression.file1[grepl(paste(unq, collapse="|"), rownames(expression.file1)),]
write.table(expression.file, file="interactors.expression.file", sep="\t")

exp.colnames <- colnames(expression.file)
exp.tumorcode <- substr(exp.colnames,1,15)
exp.tumorcode <- toupper(exp.tumorcode)
colnames(expression.file) <- exp.tumorcode

#open clinical data file
clinical.file <- read.table(args[2], header=T, fill=TRUE, sep="\t")

#possible vital status
vital_status <- c("DECEASED", "LIVING")
tumor_status <- c("WITH TUMOR", "TUMOR FREE")

for (i in 1:length(vital_status)) {

	tumor.stat <- unlist(strsplit(tumor_status[i],"\\s"))	
	tumor.stat.name <- paste(tumor.stat[-1],tumor.stat[1], sep="_")
	deceased <- extract.tumors(vital_status[i], tumor_status[i])
	write.table(deceased, file= paste(vital_status[i], tumor.stat.name , sep="_"), sep="\t")
	cat ("Computing overlap",vital_status[i],tumor.stat,"\n")
	overlapped.tumor.names <- calculate.overlap(expression.file, deceased, vital_status[i], tumor.stat.name )
	cat ("Computing correlation for",vital_status[i],tumor.stat,"\n")
	gene.corr <- calculate.correlation(overlapped.tumor.names, vital_status[i], tumor.stat.name )

}
