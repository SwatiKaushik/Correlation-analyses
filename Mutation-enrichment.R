# Enrichment analysis of interactors 
# Swati Kaushik Nov 7 2016

rm(list=ls())
args <- commandArgs(TRUE)

#To get non intersecting pairs
outersect <- function(x, y) {
     sort(c(setdiff(x, y),
            setdiff(y, x)))
 }
 
#Mutated gene sets from TCGA
mutation.file <- read.table(args[1], header= F, sep="\t", fill=TRUE)
HS.mutation.file <- levels(droplevels(mutation.fileM$V1))

#interaction file
interaction.file <- read.table(args[2], header=T)
subset.ints <- subset(interaction.file, interaction.file$zscore>2, drop=FALSE)
preys <- levels(droplevels(subset.ints$Prey))

#cancer gene census
cancer.census <- read.table(args[3], header =T, sep="\t", fill=TRUE, quote="")
cancer.census.genes <- levels(droplevels(cancer.census$Gene.Symbol))

#calculate overlap between args[1] and args[2]
overlap <- intersect(preys, HS.mutation.file)
overlap.length <- length(overlap)
preys.length <- length(preys)
ccg.length <- length(HS.mutation.file)

#calculate significance
significance.test <- rbind(c(overlap.length, preys.length), c(ccg.length, 20000))
significance.pvalue <- fisher.test(significance.test, alternative = "greater")$p.value
significance.pvalue

#extract kinases of the overlapped preys
overlpped.data.kinase <- interaction.file[grepl(paste(overlap, collapse ="|"), interaction.file$Prey),]
subset.kinase.ints <- subset(overlpped.data, overlpped.data.kinase$zscore>2, drop=FALSE)
sorted.subset.kinase.ints <- subset.kinase.ints[order(subset.kinase.ints$Prey),]
write.table(sorted.subset.kinase.ints, file="sorted.subset.kinase.ints.out", sep="\t", row.names=FALSE, quote = FALSE)

#overlap with third geneset
cancer.censusoverlap <- intersect(overlap, cancer.census.genes)
not.covered.genes <- outersect(cancer.censusoverlap, overlap)


