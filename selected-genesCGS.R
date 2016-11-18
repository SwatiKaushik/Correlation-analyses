# Calculate overlap of two genesets and plots venn diagram of the overlap
# Uses Fischer exact test to calculate significance
# Nov 18 2016

rm(list=ls())

library(VennDiagram)

args <- commandArgs(trailingOnly=TRUE)
nargs = length(args)

#Print usage
print.usage <- function(){
	 cat("Rscript selected-genesCGS.R [file1][file2]\n")
}

#check # of arguments
if (nargs < 2) {  # no. of arguments
  print.usage()
}

# check existence of file 1:cancer gene file
cancer.file <- args[[1]] 
if (! file.exists(cancer.file)) {
  cat("cancer genes file ", cancer.file,"does not exist\n")
}

# check existence of file 2:AP-MS
ap.ms.file <- args[[2]] 
if (! file.exists(ap.ms.file)) {
  cat("AP-MS genes file ", ap.ms.file ,"does not exist\n")
}

#cancer gene census
cancer.census <- read.table(cancer.file, header =T, sep="\t", fill=TRUE, quote="")
cancer.census.genes <- levels(droplevels(cancer.census$Gene.Symbol))

#data APMS
interaction.file <- read.table(ap.ms.file, header=T)
bait.unq <- unique(interaction.file$Bait)
bait.unq.d <- levels(droplevels(bait.unq))
bait.subset<- subset(interaction.file, interaction.file$zscore>2.5, drop=FALSE)
preys <- levels(droplevels(bait.subset$Prey))
preys.noTK <- preys[grep(paste(bait.unq.d, collapse="|"), preys, invert=TRUE)]

#calculate overlap
overlap <- intersect(preys.noTK, cancer.census.genes)
overlap.length <- length(overlap)
preys.length <- length(preys.noTK)
ccg.length <- length(cancer.census.genes)

#calculate significance
significance.test <- rbind(c(overlap.length, preys.length), c(ccg.length, 20000))
significance.pvalue <- fisher.test(significance.test, alternative = "greater")$p.value
#data = rbind(overlap, significance.pvalue)
write.table(overlap, file="overlap2.5.out", sep="\t")
write.table(significance.pvalue, file = "overlap2.out", append=TRUE)
cat ("significance.pvalue", significance.pvalue)

#plot venn diagram

grid.newpage()
venn.plot <- draw.pairwise.venn( area1 = preys.length, 
					area2 = ccg.length, 
					cross.area = overlap.length, 
					category =c ("Cancer gene census", "TK-interactome"),
					lty = rep("blank",2), 
					fill = c("light blue", "pink"), 
					alpha = rep(0.5, 2), 
					cat.pos = c(0, 0), 
					cat.dist = rep(0.025, 2), 
					scaled = TRUE
				 )

#require(gridExtra)
#grid.arrange(gTree(children=venn.plot), top= A , bottom="subtitle")

#writing to file
tiff(filename = "Venn_diagram2.5.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()