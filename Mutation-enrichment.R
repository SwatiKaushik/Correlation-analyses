# Enrichment analysis of interactors 
# Uses two data files - significantly mutated genes and a list of interactors in the A_B_C format
# Swati Kaushik Nov 21 2016

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
bait.unq <- unique(interaction.file$Bait)
bait.unq.d <- levels(droplevels(bait.unq))
bait.subset<- subset(interaction.file, interaction.file$zscore>2, drop=FALSE)
preys <- levels(droplevels(bait.subset$Prey))
preys.noTK <- preys[grep(paste(bait.unq.d, collapse="|"), preys, invert=TRUE)]
#preys.noTK

#cancer gene census
cancer.census <- read.table(args[3], header =T, sep="\t", fill=TRUE, quote="")
cancer.census.genes <- levels(droplevels(cancer.census$Gene.Symbol))

#calculate overlap between args[1] and args[2]
overlap <- intersect(preys.noTK, HS.mutation.file)
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

#plot venn diagram

grid.newpage()
venn.plot <- draw.pairwise.venn( area1 = preys.length, 
					area2 = ccg.length, 
					cross.area = overlap.length, 
					category =c ("SMG", "TK-interactome"),
					lty = rep("blank",2), 
					fill = c("light blue", "pink"), 
					alpha = rep(0.5, 2), 
					cat.pos = c(0, 0), 
					cat.dist = rep(0.025, 2), 
					scaled = TRUE,
					lwd=4,
					cex = 2,
					col = rep("black",2)
				 )

#require(gridExtra)
#grid.arrange(gTree(children=venn.plot), top= A , bottom="subtitle")

#writing to file
tiff(filename = "Venn_diagram2.tiff", compression = "lzw")
grid.draw(venn.plot)
dev.off()


