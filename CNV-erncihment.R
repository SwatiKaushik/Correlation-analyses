# CNV Enrichment analysis of interactors 
# Swati Kaushik Nov 21 2016

rm(list=ls())
library(VennDiagram)

args <- commandArgs(TRUE)

#Amplified/deleted cnvs from TCGA
cnv.file <- read.table(args[1], header= F, sep="\t", fill=TRUE)
HS.cnv.file <- levels(droplevels(cnv.file$V1))

#interaction file
interaction.file <- read.table(args[2], header=T)
bait.unq <- unique(interaction.file$Bait)
bait.unq.d <- levels(droplevels(bait.unq))
bait.subset<- subset(interaction.file, interaction.file$zscore>2, drop=FALSE)
preys <- levels(droplevels(bait.subset$Prey))
preys.noTK <- preys[grep(paste(bait.unq.d, collapse="|"), preys, invert=TRUE)]
#preys.noTK

#calculate overlap between args[1] and args[2]
overlap <- intersect(preys.noTK, HS.cnv.file)
overlap.length <- length(overlap)
preys.length <- length(preys.noTK)
ccg.length <- length(HS.cnv.file)

#calculate significance
significance.test <- rbind(c(overlap.length, preys.length), c(ccg.length, 20000))
significance.pvalue <- fisher.test(significance.test, alternative = "greater")$p.value
significance.pvalue

write.table(overlap, file="overlapped-cnvs.txt", sep="\t")
write.table(preys.length, file="overlapped-cnvs.txt", sep="\t", append = TRUE)
write.table(ccg.length, file="overlapped-cnvs.txt", sep="\t", append = TRUE)

#plot venn diagram

grid.newpage()
venn.plot <- draw.pairwise.venn( area1 = preys.length, 
					area2 = ccg.length, 
					cross.area = overlap.length, 
					category =c ("interactome","cnvs-PANCAN12"),
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