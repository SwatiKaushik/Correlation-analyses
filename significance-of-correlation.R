# calculate the significance of difference in correlation. It performs fischer z tranformation of the correlation data

rm(list=ls())
args <- commandArgs(TRUE)

diff.corr <- function( r1, n1, r2, n2, gene1.gene2 ){

    Z1 <- 0.5 * log( (1+r1)/(1-r1) )
    Z2 <- 0.5 * log( (1+r2)/(1-r2) )

    diff   <- Z1 - Z2
    SEdiff <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) )
    diff.Z  <- diff/SEdiff

    p.value <- 2*pnorm( abs(diff.Z), lower=F)
    cat( gene1.gene2, p.value , "\n" )

  }


data.file = read.table(args[1] , header=T, sep="\t")

for (i in 1:nrow(data.file)){
	
	gene1 <- levels(droplevels(data.file[i,2]))
	gene2 <- levels(droplevels(data.file[i,3]))
	gene1.gene2 <- (paste(gene1, gene2))
	p.value <-	diff.corr(r1=data.file[i,4], n1=1251, r2 = data.file[i,6], n2=109, gene1.gene2)

}	

