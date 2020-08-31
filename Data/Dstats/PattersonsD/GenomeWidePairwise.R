library(gplots)
library(RColorBrewer)

pairwise <- "/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/introgression/PattersonsD/summaries_20190920/thin_100bp_5_pairwise.txt"
pairwise <- read.table(pairwise,header=F, col.names=c("P1","CARD","LEW","P0","D"),stringsAsFactors = F)

cards <- unique(pairwise$CARD)
lews <- unique(pairwise$LEW)

pw.mat <- matrix(NA,nrow = length(cards), ncol=length(lews))

for (i in 1:length(cards)) {
    for (j in 1:length(lews)) {
        pw.mat[i,j] <- pairwise$D[pairwise$CARD == cards[i] & pairwise$LEW == lews[j]]
    }
}
colnames(pw.mat) <- lews
rownames(pw.mat) <- cards
YlOrRd <- brewer.pal(9,"YlOrRd")
colfunc <- colorRampPalette(YlOrRd)
col <- colfunc(16)
heatmap.2(t(pw.mat), col=col,trace="none")


### hists of chromosomes
chr2 <- "/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/abbababa/Chromosome2Pairwise.tsv"
    chr2 <- read.table(chr2,header=F, col.names=c("CARD","LEW","D"),stringsAsFactors = F)
chr5 <- "/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/abbababa/Chromosome5Pairwise.tsv"
    chr5 <- read.table(chr5,header=F, col.names=c("CARD","LEW","D"),stringsAsFactors = F)
    
chr2.mat <- matrix(NA,nrow = length(cards), ncol=length(lews))
    for (i in 1:length(cards)) {
        for (j in 1:length(lews)) {
            chr2.mat[i,j] <- chr2$D[chr2$CARD == cards[i] & chr2$LEW == lews[j]]
        }
    }
    colnames(chr2.mat) <- lews
    rownames(chr2.mat) <- cards

chr5.mat <- matrix(NA,nrow = length(cards), ncol=length(lews))
    for (i in 1:length(cards)) {
        for (j in 1:length(lews)) {
            chr5.mat[i,j] <- chr5$D[chr2$CARD == cards[i] & chr5$LEW == lews[j]]
        }
    }
    colnames(chr5.mat) <- lews
    rownames(chr5.mat) <- cards

chr2.dens <- density(chr2.mat)
chr5.dens <- density(chr5.mat)
gw.dens   <- density(pw.mat)
plot(x=0,y=0,type="n",xlab="Patterson's D",ylab="density",
     xlim = c(0,0.3), ylim = c(0,max(c(chr2.dens$y,chr5.dens$y,gw.dens$y))))
    polygon(c(gw.dens$x,rev(gw.dens$x)), 
            c(gw.dens$y,rep(0,length(gw.dens$y))), 
            lwd=2, col=adjustcolor('gray25',alpha=0.8),border="gray25")
    polygon(c(chr2.dens$x,rev(chr2.dens$x)), 
            c(chr2.dens$y,rep(0,length(chr2.dens$y))), 
            lwd=2, col=adjustcolor('darkcyan',alpha=0.8),border="darkcyan")
    polygon(c(chr5.dens$x,rev(chr5.dens$x)), 
            c(chr5.dens$y,rep(0,length(chr5.dens$y))), 
            lwd=2, col=adjustcolor('firebrick',alpha=0.8),border="firebrick")
    legend(x="topright", lwd=2,col=c("gray25","darkcyan","firebrick"),
           legend=c("genome-wide","chromosome 2", "YUP chromosome"))

