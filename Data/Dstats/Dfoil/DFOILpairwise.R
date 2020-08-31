library(gplots)
library(RColorBrewer)
par(mar = c(4,4,1,1))
setwd("/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/introgression/Dfoil/")
dfoil <- read.table("Dfoil.tsv",header=T, stringsAsFactors = F)
dfoil <- read.table("DFOILpairwise.tsv",header=T, stringsAsFactors = F)
ncomparisons <- dim(dfoil)[1]
DFO <- dfoil$DFO
DIL <- dfoil$DIL
DFI <- dfoil$DFI
DOL <- dfoil$DOL

stat <- factor(c(rep("DFO",ncomparisons),
          rep("DIL",ncomparisons),
          rep("DFI",ncomparisons),
          rep("DOL",ncomparisons)), ordered=T, levels=c("DFO","DIL","DFI","DOL"))
value <- c(DFO,DIL,DFI,DOL)
forboxplot <- data.frame(stat,value)
boxplot(value ~ stat, forboxplot, ylim = c(-1,1))
    abline(h=0,lty="dotted")

cards <- unique(dfoil$P3)
lews <- unique(dfoil$P2)

dil.mat <- matrix(NA,nrow = length(cards), ncol=length(lews))
for (i in 1:length(cards)) {
    for (j in 1:length(lews)) {
        dil.mat[i,j] <- dfoil$DIL[dfoil$P3 == cards[i] & dfoil$P2 == lews[j]]
    }
}
colnames(dil.mat) <- lews
rownames(dil.mat) <- cards

dfo.mat <- matrix(NA,nrow = length(cards), ncol=length(lews))
for (i in 1:length(cards)) {
    for (j in 1:length(lews)) {
        dfo.mat[i,j] <- dfoil$DFO[dfoil$P3 == cards[i] & dfoil$P2 == lews[j]]
    }
}
colnames(dfo.mat) <- lews
rownames(dfo.mat) <- cards

YlOrRd <- brewer.pal(9,"YlOrRd")
colfunc <- colorRampPalette(YlOrRd)
col <- colfunc(16)

heatmap.2(t(dil.mat), col=col,trace="none")
heatmap.2(t(dfo.mat), col=col,trace="none")


### explore covariates for DIL

MBpersample <- read.table("MBperSample.txt",header=T,stringsAsFactors = F)
het         <- read.table("5.het",header=T,stringsAsFactors = F)

layout(matrix(1:4,2,2,byrow=T))
par(mar = c(4,4,0,1))
plot(DFO ~ allsites, dfoil, pch = 20)
plot(DIL ~ allsites, dfoil, pch = 20)
plot(DFI ~ allsites, dfoil, pch = 20)
plot(DOL ~ allsites, dfoil, pch = 20)

### card:
cardcovar <- matrix(nrow=length(cards),ncol=3)
for (c in 1:length(cards)) {
    card <- cards[c]
    MB  <- MBpersample$MB[MBpersample$sample==card]
    hom <- het$O.HOM.[het$INDV==card] / het$N_SITES[het$INDV==card]
    avg <- mean(dfoil$DIL[dfoil$P3 == card])
    cardcovar[c,] <- c(MB,avg,hom)
}
cardcovar <- data.frame(cardcovar)
names(cardcovar) <- c("MB","DILavg","Ohomo")
### lew
lewcovar <- matrix(nrow=length(lews),ncol=3)
for (c in 1:length(lews)) {
    lew <- lews[c]
    MB  <- MBpersample$MB[MBpersample$sample==lew]
    hom <- het$O.HOM.[het$INDV==lew] / het$N_SITES[het$INDV==lew]
    avg <- mean(dfoil$DIL[dfoil$P2 == lew])
    lewcovar[c,] <- c(MB,avg,hom)
}
lewcovar <- data.frame(lewcovar)
names(lewcovar) <- c("MB","DILavg","Ohomo")

### plot covariates
plot(DILavg ~ MB,cardcovar, pch=20)
    plot(DILavg ~ Ohomo,cardcovar, pch=20)
    plot(DILavg ~ MB,lewcovar, pch=20)
    plot(DILavg ~ Ohomo,lewcovar, pch=20)
layout(matrix(1:1,1,1))
