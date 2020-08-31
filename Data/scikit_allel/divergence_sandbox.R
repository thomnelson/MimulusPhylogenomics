maindir <- "/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/scikit_allel/"

chrs <- 1:8

dxy <- NULL

for (c in chrs) {
    d <- read.table(paste0(maindir,as.character(c),".dxy.tsv"), sep="\t",header=T)
    dxy <- rbind(dxy,d)
}
dxy$span <- dxy$stop - dxy$start

hist(dxy$dxy.ccarbic - dxy$dxy.slewccar,
     breaks = seq(-0.05,0.10,by=0.0025))

hist(dxy$dxy.ccarbic, xlim=c(0,0.08),
     breaks = seq(0.0,0.30,by=0.001),
     main="",xlab="dXY, card to bicolor")

hist(dxy$dxy.ccarbic[dxy$span > 1500], xlim=c(0,0.08),
     breaks = seq(0.0,0.30,by=0.001),
     main="",xlab="dXY, card to bicolor")

### RELATEDNESS

rel <- read.table(paste0(maindir,"diversity_GATK_CE10chromonome.8.SNPs.relatedness"),
                  sep="\t",header=T)
rel <- rel[rel$INDV1 != rel$INDV2,]
rel <- rel[order(rel$RELATEDNESS_AJK,decreasing=T),]

hist(rel$RELATEDNESS_AJK, breaks=100, xlim=c(-1,1))



