### BLOCK JACKNIFE PATTERSON'S D

siteD <- "/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/introgression/PattersonsD/PatsD_CAR-C27_SNO8.tsv"
siteD <- read.table(siteD, header=T)
blocksize <- 1000000

pattersonD <- function (countABBA,countBABA) {
    return((countABBA- countBABA) / (countABBA + countBABA))
}

### get max values for each unique chromosome
chroms <- unique(siteD$CHR)
maxbps <- rep(0, length(chroms))
for (chrom in 1:length(chroms)) {
    bpmax <- max(siteD$POS[siteD$CHR == chroms[chrom]], na.rm=T)
    maxbps[chrom] <- bpmax
}

Djk <- NULL
for (i in 1:length(chroms)) {
    chrom <- chroms[i]
    blockstarts <- seq(1, maxbps[i], by=blocksize)
    blockends   <- blockstarts + (blocksize - 1)
    for (j in 1:length(blockstarts)) {
        blockindex <- siteD$CHR == chrom & 
                      siteD$POS >= blockstarts[j] & 
                      siteD$POS <= blockends[j]
        Dnoblock <- siteD[!(blockindex),]
        abba = sum(Dnoblock$ABBA)
        baba = sum(Dnoblock$BABA)
        Djk <- append(Djk, pattersonD(abba,baba))
    }
}
