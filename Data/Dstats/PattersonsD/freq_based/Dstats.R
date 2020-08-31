###
### FUNCTIONS FOR COMPUTING ADMIXTURE STATS
###   IMPLEMENTED FROM PYTHON FUNCTIONS IN 
###   SIMON MARTIN'S GENOMICS.PY MODULE
###


    

    
f4 <- function (p1,p2,p3,p4) {
    f = (1 - p1)*p2*p3*(1-p4) - p1 * (1-p2)*p3*(1-p4)
    return(f)
}

f4_c <- function (p1,p2,p3,p4) {
    f = f4(p1,p2,p3,p4) + f4(1-p1,1-p2,1-p3,1-p4)
    return(f)
}

D_new <- function (p1,p2,p3,p4) {
    D = sum(f4_c(p1,p2,p3,p4)) / 
        (sum(
            (1 - p1)*p2*p3*(1-p4) + 
                p1 * (1-p2)*p3*(1-p4)+p1*(1-p2)*(1-p3)*p4 + 
                (1-p1) * p2*(1-p3)*p4))
    return(D)
}

fd_new <- function (p1,p2,p3,p4) {
    # define donor population frequency
    pd = p2* (p2>p3) + p3*(p3>=p2)
    fd = sum(f4_c(p1,p2,p3,p4)) / sum(f4_c(p1,pd,pd,p4))
    return(fd)
}

freqs <- read.table("/Users/thom/Dropbox/2.Fishman_Lab/Data/LewCard/popgenetics/introgression/PattersonsD/freq_based/freq_REF.5.tsv",
                    header =T)
freqs[freqs$bic == 1,] <- 1 - freqs[freqs$bic == 1,]
freqs$pos <- abs(freqs$pos)

windsize <- 100
stepsize <- 10
winstart <- 1

snpindex <- NULL
fdwindow <- NULL
fdwindow.lc <- NULL # use prev species tree
bpmid    <- NULL

pb <- txtProgressBar(1,dim(freqs)[1],char=" o_O ", style=3)
while (winstart < dim(freqs)[1]) {
    snpindex <- append(snpindex, winstart)
    snps     <- freqs[winstart:(winstart+windsize-1)  ,]
    bpmid    <- append(bpmid, mean(snps$pos))
    fdwind   <- fd_new(snps$verb,snps$ccar,snps$slew,snps$bic)
    fdwind.lc<- fd_new(snps$ccar,snps$slew,snps$verb,snps$bic)
    fdwindow <- append(fdwindow,fdwind)
    fdwindow.lc <- append(fdwindow.lc,fdwind.lc)
    winstart <- winstart + stepsize
    setTxtProgressBar(pb, winstart)
}

layout(matrix(1:2,2,1,byrow=T))
par(mar = c(4,4,1,1))
plot(snpindex,fdwindow.lc,type='l', ylim = c(-1,1))
plot(bpmid,fdwindow.lc,type='l',    ylim = c(-1,1))
layout(matrix(1:1,1,1,byrow=T))
