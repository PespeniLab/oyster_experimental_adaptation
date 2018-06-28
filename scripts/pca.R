
library(pcadapt)
library(dplyr)
library(stringr)

mydata <- read.table("~/oyster/analysis/oyster.snp_rc", header=FALSE)
colnames(mydata) <- c("chr","pos","rc","allele_count","allele_states","deletion_sum","  snp_type"," major_alleles(maa)","minor_alleles(mia)","maa_1","maa_2","maa_3","maa_4","maa_5","maa_6","mia_1","mia_2","mia_3","mia_4","mia_5","mia_6")

nm <- c("maa_1","maa_2","maa_3","maa_4","maa_5","maa_6")
out <- as.data.frame(matrix(ncol=6, nrow=nrow(mydata)))
colnames(out) <- nm
for (i in nm){
    dat <- mydata[,grep(i, colnames(mydata))]
    num <- as.numeric(str_split_fixed(dat, "/", 2)[,1])
    denom <- as.numeric(str_split_fixed(dat, "/", 2)[,2])
    out[,i] <- num/denom
}

dat_t <- t(out)

filename <- read.pcadapt(dat_t,type="pool")
x <- pcadapt(filename,K=6) # calc actual pca
plot(x, option = "screeplot", pop = poplist.names)

poplist.names <- c("LL", "LL selected", "OP", "OP selected", "TB", "TB selected")

plot(x, option = "scores", i = 1, j = 2, pop = poplist.names)



# get proportion of total variance explained:
x$singular.values[1]/sum(x$singular.values)*100
x$singular.values[2]/sum(x$singular.values)*100


# jitter 

#x$scores[6,1] <- x$scores[6,1] - 0.015
#x$scores[2,2] <- x$scores[2,2] - 0.015
#x$scores[1,3] <- x$scores[1,3] - 0.015


sp <- c(21,21, 23,23,22,22)
all.col <- c("firebrick3", "firebrick3", "royalblue3", "royalblue3", "springgreen4", "springgreen4")
bg.col <-  c("firebrick3", NA,"royalblue3", NA,"springgreen4", NA)

png("~/oyster/figures/pca.png", res=300, height=85, width=85, units="mm")

#dev.new(width=3.35, height=3.35, units="mm")
par(mfrow = c(1, 1), mar=c(3, 3, 1, 1), mgp=c(3, 1, 0), las=0)
plot(y=x$scores[,1], x=x$scores[,2],
    pch=sp[as.factor(poplist.names)],
    cex=1.6,
    bg=bg.col[as.factor(poplist.names)],
    col=all.col[as.factor(poplist.names)],
    lwd=1.5,
    ylab="",
    xlab="",
    cex.lab=1.1, cex.axis=1,
    xaxt="n",yaxt="n"
)

axis(1, mgp=c(2, .4, 0), cex.axis=0.6) # second is tick mark labels
axis(2, mgp=c(2, .4, 0), cex.axis=0.6)

title(xlab="PC1: 33.0%", line=2, cex.lab=0.8)
title(ylab="PC2: 26.5%", line=2, cex.lab=0.8)

legend("topleft", pch=c(21,21, 23,23,22,22),
    pt.bg=c("firebrick3", NA,"royalblue3", NA,"springgreen4", NA),
    col=c("firebrick3","firebrick3", "royalblue3", "royalblue3", "springgreen4", "springgreen4"),
    legend=c("LL", "LL selected", "OP", "OP selected", "TB", "TB selected"),
    pt.cex=1.1, cex=0.6, pt.lwd=1.5)


dev.off()

