
library(poolfstat)
setwd("~/oyster/analysis")
pdat <- popsync2pooldata(sync.file="variants.sync", poolsizes = c(40,40,30),
    poolnames = c("LL_33", "OP_33","TB_33"), min.maf=0, min.rc=0,
     max.cov.per.pool = 1e+07)

# Data consists of 136448 SNPs for 6 Pools

pooldata2genobaypass(pdat, writing.dir = "~/oyster/analysis", prefix = "",
  subsamplesize = -1)# note that -1 means don't subsample

fst <- computePairwiseFSTmatrix(pdat, method = "Anova",
  min.cov.per.pool = -1, max.cov.per.pool = 1e+07, min.maf = -1,
  output.snp.values = FALSE)

fst.all <- computePairwiseFSTmatrix(pdat, method = "Anova",
  min.cov.per.pool = -1, max.cov.per.pool = 1e+07, min.maf = -1,
  output.snp.values = TRUE)

head(fst.all[[3]])

pairwise <- as.data.frame(fst.all[[3]])

pairwise[is.na(pairwise)] <- 0
pairwise[pairwise < 0] <- 0

colMeans(pairwise)

write.table(file="~/oyster/analysis/fst_pairwise.txt", pairwise, sep="\t",
              row.names=FALSE, quote=FALSE)

#### fst corr plot

#fst_corrplot.R
# plot fst in corr plot. for visualization


library(stringr)
library(RColorBrewer)
#library(qqman)
library(data.table)
library(corrplot)

pairwise <- read.table("~/oyster/analysis/fst_pairwise.txt", header=TRUE)

#calculate mean fst
fst_mean <- as.data.frame(matrix(nrow=1, ncol=3))
colnames(fst_mean) <- colnames(pairwise)
fst_mean <- colMeans(pairwise)

#need to make matrix of fst values for corrplot
name <- unique(substr(names(fst_mean), 1,2))
name <- c(name, "TB")
values <- as.data.frame(matrix(nrow=3, ncol=3))
rownames(values) <- name
colnames(values) <- name

for (i in names(fst_mean)){
    pop1 <- substr(i, 1, 2)
    pop2 <- substr(i, 10, 11)
    col <- which( colnames(values) == pop1 )
    row <- which( rownames(values) == pop2 )
    values[row,col]<- fst_mean[i]
    col <- which( colnames(values) == pop2 )
    row <- which( rownames(values) == pop1 )
    values[row,col]<- fst_mean[i]
}

values <- sapply(values, round, 3)
values[is.na(values)] <- 0
rownames(values) <- colnames(values)


pdf("~/oyster/figures/fst.corr.pdf",height=3.9, width=3.9)
par(mfrow = c(1, 1), mar=c(3, 3, 1, 1), mgp=c(3, 1, 0), las=0)

corrplot(values, method="color", type="upper",cl.lim=c(0,max(values)),
    col=colorRampPalette(c("red","white","blue"))(200),
    is.corr=FALSE, addCoef.col="black",
    tl.col="black", cl.pos ="n",
    diag=FALSE, tl.cex = 1,tl.srt=45,
    number.digits = 3)

dev.off()

