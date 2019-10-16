# pcadapt


library(pcadapt)
library(dplyr)
library(stringr)
library(scales)
library(qvalue)

af <- read.table("~/oyster/analysis/filtered_allele_freqs.txt", header=TRUE)
dat_t <- t(af[,2:7])


filename <- read.pcadapt(dat_t,type="pool")

res <- pcadapt(filename,min.maf=0.001, K=6) # calc actual pca

plot(res, option = "screeplot")

plot(res, option = "scores", pop = colnames(af[2:length(colnames(af))]))
plot(res, option = "scores", pop = colnames(af[2:length(colnames(af))]), i = 1, j = 3,)


# first, what is the signal from population structure:

res <- pcadapt(filename,min.maf=0.025, K = 2) # calc actual pca

summary(res)

plot(res , option = "manhattan")
plot(res, option = "qqplot", threshold = 0.05)
hist(res$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(res, option = "stat.distribution")


# outliers 

padj <- qvalue(res$pvalues)$qvalue
alpha <- 0.05
outlier.pval <- padj[which(padj < alpha)]
outliers <- which(padj < alpha)
length(outliers)
#[1] 14619


snp_pc <- get.pc(res, outliers)

length(which(snp_pc$PC == "1"))
# [1] 8173
length(which(snp_pc$PC == "2"))
# [1] 7592

outliers_1 <- snp_pc[which(snp_pc$PC == "1"),]
outliers_2 <- snp_pc[which(snp_pc$PC == "2"),]


# write output file of results

chr <- sapply(strsplit(as.character(af$SNP), ":", fixed=TRUE), "[", 1)
pos <- sapply(strsplit(as.character(af$SNP), ":", fixed=TRUE), "[", 2)


out_txt <- data.frame(CHR = chr, POS= pos, SNP= af$SNP, pval=res$pvalues, qval = padj)
out_txt$PC <- NA

out_txt$PC <- snp_pc$PC[match(row.names(out_txt), snp_pc$SNP)]

head(out_txt)

out_txt[which(out_txt$qval < 0.05),]


write.table(file="~/oyster/analysis/pcadapt_pc1_2.txt", out_txt, quote=FALSE, row.names=FALSE)



##################
#
# with PC 3
#
##################


res_3 <- pcadapt(filename,min.maf=0.025, K = 3) # calc actual pca

summary(res_3)

padj <- qvalue(res_3$pvalues)$qvalue
alpha <- 0.05
outlier.pval <- padj[which(padj < alpha)]
outliers <- which(padj < alpha)
length(outliers)
# [1] 15765

sum(!is.na(res_3$pvalues))

head(res_3$maf)
head(res_3$pvalues)


snp_pc <- get.pc(res_3, outliers)

length(which(snp_pc$PC == "1"))
# 6189
length(which(snp_pc$PC == "2"))
# [1] 5488
length(which(snp_pc$PC == "3"))
# [1] 4088

length(which(snp_pc$PC == "1" | snp_pc$PC == "2"  ))

# write output file of results

chr <- sapply(strsplit(as.character(af$SNP), ":", fixed=TRUE), "[", 1)
pos <- sapply(strsplit(as.character(af$SNP), ":", fixed=TRUE), "[", 2)


out_txt <- data.frame(CHR = chr, POS= pos, SNP= af$SNP, pval=res_3$pvalues, qval = padj)
out_txt$PC <- NA

out_txt$PC <- snp_pc$PC[match(row.names(out_txt), snp_pc$SNP)]

head(out_txt)

out_txt[which(out_txt$qval < 0.05),]


write.table(file="~/oyster/analysis/pcadapt_pc1_2_3.txt", out_txt, quote=FALSE, row.names=FALSE)
