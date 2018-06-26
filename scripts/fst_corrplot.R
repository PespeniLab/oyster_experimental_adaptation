#fst_corrplot.R
# plot fst in corr plot. for visualization


library(stringr)
library(RColorBrewer)
#library(qqman)
library(data.table)

fst <- as.data.frame(fread("~/oyster/analysis/oyster.fst", header=FALSE))
#fst.win <- as.data.frame(fread("~/oyster/analysis/oyster.window.fst", header=FALSE))


#name columns
p.names <- c("CHR", "BP", "snp_num", "frac_cov", "min_cov",
    "LL_33_vs_LL_05", "LL_33_vs_OP_33", "LL_33_vs_OP_05", "LL_33_vs_TB_33", "LL_33_vs_TB_05",
    "LL_05_vs_OP_33", "LL_05_vs_OP_05", "LL_05_vs_TB_33", "LL_05_vs_TB_05",
    "OP_33_vs_OP_05", "OP_33_vs_TB_33", "OP_33_vs_TB_05",
    "OP_05_vs_TB_33", "OP_05_vs_TB_05",
    "TB_33_vs_TB_05")
colnames(fst) <- p.names

#change chr to numbers
fst$CHR <- as.numeric(gsub("Contig_", "", fst$CHR))

# convert to numeric. remove formatting and pop labs, etc.

for(i in 6:length(p.names)){
    fst[,p.names[i]] <- as.numeric(str_split_fixed(fst[,p.names[i]], "=", 2)[,2])
}


#calculate mean fst
fst_mean <- as.data.frame(matrix(nrow=1, ncol=15))
colnames(fst_mean) <- p.names[6:20]


for(i in 6:length(p.names)){
    fst_mean[,p.names[i]] <- mean(fst[,p.names[i]])
}

#need to make matrix of fst values for corrplot
name <- unique(substr(colnames(fst_mean), 1,5))
name <- c(name, "TB_05")
values <- as.data.frame(matrix(nrow=6, ncol=6))
rownames(values) <- name
colnames(values) <- name

for (i in colnames(fst_mean)){
    pop1 <- substr(i, 1, 5)
    pop2 <- substr(i, 10, 14)
    col <- which( colnames(values) == pop1 )
    row <- which( rownames(values) == pop2 )
    values[row,col]<- fst_mean[1,i]
    col <- which( colnames(values) == pop2 )
    row <- which( rownames(values) == pop1 )
    values[row,col]<- fst_mean[1,i]
}

values <- sapply(values, round, 3)
values[is.na(values)] <- 0
rownames(values) <- colnames(values)


png("~/oyster/figures/fst.corr.png", , res=300, height=100, width=100, units="mm")
par(mfrow = c(1, 1), mar=c(3, 3, 1, 1), mgp=c(3, 1, 0), las=0)

corrplot(values, method="color", type="upper",cl.lim=c(0,max(values)),
    col=colorRampPalette(c("red","white","blue"))(200),
    is.corr=FALSE, addCoef.col="black",
    tl.col="black", cl.pos ="n",
    diag=FALSE, tl.cex = 1,tl.srt=45,
    number.digits = 3)

dev.off()

