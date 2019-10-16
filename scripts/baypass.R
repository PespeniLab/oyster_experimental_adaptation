library(scales)


#######################################################
# identify loci associated with low salinity
#######################################################

#Signifi-cance was assessed based on the Bayes Factor (BF) between modelsand according to Jeffrey’s rule (1961); 
    #that is, markers with mdecisive evidence (BF > 20) were retained

covaux.snp.res=read.table("~/oyster/analysis/baypass/all_aux_summary_betai.out",h=T)

# db are in units of 10 × log10(BF)
length(which(covaux.snp.res$BF.dB > 20))
# [1] 847

# read in snp name
snp_in <- read.table("~/oyster/analysis/snpdet")

# read in outliers from PCAdapt:
pcad <- read.table("~/oyster/analysis/pcadapt_pc1_2.txt", header=TRUE)

# make data frame
df <- data.frame(CHR=snp_in$V1, POS=snp_in$V2, 
    BF = covaux.snp.res$BF.dB., snp_num=seq(1, nrow(snp_in)))
df$SNP <- paste(df$CHR, df$POS, sep=":")
df <- merge(df, pcad, by="SNP")

# pull out only significant
sig <- df[(which(df$BF > 20)),]

## save output

write.table(file="~/oyster/analysis/go_enrich/baypass_only.txt",
          sig,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=",")



# pull out those sig from baypass and sig from PCAdapt
sig <- df[(which(df$BF > 20  & df$qval < 0.05)),]

## save output

write.table(file="~/oyster/analysis/go_enrich/baypass_salinityresp.txt",
          sig,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=",")


##############################################################
############
############ pc 3
############
##############################################################

pc3 <- read.table("~/oyster/analysis/pcadapt_pc1_2_3.txt", header=TRUE)

colnames(pc3) <- c("CHR", "POS", "SNP", "pval_pc3", "qval_pc3", "PC3")

df <- merge(df, pc3, by="SNP")

length(which(df$BF > 20  & df$pval_pc3 != "NA"))
# 718
length(which(df$BF > 20))
# 847

sig <- df[(which(df$BF > 20  & df$qval_pc3 < 0.05 & df$PC3 == 3)),]

nrow(sig)
# 333, or 46%

write.table(file="~/oyster/analysis/go_enrich/baypass_pc3.txt",
          sig,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=",")


sig <- sig[(which(sig$qval < 0.05)),]

nrow(sig)
# 97


write.table(file="~/oyster/analysis/go_enrich/baypass_pc3_pc12.txt",
          sig,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=",")
