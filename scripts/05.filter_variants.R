# filter variants

## filtering raw variants

library(stringr)

dat <- read.table("~/oyster/analysis/varscan_out.txt", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/oyster/analysis/varscan_out.txt", stringsAsFactors=FALSE, nrows=1)

pops <- c("LL_33", "OP_33", "TB_33")

colnames(dat) <- c(datnm[1,1:10], pops)

nrow(dat)
# 558530
# first, remove all where the number of samples not covered/called is 0
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1)
#[1] 204021
# now, filter to include only variable sites
#dat2 <- dat1[which(dat1$SamplesRef < 6 & dat1$SamplesHom < 6 ),]

dat2 <- dat1

# filter by coverage:

depthout <- as.data.frame(matrix(nrow=nrow(dat2), ncol=length(pops)))
colnames(depthout) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat2[,grep(i_pop, colnames(dat2))]
        depth <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,2])
        # sum up reads
        depthout[,grep(i_pop, colnames(depthout))] <- depth

    }

head(depthout)
colMeans(depthout)
hist(as.matrix(depthout))
quantile(as.matrix(depthout), c(0.5, 0.975, 0.99))

hi_cv <- (apply(depthout, 1, function(x) {ifelse((length(which(x > 977)) > 0), FALSE, TRUE)}))
sum(hi_cv)
#[1] 196597
dat3 <- dat2[hi_cv,]

#many of these are skewed by indels. only keep reads where depth of actual bialleleic snps > 30
# from the manual: Also, VarScan reports variants on a biallelic basis.
    #That is, for a given SNP call, the "reads1" column is the number of
    #reference-supporting reads (RD), and the "reads2" column is the number of
    #variant-supporting reads (AD).
    #There may be additional reads at that position showing other bases (SNP or indel variants).
    #If these other variants meet the calling criteria, they will be reported in
    #their own line. If not, it may look like you have "missing" reads.
# columns for each call are:
    #consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
keep <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(keep) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])
        # sum up reads
        keep[,grep(i_pop, colnames(keep))] <- (maj+ minor)

    }

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 40)) > 0), FALSE, TRUE)}))

sum(low_cv)
#[1] 150435

dat4 <- dat3[low_cv,]
nrow(dat4)

dat3 <- dat4

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]

# save filtered genotypes

dat3 <- dat4
nrow(dat3)
#[1] 142369

# here calculate allele freqs
# columns for each call are: consensus genotype, total depth, num of read 1, num of read 2, allele freq, pvalue of snp.
af <- as.data.frame(matrix(nrow=nrow(dat3), ncol=length(pops)))
colnames(af) <- pops
    # cycle through each population
    for(i_pop in pops){

        tmp_pop <- dat3[,grep(i_pop, colnames(dat3))]
        maj <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,3])
        minor <- as.numeric(str_split_fixed(tmp_pop, ":", n=6)[,4])

        # calculate af
        af[,grep(i_pop, colnames(af))] <- maj/(maj+ minor)

    }


head(dat3[(which(rowSums(af) == 0)),])

# get rid of invariant sites. these are site different from the reference, but not variable in our pops
dat4 <- dat3[(which(rowSums(af) > 0)),]
af <- af[(which(rowSums(af) > 0)),]
nrow(dat4)
#142073

head(af)

# low maf cut off of < 0.05 in at least 4 groups.
## note that this corresponds to 2 reads at 40x, which seems reasonable.

low_maf <- c()
for(i in 1:nrow(af)){
    low_maf[i] <- ifelse(sum(af[i,] > 0.05 | af[i,] < 0.95) < 1, FALSE, TRUE)
}
sum(low_maf)

af.out <- (cbind(paste(dat4$Chrom, dat4$Position, sep=":"),af))
colnames(af.out) <- c("SNP", colnames(af))

dat4 <- dat4[low_maf,]
af.out <- af.out[low_maf,]
nrow(dat4)


write.table(file="~/oyster/analysis/filtered_variants.txt", dat4, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/oyster/analysis/filtered_allele_freqs.txt",
     af.out, sep="\t", row.names=FALSE, quote=FALSE)

# write chr and position to keep

write.table(file="~/oyster/analysis/snps.keep",
    cbind(dat4$Chrom, dat4$Position), sep="\t",
    row.names=FALSE, quote=FALSE, col.names=FALSE)
