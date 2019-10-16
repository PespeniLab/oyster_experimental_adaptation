# filter variants

## filtering raw variants

library(stringr)

dat <- read.table("~/oyster/analysis/varscan_out.txt", stringsAsFactors=FALSE, skip=1)

datnm <- read.table("~/oyster/analysis/varscan_out.txt", stringsAsFactors=FALSE, nrows=1)

pops <- c("LL_33", "LL_05", "OP_33", "OP_05", "TB_33", 'TB_05')

colnames(dat) <- c(datnm[1,1:10], pops)

nrow(dat)
# 906823
# first, remove all where the number of samples not covered/called is 0
dat1 <- dat[which(dat$SamplesNC == 0),]
nrow(dat1)
#[1] 178388
# now, filter to include only variable sites
#dat2 <- dat1[which(dat1$SamplesRef < 6 & dat1$SamplesHom < 6 ),]

dat2 <- dat1


# filter by coverage:
filt <- (median(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2])))*3
# median coverage is 593
dat3 <- dat2[(which(as.numeric((str_split_fixed(dat2[,5] , ":", n=6))[,2]) < filt)),]
nrow(dat3)
#[1] 157665


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

low_cv <- (apply(keep, 1, function(x) {ifelse((length(which(x < 30)) > 0), FALSE, TRUE)}))

sum(low_cv)
#[1] 155542

dat4 <- dat3[low_cv,]
nrow(dat4)

dat3 <- dat4
nrow(dat3)
#[1] 155542

# to get rid of multiallele sites or deletions/insertions
dat4 <- dat3[(which(as.numeric(lapply(strsplit(dat3[,4],","), length)) == 1)),]

# save filtered genotypes

dat3 <- dat4
nrow(dat3)
#[1] 145671

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

# get rid of invariant sites
dat4 <- dat3[(which(rowSums(af) > 0)),]
dat3 <- dat4
af <- af[(which(rowSums(af) > 0)),]
nrow(dat3)
#145396

# get ride of invariant sites in 33 treatments

keep.out <- (which(af$LL_33 < 1 & af$LL_33 > 0 | af$OP_33 < 1 & af$OP_33 > 0 | af$TB_33 < 1 & af$TB_33 > 0 ))
#keep.out <- (which(af$LL_33 < 1 & af$LL_33 > 0.025 | af$OP_33 < 1 & af$OP_33 > 0.025 | af$TB_33 < 1 & af$TB_33 > 0.025 ))

af <- af[keep.out,]
dat3 <- dat3[keep.out,]


af.out <- (cbind(paste(dat3$Chrom, dat3$Position, sep=":"),af))

colnames(af.out) <- c("SNP", colnames(af))

nrow(af.out)
nrow(dat3)
#[1] 136448
#[1] 136314

write.table(file="~/oyster/analysis/filtered_variants.txt", dat3, sep="\t",
              row.names=FALSE, quote=FALSE)

write.table(file="~/oyster/analysis/filtered_allele_freqs.txt",
     af.out, sep="\t", row.names=FALSE, quote=FALSE)

# write chr and position to keep

write.table(file="~/oyster/analysis/snps.keep",
    cbind(dat3$Chrom, dat3$Position), sep="\t",
    row.names=FALSE, quote=FALSE, col.names=FALSE)
