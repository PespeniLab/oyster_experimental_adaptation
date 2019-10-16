

use the dadi pipeline from Portik: https://github.com/dportik/dadi_pipeline

Cite the following pub:

    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O., Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017. Evaluating mechanisms of diversification in a Guineo-Congolian forest frog using demographic model selection. Molecular Ecology 26: 5245-5263. https://doi.org/10.1111/mec.14266


Note that for all of this it was right as dadi was being updated to python3 and the pipeline was largely python2 based. As a result, these are all python2. The pipeline was downloaded from github June 10, 2019. 


```r

#thin ~/oyster/analysis/baypass.33.geno.in so only has 1 snp per contig

bp.in <- read.table("~/oyster/analysis/baypass.33.geno.in")
snp_in <- read.table("~/oyster/analysis/snpdet")

head(snp_in)

contig <- (unique(snp_in$V1))

df.out <- as.data.frame(matrix(ncol=12, nrow=length(contig)))

colnames(df.out) <- c("Ingroup", "Outgroup",
                    "Allele1", "LL1", "OP1", "TB1",
                    "Allele2", "LL2", "OP2", "TB2",
                    "GeneID", "Position" )

ct.1 <- 0
for(i in 1:length(contig)){
  tmp_idx <- which(snp_in$V1 == contig[i])
  if(length(tmp_idx) > 1){
      i.med <- round(median(as.numeric(tmp_idx)))
      df.out[i,c(4,5,6)] <- bp.in[i.med,c(1,3,5)]
      df.out[i,c(8,9,10)] <- bp.in[i.med,c(2,4,6)]
      df.out$Allele1[i] <- as.character(snp_in$V4[i.med])
      df.out$Allele2[i] <- as.character(snp_in$V5[i.med])
      df.out$GeneID[i] <- as.character(snp_in$V1[i.med])
      df.out$Position[i] <- snp_in$V2[i.med]

  }

  else{
    i.med <- tmp_idx
      df.out[i,c(4,5,6)] <- bp.in[i.med,c(1,3,5)]
      df.out[i,c(8,9,10)] <- bp.in[i.med,c(2,4,6)]
      df.out$Allele1[i] <- as.character(snp_in$V4[i.med])
      df.out$Allele2[i] <- as.character(snp_in$V5[i.med])
      df.out$GeneID[i] <- as.character(snp_in$V1[i.med])
      df.out$Position[i] <- snp_in$V2[i.med]
      df.out$Ingroup[i] <- c("---")
      df.out$Outgroup[i] <- c("---")

    print(paste("only 1!: total is ", ct.1, " out of total snps: ", i))
    ct.1 <- ct.1 + 1
  }

  if(i%%1000 == 0){print(i)}
}


write.table(file="~/oyster/analysis/dadi/dadi.thinned.in",
          df.out,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=" ")


df.out <- as.data.frame(matrix(ncol=12, nrow=nrow(bp.in)))

colnames(df.out) <- c("Ingroup", "Outgroup",
                    "Allele1", "LL1", "OP1", "TB1",
                    "Allele2", "LL2", "OP2", "TB2",
                    "GeneID", "Position" )

df.out[,c(4,5,6)] <- bp.in[,c(1,3,5)]
df.out[,c(8,9,10)] <- bp.in[,c(2,4,6)]
df.out$Allele1 <- as.character(snp_in$V4)
df.out$Allele2 <- as.character(snp_in$V5)
df.out$GeneID <- as.character(snp_in$V1)
df.out$Position <- snp_in$V2


write.table(file="~/oyster/analysis/dadi/dadi.all.in",
          df.out,
          col.names=TRUE, row.names=FALSE, quote=FALSE,sep=" ")


```

```bash

# clean up the input file
sed -i 's/NA/---/g' ~/oyster/analysis/dadi/dadi.thinned.in

sed -i 's/LL2/LL/g' ~/oyster/analysis/dadi/dadi.thinned.in
sed -i 's/OP2/OP/g' ~/oyster/analysis/dadi/dadi.thinned.in
sed -i 's/TB2/TB/g' ~/oyster/analysis/dadi/dadi.thinned.in

sed -i 's/LL1/LL/g' ~/oyster/analysis/dadi/dadi.thinned.in
sed -i 's/OP1/OP/g' ~/oyster/analysis/dadi/dadi.thinned.in
sed -i 's/TB1/TB/g' ~/oyster/analysis/dadi/dadi.thinned.in

# parse down to only LL and OP
cut -f 1,2,3,4,5,7,8,9,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.thinned.in > ~/oyster/analysis/dadi/dadi.llop.in
cut -f 1,2,3,4,6,7,8,10,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.thinned.in > ~/oyster/analysis/dadi/dadi.lltb.in
cut -f 1,2,3,5,6,7,9,10,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.thinned.in > ~/oyster/analysis/dadi/dadi.optb.in

# repeat for all

sed -i 's/NA/---/g' ~/oyster/analysis/dadi/dadi.all.in

sed -i 's/LL2/LL/g' ~/oyster/analysis/dadi/dadi.all.in
sed -i 's/OP2/OP/g' ~/oyster/analysis/dadi/dadi.all.in
sed -i 's/TB2/TB/g' ~/oyster/analysis/dadi/dadi.all.in

sed -i 's/LL1/LL/g' ~/oyster/analysis/dadi/dadi.all.in
sed -i 's/OP1/OP/g' ~/oyster/analysis/dadi/dadi.all.in
sed -i 's/TB1/TB/g' ~/oyster/analysis/dadi/dadi.all.in

# parse down to only LL and OP

#cut -f 1,2,3,4,5,7,8,9,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.all.in > ~/oyster/analysis/dadi/dadi.all.llop.in
#cut -f 1,2,3,4,6,7,8,10,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.all.in > ~/oyster/analysis/dadi/dadi.all.lltb.in
#cut -f 1,2,3,5,6,7,9,10,11,12 -d ' ' ~/oyster/analysis/dadi/dadi.all.in > ~/oyster/analysis/dadi/dadi.all.optb.in



```


## run dadi to get sfs

just checking things out here.

```python
import numpy
import dadi
#import matplotlib
import matplotlib.pyplot as pyplot
# convert snp data to FS

## first, make the snp data into python dictionary
filepath = '/users/r/b/rbrennan/oyster/analysis/dadi/dadi.llop.in'

dd = dadi.Misc.make_data_dict(filepath)

## then make into FS
### projects is the population sample sizes for the FS- remember 2 per iniv (because diploid)
##### this will down sample reads for each snp to the sample size given here.
### polarized means unfolded or folded. False gives folded.
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['LL', 'OP'],
                            projections=[40, 40],
                            polarized=False)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['LL', 'TB'],
                            projections=[40, 40],
                            polarized=False)

# Frequency spectra are stored in dadi.Spectrum object
## this is similar to numpy.masked_array
pyplot.figure() # to generate empty fig
dadi.Plotting.plot_single_2d_sfs(fs)
pyplot.show()


fsmarg = fs.marginalize([1]) # remove pop 2 keep pop1

pyplot.figure() # to generate empty fig
dadi.Plotting.plot_1d_fs(fsmarg)
pyplot.show()

filepath = '/users/r/b/rbrennan/oyster/analysis/dadi/dadi.optb.in'

dd = dadi.Misc.make_data_dict(filepath)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['OP', 'TB'],
                            projections=[40, 40],
                            polarized=False)

# Frequency spectra are stored in dadi.Spectrum object
## this is similar to numpy.masked_array
pyplot.figure() # to generate empty fig
dadi.Plotting.plot_single_2d_sfs(fs)
pyplot.show()


fsmarg = fs.marginalize([1]) # remove pop 2 keep pop1

pyplot.figure() # to generate empty fig
dadi.Plotting.plot_1d_fs(fsmarg)
pyplot.show()


filepath = '/users/r/b/rbrennan/oyster/analysis/dadi/dadi.lltb.in'

dd = dadi.Misc.make_data_dict(filepath)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['OP', 'TB'],
                            projections=[40, 40],
                            polarized=False)

# Frequency spectra are stored in dadi.Spectrum object
## this is similar to numpy.masked_array
pyplot.figure() # to generate empty fig
dadi.Plotting.plot_single_2d_sfs(fs)
pyplot.show()


fsmarg = fs.marginalize([1]) # remove pop 2 keep pop1

pyplot.figure() # to generate empty fig
dadi.Plotting.plot_1d_fs(fsmarg)
pyplot.show()


```

## run optimization

These are broken up just to speed things up

```bash

nohup python2 ~/oyster/analysis/dadi/dadi_pipeline/optimize_1.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_1.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_1.stdout_$(date +"%F_%R").txt &

echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_1.pid

nohup python2 ~/oyster/analysis/dadi/dadi_pipeline/optimize_2.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_2.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_2.stdout_$(date +"%F_%R").txt &
echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_2.pid


nohup python2 ~/oyster/analysis/dadi/dadi_pipeline/optimize_3.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_3.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_3.stdout_$(date +"%F_%R").txt &
echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_3.pid



##### for LL and TB

#LL vs TB

nohup python ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_1.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_1.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_1.stdout_$(date +"%F_%R").txt &
echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBOP_1.pid

nohup python ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_2.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_2.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_2.stdout_$(date +"%F_%R").txt &
echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_2.pid

nohup python ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_3.py 2> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_3.stderr_$(date +"%F_%R").txt 1> ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_3.stdout_$(date +"%F_%R").txt &
echo $! > ~/oyster/analysis/dadi/dadi_pipeline/optimize_TBLL_3.pid

```


## check model fits

get and save the best model parameters for each analysis:

```r

nomig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_OP.no_mig.optimized.txt", header=T, sep="\t")
asym_mig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_OP.asym_mig.optimized.txt", header=T, sep="\t")
sym_mig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_OP.sym_mig.optimized.txt", header=T, sep="\t")

nomig[which.max(nomig$log.likelihood),]
nomig[which.min(nomig$AIC),]
# 0.142,0.1736,0.0142
# log.likelihood: -1583.95  AIC: 3173.9

asym_mig[which.max(asym_mig$log.likelihood),]
asym_mig[which.min(asym_mig$AIC),]
# log.likelihood: -1515.31  AIC: 3040.62
# 0.38,0.4034,4.1025,6.0181,0.1147
# nu1, nu2, m12, m21, T

# pop1 = LL, pop2= OP

sym_mig[which.max(sym_mig$log.likelihood),]
sym_mig[which.min(sym_mig$AIC),]
# log.likelihood: -1518.32  AIC: 3044.64
# 0.8941,1.0828,2.1383,4.4325
# nu1, nu2, m, T

nomig$group <- substr(nomig$Replicate, 1, 7)

nomig <- nomig[-c(13,14,15,16),]

boxplot(as.numeric(as.character(nomig$log.likelihood)) ~ as.factor(nomig$group), col="grey", ylim=c(-3000, -1500))


### then do ll and tb
nomig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_TB.no_mig.optimized.txt", header=T, sep="\t")
asym_mig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_TB.asym_mig.optimized.txt", header=T, sep="\t")
sym_mig <- read.csv("~/oyster/analysis/dadi/dadi_pipeline/LL_TB.sym_mig.optimized.txt", header=T, sep="\t")


nomig[which.max(nomig$log.likelihood),]
nomig[which.min(nomig$AIC),]
# 0.221,0.2503,0.0371
# log.likelihood:  -1985.41  AIC:3976.82

asym_mig[which.max(asym_mig$log.likelihood),]
asym_mig[which.min(asym_mig$AIC),]
# log.likelihood: -1805.33  AIC: 3620.66
# 1.5248,2.3771,0.8413,0.3438,10.5383
# nu1, nu2, m12, m21, T

# pop1 = LL, pop2= OP

sym_mig[which.max(sym_mig$log.likelihood),]
sym_mig[which.min(sym_mig$AIC),]
# log.likelihood: -1818.51  AIC: 3645.02
# 3.0209,3.3259,0.3463,15.2636
# nu1, nu2, m, T


```

### look at the actual model fits. residuals, etc.

```python
#import copy_reg as copyreg
import numpy
import dadi
import matplotlib
import matplotlib.pyplot as pyplot
import Plotting_Functions
import Optimize_Functions
import Models_2D
# convert snp data to FS
## first, make the snp data into python dictionary
filepath = '/Users/reidbrennan/Documents/UVM/Oyster/dadi_analysis/dadi.llop.in'

dd = dadi.Misc.make_data_dict(filepath)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['LL', 'OP'],
                            projections=[40, 40],
                            polarized=False)


# Frequency spectra are stored in dadi.Spectrum object
## this is similar to numpy.masked_array
pyplot.figure() # to generate empty fig
dadi.Plotting.plot_single_2d_sfs(fs)
pyplot.show()


# plotting using the dadi pipeline

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
pop_ids=['LL', 'OP']
prefix = "_".join(pop_ids)

#Make sure to define your extrapolation grid size.
pts = [200,220,240]


### NO MIGRATION

#Provide best optimized parameter set for empirical data.
#These will come from previous analyses you have already completed

#Fit the model using these parameters and return the model SFS.
#Here, you will want to change the "sym_mig" and sym_mig arguments to match your model function,
#but everything else can stay as it is. See above for argument explanations.

# fs: spectrum object name
# pts: grid size for extrapolation, list of three values
# outfile: prefix for output naming
# model_name: a label help name the output files; ex. "sym_mig"
# func: access the model function from within this script or from a separate model script
# in_params: the previously optimized parameter values to use
# fs_folded: A Boolean value indicating whether the empirical fs is folded (True) or not (False)
emp_params = [0.142,0.1736,0.0142]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "no_mig", Models_2D.no_mig, emp_params, fs_folded=True)

vmin_val = float(0.01)
Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig", vmin_val = vmin_val)


#### SYM MIG

emp_params = [0.8941,1.0828,2.1383,4.4325]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, emp_params, fs_folded=True)

Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig",vmin_val = vmin_val)

#### ASYM MIG

emp_params = [0.38,0.4034,4.1025,6.0181,0.1147]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, emp_params, fs_folded=True)

Plotting_Functions.Plot_2D(fs, model_fit, prefix, "asym_mig",vmin_val = vmin_val)

#############################
#############################
# rerun for ll TB
#############################
#############################
filepath = '/Users/reidbrennan/Documents/UVM/Oyster/dadi_analysis/dadi.lltb.in'

dd = dadi.Misc.make_data_dict(filepath)

fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['LL', 'TB'],
                            projections=[40, 40],
                            polarized=False)


# Frequency spectra are stored in dadi.Spectrum object
## this is similar to numpy.masked_array
pyplot.figure() # to generate empty fig
dadi.Plotting.plot_single_2d_sfs(fs)
pyplot.show()

# plotting using the dadi pipeline

#create a prefix based on the population names to label the output files
#ex. Pop1_Pop2
pop_ids=['LL', 'OP']
prefix = "_".join(pop_ids)

#Make sure to define your extrapolation grid size.
pts = [200,220,240]


### NO MIGRATION

emp_params = [0.221,0.2503,0.0371]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "no_mig", Models_2D.no_mig, emp_params, fs_folded=True)

vmin_val = float(0.01)
Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig", vmin_val = vmin_val)


#### SYM MIG

emp_params = [3.0209,3.3259,0.3463,15.2636]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "sym_mig", Models_2D.sym_mig, emp_params, fs_folded=True)

Plotting_Functions.Plot_2D(fs, model_fit, prefix, "sym_mig",vmin_val = vmin_val)

#### ASYM MIG

emp_params = [1.5248,2.3771,0.8413,0.3438,10.5383]

model_fit = Plotting_Functions.Fit_Empirical(fs, pts, prefix, "asym_mig", Models_2D.asym_mig, emp_params, fs_folded=True)

Plotting_Functions.Plot_2D(fs, model_fit, prefix, "asym_mig",vmin_val = vmin_val)



```


