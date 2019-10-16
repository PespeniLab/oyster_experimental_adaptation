#Rapid adaptation to experimental salinity

This repository hold scripts to run all analyses for experimental responses to salinity selection in the oyster Ostrea lurida.

To fill in:

Date of experiment:
Date of sequencing:

## Data availability

All raw sequence data can be accessed on NCBI, bioproject:

Other relevant files are supplement with the manuscript. 

## Scripts

Below are all scripts necessary to run the full analysis for the manuscript. A short description accompanies each.

align:`align_star.sh`  
downsample bams: `downsamp.sh`  
convert to pileup: `to_pileup.sh`  
call variants: `varscan.sh`  
Filter variants: `filter_variants.R`  
convert to sync: `to_sync.R`  


calc Fst: `fstats.R`  
PCA: `pca.R`  
PCAdapt: `PCAdapt.R`  
Run baypass, no downstream analyses: `baypass_run.sh`  
pull out significant baypass snp sets: `baypass.R`  
make go annotations: `go_annotations.py`  
convert to topgo format: `to_topgo.R`  
Run topgo: `TopGO.R`  

dadi analysis: in directory: dadi, see `run_dadi.md` for details.  

Fig1: `fstats.R` and `heir_clust.R` and `map.R`  
Fig2: `pca.R`  
Fig3: `snp_plot.R`  
Fig4: `venn.R`  
Fig5: in `run_dadi.md`  

