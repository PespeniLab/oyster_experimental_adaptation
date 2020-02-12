# Rapid adaptation to experimental salinity

This repository hold scripts to run all analyses for experimental responses to salinity selection in the oyster Ostrea lurida.

To fill in:

Date of experiment:
Date of sequencing:

## Data availability

All raw sequence data can be accessed on NCBI, bioproject:

GO annotations and reference transcriptome were pulled from:

    Maynard, A., Bible, J.M., Pespeni, M.H., Sanford, E. and Evans, T.G., 2018. Transcriptomic responses to extreme low salinity among locally adapted populations of Olympia oyster (Ostrea lurida). Molecular ecology, 27(21), pp.4225-4240.

The transcriptome can be found at NCBI GEO GSE98355

Other relevant files are supplement with the manuscript.

## Scripts

Below are all scripts necessary to run the full analysis for the manuscript. A short description accompanies each.

### data processing
- align to reference transcriptome:`align_star.sh`
- downsample bams to get pooled counts: `downsamp.sh`
- convert to pileup: `to_pileup.sh`
- call variants: `varscan.sh`
- Filter variants: `filter_variants.R`
- convert to sync for downstream analyses: `to_sync.R`

### Analyses:

- calc Fst: `fstats.R`
- run PCA: `pca.R`
- identify local adaptation loci with PCAdapt: `PCAdapt.R`
- Run baypass to identify salinity responsive loci; no downstream analyses: `baypass_run.sh`
- pull out significant baypass snp sets: `baypass.R`
- make GO annotations: `go_annotations.py`
- convert to topgo format: `to_topgo.R`
- Run topgo: `TopGO.R`
- dadi analysis: in directory: dadi, see `run_dadi.md` for details.

### Figures

- Fig1: `fstats.R` and `heir_clust.R` and `map.R`
- Fig2: `pca.R`
- Fig3: `snp_plot.R`
- Fig4: `venn.R`
- Fig5: in `run_dadi.md`
