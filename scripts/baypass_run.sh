
# running baypass
# for these, first need to run the core model. After this, run the aux model to identify shifts in allele freq associated with the covariate.
#lets subset to only the control pops:

cat ~/oyster/analysis/genobaypass| cut -f 1,2,5,6,9,10 -d " " > ~/oyster/analysis/baypass.33.geno.in


#haploid pool size of each population
echo "40 40 40 40 30 32" > ~/oyster/analysis/baypass.psize.in
echo "40 40 30" > ~/oyster/analysis/baypass.33.psize.in

# and make a covariate file

echo "0 1 0 1 0 1" > ~/oyster/analysis/baypass.covariate.in
echo "0 1 0 1 0 0" > ~/oyster/analysis/baypass.covariate_pcadapt.in

cd ~/oyster/analysis/baypass

# run core model with only 33 pops, to get heirarchical clust for fig 1. Show pop structure

# 33 only, core
~/bin/baypass_2.1/sources/g_baypass -npop 3 -gfile ~/oyster/analysis/baypass.33.geno.in \
        -poolsizefile ~/oyster/analysis/baypass.33.psize.in  \
        -outprefix 33_core \
        -burnin 5000 \
        -pilotlength 1500 \
        -nval 2000 \
        -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_33.stdout_$(date +"%F_%R").txt

## next run all pops to identify shifts in allele freq

# run all pops, core model:

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/genobaypass \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -outprefix all_core \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all.stdout_$(date +"%F_%R").txt

# aux model to identify shifts in response to salinity in all the populations

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/genobaypass \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -efile ~/oyster/analysis/baypass.covariate.in \
    -outprefix all_aux \
    -omegafile all_core_mat_omega.out \
    -scalecov \
    -auxmodel \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all_aux.stdout_$(date +"%F_%R").txt

# shifts in allele frequency in only LL and OP, not TB.

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/genobaypass \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -efile ~/oyster/analysis/baypass.covariate_pcadapt.in \
    -outprefix all_aux_pcadapt \
    -omegafile all_core_mat_omega.out \
    -scalecov \
    -auxmodel \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_allpcadapt_aux.stdout_$(date +"%F_%R").txt
