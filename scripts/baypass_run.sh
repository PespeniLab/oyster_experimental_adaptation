
#lets subset to only the control pops:

cat ~/oyster/analysis/baypass.geno.in | cut -f 1,2,5,6,9,10 -d " " > ~/oyster/analysis/baypass.33.geno.in

#haploid pool size of each population
echo "40 40 40 40 30 32" > ~/oyster/analysis/baypass.psize.in
echo "40 40 30" > ~/oyster/analysis/baypass.33.psize.in

# and make a covariate file. specified by -efile
echo "1 -1 -1" > ~/oyster/analysis/baypass.33.covariate.in

echo "-1 1 -1 1 -1 1" > ~/oyster/analysis/baypass.covariate.in

cd ~/oyster/analysis/baypass

# 33 only, core
~/bin/baypass_2.1/sources/g_baypass -npop 3 -gfile ~/oyster/analysis/baypass.33.geno.in \
    -poolsizefile ~/oyster/analysis/baypass.33.psize.in  \
    -outprefix 33_core \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_33.stdout_$(date +"%F_%R").txt

# 33 only, aux to identify salinity associations

~/bin/baypass_2.1/sources/g_baypass -npop 3 -gfile ~/oyster/analysis/baypass.33.geno.in \
    -poolsizefile ~/oyster/analysis/baypass.33.psize.in  \
    -efile ~/oyster/analysis/baypass.33.covariate.in \
    -outprefix 33_aux \
    -omegafile 33_core_mat_omega.out \
    -auxmodel \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_33_aux.stdout_$(date +"%F_%R").txt

# run all pops:

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/baypass.geno.in \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -outprefix all_core \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all.stdout_$(date +"%F_%R").txt

# all, aux to identify shifts in response to salinity

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/baypass.geno.in \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -efile ~/oyster/analysis/baypass.covariate.in \
    -outprefix all_aux \
    -omegafile all_core_mat_omega.out \
    -auxmodel \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all_aux.stdout_$(date +"%F_%R").txt
