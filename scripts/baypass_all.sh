
# run all pops:

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/genobaypass \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -outprefix all_core \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all.stdout_$(date +"%F_%R").txt

# all, aux to identify shifts in response to salinity

~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/genobaypass \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -efile ~/oyster/analysis/baypass.covariate.in \
    -outprefix all_aux \
    -omegafile all_core_mat_omega.out \
    -scalecov \
    -auxmodel \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_all_aux.stdout_$(date +"%F_%R").txt
