~/bin/baypass_2.1/sources/g_baypass -npop 6 -gfile ~/oyster/analysis/baypass/G.all_pods \
    -poolsizefile ~/oyster/analysis/baypass.psize.in  \
    -outprefix all_core_pods \
    -nthreads 4 2>&1 | tee ~/oyster/log_out/baypass_allpod.stdout_$(date +"%F_%R").txt
