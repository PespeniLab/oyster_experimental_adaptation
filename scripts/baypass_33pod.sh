~/bin/baypass_2.1/sources/g_baypass -npop 3 -gfile ~/oyster/analysis/baypass/G.33_pods \
    -poolsizefile ~/oyster/analysis/baypass.33.psize.in  \
    -outprefix 33_core_pods \
    -nthreads 8 2>&1 | tee ~/oyster/log_out/baypass_33pod.stdout_$(date +"%F_%R").txt
