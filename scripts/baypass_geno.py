import numpy as np

snp_path = '/users/r/b/rbrennan/oyster/analysis/oyster.snp_rc'

#make empty array to save output
snp_out = np.empty(shape=(76633,12), dtype = object)# row, col

with open(snp_path) as master_file:
    header_line = next(master_file) # skip header row
    for idx, line in enumerate(master_file):
        tmp_snp = line.split("\t")[9:]

        tmp_idx = 0
        for af in tmp_snp[0:6]:
            snp_out[idx, tmp_idx] = af.split("/")[0]
            snp_out[idx,tmp_idx+1] = int(af.split("/")[1]) - int(af.split("/")[0])
            tmp_idx = tmp_idx+2

np.savetxt('/users/r/b/rbrennan/oyster/analysis/baypass.geno.in', snp_out, fmt='%s',delimiter=' ')
