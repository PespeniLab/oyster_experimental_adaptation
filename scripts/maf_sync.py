# filter sync file by maf

from __future__ import division

l = [3,4,5,6,7,8]

with open('/users/r/b/rbrennan/oyster/data/oyster_initial.sync', "r") as master_file, open('/users/r/b/rbrennan/oyster/data/oyster.sync', 'w') as myfile:
    for line in master_file:
        maf_count = 0
        del_count = 0
        for i in l:
            tmp = line.split("\t")[i].split("\n")[0].split(":")
            tmp_num = map(int, tmp)
            tmp_base = tmp_num[0:4]
            del_count = del_count + tmp_num[5]
            tmp_grt = [i for i in tmp_base if i > 0]
            if len(tmp_grt) > 1: #check if the pop is variable for this position
                tmp_af = min(tmp_grt)/(min(tmp_grt)+max(tmp_grt))
            else:
                tmp_af = 0
            if tmp_af >= 0.05: maf_count = maf_count + 1
        if (maf_count > 0 and del_count == 0):
            myfile.write(line)
