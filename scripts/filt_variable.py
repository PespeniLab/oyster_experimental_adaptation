
import re

a_base=['A','a']
c_base=['C','c']
t_base=['T','t']
g_base=['G','g']

l = [4,7,10,13,16,19]

ix = 0
good_snp = 0
with open('/users/r/b/rbrennan/oyster/data/oyster.filter.mpileup', "r") as master_file, open('oyster.final.mpileup', 'w') as myfile:
    for line in master_file:
        a_count = 0
        c_count = 0
        t_count = 0
        g_count = 0
        dep = 0
        for i in l:
            tmp = line.split("\t")[i]
            a_count = a_count + len([b for b in tmp if b in a_base])
            c_count = c_count + len([b for b in tmp if b in c_base])
            t_count = t_count + len([b for b in tmp if b in t_base])
            g_count = g_count + len([b for b in tmp if b in g_base])
            # calc maf > 0.01
            dep = dep + int(line.split("\t")[i-1])
        bi_allele = sum([a_count > 0, c_count > 0, t_count > 0, g_count > 0])
        # removing indels with the regular expression below
        if (a_count > 1 or c_count > 1 or t_count > 1 or g_count > 1) and bool(re.search(r"[+-](\d+)", tmp)) is False and (float((a_count + c_count + t_count + g_count)/float(dep)) >= 0.01) and bi_allele == 1:
            myfile.write(line)
            good_snp = good_snp + 1
            if good_snp % 10000 == 0: print "number of good snps:", good_snp
        if ix % 50000 == 0: print "iteration", ix
        ix = ix + 1

