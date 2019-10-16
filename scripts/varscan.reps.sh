cd ~/oyster/data

#awk '{ if ($4 >= 20 && $4 <= 558 && $7 >= 20 && $7 <= 558 && $10 >= 20 && $10 <= 558 && $13 >= 20 && $13 <= 558 && $16 >= 20 && $16 <= 558 && $19 >= 20 && $19 <= 558 ) print $0 }'  oyster.mpileup > oyster.filter.mpileup

#cat oyster.filter.mpileup | wc -l

java -jar ~/bin/VarScan.v2.4.3.jar mpileup2snp ~/oyster/data/oyster.reps.mpileup --min-coverage 10 --min-reads 2 --min-avg-qual 20 --min-var-freq 0.025 --variants --p-value 0.1 > ~/oyster/analysis/varscan_reps.txt

