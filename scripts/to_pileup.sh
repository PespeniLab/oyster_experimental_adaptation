#remove previously marked duplicates:
cd ~/oyster/data

for i in LL_33 LL_5 OP_33 OP_5 TB_33 TB_5; do

    echo "starting ${i}"

    samtools view -h ${i}.merged.bam | awk 'BEGIN{OFS="\t"}{if(NF>5) {if(and($2,1024)) {$2-=1024}} print $0}' | samtools view -Sbo ${i}.unmarked.bam -

    echo "${i} done"

done

# -A will include any dups.
samtools mpileup -E -Q 20 --skip-indels -d 2000 -f /data/oyster/transcriptome/oyster_assembly.fa -o oyster.mpileup LL_33.unmarked.bam LL_5.unmarked.bam OP_33.unmarked.bam  OP_5.unmarked.bam  TB_33.unmarked.bam  TB_5.unmarked.bam

# filter pileup
# include sites with > 10x, less than 150x
# also remove ^ and $ start and ends of reads. remove them. ^ has qual score follwoing too

awk '{ if ($4 >= 20 && $4 <= 350 && $7 >= 20 && $7 <= 350 && $10 >= 20 && $10 <= 350 && $13 >= 20 && $13 <= 350 && $16 >= 20 && $16 <= 350 && $19 >= 20 && $19 <= 350 ) print $0 }'  oyster.mpileup | sed 's/\^.//g' | sed 's/\$//g'> oyster.filter.mpileup
