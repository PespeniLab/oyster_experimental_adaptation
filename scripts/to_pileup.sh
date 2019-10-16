#remove previously marked duplicates:
cd ~/oyster/data

for i in LL_33 LL_5 OP_33 OP_5 TB_33 TB_5; do

    echo "starting ${i}"

    samtools view -h ${i}.merged.bam | awk 'BEGIN{OFS="\t"}{if(NF>5) {if(and($2,1024)) {$2-=1024}} print $0}' | samtools view -Sbo ${i}.unmarked.bam -

    echo "${i} done"

done

# -A will include any dups.
samtools mpileup -E -Q 20 --skip-indels -d 2000 -f /data/oyster/transcriptome/oyster_assembly.fa -o oyster.mpileup LL_33.unmarked.bam LL_5.unmarked.bam OP_33.unmarked.bam  OP_5.unmarked.bam  TB_33.unmarked.bam  TB_5.unmarked.bam

