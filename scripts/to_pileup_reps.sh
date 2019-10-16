
#remove previously marked duplicates:
cd ~/oyster/data

for i in `ls /data/oyster/aligned/STAR/*.bam | cut -f1 -d'.' | sort -u` ; do

    echo "starting ${i}"

    samtools view -h /data/oyster/aligned/${i}.bam | awk 'BEGIN{OFS="\t"}{if(NF>5) {if(and($2,1024)) {$2-=1024}} print $0}' | samtools view -Sbo ${i}.unmarked.bam -

    echo "${i} done"

done


ls ~/oyster/data | grep 'bam'| sort -u | grep -E '_1.bam|_2.bam|_3.bam|_4.bam' > bam.file.in

samtools mpileup -E -Q 20 --skip-indels -d 8000 -f /data/oyster/transcriptome/oyster_assembly.fa -o oyster.reps.mpileup --bam-list bam.file.in
