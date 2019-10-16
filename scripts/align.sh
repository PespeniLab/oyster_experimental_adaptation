#!/bin/bash -l

my_bwa=~/bin/bwa/bwa
my_samtools=~/bin/samtools-1.6/samtools
bwagenind=/data/oyster/transcriptome/oyster_assembly.fa
my_samblstr=~/bin/samblaster/samblaster


# index reference
#$my_bwa index oyster_assembly.fa

cd /data/oyster/fastq/cleanreadsNHC/


for sample in `ls /data/oyster/fastq/cleanreadsNHC/*.fq.gz |cut -f 6- -d "/" | cut -f1 -d'.' | sort -u`
do
    samp=$(echo $sample | cut -f 1-2 -d "_")
    lib=$(echo $sample | cut -f 1-3 -d "_")
    echo $sample
    echo $samp
    echo $lib
    rg=$(echo \@RG\\tID:$lib\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$lib)
    #echo $rg

    $my_bwa mem -t 6 -R $rg $bwagenind ${sample}.nhc.cl.fq.gz | \
    $my_samtools view -h -u - | \
    $my_samtools sort - -O bam -o /data/oyster/aligned/$lib.bam
done

