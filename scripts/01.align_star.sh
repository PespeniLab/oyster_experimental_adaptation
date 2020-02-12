!/bin/bash -l

# do only once:
#~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR  --runMode genomeGenerate --genomeChrBinNbits 16 limitGenomeGenerateRAM=36052865749 --runThreadN 10 sjdbGTFfile /data/oyster/transcriptome/oyster_assembly.gtf --genomeDir /data/oyster/transcriptome/star_genome/ --genomeFastaFiles /data/oyster/transcriptome/oyster_assembly.fa

cd /data/oyster/fastq/cleanreadsNHC/


for sample in `ls /data/oyster/fastq/cleanreadsNHC/*.fq.gz |cut -f 6- -d "/" | cut -f1 -d'.' | sort -u`
do

samp=$(echo $sample | cut -f 1-2 -d "_")
lib=$(echo $sample | cut -f 1-3 -d "_")

echo ${lib}

~/bin/STAR-2.7.0e/bin/Linux_x86_64/STAR --genomeDir /data/oyster/transcriptome/star_genome --runThreadN 10 --readFilesCommand zcat --readFilesIn ${sample}.nhc.cl.fq.gz --outFileNamePrefix /data/oyster/aligned/STAR/${lib}_ --outSAMattrRGline ID:${lib} SM:${lib} --outSAMtype BAM Unsorted

done
