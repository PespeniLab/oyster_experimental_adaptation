
cd ~/oyster/data

# -A will include any dups.
samtools mpileup -E -Q 20 --skip-indels -d 2000 -f /data/oyster/transcriptome/oyster_assembly.fa -o oyster.mpileup LL_33.merged.bam OP_33.merged.bam TB_33.merged.bam

