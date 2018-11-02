#convert to pileup:
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/LL_33.merged.bam > ~/oyster/data/LL_33.pileup
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/LL_5.merged.bam > ~/oyster/data/LL_5.pileup
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/OP_33.merged.bam > ~/oyster/data/OP_33.pileup
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/OP_5.merged.bam > ~/oyster/data/OP_5.pileup
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/TB_33.merged.bam > ~/oyster/data/TB_33.pileup
#~/bin/samtools-1.6/samtools mpileup ~/oyster/data/TB_5.merged.bam > ~/oyster/data/TB_5.pileup

# calc pi

 for i in $(ls ~/oyster/data/ | grep '\.pileup' | cut -f 1 -d '.'); do

    echo "Status: starting $i"

    perl ~/bin/popoolation_1.2.2/Variance-at-position.pl --pool-size 20 --min-qual 20 \
    --min-coverage 10 --min-count 1 --max-coverage 564 --min-covered-fraction 0.3 \
    --pileup ~/oyster/data/${i}.pileup \
    --gtf /data/oyster/transcriptome/oyster_assembly.gtf \
    --output ~/oyster/analysis/${i}.pi --measure pi

    echo "Status: $i popoolation done"

done

