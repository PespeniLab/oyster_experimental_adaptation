#convert to pileup:
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/LL_33.merged.bam > ~/oyster/data/LL_33.pileup
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/LL_5.merged.bam > ~/oyster/data/LL_5.pileup
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/OP_33.merged.bam > ~/oyster/data/OP_33.pileup
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/OP_5.merged.bam > ~/oyster/data/OP_5.pileup
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/TB_33.merged.bam > ~/oyster/data/TB_33.pileup
~/bin/samtools-1.6/samtools mpileup ~/oyster/data/TB_5.merged.bam > ~/oyster/data/TB_5.pileup

# calc pi
perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/LL_33.pileup --output ~/oyster/analysis/LL_33.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/LL_33.pi | grep -v 'na' > ~/oyster/analysis/LL_33.pi1
mv ~/oyster/analysis/LL_33.pi1 ~/oyster/analysis/LL_33.pi

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/LL_5.pileup --output ~/oyster/analysis/LL_5.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/LL_5.pi | grep -v 'na' > ~/oyster/analysis/LL_5.pi1
mv ~/oyster/analysis/LL_5.pi1 ~/oyster/analysis/LL_5.pi

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/OP_33.pileup --output ~/oyster/analysis/OP_33.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/OP_33.pi | grep -v 'na' > ~/oyster/analysis/OP_33.pi1
mv ~/oyster/analysis/OP_33.pi1 ~/oyster/analysis/OP_33.pi

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/OP_5.pileup --output ~/oyster/analysis/OP_5.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/OP_5.pi | grep -v 'na' > ~/oyster/analysis/OP_5.pi1
mv ~/oyster/analysis/OP_5.pi1 ~/oyster/analysis/OP_5.pi

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/TB_33.pileup --output ~/oyster/analysis/TB_33.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/TB_33.pi | grep -v 'na' > ~/oyster/analysis/TB_33.pi1
mv ~/oyster/analysis/TB_33.pi1 ~/oyster/analysis/TB_33.pi

perl ~/bin/popoolation_1.2.2/Variance-sliding.pl --input ~/oyster/data/TB_5.pileup --output ~/oyster/analysis/TB_5.pi --measure pi --window-size 500 --step-size 500 --min-count 2 --min-coverage 4 --max-coverage 2000 --min-qual 10  --min-covered-fraction 0.5 --pool-size 40

cat ~/oyster/analysis/TB_5.pi | grep -v 'na' > ~/oyster/analysis/TB_5.pi1
mv ~/oyster/analysis/TB_5.pi1 ~/oyster/analysis/TB_5.pi


