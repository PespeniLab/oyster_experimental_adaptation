
#count up reads:

cd /data/oyster/aligned

for i in $(ls *.bam | cut -f -1 -d "." | uniq )

do {
ONE=$(samtools view -F 4 ${i}.bam | wc -l) #mapped reads
TWO=$(samtools view -F 4 -q 20 ${i}.bam | wc -l) #mapped reads with quality of 20
THREE=$(samtools view -f 1024 ${i}.bam | wc -l) #pcr duplicates
FOUR=$(samtools view  ${i}.bam | wc -l) #total
echo ${i},$ONE,$TWO,$THREE,$FOUR
 } >> ~/oyster/count.aligned_qual.txt

done


#need to downsample to the lowest number of mapped reads in a group.

cd /data/oyster/aligned


for i in LL_33 LL_5 OP_33 OP_5 TB_33 TB_5; do

    cat ~/oyster/count.aligned_qual.txt | grep $i > tmp.ct

    low=$(sort -t, -nk2 tmp.ct | cut -f 3 -d ',' | head -n 1)
    low1=$(sort -t, -nk2 tmp.ct)
    echo $low1
    echo $i
    echo $low
    for samp in "_1" "_2" "_3" "_4"; do
        id=$(echo ${i}${samp})
        spct=$(cat tmp.ct | grep ${id} | cut -f 3 -d ',')
        prop1=$(echo "scale=5 ; $low / $spct"|bc)
        echo $samp
        echo $prop1
        if (( prop1 < 1)); then
            echo "SUBSAMPLE ${i}${samp}"
            samtools view -hb -F 4 -q 20 ${i}${samp}.bam | \
            /data/programs/sambamba_v0.6.0 view -f bam -t 4 --subsampling-seed=3 -s ${prop1} /dev/stdin | samtools sort > /data/oyster/subsample/${i}${samp}.subsample.bam
        else
            echo "LOWEST ${i}${samp}"
            ~/bin/samtools-1.6/samtools view -hb -F 4 -q 20 ${i}${samp}.bam | ~/bin/samtools-1.6/samtools sort > /data/oyster/subsample/${i}${samp}.subsample.bam
        fi
    echo "DONE with 1 LOOP"
    done

done

# one sample in TB_33 has very low reads. remove this one.

for i in TB_33 ; do

    cat ~/oyster/count.aligned_qual.txt | grep $i | grep -v 'TB_33_4' > tmp.ct

    low=$(sort -t, -nk2 tmp.ct | cut -f 3 -d ',' | head -n 1)
    low1=$(sort -t, -nk2 tmp.ct)
    echo $low1
    echo $i
    echo $low
    for samp in "_1" "_2" "_3"; do
        id=$(echo ${i}${samp})
        spct=$(cat tmp.ct | grep ${id} | cut -f 3 -d ',')
        prop1=$(echo "scale=5 ; $low / $spct"|bc)
        echo $samp
        echo $prop1
        if (( prop1 < 1)); then
            echo "SUBSAMPLE ${i}${samp}"
            samtools view -hb -F 4 -q 20 ${i}${samp}.bam | \
            /data/programs/sambamba_v0.6.0 view -f bam -t 4 --subsampling-seed=3 -s ${prop1} /dev/stdin | samtools sort > /data/oyster/subsample/${i}${samp}.subsample.bam
        else
            echo "LOWEST ${i}${samp}"
            ~/bin/samtools-1.6/samtools view -hb -F 4 -q 20 ${i}${samp}.bam | ~/bin/samtools-1.6/samtools sort > /data/oyster/subsample/${i}${samp}.subsample.bam
        fi
    echo "DONE with 1 LOOP"
    done

done

rm ~/oyster/count.subsample.txt

cd /data/oyster/subsample

for i in $(ls *.bam | cut -f -1 -d "." | uniq )

do {
ONE=$(samtools view  ${i}.subsample.bam | wc -l) #total
echo ${i},$ONE
 } >> ~/oyster/count.subsample.txt

done


# merge by pop

cd /data/oyster/subsample

for i in $(ls *.bam | cut -f 1-2 -d "_" | uniq ); do {

    ls -d $PWD/*.bam | grep ${i} > bam_list.txt

    ~/bin/bamtools/bamtools merge -list bam_list.txt | ~/bin/samtools-1.6/samtools sort > ~/oyster/data/${i}.merged.bam

    echo "${i} done"

}

done

