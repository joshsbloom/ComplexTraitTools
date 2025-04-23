#!/bin/bash
#recommend debug running this interactively
#root_dir should be somewhere on /u/project/kruglyak
root_dir='/data2/xQTL/FluconazoleV1_041825/'
bam_dir="${root_dir}/bam_files/"
count_dir="${root_dir}/count_files/"
fastq_dir="${root_dir}/fastq/"

#remove the asterisks from your column names, then extract column 4 and remove first line  
#rename BY to S288C
keyfile='/home/jbloom/Dropbox/code/ComplexTraitTools/xQTLStats/examples/FluconazoleV1_041825/flatkey.csv'

#on hoffman2
vcffile='/u/project/kruglyak/jsbloom/PUBLIC_SHARED/yeast/1011/mega_filtered.vcf.gz'
#vcffile='/media/hoffman2/PUBLIC_SHARED/yeast/1011/mega_filtered.vcf.gz'

#indexed sacCer3 genome
genome='/u/project/kruglyak/jsbloom/PUBLIC_SHARED/yeast/1011/sacCer3.fasta'
#genome='/media/hoffman2/PUBLIC_SHARED/yeast/1011/sacCer3.fasta'

nthreads=8

#extract sample name from key file
sname=($(tail -n +2 ${keyfile} | cut -d ',' -f 4 ))

mkdir ${bam_dir}
mkdir ${count_dir}
mkdir ${fastq_dir}
var=1
for t in ${sname[@]}
do 
    echo $t
    echo $var
    # note, should cat together files from multiple lanes if it was run on an instrument with multiple lanes
    # note, use of star expansion for basespace downloader adding additional crap to folder name
    bwa mem -t $nthreads -R '@RG\tID:${t}\tSM:{$t}' ${genome} ${fastq_dir}${t}*/${t}_S${var}_L001_R1_001.fastq.gz ${fastq_dir}${t}*/${t}_S${var}_L001_R2_001.fastq.gz | \
    sambamba view --nthreads=${nthreads} --show-progress --sam-input --format=bam --with-header /dev/stdin | \
    sambamba sort --nthreads=${nthreads} --show-progress --tmpdir=/tmp --out=${bam_dir}${t}.bam /dev/stdin
    gatk ASEReadCounter --min-mapping-quality 20 -I ${bam_dir}${t}.bam -V ${vcffile} -R $genome -O ${count_dir}${t}.txt
    var=$((var+1))
done
