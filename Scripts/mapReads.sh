#!/bin/bash
# this script bwa mapped PE fastq file to hg38 with Zhang...Ren paper parameters
# wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# gunzip hg38.fa.gz
# bwa index hg38.fa
ref=$1
fastq=$(echo $2 | awk -F "." '{print $1 ".demultiplexed"}')
donor=$(echo $2 | awk -F "_" '{print $1}' | awk -F "-" '{print $5}')
organ=$(echo $2 | awk -F "_" '{print $2}' | awk -F "." '{print $1}')
out=$3
dir=$4
cd $dir
gunzip $fastq.R1.fastq.gz
gunzip $fastq.R2.fastq.gz
# report read length
echo "number of reads and read length in fastq files"
R1=($(awk '{if(NR%4==2) print length($1)}' $fastq.R1.fastq | sort -n | uniq -c))
R2=($(awk '{if(NR%4==2) print length($1)}' $fastq.R2.fastq | sort -n | uniq -c))
echo ${R1[0]}
echo ${R2[0]}
# start genome mapping
echo "start genome mapping"
bwa mem -M -k 32 -t 8 $ref $fastq.R1.fastq $fastq.R2.fastq > $donor.sam
# report SAMstats output to file
SAMstats --sorted_sam_file $donor.sam --outf $donor.SAMstats.txt
echo "finish mapping"
# start filtering
echo "start reads filtering"
samtools view -bS -f 2 -F 0x70c -q 30 $donor.sam | samtools sort - $donor.filt.sorted
echo "# reads uniquely mapped to genome"
unique=$(samtools view -@ 6 -q 1 -F 4 -F 256 -h $donor.filt.sorted.bam | grep -v -E -e '\bXA:Z:' -e '\bSA:Z:' | wc -l)
echo "# reads multi mapped to genome"
multi=$(samtools view -@ 6 -b -f 256 -F 4 -h $donor.filt.sorted.bam | wc -l)
total=$((${R1[0]} + ${R2[0]}))
length=$(awk -v r1=${R1[0]} -v r2=${R2[0]} '{print (r1+r2)/2;}')
uniquepct= $(awk -v unique=$unique -v total=$total '{print (unique/total);}')
multipct=$(awk -v multi=$multi -v total=$total '{print (multi/total); }')
# report filtered mapping stats
samtools flagstat $donor.filt.sorted.bam > $donor.filter.flagstat.txt
samtools index $donor.filt.sorted.bam
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$organ-$donor" "$total" "$length" "$unique" "$uniquepct" "$multi" "$multipct" "$fastq" >> $out
