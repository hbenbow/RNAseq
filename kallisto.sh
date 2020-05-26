#!/bin/bash

for i in  $(ls ../../clean_fastq/* | cut -c19-22 | uniq)
do
	kallisto quant -i index -o  ${i}   ../../clean_fastq/${i}_1.fq.gz ../../clean_fastq/${i}_2.fq.gz
done
