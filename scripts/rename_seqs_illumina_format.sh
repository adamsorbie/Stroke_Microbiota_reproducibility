#!/bin/bash
# Author: Adam Sorbie 
# Date: 15/04/21
# Version 0.6.0

# $1 = delim
# $2 = last field
# $3 = fwd pattern e.g. R1.fastq.gz
# $4 = rev pattern e.g. R1.fastq.gz
# rename forward reads
for i in *${3};
do
  name=$(echo "$i" | cut -f1-"$2" -d "$1")
	mv -n "$i" "${name}_S1_L001_R1_001.fastq.gz"
done
# rename reverse reads
for i in *${4};
do
  name=$(echo "$i" | cut -f1-"$2" -d "$1")
	mv -n "$i" "${name}_S1_L001_R2_001.fastq.gz"
done
