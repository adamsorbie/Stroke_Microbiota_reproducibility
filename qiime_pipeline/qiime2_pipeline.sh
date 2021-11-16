#!/bin/bash 

# Author: Adam Sorbie
# Date: 17/09/21
# Version: 0.0.1 

 while getopts g:G:m:M:p: flag
 do
   case "${flag}" in
     g) f_primer=${OPTARG};;
     G) r_primer=${OPTARG};;
     m) min_len=${OPTARG};;
     f) trunc_f=${OPTARG};;
     r) trunc_r=${OPTARG};;
     n) maxEE_f=${OPTARG};;
     N) maxEE_r=${OPTARG};;
     *) echo "usage: $0 [-g] [-G] [-m] [-f] [-r] [-n] [-N]" >&2
        exit 1 ;;
   esac
done

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi

cd qiime2_pipeline || exit

# run cutadapt 
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza
  --p-front-f $f_primer
  --p-front-r $r_primer
  --p-minimum-length $min_len


qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trunc-len-f $trunc_f \
  --p-trunc-len-r $trunc_r \
  --p-max-ee-f $maxEE_f
  --p-max-ee-r $maxEE_r
  --p-trunc-q 3
  --p-chimera-method 'consensus'
  --p-n-reads-learn 100000000
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza


## Write out abundance table and filter by 0.25% 

qiime tools export table.qza --output-dir .
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv

qiime tools export rep-seqs.qza --output-dir .

# Filter 
python filter_low_abundant.py
bash filter_fasta.sh 

## Re-import filtered table 
qiime tools import 


## Taxonomy assignment 


