#!/bin/bash 

# Author: Adam Sorbie
# Date: 17/09/21
# Version: 0.0.1 

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

## Re-import filtered table 



## Taxonomy assignment 


