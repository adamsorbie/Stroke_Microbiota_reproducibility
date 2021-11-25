#!/bin/bash 

# Author: Adam Sorbie
# Date: 25/11/21
# Version: 0.5.1 

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 
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

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi

mkdir qiime2_pipeline && cd "$_"

# import fastq files 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $path \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime demux summarize \
    --i-data demux-paired-end.qza \
    --o-visualization ./demux-paired-end.qzv


