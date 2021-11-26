#!/bin/bash 

# Author: Adam Sorbie
# Date: 25/11/21
# Version: 0.5.3 

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 
while getopts p:o: flag
 do
   case "${flag}" in
     p) path=${OPTARG};;
     o) output=${OPTARG};; 
     *) echo "usage: $0 [-p] [-o]" >&2
        exit 1 ;;
   esac
done

# check qiime2 environment activated and exit if not 
env=$CONDA_DEFAULT_ENV
if [[ $env != *"qiime2"* ]]; then
  echo "Please activate qiime2 environment"
  exit
fi

mkdir -p $output

# create manifest file if required 
ls ${path}/*.gz | cut -f9 -d"/" > q2_manifest_sn
ls ${path}/*.gz > q2_manifest_fp 
paste -d, q2_manifest_sn q2_manifest_fp > q2_manifest_tmp.csv
rm -rf q2_manifest_fp q2_manifest_sn 

python create_manifest.py q2_manifest_tmp.csv ${output}/q2_manifest.csv 

rm -rf q2_manifest_tmp.csv 

cd $output || exit 
# import fastq files 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $path \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux-paired-end.qzv


