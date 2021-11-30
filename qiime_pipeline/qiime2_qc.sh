#!/bin/bash 

# Author: Adam Sorbie
# Date: 25/11/21
# Version: 0.5.3 

# This script assumes you have installed the latest version of QIIME2 and activated your conda environment 
phred=33
while getopts p:o: flag
 do
   case "${flag}" in
     p) path=${OPTARG};;
     o) output=${OPTARG};; 
     e) encoding=${OPTARG};;
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
ls ${path}/*R1_001.fastq.gz | xargs -n 1 basename | cut -f1 -d"_" > q2_manifest_sn
ls ${path}/*R1_001.fastq.gz > q2_manifest_ffp
ls ${path}/*R2_001.fastq.gz > q2_manifest_rfp
paste q2_manifest_sn q2_manifest_ffp q2_manifest_rfp > q2_manifest_tmp.tsv
rm -rf q2_manifest_sn q2_manifest_ffp q2_manifest_rfp
python create_manifest.py q2_manifest_tmp.tsv ${output}/q2_manifest
rm -rf q2_manifest_tmp.tsv 

cd $output || exit 
if [[ "$phred" == 64 ]]; 
then
  # import fastq files 
  qiime tools import \
    --type SampleData[PairedEndSequencesWithQuality] \
    --input-path ${output}/q2_manifest \
    --output-path ${output}/demux-paired-end.qza \
    --input-format PairedEndFastqManifestPhred64
elif [[ "$phred" == 33 ]];
then
  cd $output || exit 
  # import fastq files 
  qiime tools import \
    --type SampleData[PairedEndSequencesWithQuality] \
    --input-path ${output}/q2_manifest \
    --output-path ${output}/demux-paired-end.qza \
    --input-format PairedEndFastqManifestPhred33V2
else 
  echo "Not a supported encoding"
fi 

qiime demux summarize \
  --i-data ${output}/demux-paired-end.qza \
  --o-visualization ${output}/demux-paired-end.qzv


