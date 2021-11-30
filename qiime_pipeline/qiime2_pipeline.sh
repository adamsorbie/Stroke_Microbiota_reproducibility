#!/bin/bash 

# Author: Adam Sorbie
# Date: 25/11/21
# Version: 0.5.1 


 while getopts p:g:G:m:f:r:n:N flag
 do
   case "${flag}" in
     p) path=${OPTARG};; 
     g) f_primer=${OPTARG};;
     G) r_primer=${OPTARG};;
     m) min_len=${OPTARG};;
     f) trunc_f=${OPTARG};;
     r) trunc_r=${OPTARG};;
     n) maxEE_f=${OPTARG};;
     N) maxEE_r=${OPTARG};;
    
     *) echo "usage: $0 [-p] [-g] [-G] [-m] [-f] [-r] [-n] [-N]" >&2
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

cd $path || exit

# run cutadapt 
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza
  --p-front-f $f_primer
  --p-front-r $r_primer
  --p-minimum-length $min_len
  --o-trimmed-sequences demux-paired-end-trimmed.qza 

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end-trimmed.qza \
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

qiime tools export --input-path rep-seqs.qza --output-path rep_seqs.fasta 

# Filter 
python filter_low_abundant.py
bash filter_fasta.sh -i rep_seqs.fasta -p keep_asvs.tsv -o rep_seqs_filt.fasta
# Re-convert to biom 
biom convert -i feature-table_filt.tsv -o feature-table_filt.biom --table-type="OTU table" --to-hdf5

## Re-import filtered table and fasta file
qiime tools import \
  --input-path feature-table_filt.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \ 
  --output-path feature-table_filt.qza

qiime tools import \
  --input-path rep_seqs_filt.fasta \
  --output-path rep_seqs_filt.qza \
  --type 'FeatureData[Sequence]'


## Taxonomy assignment 
wget \
  -O "silva-138-99-nb-classifier.qza" \
  "https://data.qiime2.org/2021.11/common/silva-138-99-nb-classifier.qza"


qiime feature-classifier classify-sklearn \
  --i-reads rep_seqs_filt.qza \
  --i-classifier silva-138-99-nb-classifier.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime feature-table tabulate-seqs \
  --i-data rep_seqs_filt.qza \
  --o-visualization rep_seqs_filt.qzv

## Tree generation 
qiime phylogeny align-to-tree-mafft-fasttree \ 
  --i-sequenes rep_seqs_filt.qza
  --o-alignment rep_seqs_filt_aln.qza
  --o-masked-alignment rep_seqs_filt_aln_mask.qza
  --o-tree unrooted_tree.qza
  --o-rooted-tree rooted_tree.qza
