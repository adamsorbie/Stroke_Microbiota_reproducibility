Microbiota\_Analysis
================
Adam Sorbie
26/10/2021

# Analysis of Microbiota data

This document contains code to analyse the microbiota data generated in
Sorbie, Jimenez and Benakis 2021 and is intended to act as template for
other stroke researchers to conduct their own analyses.

## Set-up 

To be able to perform these analyses you will first need to clone this repository: 
(type ```git clone git@github.com:adamsorbie/Stroke_Microbiota_reproducibility.git```). 
You will also need a working installation of Qiime2 to follow this tutorial. Make sure you have followed the first section on the main page [here](https://github.com/adamsorbie/Stroke_Microbiota_reproducibility) (up to Data processing pipeline heading). 

### Installation 

If you have successfully installed conda, the installation of Qiime2 should be a breeze. 

We can run the following commands to instally everything we need and make it easier to activate our Qiime2 conda env. 
```
wget https://data.qiime2.org/distro/core/qiime2-2022.2-py38-linux-conda.yml
conda env create -n qiime2-2022.2 --file qiime2-2022.2-py38-linux-conda.yml
rm qiime2-2022.2-py38-linux-conda.yml

# add short alias to activate qiime environment quickly
echo "alias q2='conda activate qiime2-2022.2'" >> ~/.bashrc && source ~/.bashrc 
```
## Data Processing 


The Qiime2 processing pipeline mirrors that of the bash/R based pipeline, using the same tools under the Qiime2 wrapper. 

Before doing anything we need to make sure to activate our Qiime2 conda env. Thanks to the handy alias we defined above we can do this like so: 

```
q2
```

To perform sample QC etc, we first run the bash script ```qiime2_qc.sh```, providing the path of fastq files and an output folder for the qiime2 artifact and visualization. 

Note that if data is not in phred33 fastq format it must be converted before import. Qiime2 can do this internally, however it is extremely slow. If you are unsure whether your data are phred33 or phred64, you can ask your sequencing provider or alternatively, use vsearch to estimate the type. Vsearch is already installed within your qiime environment, thus you can just run the command ```vsearch --fastqchars``` on a single unzipped fastq file to check. e.g. 

```
vsearch --fastq_chars RD118_S1_L001_R1_001.fastq
``` 

If you need to convert to phred33 a bash script to perform this using the very fast bbmap package is provided. 

First install bbmap: 

```
conda install -c bioconda bbmap
```

Then run the script ```convert_phred.sh``` like so, supplying the file path to your fastq files 

```
bash convert_phred.sh -p ~/data_win/16S/iScience-review/data/fastq
```

Before proceeding to Qiime2, it's important to check that the following python scripts are in the working directory or in your path: 

```create_manifest.py``` and ```filter_low_abundant.py```. Without these the bash scripts below will not function correctly. 

The QC script can be ran like so: 
```
bash qiime2_qc.sh -p ~/data_win/16S/iScience-review/data/fastq/phred33/ -o ~/data_win/16S/iScience-review/data/q2_out
```

"qiime2_qc.sh" script has the following required flags: ```-p``` = full path to read files. String. ```-o``` = desired output folder (will be created if not existing). String. 


This will output a file called "demux-paired-end.qza" and additionally "demux-paired-end.qzv". The second file can be used to visualise sequence quality, either by calling ```qiime tools view``` or by using the online tool available here: https://view.qiime2.org/. 

Once you have inspected the read quality and determined the optimal values for trimming, you can proceed to run the full pipeline below: 
```
bash qiime2_pipeline.sh -p ~/data_win/16S/iScience-review/data/q2_out/ -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -m 250 -f 240 -r 220 -n 2 -N 5 -t 4
```

"qiime2_pipeline.sh" script has the following required flags: ```-p``` = full path to output folder created above. String. ```-g``` = Sequence of forward primer. String. ```-G``` = Sequence of reverse primer. String. ```-m``` = Minimum sequence length for trimming. Integer. ```-f``` = Number of bp to truncate forward reads to. Integer. ```-r``` = Number of bp to truncate reverse reads to. Integer. ```-n``` = Expected errors, forward reads. Integer. ```-N``` = Expected errors, forward reads. Integer.

The ```-t``` is optional and sets the number of threads. Integer (default:1). 

This script will import your raw reads into Qiime2, perform adapter and primer trimming and then denoise with DADA2. Since Qiime2 does not implement abundance filtering, the resulting tables are exported and filtered to remove those that do not reach 0.25% relative abundance in at least one sample. Although this is believed to not be necessary for denoising methods, recent work has shown that this reduces the impact of spurious sequences on results (Reitmeier et al., 2021, ISME Communications). Filtered tables are then re-imported into Qiime2 to add taxonomy and generate a phylogenetic tree. 

Note that compared to the R-based pipeline, the qiime2 pipeline runs very slow, and make take up to a day to process moderately sized datasets. 

## Data Analysis 

UNDER CONSTRUCTION - We are currently working on implementing a tutorial on an analysis pipeline within Qiime2 to mirror that of the R-based pipeline. For now please use the Qiime2 tutorials available here: https://docs.qiime2.org/2022.2/tutorials/



## Data Normalisation



## Alpha and Beta Diversity

```
qiime diversity core-metrics-phylogenetic --i-table q2_out/feature-table-filt.qza --i-phylogeny q2_out/rooted-tree.qza --m-metadata-file metadata_subset.tsv --p-sampling-depth 120000 --output-dir ./core-metrics-phylo-results 
```

```
qiime diversity core-metrics --i-table q2_out/feature-table-filt.qza --m-metadata-file metadata_subset.tsv --p-sampling-depth 120000 --output-dir ./core-metrics-results
```
### Plotting and Statistics 




## Differentially abundant taxa
