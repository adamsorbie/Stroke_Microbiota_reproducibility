# Stroke_Microbiota_reproducibility
Repository for Sorbie, Jimenez and Benakis, iScience submitted

<figure>
<img src="https://user-images.githubusercontent.com/23216921/141769492-48fb6c60-d466-4242-a270-ad910ed666fc.png" width="350" height="550">
<figcaption>Bi-directional interaction of Microbiota and Stroke.</figcaption>
</figure>

Recently, the interplay between the microbiota and ischemic stroke is increasingly being recognised, however differences in methodology, study populations and data analysis have led to a lack of consistent findings between different studies. In our manuscript we aim to provide a framework to conduct stroke-microbiota studies in a rigorous and reproducible manner. This github repository accompanies the manuscript, providing code and a tutorial for processing and analysis of 16S sequencing data.  

## Data processing and Analysis 

Here we provide two pipeline scripts to process 16S microbiota data. To reach a wider audience, one is written in R and the other in shell/python via QIIME2 (currently only alpha version) so researchers can choose whichever they are most comfortable with. 


### Installation 

Although this analysis pipeline can in theory be ran on Windows 10/11, it works better on unix-based systems. If you are a Windows user, we recommend installing the Windows subystem for Linux. 

To install on Windows 10:

Open a powershell (search for powershell in the search bar) **as an administrator** (right click on the program) and copy and paste the following command:

```
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
```

You will then be prompted to restart your computer. After booting up again, open the Microsoft store app and search for your preferred linux distribution, for example Ubuntu or Debian. 

Once installed, open the terminal and provide a username and password. From here on out, all commands are entered into the __linux__ terminal/shell. To download the analysis scripts, type or copy paste the following into your terminal:

```
git clone https://github.com/adamsorbie/Stroke_Microbiota_reproducibility.git
```

if that doesn't work then you may need to install git first, you can do this by typing:

```
sudo apt update
sudo apt upgrade
```

followed by:

```
sudo apt install git
```
then repeat the above command. Note that this command will only work if you installed the Ubuntu or Debian distributions. 


### Dependencies 

For the R version the following dependencies must be installed and available in your path: 

FastTree
MAFFT 

FASTQC
MultiQC

The first two can be installed like so: 

```
sudo apt install fasttree
sudo apt install mafft
```

The second two are a little more complicated and require a working conda installation first: 

First download the conda installer
```
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh
```

Run this like so and follow the on-screen instructions: 
```
bash Miniconda3-py38_4.10.3-Linux-x86_64.sh
```
Restart the shell and finally install FASTQC and Multiqc: 

```
conda install -c bioconda fastqc 
conda install -c multiqc=
```



### Running the data-processing pipeline 

Here we provide a tutorial of the R-based and QIIME2 based pipelines. 

#### Before starting 

* Ensure reads are named according to illumina naming format: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm. If your filenames do not follow this convention, the script "rename_seqs_illumina_format.sh" is provided to do this for you. 
* Note that this script uses positional instead of named arguments. 
To run: ```bash rename_seqs_illumina_format.sh "-" 3 R1.fastq R2.fastq ~/data_win/16S/iScience-review/data/fastq```

The first argument is the delimeter separating sample name from the rest of the file name while the number represents the last field to split. Then we have the consistent pattern in the forward and reverse reads, usually R1* and R2* and finally the path where the files are located. Note that this script assumes your files are gzipped. 

Using this naming scheme as an example ```Sample1-Mouse-Stool-illumina-R1.fastq.gz```, running the above script would yield ```Sample1-Mouse-Stool_S1_L001_R1_001.fastq.gz```

#### QC of FASTQ reads 

Raw reads should be checked for any quality issues before any further processing. The script "qc.R" can be used to analyze read quality via DADA2. 

The "qc.R" script has the following __required__ flags: 
    ```-p``` = full path to read files. String
    ```-o``` = full path to output directory. If not existing will be created. String

Example used in our study.  

```
Rscript dada2_qc.R -p ~/data_win/16S/iScience-review/data/fastq/ -o ~/data_win/16S/iScience-review/data/fastq/qc_out/
```

You will of course have to modify the input and output paths to where your files are located in the filesystem. 

Output: 

![quality_profile_f](https://user-images.githubusercontent.com/23216921/137965016-32fb18a2-da10-4167-84bb-92045dd96675.jpg)

![quality_profile_r](https://user-images.githubusercontent.com/23216921/137965029-660d5a34-30a9-4d1b-b8f5-4405e81c47a1.jpg)

Quality is good in both the forward and reverse reads here but starts to drop slightly around 250bp on the forward and 220bp on the reverse.  

#### Trimming adapter and primer sequences 

Removal of non-biological sequences is an important step before denoising reads and calling Amplicon sequence variants (ASVs)as their presence can inflate samples diversity and lead to many spurious ASVs. Here we use cutadapt to remove these sequences. A bash script "run_cutdapt.sh" is provided to automate this process. 

"run_cutadapt.sh" script has the following __required__ flags: 
    ```-g``` = Forward primer sequence. String
    ```-G``` = Reverse primer sequence. String
    ```-m``` = Minimum length of sequences - larger sequences will be removed (Set this ~10% lower than expected length). String.
    ```-M``` = Maximum length of sequences - larger sequences will be removed (Set this ~10% higher than expected length). String.
    ```-p``` = full path to read files. String.

Example used in our study.  
```
bash run_cutadapt.sh -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -m 250 -M 300 -p ~/data_win/16S/iScience-review/data/fastq/
```

Cutadapt outputs a log file found in: "fastq_directors/trimmed_primer". We recommend double-checking this to ensure primer removal went smoothly. The output should look something like this: 

```
This is cutadapt 3.4 with Python 3.6.13
Command line parameters: -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT -m 250 -M 300 -j 0 -o RD118_trimmed_primer_R1_001.fastq.gz -p RD118_trimmed_primer_R2_001.fastq.gz RD118_S1_L001_R1_001.fastq.gz RD118_S1_L001_R2_001.fastq.gz
Processing reads on 12 cores in paired-end mode ...
Finished in 7.04 s (25 us/read; 2.43 M reads/minute).

=== Summary ===

Total read pairs processed:            285,040
  Read 1 with adapter:                 284,274 (99.7%)
  Read 2 with adapter:                 284,157 (99.7%)
Pairs that were too short:                  12 (0.0%)
Pairs that were too long:                    0 (0.0%)
Pairs written (passing filters):       285,028 (100.0%)

Total basepairs processed:   168,743,680 bp
  Read 1:    84,371,840 bp
  Read 2:    84,371,840 bp
Total written (filtered):    157,669,485 bp (93.4%)
  Read 1:    78,981,670 bp
  Read 2:    78,687,815 bp

=== First read: Adapter 1 ===

Sequence: GTGCCAGCMGCCGCGGTAA; Type: regular 5'; Length: 19; Trimmed: 284274 times

No. of allowed errors:
1-9 bp: 0; 10-19 bp: 1

Overview of removed sequences
length	count	expect	max.err	error counts
13	2	0.0	1	2
14	4	0.0	1	2 2
15	11	0.0	1	6 5
16	59	0.0	1	12 47
17	744	0.0	1	258 486
18	13228	0.0	1	4747 8481
19	269610	0.0	1	264518 5092
20	601	0.0	1	282 319
21	1	0.0	1	1
22	2	0.0	1	2
171	10	0.0	1	10
176	1	0.0	1	1
196	1	0.0	1	1
```

#### Running DADA2 

Now that sample quality has been checked and non-biological sequences removed we can proceed with DADA2. 


An R script "run_dada2.R" is provided to automate each step in the process. 

"run_dada2.R" Rscript has the following __required__ flags: 
    ```-p``` = full path to trimmed read files - (trimmed_primer directory generated in last step). String
    ```-f``` = Length to truncate forward reads to. Integer.
    ```-r``` = Length to truncate reverse reads to. Integer.
    ```-o``` = output folder. String.
    
  and additionally these optional flags: 
    ```-q``` = whether to output dada2 qc metrics/plots. Boolean. (defaults to True) 
    ```-n``` = expected errors forward reads. Integer. (default: 2)
    ```-N``` = expected errors forward reads. Integer. (default: 2)
    ```-t``` = Number of threads to use in parallel. String (default: 2)
  
  
```
Rscript run_dada2.R -p ~/data_win/16S/iScience-review/data/fastq/trimmed_primer/ -f 240 -r 220 -o ~/data_win/16S/iScience-review/data/dada2_out/ -n 2 -N 5 -t 6
```

This should create a folder called dada2_out (or whatever you chose to name) which contains the results. In the folder you should find the following files: 

* "aligned.fasta" - Aligned ASV sequences (used as input to generate tree)
* "ASV_seqs.fasta" - Unique sequences for each ASV 
* "ASV_seqtab.tab" - Feature count table
* "ASV_seqtab_tax.tab" - Feature count table with additional taxonomy column
* "ASV_tree.tre" - Phylogenetic tree of ASVs in newick format
* "error_plot_f.jpg/error_plot_r.jpg" - Plot of estimated error rates from dada2 "learnErrors" function 
* "parameter_log.txt" - log of input parameters, i.e. commands provided to run_dada2.R
* "sequence_len_distr.png" - plot of the sequence length distribution 
* "study_stats.txt" - Filtering and merging statistics 

Before proceeding with the results it's important to take a look at some of these files to check everything went ok with the pipeline: 

The ```study_stats.txt``` is one of the most important. In any case, and indeed at any step, the majority of reads should be retained. You can check the percentage of reads which made it through the entire pipeline in the "percent_reads_retained" column. Generally, around 60% is good and indicates that these steps likely worked correctly. As we implement an additional filtering step based on the findings from [Reitmeier et al 2021, ISME Communications](https://www.nature.com/articles/s43705-021-00033-z?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+ismecomms%2Frss%2Fcurrent+%28ISME+Communications+-+Latest%29) where only ASVs with a relative abundance of greater than 0.25% are retained. 

Users should additionally check the error rate plots: in which the black line should fit fairly well with the black points, the sequence length distribution plot: which should generally show one peak and finally the feature table with taxonomy. If many taxa are unclassified and your samples were not from an uncommonly profiled environment, then this may indicate that an error occured with taxonomic classification. 

Provided everything looks ok, we can proceed with the analysis of the resulting data. 

### Data-Analysis 

in addition to our data processing pipeline, we have included a set of R functions, which perform many common tasks (e.g. Data normalisation, calculation of alpha and beta diversity metrics). A tutorial on how to use these is provided [here](https://github.com/adamsorbie/Stroke_Microbiota_reproducibility/blob/main/R_based_pipeline/analysis/stroke_microbiota_analysis.md) to help users learn 16S data analysis or indeed analyse their own data. 

For those who prefer to use the Qiime2 based pipeline, the associated tutorial can be found [here](under construction) 
 
