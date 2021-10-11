# Stroke_Microbiota_reproducibility
Repository for Sorbie, Jimenez and Benakis, iScience 2021 


## Data processing and Analysis 

Here we provide two pipeline scripts to process 16S microbiota data. To reach a wider audience, one is written in R and the other in shell/python via QIIME2 so researchers can choose whichever they are most comfortable with. 


### Installation 

Although this analysis pipeline can in theory be ran on Windows 10/11, it works better on unix-based systems. If you are a Windows user, we recommend installing the Windows subystem for Linux. 

To install on Windows 10:

Open a powershell (search for powershell in the search bar) **as an administrator** (right click on the program) and copy and paste the following command:

```
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
```

You will then be prompted to restart your computer. After booting up again, open the Microsoft store app and search for your preferred linux distribution, for example Ubuntu or Debian. 

Once installed, open the terminal and provide a username and password. To download the analysis scripts, type or copy paste the following into your terminal:

```
git clone https://github.com/adamsorbie/Stroke_Microbiota_reproducibility.git
```

if that doesn't work then you may need to install git first, you can do this by typing:

```
sudo apt update
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

They can be installed like so: 

```
sudo apt install fasttree
sudo apt install mafft
```

### Running the analysis pipeline 

Here we provide a tutorial of the R-based and QIIME2 based pipelines. 

#### Before starting 

* Ensure reads are named according to illumina naming format: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm. If your filenames do not follow this convention, the script "rename_seqs_illumina_format.sh" is provided to do this for you. 
* 

#### QC of FASTQ reads 

Raw reads should be checked for any quality issues before any further processing. The script "qc.sh" can be used to run FASTQC/MultiQC, providing an overview of read quality and also the FIGARO tool from Zymo which suggests trimming and filtering parameters for DADA2 downstream. Unfortunately, reads with different lengths are not currently supported by FIGARO so if this is the case with your data, you will need to determine these parameters yourself. 

The "qc.sh" script has the following __required__ flags: 
    ```-a``` = Expected amplicon length - The length of the (merged) amplicon not including primers. Integer.
    ```-f``` = Length of the forward primer (number of bases). Integer.
    ```-r``` = Length of the reverse primer (number of bases). Integer.
    ```-p``` = full path to read files. String
    ```-o``` = full path to output directory. If not existing will be created. String

An optional ```-m ``` flag is also available. This sets the minium overlap between paired reads and is useful when the overlap between the forward and reverse reads is less than ~50. 

Example used in our study.  

```
bash qc.sh -a 291 -f 19 -r 20 -p ~/data_win/16S/iScience-review/data/fastq/ -o ~/data_win/16S/iScience-review/data/fastq/qc_out/
```

Output: 

![fastqc_per_base_sequence_quality_plot](https://user-images.githubusercontent.com/23216921/136770720-636ec9f8-9ccf-446b-9690-bcf434b402b8.png)

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
