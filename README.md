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

#### QC of FASTQ reads 

Raw reads should be checked for any quality issues before any further processing. The script "qc.sh" can be used to run FASTQC/MultiQC, providing an overview of read quality and also the FIGARO tool from Zymo which suggests trimming and filtering parameters for DADA2 downstream. Unfortunately, reads with different lengths are not currently supported by FIGARO so if this is the case with your data, you will need to determine these parameters yourself. 

The "qc.sh" script has the following __required__ flags: 
    ```-a``` = Expected amplicon length - The length of the amplicon not including primers  
    ```-f``` = Length of the forward primer (number of bases)
    ```-r``` = Length of the reverse primer (number of bases)
    ```-p``` = full path to read files 
    ```-o``` = full path to output directory. If not existing will be created 

An additional ```-m ``` flag is also available. This sets the minium overlap between paired reads and is useful when the overlap between the forward and reverse reads is less than ~50. 

