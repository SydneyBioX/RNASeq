# RNASeq
This repo contains Snakemake scripts for processing RNAseq (both bulk and single cell). These scripts can be used on Savona/Verona or the HPC.  

## Installing Snakemake
Snakemake is a pipeline language written in Python. To use it, you will need to make sure it is installed on the server you are using. On savona, you can install using the command. 

```
virtualenv -p python3 .venv
source .venv/bin/activate
pip install snakemake
```
If you do it this way, in a virualenv, next time you want to use snakemake you will have to start the virtualenv again. 
To test that the install worked, you should be able to run `snakemake` on the command line. 

## Directory structure
In this first vanilla version of the pipeline, I am assuming you have paired-end gzipped fastq files. 
Move your fastq files to a new directory called `fastq` and make sure the files are have the names `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`. Now you need the `Snakefile` and the `config.yaml` file in the same directory. The directory structure of your project should look like this: 

```
+-- config.yaml
+-- Snakefile
+-- fastq
|   +-- sample1_1.fastq.gz
|   +-- sample1_2.fastq.gz
+-- processed
|   +-- bam
|   |   +-- sample1_sorted.bam
|   |   +-- sample1_sorted.bam.bai 
```

## Config.yaml 
Change this file to reflect the species you want to align, and the names of the sample you have in your dataset. 

## Running the pipeline 

To run the pipeline is just one line of code. 

```
snakemake --cores 8 
```
