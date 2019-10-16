# Yunnan Metabarcoding School : Day 2 - from sequences to OTU table
This repository contains the data, pre-processing steps, techniques and software that we are going to use on day 2 of YMBS. If you cannot get the data from the repository, then please get them from me (USB stick). 

The major steps that we are going to do today. 
1. Pre-processing steps: Quality checks, adapter removal, quality trimming, error correction, merging of reads.
2. Sorting and Filtering: Demultiplexing, filtering for features
3. Clustering: From sequences to OTU table. 

## Data 
Before we get started on the exercises, let us make sure that we have the sequencing data. Let us store the raw sequencing data and move it to a directory called rawdata.
```bash
mkdir rawdata
mv TOG*gz rawdata
```
There are 3 paired end sequencing files, each belonging to a different pool. So, we will start with making 3 directories, called pool1, pool2 and pool3. 
```bash
mkdir pool1 pool2 pool3
```
Now we are ready to start processing. 

## Pre-processing steps
Here we outline the steps we will take before we move on to treating our DNA as metabarcoding samples.

### Quality Checks using FastQC
First, we want to make sure that the read sequences that we got from our sequencing facility look good. So, we will run FastQC for a quick quality check.
```bash
cd rawdata
fastqc TOG-course-pool1-R1.fastq.bz2
```
Let us start with examining the results from fastqc. Now, run fastqc for the other sequencing files. 

### Adapter removal
Remember what the sequenced fragment looked like? Before we jump headfirst into analysing the data and getting OTUs, we need to strip away all the non-biological parts of the sequenced fragment. The first step in that is to removed the Illumina sequencing adapters. To do this, we will use Adapter Removal. Make sure you are in your project directory before you start.
```bash
AdapterRemoval --file1 TOG-course-pool1-R1.fastq.bz2 --file2 TOG-course-pool1-R2.fastq.bz2 --basename pool1/pool1 --gzip
```
This takes in the sequencing file for read1 and read2, and the output prefix. Since we want to be nice to our computers and not take up too much space, so let us zip our output files as well. What do the output files look like? Which files would you continue working with?

### 
