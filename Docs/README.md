# reStrainingOrder - Mouse Strain Identification

## User Guide - v0.1.0
#### 24 July, 2019

This User Guide outlines how reStrainingOrder tools and gives more details for each individual step.


# 1) Quick Reference


Bismark needs a working version of Perl and it is run from the command line. Furthermore, [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2) or [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) needs to be installed on your computer. For more information on how to run Bismark with Bowtie 2 please go to the end of this manual.

First you need to download a reference genome and place it in a genome folder. Genomes can be obtained e.g. from the [Ensembl](http://www.ensembl.org/info/data/ftp/index.html/) or [NCBI](ftp://ftp.ncbi.nih.gov/genomes/) websites. For the example below you would need to download the _Homo sapiens_ genome. Bismark supports reference genome sequence files in `FastA` format, allowed file extensions are either either `.fa` or `.fasta`. Both single-entry and multiple-entry `FastA` files are supported.

We would like to hear your comments or suggestions! Please e-mail [felix.krueger@babraham.ac.uk](mailto:felix.krueger@babraham.ac.uk)

## Installation notes


```
tar xzf reStrainingOrder_v0.X.Y.tar.gz
```

## Dependencies
Bismark requires a working of Perl and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2) (or [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)) to be installed on your machine. Bismark will assume that the Bowtie 2/ HISAT2 executable is in your path unless the path to Bowtie/ HISAT2 is specified manually with:
```
--path_to_bowtie2 </../../bowtie2> or 
--path_to_hisat2 </../../hisat2>
```

## Hardware requirements


## test data set
A test data set is available for download ...


## Which kind of files are supported?


    
# Credits
reStrainingOrder was written by Felix Krueger at the [Babraham Bioinformatics Group](http://www.bioinformatics.babraham.ac.uk/).

![Babraham Bioinformatics](Images/bioinformatics_logo.png)
