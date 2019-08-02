# reStrainingOrder - Mouse Strain Identification

## User Guide - v0.1.0
#### 02 August, 2019

This User Guide outlines how reStrainingOrder tools work and gives more details for each individual step.

[<img title="Babraham Bioinformatics" style="float:right;margin:20px 20 20 600px" id="Babraham Bioinformatics" src="Images/logo.png" height="88" >](http://www.bioinformatics.babraham.ac.uk/index.html)


Last update: 02/07/2019

#### Table of Contents
* [Introduction](#version-060)
* [Methodology](#adaptive-quality-and-adapter-trimming-with-trim-galore)
  1. [Quality Trimming](#step-1-quality-trimming)
  2. [Adapter Trimming](#step-2-adapter-trimming)
    - [Auto-detection](#adapter-auto-detection)
    - [Manual adapter sequence specification](#manual-adapter-sequence-specification)
  3. [Removing Short Sequences](#step-3-removing-short-sequences)
  4. [Specialised Trimming - hard- and Epigenetic Clock Trimming](#step-4-specialised-trimming)
* [Full list of options for Trim Galore!](#full-list-of-options-for-trim-galore)
  * [RRBS-specific options](#rrbs-specific-options-mspi-digested-material)
  * [Paired-end specific options](#paired-end-specific-options)

# 1) Quick Reference


Bismark needs a working version of Perl and it is run from the command line. Furthermore, [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2) or [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) needs to be installed on your computer. For more information on how to run Bismark with Bowtie 2 please go to the end of this manual.

First you need to download a reference genome and place it in a genome folder. Genomes can be obtained e.g. from the [Ensembl](http://www.ensembl.org/info/data/ftp/index.html/) or [NCBI](ftp://ftp.ncbi.nih.gov/genomes/) websites. For the example below you would need to download the _Homo sapiens_ genome. Bismark supports reference genome sequence files in `FastA` format, allowed file extensions are either either `.fa` or `.fasta`. Both single-entry and multiple-entry `FastA` files are supported.

We would like to hear your comments or suggestions! Please e-mail [felix.krueger@babraham.ac.uk](mailto:felix.krueger@babraham.ac.uk)


## Which kind of files are supported?

reStrainingOrder should work with 

## Installation notes


```
tar xzf reStrainingOrder_v0.X.Y.tar.gz
```

## Dependencies
reStrainingOrder requires a working version of Perl and [Samtools](http://samtools.sourceforge.net/) to be installed on your machine. It is assumed that the samtools executable is in your `PATH` unless the path is specified manually with:
```
--samtools_path </../../samtools>

```

## Hardware requirements


# The reStrainingOrder workflow in more detail

## Specific considerations for more specialised applications or software
### Paired-end:  

In paired-end mode, both reads are used for the classification. Read pairs with conflicting reads (tag CF) or pairs containing both tags G1 and G2 are considered conflicting and are not reported by default. Reporting of these reads can be enabled using the option `--conflicting`.

Singleton alignments in the allele-tagged paired-end file (which is the default for e.g. TopHat) are also sorted into the above four files. Specifying `--singletons` will write these alignments to special singleton files instead (ending in `*_st.bam`).


### Hi-C data: 

Assumes data processed with [HiCUP](www.bioinformatics.babraham.ac.uk/projects/hicup/ "HiCUP on the Babraham Bioinformatics website") as input, i.e. the input BAM files are by definition paired-end and Reads 1 and 2 follow each other. Hi-C sorting discriminates several more possible read combinations:

```G1-G1
G2-G2
G1-UA
G2-UA
G1-G2
UA-UA
```


Again, read pairs containing a conflicting reads (tag `XX:Z:CF`) are not printed out by default, but this may be enabled using the option `--conflicting`. For an example report please see below.


### RNA-Seq alignments with STAR: 

Alignment files produced by the Spliced Transcripts Alignment to a Reference [(STAR) aligner](https://github.com/alexdobin/STAR/ "STAR") also work well with SNPsplit, however a few steps need to be adhered to make this work:

**1)** Since SNPsplit only recognises the CIGAR operations M, I, D and N (see above), alignments need to be run in end-to-end mode and not using local alignments (which may result in soft-clipping). This can be accomplished using the option: `--alignEndsType EndToEnd`

**2)** SNPsplit requires the `MD:Z:` field of the BAM alignment to work out mismatches involving masked N positions. Since STAR doesn’t report the `MD:Z:` field by default it needs to be instructed to do so, e.g.: `--outSAMattributes NH HI NM MD`

**3)** To save some time and avoid having to sort the reads by name, STAR can be told to leave R1 and R2 following each other in the BAM file using the option: `--outSAMtype BAM Unsorted`

### Alignments with HISAT2:

DNA or RNA alignment files produced by [HISAT2](https://github.com/infphilo/hisat2 "HISAT2") also work well with SNPsplit if you make sure that HISAT2 doesn’t perform soft-clipping. At the time of writing the current version of HISAT2 (2.0.3-beta) does perform soft-clipping (CIGAR operation: S) even though this is not well documented (or in fact the documentation on Github suggests that the default mode is end-to-end which should not perform any soft-clipping whatsoever). Until the end-to-end mode works as expected users will have to set the penalty for soft-clipping so high that it is effectively not performed (`--sp` is the option governing the soft-clipping penalty). We suggest adding the following option to the HISAT2 command: `--sp 1000,1000`

**EDIT:** HISAT2 does now also have an option `--no-softclip` which should have the same effect.

### Alignments with BWA

The other very popular Burrows-Wheeler Aligner ([BWA](http://bio-bwa.sourceforge.net/ "BWA homepage at SourceForge")) is unfortunately **not compatible** with SNPsplit alignment sorting as was [disscussed in more detail here](https://github.com/FelixKrueger/SNPsplit/issues/19 "BWA alignments are not compatible with SNPsplit"). Briefly, the reason for this is that BWA randomly replaces the ambiguity nucleobase `N` in the reference by either `C`, `A`, `T` or `G`, thereby rendering an N-masked allele-sorting process impossible.

### Bisulfite-Seq data: 

This mode assumes input data has been processed with the bisulfite mapping tool [Bismark](https://github.com/FelixKrueger/Bismark "Bismark project page on Github"). SNPsplit will run a quick check at the start of a run to see if the file provided appears to be a Bismark file, and set the flags `--bisulfite` and/or `--paired` automatically. In addition it will perform a quick check to see if a paired-end file appears to have been positionally sorted, and if not will set the `--no_sort` flag (this data is extracted from the `@PG` header line). Paired-end (`--paired`) or bisulfite (`--bisulfite`) mode can also be set manually. Paired-end mode requires Read 1 and Read 2 of a pair to follow each other in consecutive lines. Please be aware that files should not be name-sorted using [`Sambamba`](https://github.com/FelixKrueger/Bismark/issues/170).

##### Utilisation of SNP positions and allele assignment of bisulfite reads
In contrast to the standard mode of using all known SNP positions, SNPs involving C to T transitions may not be used for allele-specific sorting since they might reflect either a SNP or a methylation state. This includes all of the following Reference/SNP combinations: 

C/T or T/C for forward strand alignments, and G/A or A/G for reverse strand alignments. 

The number of SNP positions that have been skipped because of this bisulfite ambiguity is reported in the report file. These positions may however be used to assign opposing strand alignments since they do not involve C to T transitions directly. For that reason, the bisulfite call processing also extracts the bisulfite strand information from the alignments in addition to the base call at the position involved. For any SNPs involving C positions that are not C to T SNPs both methylation states, i.e. C and T, are allowed to match the C position.

For SNPs which are masked by Ns in the genome no methylation call will be performed, i.e. they will receive a `.` (dot) in the methylation call string. This means that SNP positions are used for allele-sorting only but do not participate in calling methylation. While this may reduce the number of total methylation calls somewhat it effectively eliminates the problem of using positions with potentially incorrect methylation status.


 
## SNPsplit genome preparation

`SNPsplit_genome_preparation` is designed to read in a variant call file from the Mouse Genomes Project (e.g. this [latest file](ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz)) and generate new genome versions where the strain SNPs are either incorporated into the new genome (full sequence) or masked by the ambiguity nucleobase `N` (**N-masking**).

`SNPsplit_genome_preparation` may be run in two different modes:

#### Single strain mode:
**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with   strain <name>
**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated either as N-masking (default), or full sequence (option `--full_sequence`)


## test data set
A test data set will be available for download ...




    
# Credits
reStrainingOrder was written by Felix Krueger at the [Babraham Bioinformatics Group](http://www.bioinformatics.babraham.ac.uk/).

![Babraham Bioinformatics](Images/bioinformatics_logo.png)
