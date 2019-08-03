[<img title="Babraham Bioinformatics" align="left" style="float:right;margin:50px 200px 200px 200px" id="Babraham Bioinformatics" src="Images/bioinformatics_logo.png" height="80">](http://www.bioinformatics.babraham.ac.uk/index.html)

<img title="Odd One Out" align="right" width="200" id="header_img" src="Images/mice_logo.png">

# reStrainingOrder - Mouse Strain Identification

## User Guide - v0.1.0

This User Guide outlines how reStrainingOrder works and gives details for each step.


#### Table of Contents
* [Quick Reference](#quick-reference)
  1. [Supported file types](#which-kind-of-files-are-supported)
  2. [Installation Notes](#step-1-quality-trimming)
  2. [Adapter Trimming](#step-2-adapter-trimming)
    - [Auto-detection](#adapter-auto-detection)
    - [Manual adapter sequence specification](#manual-adapter-sequence-specification)
 * [reStrainingOrder Workflow](#the-reStrainingOrder-workflow-in-more-detail)
  3. [Removing Short Sequences](#step-3-removing-short-sequences)
  4. [Specialised Trimming - hard- and Epigenetic Clock Trimming](#step-4-specialised-trimming)
* [Full list of options for Trim Galore!](#full-list-of-options-for-trim-galore)
  * [RRBS-specific options](#rrbs-specific-options-mspi-digested-material)
  * [Paired-end specific options](#paired-end-specific-options)

## Quick Reference


### Which kind of files are supported?

reStrainingOrder should work with most types of Illumina sequencing reads. More specifically, we have tested it with ChIP- and Input-seq, RNA-seq as well as different tpyes of Bisulfite-seq (WGBS, PBAT). Aligners that were shown to work well with the N-masked genome approach inlcude `Bowtie2`, `HISAT2`, `STAR` and `Bismark`.


#### Feedback
We would like to hear your comments or suggestions! Please e-mail [me here](mailto:felix.krueger@babraham.ac.uk)!


### Installation notes

Just download the latest version under releases, and extract the tar archive into a folder. Done.
```
tar xzf reStrainingOrder_v0.X.Y.tar.gz
```

### Dependencies
reStrainingOrder requires a working version of Perl and [Samtools](http://samtools.sourceforge.net/) to be installed on your machine. It is assumed that the samtools executable is in your `PATH` unless the path is specified manually with:
```
--samtools_path </../../samtools>
```

### Hardware requirements

While the genome preparation is not very resource hungry, the alignment and scoring part are a bit more demanding. Assuming single core operation, alignments to the the multi-strain genome typically take up to 5GB of RAM for Bowtie2 or HISAT2, Bismark alignments may need between 12 and 18GB (directional/PBAT or non directional alignments). The scoring part `reStrainingOrder` currently takes ~20GB of RAM. For future versions we may look into reducing this memory footprint somewhat. 

# The reStrainingOrder workflow in more detail

## Step I - Genome preparation 

This is a one-off process, performed by `reStraining`.

`reStraining` is designed to read in a variant call file from the Mouse Genomes Project (download e.g. from this location: ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz) and generate a new genome version where all positions found as a SNP in any of the strains (currently 35 strains) are masked by the ambiguity nucleobase `N` (**N-masking**).

Here is a sample command for this step:

```
reStraining --vcf mgp.v5.merged.snps_all.dbSNP142.vcf.gz --reference /bi/scratch/Genomes/Mouse/GRCm38/
```

This creates a folder (called `SNPs_directory`) to store the SNPs per chromosome. These files look like this

```
>11
11	3100106	G	C	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
11	3100127	A	C	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
11	3100380	C	A	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
```
And may be deleted afterwards to save disk space if you have no further use for them.

It will also write out matrix file for chromosome 1, which is in a similar format:

```
Chromosome	Position	REF	ALT	129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	AKR_J	A_J	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	KK_HiJ	LEWES_EiJ	LP_J	MOLF_EiJ	NOD_ShiLtJ	NZB_B1NJ	NZO_HlLtJ	NZW_LacJ	PWK_PhJ	RF_SEA_GnJ	SPRET_EiJ	ST_bJ	WSB_EiJ	ZALENDE_EiJ
1	3000023	C	A	1	1	0	0	0	0	1	1	0	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0
1	3000126	G	T	1	1	0	0	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	0	1	1	0	1
1	3000185	G	T	1	1	1	1	0	0	1	1	0	0	0	0	1	1	0	1	0	1	1	1	1	0	0	1
1	3000234	G	A	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
```

This matrix file is used as input for the SNP scoring (Step III: reStrainingOrder, see below). 

#### Indexing the MGP genome

**Bowtie2:**
```
bowtie2-build --threads 8 chr10.N-masked.fa,chr11.N-masked.fa,chr12.N-masked.fa,chr13.N-masked.fa,chr14.N-masked.fa,chr15.N-masked.fa,chr16.N-masked.fa,chr17.N-masked.fa,chr18.N-masked.fa,chr19.N-masked.fa,chr1.N-masked.fa,chr2.N-masked.fa,chr3.N-masked.fa,chr4.N-masked.fa,chr5.N-masked.fa,chr6.N-masked.fa,chr7.N-masked.fa,chr8.N-masked.fa,chr9.N-masked.fa,chrMT.N-masked.fa,chrX.N-masked.fa,chrY.N-masked.fa MGP.N-masked
```
**HISAT2:**
```
hisat2-build --threads 8 chr10.N-masked.fa,chr11.N-masked.fa,chr12.N-masked.fa,chr13.N-masked.fa,chr14.N-masked.fa,chr15.N-masked.fa,chr16.N-masked.fa,chr17.N-masked.fa,chr18.N-masked.fa,chr19.N-masked.fa,chr1.N-masked.fa,chr2.N-masked.fa,chr3.N-masked.fa,chr4.N-masked.fa,chr5.N-masked.fa,chr6.N-masked.fa,chr7.N-masked.fa,chr8.N-masked.fa,chr9.N-masked.fa,chrMT.N-masked.fa,chrX.N-masked.fa,chrY.N-masked.fa MGP.N-masked
```

**Bismark:**
```
bismark_genome_preparation --verbose --parallel 2 .
```

## Step II - Alignments to the MGP genome



### Specific considerations for more specialised applications or software

#### RNA-Seq alignments with STAR: 

Alignment files produced by the Spliced Transcripts Alignment to a Reference [(STAR)](https://github.com/alexdobin/STAR/ "STAR")aligner also should also work well, however a few steps need to be adhered to make this work:

**1)** Since reStrainingOrder only recognises the CIGAR operations M, I, D and N, alignments need to be run in end-to-end mode and not using local alignments (which may result in soft-clipping). This can be accomplished using the option: `--alignEndsType EndToEnd`

**2)** reStrainingOrder requires the `MD:Z:` field of the BAM alignment to work out mismatches involving masked N positions. Since STAR doesn’t report the `MD:Z:` field by default it needs to be instructed to do so, e.g.: `--outSAMattributes NH HI NM MD`

**3)** To save some time and avoid having to sort the reads by name, STAR can be told to leave R1 and R2 following each other in the BAM file using the option: `--outSAMtype BAM Unsorted`

### Alignments with HISAT2:

DNA or RNA alignment files produced by [HISAT2](https://github.com/infphilo/hisat2 "HISAT2") also work well if you make sure that HISAT2 doesn’t perform soft-clipping. At the time of writing HISAT2 does perform soft-clipping (CIGAR operation: `S`) by default, so you need to specify the option `--no-softclip`.

### Alignments with BWA

The other very popular Burrows-Wheeler Aligner ([BWA](http://bio-bwa.sourceforge.net/ "BWA homepage at SourceForge")) is unfortunately **not compatible** with reStrainingOrder processing as was [disscussed in more detail here](https://github.com/FelixKrueger/SNPsplit/issues/19 "BWA alignments are not compatible with SNPsplit"). Briefly, the reason for this is that BWA randomly replaces the ambiguity nucleobase `N` in the reference by either `C`, `A`, `T` or `G`, thereby rendering an N-masked allele-sorting process impossible.

### Bisulfite-Seq data: 

This mode assumes input data has been processed with the bisulfite mapping tool [Bismark](https://github.com/FelixKrueger/Bismark "Bismark project page on Github"). SNPsplit will run a quick check at the start of a run to see if the file provided appears to be a Bismark file, and set the flags `--bisulfite` and/or `--paired` automatically. In addition it will perform a quick check to see if a paired-end file appears to have been positionally sorted, and if not will set the `--no_sort` flag (this data is extracted from the `@PG` header line). Paired-end (`--paired`) or bisulfite (`--bisulfite`) mode can also be set manually. Paired-end mode requires Read 1 and Read 2 of a pair to follow each other in consecutive lines. Please be aware that files should not be name-sorted using [`Sambamba`](https://github.com/FelixKrueger/Bismark/issues/170).

##### Utilisation of SNP positions and allele assignment of bisulfite reads
In contrast to the standard mode of using all known SNP positions, SNPs involving C to T transitions may not be used for allele-specific sorting since they might reflect either a SNP or a methylation state. This includes all of the following Reference/SNP combinations: 

C/T or T/C for forward strand alignments, and G/A or A/G for reverse strand alignments. 

The number of SNP positions that have been skipped because of this bisulfite ambiguity is reported in the report file. These positions may however be used to assign opposing strand alignments since they do not involve C to T transitions directly. For that reason, the bisulfite call processing also extracts the bisulfite strand information from the alignments in addition to the base call at the position involved. For any SNPs involving C positions that are not C to T SNPs both methylation states, i.e. C and T, are allowed to match the C position.

For SNPs which are masked by Ns in the genome no methylation call will be performed, i.e. they will receive a `.` (dot) in the methylation call string. This means that SNP positions are used for allele-sorting only but do not participate in calling methylation. While this may reduce the number of total methylation calls somewhat it effectively eliminates the problem of using positions with potentially incorrect methylation status.


 


#### Single strain mode:
**1)** The VCF file is read and filtered for high-confidence SNPs in the strain specified with   strain <name>
**2)** The reference genome (given with `--reference_genome <genome>`) is read into memory, and the filtered high-confidence SNP positions are incorporated either as N-masking (default), or full sequence (option `--full_sequence`)


## test data set
A test data set will be available for download ...




    
# Credits
reStrainingOrder was written by Felix Krueger at the [Babraham Bioinformatics Group](http://www.bioinformatics.babraham.ac.uk/).

![Babraham Bioinformatics](Images/bioinformatics_logo.png)
