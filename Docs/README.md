[<img title="Babraham Bioinformatics" align="left" style="float:right;margin:50px 200px 200px 200px" id="Babraham Bioinformatics" src="Images/bioinformatics_logo.png" height="80">](http://www.bioinformatics.babraham.ac.uk/index.html)

<img title="Odd One Out" align="right" width="200" id="header_img" src="Images/mice_logo.png">

# reStrainingOrder - Mouse Strain Identification

## User Guide - v0.1.0

This User Guide outlines how reStrainingOrder works and gives details for each step.


#### Table of Contents
- [Quick Reference](#quick-reference)
- [reStrainingOrder Workflow](#the-reStrainingOrder-workflow-in-more-detail)
  1. [Step I: Genome preparation](#Step-I---Genome-preparation)
    - [Running reStraining](#a\)-running-reStraining)
       1. [Supported file types](#which-kind-of-files-are-supported)
       2. [Installation Notes](#installation-notes)
    - [Indexing the MGP genome](#b\)-indexing-the-MGP-genome)
  2. [Step II: Alignments to the MGP N-masked genome](#Step-II---Alignments-to-the-MGP-genome)
  3. [Step III: Scoring SNPs](#Step-III---Scoring-SNPs)


## Quick Reference


#### Which kind of files are supported?

reStrainingOrder should work with most types of Illumina sequencing reads. More specifically, we have tested it with ChIP- and Input-seq, RNA-seq as well as different tpyes of Bisulfite-seq (WGBS, PBAT). Aligners that were shown to work well with the N-masked genome approach include `Bowtie2`, `HISAT2`, `STAR` and `Bismark`.




#### Installation notes

Just download the latest version under [releases](https://github.com/FelixKrueger/reStrainingOrder/releases), and extract the tar archive into a folder. Done.
```
tar xzf reStrainingOrder_v0.X.Y.tar.gz
```

#### Dependencies

reStrainingOrder requires a working version of Perl and [Samtools](http://samtools.sourceforge.net/) to be installed on your machine. It is assumed that the samtools executable is in your `PATH` unless the path is specified manually with:
```
--samtools_path </../../samtools>
```

#### Hardware requirements

While the genome preparation is not very resource hungry, the alignment and scoring part are a bit more demanding. Assuming single core operation, alignments to the multi-strain genome typically take up to 5GB of RAM for Bowtie2 or HISAT2, Bismark alignments may need between 12 and 18GB (directional/PBAT or non-directional alignments). The scoring part `reStrainingOrder` currently takes ~20GB of RAM. For future versions we may look into reducing this memory footprint somewhat. 

#### Feedback

We would like to hear your comments or suggestions! Please e-mail [me here](mailto:felix.krueger@babraham.ac.uk)!


# The reStrainingOrder workflow in more detail

## Step I - Genome preparation 

### a) Running `reStraining`

This is a one-off process.

`reStraining` is designed to read in a variant call file from the Mouse Genomes Project (download e.g. from this location: ftp://ftp-mouse.sanger.ac.uk/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz (FTP links are not rendered nicely in Github markdown)) and generate a new genome version where all positions found as a SNP in any of the strains (currently 35 different ones) are masked by the ambiguity nucleobase `N` (**N-masking**). The entire process of filtering through ~80 million SNP positions and preparing the N-masked genome typically takes four hours on our server and requires some 6GB of memory.

Here is a sample command for this step:

```
reStraining --vcf mgp.v5.merged.snps_all.dbSNP142.vcf.gz --reference /bi/scratch/Genomes/Mouse/GRCm38/
```

This command:
 * creates a folder for the new N-masked genome
 * produces a high confidence SNP matrix for chromosome 1
 * creates a folder to store the SNP information per chromosome
 * generates a SNP filtering and genome preparation report

**N-masked genome folder**

This folder (called `MGP_strains_N-masked`) and its FastA contents are vital for subsequent steps. For sample commands to index the new N-masked sequence files please [see below](#b\)-indexing-the-mgp-genome). 


**Chromosome 1 matrix file**

The genome preparation command writes out a matrix file for chromosome 1 only (called `MGPv5_SNP_matrix_chr1.txt.gz`), which is in the following format:

```
Chromosome	Position	REF	ALT	129P2_OlaHsd	129S1_SvImJ	129S5SvEvBrd	AKR_J	A_J	BALB_cJ	BTBR_T+_Itpr3tf_J	BUB_BnJ	C3H_HeH	C3H_HeJ	C57BL_10J	C57BL_6NJ	C57BR_cdJ	C57L_J	C58_J	CAST_EiJ	CBA_J	DBA_1J	DBA_2J	FVB_NJ	I_LnJ	KK_HiJ	LEWES_EiJ	LP_J	MOLF_EiJ	NOD_ShiLtJ	NZB_B1NJ	NZO_HlLtJ	NZW_LacJ	PWK_PhJ	RF_SEA_GnJ	SPRET_EiJ	ST_bJ	WSB_EiJ	ZALENDE_EiJ
1	3000023	C	A	1	1	0	0	0	0	1	1	0	1	0	0	1	0	0	0	1	0	0	0	0	0	0	0
1	3000126	G	T	1	1	0	0	0	1	1	1	0	1	1	1	1	1	1	0	1	1	0	0	1	1	0	1
1	3000185	G	T	1	1	1	1	0	0	1	1	0	0	0	0	1	1	0	1	0	1	1	1	1	0	0	1
1	3000234	G	A	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
```
A score of 0 for a strain indicates that a given strain has the `REF` base at this position, a call of 1 means that it contains the `ALT` base with high confidence. This matrix file is used as input for the SNP scoring process (reStrainingOrder, [see below](#step-III---scoring-snps) ). 

The matrix is written out for a single chromosome only to use less memory in the scoring process. In theory one could use any other chromosome as well (or even the whole genome, but with 70M positions this would be challenging...!). This is the SNP filtering summary:

```
SNP position summary for all MGP strains (based on mouse genome build GRCm38)
===========================================================================

Positions read in total:	78,772,544
Positions skipped because the REF/ALT bases were not well defined:	960,167
Positions discarded as no strain had a high confidence call:	7,936,128

Positions printed to THE CHR1 MATRIX in total:	5,506,653
```

**Please note:** that only positions that have a single `REF/ALT` genotype were considered (i.e. positions with several ALT positions for different strains (e.g. `REF: A`, `ALT: C,T`) were skipped for simplicity. Also, positions where the reference sequence did not have a DNA base or positions with no high confidence SNP call in any of the strains were skipped entirely. 

In total, the chr1 matrix file contains ~5.5 million positions that were of high quality in one or more strains.

**SNP folder**

The folder `SNPs_directory` stores files containing SNPs per chromosome. The files are in this format (similar to the chr1 matrix file but without header, shown here for chr11):

```
>11
11	3100106	G	C	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
11	3100127	A	C	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
11	3100380	C	A	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
```

If you have no further use for these files they may be deleted afterwards to save disk space (they were used for the N-masking process and have already served their purpose).


**SNP filtering and genome preparation reports**

These files are intended for record keeping purposes.


### b) Indexing the MGP genome

This is again a one-off process.

Here are sample commands for some popular aligners using 4 cores each (just copy-paste). Depending on your resources this process may take up to a few hours.

**Bowtie2:**
```
bowtie2-build --threads 4 chr10.N-masked.fa,chr11.N-masked.fa,chr12.N-masked.fa,chr13.N-masked.fa,chr14.N-masked.fa,chr15.N-masked.fa,chr16.N-masked.fa,chr17.N-masked.fa,chr18.N-masked.fa,chr19.N-masked.fa,chr1.N-masked.fa,chr2.N-masked.fa,chr3.N-masked.fa,chr4.N-masked.fa,chr5.N-masked.fa,chr6.N-masked.fa,chr7.N-masked.fa,chr8.N-masked.fa,chr9.N-masked.fa,chrMT.N-masked.fa,chrX.N-masked.fa,chrY.N-masked.fa MGP.N-masked
```
**HISAT2:**
```
hisat2-build --threads 4 chr10.N-masked.fa,chr11.N-masked.fa,chr12.N-masked.fa,chr13.N-masked.fa,chr14.N-masked.fa,chr15.N-masked.fa,chr16.N-masked.fa,chr17.N-masked.fa,chr18.N-masked.fa,chr19.N-masked.fa,chr1.N-masked.fa,chr2.N-masked.fa,chr3.N-masked.fa,chr4.N-masked.fa,chr5.N-masked.fa,chr6.N-masked.fa,chr7.N-masked.fa,chr8.N-masked.fa,chr9.N-masked.fa,chrMT.N-masked.fa,chrX.N-masked.fa,chrY.N-masked.fa MGP.N-masked
```

**Bismark:**
```
bismark_genome_preparation --verbose --parallel 2 .
```

## Step II - Alignments to the MGP genome

Here I will only briefly mention a few aspects that affect the alignment step (in no particular order):

- Reads should be adapter- and quality-trimmed before aligning them to the `MGP.N-masked` genome. As we are scoring single bp matches/mismatches later on, sequences should be of good quality. To avoid a "garbage-in-garbage-out" scenario, a single run through [Trim Galore](https://github.com/FelixKrueger/TrimGalore) or similar tool should suffice.

- Alignments to N-masked genomes take (considerably) longer than to their non-N-masked counterparts. In the interest of time, it might thus be advisable to use only a subset of the original input FastQ file(s). We have run tests files down-sampled to 10, 5 or 2 million reads, which all arrived at the same answer. Keep in mind though that reStrainingOrder only assays chromosome 1 which is ~7% of the entire mouse genome. Example: if you started from say merely 1 million reads, only a few (ten-)thousand of them might eventually end up aligning to chromosome 1, and out of those you would then only look at the ones covering N-masked positions... I think our suggestion would be: down-sampling: yes, but not not too far down (maybe 5-20M as a rough guideline?).

- The aligner you wish to employ needs to be capable of aligning to genomes containing `N` nucleotides. This has already been documented for [SNPsplit here](https://github.com/FelixKrueger/SNPsplit/blob/master/SNPsplit_User_Guide.md#specific-considerations-for-more-specialised-applications-or-software), but I will quickly list the most important points again here.


### Specific considerations for more specialised applications or software

### Alignments with Bowtie2

The output format of Bowtie2 is compatible with reStrainingOrder in its default end-to-end mode. Alignments in `--local` mode are not supported.

### Alignments with Bismark

The output format of Bismark files with either Bowtie2 or HISAT2 are generally compatible with reStrainingOrder in their default mode (end-to-end mapping). Alignments in `--local` mode are not supported.

### Alignments with HISAT2:

DNA or RNA alignment files produced by [HISAT2](https://github.com/infphilo/hisat2 "HISAT2") also work well if you make sure that HISAT2 doesn’t perform soft-clipping. At the time of writing HISAT2 does perform soft-clipping (CIGAR operation: `S`) by default, so you need to specify the option `--no-softclip`.

#### RNA-Seq alignments with STAR: 

Alignment files produced by the Spliced Transcripts Alignment to a Reference [(STAR)](https://github.com/alexdobin/STAR/ "STAR") aligner also should also work well, however a few steps need to be adhered to make this work:

**1)** Since reStrainingOrder only recognises the CIGAR operations M, I, D and N, alignments need to be run in end-to-end mode and not using local alignments (which may result in soft-clipping). This can be accomplished using the option: `--alignEndsType EndToEnd`

**2)** reStrainingOrder requires the `MD:Z:` field of the BAM alignment to work out mismatches involving masked N positions. Since STAR doesn’t report the `MD:Z:` field by default it needs to be instructed to do so, e.g.: `--outSAMattributes NH HI NM MD`

**3)** To save some time and avoid having to sort the reads by name, STAR can be told to leave R1 and R2 following each other in the BAM file using the option: `--outSAMtype BAM Unsorted`

### Alignments with BWA

The other very popular Burrows-Wheeler Aligner ([BWA](http://bio-bwa.sourceforge.net/ "BWA homepage at SourceForge")) is unfortunately **not compatible** with reStrainingOrder processing as was [disscussed in more detail here](https://github.com/FelixKrueger/SNPsplit/issues/19 "BWA alignments are not compatible with SNPsplit"). Briefly, the reason for this is that BWA randomly replaces the ambiguity nucleobase `N` in the reference by either `C`, `A`, `T` or `G`, thereby rendering an N-masked allele-sorting process impossible.

### Bisulfite-Seq data: 

This mode assumes input data has been processed with the bisulfite mapping tool [Bismark](https://github.com/FelixKrueger/Bismark "Bismark project page on Github"). reStrainingOrder will run a quick check at the start of a run to see if the file provided appears to be a Bismark file, and set the flags `--bisulfite` automatically. Bisulfite (`--bisulfite`) mode can also be set manually.

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


## Step III - Scoring SNPs

    
# Credits
reStrainingOrder was written by Felix Krueger at the [Babraham Bioinformatics Group](http://www.bioinformatics.babraham.ac.uk/).

![Babraham Bioinformatics](Images/bioinformatics_logo.png)
