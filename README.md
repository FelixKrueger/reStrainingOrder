[![Build Status](https://travis-ci.org/FelixKrueger/reStrainingOrder.svg?branch=master)](https://travis-ci.org/FelixKrueger/reStrainingOrder)

<img title="Odd One Out" align="right" id="header_img" src="Docs/Images/mice_logo.png">

# reStrainingOrder - why do we need one?
reStrainingOrder is intended as QC tool that attempts to identify the genotype of pure strain or hybrid mouse samples. It can be be used to check public data as well as provide useful insight into mouse strains commonly used in your own lab.

To do this, reStrainingOrder harnesses single-nucleotide polymorphism (SNP) information collected by the [Mouse Genomes Project](https://www.mousegenomes.org/), and constructs a fully N-masked genome similar to the approach of [SNPsplit](http://felixkrueger.github.io/SNPsplit/). The project has been updated to work with the latest release of the Mouse Genomes Project, which means that it will now assume the GRCm39 mouse genome build by default, and use the latest SNP annotation file (v8: [mgp_REL2021_snps.vcf.gz](https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz)).

reStrainingOrder is intended to work with most common types of Illumina sequencing - including `RNA-seq`, `ChIP-seq`, `ATAC-seq` or any kind of `Bisulfite-seq`. Supported aligners include [`Bowtie2`](https://github.com/BenLangmead/bowtie2), [`HISAT2`](https://ccb.jhu.edu/software/hisat2/index.shtml), [`STAR`](https://github.com/alexdobin/STAR), and [`Bismark`](https://github.com/FelixKrueger/Bismark) (Oxford comma, anyone?).

#### Pure Strain example report

Here is an example of a pure [C57BL/6 strain (Black6)](https://www.bioinformatics.babraham.ac.uk/projects/reStrainingOrder/pure_strain_example.html) strain. If the single-strain compatibility score is > 99%, one doesn't really need to look any further. The allele-ratio also indicated that the sample is purely Black6.

#### Hybrid example report

Here is an example of a [129S1/CAST hybrid](https://www.bioinformatics.babraham.ac.uk/projects/reStrainingOrder/129_CAST_hybrid_example.html) strain (public data taken from this [GEO entry](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM753570), aligned with Bismark/Bowtie2). The Hybrid Strain Compatibility Scores indicate that a 129S1_SvImJ/CAST_EiJ hybrid is > 99% compatible with the data. The allele-ratios between 129 and CAST are almost 1:1.


## Installation

Most of reStrainingOrder is written in Perl and is executed from the command line. To install reStrainingOrder simply download the latest release of the code from the [Releases page](https://github.com/FelixKrueger/reStrainingOrder/releases) and extract the files into a reStrainingOrder installation folder. The plotly.js JavaScript library comes packaged together with reStrainingOrder and does not need to be installed separately.

reStrainingOrder requires the following tools installed and ideally available in the `PATH`:
- [Samtools](http://samtools.sourceforge.net/)


## Documentation
The reStrainingOrder documentation can be found here: [reStrainingOrder User Guide](./Docs/README.md)

## Licences

reStrainingOrder itself is free software, `reStrainingReport` produces HTML graphs powered by [Plot.ly](https://plot.ly/javascript/) which are also free to use and look at!

## Credits
This project was started as part of the 2018 Cambridge area bioinformatics hackathon; reStrainingOrder was written by Felix Krueger, part of [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk), now part of Altos Bioinformatics.
<p align="center"> <img title="Babraham Bioinformatics" id="logo_img" src="Docs/Images/bioinformatics_logo.png"></p>
