# reStrainingOrder Changelog

## Changelog for version 0.5.0 [unreleased]

### reStraining

- Chromosomes are now processed in the order of the VCF header, so the genome preparation report is the same on every run. Addressed in [#16](https://github.com/FelixKrueger/reStrainingOrder/pull/16).
- Chromosome 1 is now filtered under its UCSC name (`chr1`) too. Addressed in [#18](https://github.com/FelixKrueger/reStrainingOrder/pull/18).

### reStrainingOrder

- Fixed a crash at the end of every run with `--output_dir`, and the HTML report now covers only the current run. Addressed in [#13](https://github.com/FelixKrueger/reStrainingOrder/pull/13).
- Inputs that are not `.sam` or `.bam` are now refused instead of being overwritten by their own report. Addressed in [#14](https://github.com/FelixKrueger/reStrainingOrder/pull/14).
- The strain table is now sorted numerically, with ties in strain order. Addressed in [#15](https://github.com/FelixKrueger/reStrainingOrder/pull/15).
- The word `sorted.` is now dropped from report names only where it starts a word. Addressed in [#17](https://github.com/FelixKrueger/reStrainingOrder/pull/17).
- With several input files, each file is now scored on its own reads only. Addressed in [#19](https://github.com/FelixKrueger/reStrainingOrder/pull/19).
- Bismark files are now detected from any `@PG` line, not only the last one. Addressed in [#20](https://github.com/FelixKrueger/reStrainingOrder/pull/20).
- One bisulfite branch now counts an unexpected base as OTHER, like the others. Addressed in [#21](https://github.com/FelixKrueger/reStrainingOrder/pull/21).
- Reads with soft clips, splices or indels are now scored, so fewer reads are reported as not containing one of the expected bases. Addressed in [#22](https://github.com/FelixKrueger/reStrainingOrder/pull/22).
- The strain tallies and hybrid detection are much faster and use far less memory, and a malformed row in the SNP file now stops the run. Addressed in [#23](https://github.com/FelixKrueger/reStrainingOrder/pull/23).
- Blank lines in the SNP file are now skipped, and errors in it name the right line.
- Secondary and supplementary alignments are now skipped and counted in the report.
- Reads at a SNP position that have no sequence or CIGAR, or whose sequence, CIGAR and MD tag disagree, are now skipped and counted in the report.

## Changelog for version 0.4.0 [release on 12 April 2023]

Updated documentation to reflect the changes made by swtiching over to the `GRCm39` mouse genome, as well as the v8 annotation files. The v5 and v7 versions are probably no longer available for download, but we have left the option `--v7` in there for backward compatibility for the time being.

### reStraining

- Updated the genome preparation to now work with the latest (v8) genome annotation file (mgp_REL2021_snps.vcf.gz) from the [Mouse Genomes Project](https://www.mousegenomes.org/). The mgp_v5 version for the now outdated GRCm38 genome are now no longer supported (since this is primarily a screening tool anyway...).


## Changelog for version 0.3.0 [release on 27 March 2022]

### reStraining

- Added new option `--v7_VCF` to allow using the file 'mgp_REL2005_snps_indels.vcf.gz' instead of the mgp.v5 file mentioned in the main documentation. This file contains both SNP and INDEL information, but INDELs are skipped. The mgp.v5 file is still listed as Current_SNPs, so the v7 file is still considered experimental. The v7 file contains a number of additional strains, and was released in May2020. Get it here: ftp://ftp-mouse.sanger.ac.uk/REL-2004-v7-SNPs_Indels/mgp_REL2005_snps_indels.vcf.gz 

### reStrainingOrder

- Fixed an issue where allele-ratios were always assiged to the first strain in pairwise comparisons. Addressed in [#2](https://github.com/FelixKrueger/reStrainingOrder/issues/2). 


## Changelog for version 0.2.0 [release on 12 July 2021]

- Fixed the allele-ratio plot/table in the HTML report (it was picking these ratios in a rather random manner). Addressed in [#1](https://github.com/FelixKrueger/reStrainingOrder/issues/1). 

## Changelog for reStrainingOrder v0.1.0. [release on 15 August 2019]

- Initial release. All basic tools seem to be working. Documentation and help need a little more work, but this will have to wait until after the holidays.

- Travis CI also needs some TLC (read: there are currently no use useful checks being carried out).
