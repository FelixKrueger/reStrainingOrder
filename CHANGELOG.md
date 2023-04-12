# reStrainingOrder Changelog

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
