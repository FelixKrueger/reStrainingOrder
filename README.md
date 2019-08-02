<img title="Odd One Out" align="right" id="header_img" src="Docs/Images/mice_logo.png">

# Why do we need one?
reStrainingOrder is intended as QC tool that attempts to predict the genotype of pure strain or hybrid mouse samples. It can be be used to check public data as well as provide useful insight into mouse strains commonly used in your own lab.

To do this, reStrainingOrder harnesses single-nucleotide polymorphism (SNP) information collected by the Mouse Genomes Project (MGP, http://www.sanger.ac.uk/science/data/mouse-genomes-project).


## Installation

Most of reStrainingOrder is written in Perl and is executed from the command line. To install reStrainingOrder simply download the latest release of the code from the [Releases page](https://github.com/FelixKrueger/SNPsplit/releases) and extract the files into a SNPsplit installation folder. The plotly.js JavaScript library comes packaged together with reStrainingOrder and does not need to be installed separately.

reStrainingOrder requires the following tools installed and ideally available in the `PATH`:
- [Samtools](http://samtools.sourceforge.net/)


## Documentation
The reStrainingOrder documentation can be found here: [reStrainingOrder User Guide](./reStrainingOrder_User_Guide.md)


## Links

This project was started as part of the 2018 Cambridge area bioinformatics hackathon (https://www.cambiohack.uk/).

## Licences

reStrainingOrder itself is free software, `reStrainingReport` produces HTML graphs powered by [Plot.ly](https://plot.ly/javascript/) which are also free to use and look at!

## Credits

reStrainingOrder was written by Felix Krueger, part of the [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk) group.

<p align="center"> <img title="Babraham Bioinformatics" id="logo_img" src="Images/bioinformatics_logo.png" width=300></p>
