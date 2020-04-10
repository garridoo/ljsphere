## Scripts of Wippel, Tao, *et. al*, Host adaptation and invasiveness of commensals in the *Lotus* and *Arabidopsis* root microbiota.

correspondence to Ruben Garrido-Oter

garridoo@mpipz.mpg.de

These scripts are made available to facilitate the reproducibility of our research. If you re-use any or part of this code, please reference with comments and cite our paper. Raw data and intermediate results necessary to run these scripts can also be downloaded [here](http://www.mpipz.mpg.de/R_scripts).

---------------------------

### Accession numbers

Raw *16S* rRNA amplicon reads were deposited in the European Nucleotide Archive (ENA) under the accession number [PRJEB37695](XXX). Similarly, sequencing reads and genome assemblies of the *Lj*-SPHERE core collection were uploaded to the same database with the accession number [PRJEB37696](XXX).

### Scripts used for processing data, creating the figures and performing the statistical analysis reported in the manuscript.

#### Analysis of culture independent and IRL *16S* rRNA amplicon data

Raw data corresponding to the natural community greenhouse experiments with *Lotus* and *Arabidopsis* data (Zgadzaj *et al.*, 2019; Harbort *et al.*, 2020) along with intermediate results (ASV table, metadata, etc.) can be found [here](http://www.at-sphere.com/ljsphere/LjAt_NC.tar.gz).

[LjAt_NC_config.sh](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC_config.sh) and
[LjAt_NC.sh](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC.sh):
bash pipeline to pre-process raw sequencing amplicon data from the greenhouse experiments and obtain an ASV table.

[LjAt_NC.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC.R),
[LjAt_NC_load_data.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC_load_data.R),
[LjAt_NC_diversity.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC_diversity.R), and
[LjAt_NC_rank.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_NC_rank.R)
R script to process the natural community ASV table and perform downstream analysis of alpha- and beta-diversity as well as rank-abundance plots and enrichment tests (Fig. 1a-b; Fig. S1).

Raw sequencing data and intermediate results from the *At*-IRL and *Lj*-IRL as well as corresponding natural community inocula can be downloaded [here](http://www.at-sphere.com/ljsphere/LjAt_IRL.tar.gz).

[LjAt_IRL_config.sh](https://github.com/garridoo/ljsphere/blob/master/LjAt_IRL_config.sh) and
[LjAt_IRL.sh](https://github.com/garridoo/ljsphere/blob/master/LjAt_IRL.sh)
bash pipeline to pre-process raw sequencing amplicon data and obtain a combined OTU table from all IRL data.

[IRLs.R](https://github.com/garridoo/ljsphere/blob/master/IRLs.R),
[IRLs_load_data.R](https://github.com/garridoo/ljsphere/blob/master/IRLs_load_data.R), and
[IRLs_recovery_rates.R](https://github.com/garridoo/ljsphere/blob/master/IRLs_recovery_rates.R):
R scripts used to cross-reference data from shallow sequencing of IRL plates (CFUs) and culture-independent profiling of start inocula, and to estimate species recovery rates (Fig. 1c-f).

[plotting_parameters.R](https://github.com/garridoo/ljsphere/blob/master/plotting_parameters.R), and
[plotting_functions.R](https://github.com/garridoo/ljsphere/blob/master/plotting_functions.R):
R scripts containing plotting parameters such as colors, ggplot2 themes, etc. as well as auxiliary functions for generating the figures reported in the paper (some of the details may vary with respect to the published version).

#### Whole-genome assembly, quality control, and annotation of the *Lj*-SPHERE core culture collection

Raw sequencing data (FASTQ files) from the core *Lj*-SPHERE culture collection as well as assemblies (FNA files), nucleotide ORFs (FFN), amino acid ORFs (FAA), GFF files, KEGG annotations (KO), reference *16S* rRNA sequences, AMPHORA marker gene alignments, and metadata (used for Fig. 5) can be downloaded in bulk [here](http://www.at-sphere.com/ljsphere/LjAt_WGS.tar.gz).

[assembly.functions.sh](https://github.com/garridoo/ljsphere/blob/master/assembly.functions.sh): bash script containing auxiliary functions for whole-genome assembly.

[assembly.sh](https://github.com/garridoo/ljsphere/blob/master/assembly.sh): script used to assemble the genomes using SOAPdenovo and A5. It can be run in parallel using either the custom script bellow or the gnu parallel suit.

[parallel.sh](https://github.com/garridoo/lsphere/blob/master/parallel.sh): custom script to run bash functions in parallel in a multi-core machine.

[assembly_stats.R](https://github.com/garridoo/ljsphere/blob/master/assembly_stats.R): R script used to generate assembly statistics as well as GC and k-mer spectral projections. The output of this script contains clean assemblies (all contigs <1,000 bp are removed) as well as a PDF file containing a report which was used to manually inspect for likely contaminated assemblies.

[functional_MDS.R](https://github.com/garridoo/ljsphere/blob/master/functional_MDS.R): R script to perform dimensionality reduction of genome functional profiles on the genomes from the *Lj*- and *At*-SPHERE collections.

#### Reference-based analysis of *16S* rRNA amplicon data from SynCom experiments

Raw sequencing data from all SynCom experiments (FASTQ files and mapping files including barcodes, biological and technical metadata, etc.) along with intermediate results (ASV tables, etc.) can be downloaded [here](http://www.at-sphere.com/ljphere/LjAt_SC.tar.gz).

The experiments reported in the paper correspond to independent Illumina MiSeq sequencing runs as follows:

| Sequencing run | Experiment ID            |
| -------------- |--------------------------|
| AtLj_001        | SynCom experiment A      |
| AtLj_002        | SynCom experiment B      |
| AtLj_003        | SynCom experiment D      |
| AtLj_004        | SynCom experiment F      |
| AtLj_005        | SynCom experiment E      |
| AtLj_006        | SynCom experiment C      |
| MDA10           | Millifluidics experiment (G) |

[LjAt_SC_load_data.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_load_data.R): R script used to load data and metadata from the SynCom competition experiments.

[LjAt_SC_diversity.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_diversity.R) and
[cpcoa.func.R](https://github.com/garridoo/ljsphere/blob/master/cpcoa.func.R): R scripts employed for beta-diversity analysis, dimensionality reduction, and permutation analyses of variance (Figs. 3B-D, and 4B).

[LjAt_SC_host_preference.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_host_preference.R): script used for plotting host preference boxplots from competition experiments (Fig. 3E-G).

[LjAt_SC_invasion_load_data.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_invasion_load_data.R): R script used to load data and metadata from the invasion and persistence experiment.

[LjAt_SC_invasion.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_invasion.R): R script used for plotting host preference boxplots from the invasion and persistence experiment (Fig. 4D-E).

[hpi.R](https://github.com/garridoo/ljsphere/blob/master/hpi.R), and
[hpi_invasion.R](https://github.com/garridoo/ljsphere/blob/master/hpi_invasion.R): R scripts used to calculate host preference and invasiveness indices (Fig. 5).

[LjAt_MDA_diversity.R](https://github.com/garridoo/ljsphere/blob/master/LjAt_SC_diversity.R): R script used for beta-diversity analysis, dimensionality reduction, and permutation analyses of variance of *in vitro* millifluidics experiments (Fig. S4).

---------------------------

For any questions regarding these scripts, please contact

Ruben Garrido-Oter

garridoo@mpipz.mpg.de
