#!/bin/bash

# config file for 16S data analysis
# guan@mpipz.mpg.de, modified from garridoo@mpipz.mpg.de

# set python path
export PYTHONPATH=/biodata/dep_psl/grp_rgo/tools/lib/python2.7/site-packages/:$PYTHONPATH

# working_dir is the folder where (intermediate) results are stored
working_dir="/biodata/dep_psl/grp_rgo/guan/OTU_At_Lj_Roche_MiSeq/02.results/"

# data_dir is the folder where the raw data (and maping file) are deposited
data_dir="/biodata/dep_psl/grp_rgo/guan/OTU_At_Lj_Roche_MiSeq/"

# paths to relevant reference database files
refdata_dir="/biodata/dep_psl/grp_rgo/tools/16s/ref_data"

# path to folder containing usearch binary and python scripts
usearch_dir="/biodata/dep_psl/grp_rgo/tools/usearch/"


gold_db=$refdata_dir"/cs_gold.fa"
gg_core_aligned_db=$refdata_dir"/gg_13_8_otus/rep_set_aligned/97_otus.fasta"

silva_core_db="/biodata/dep_psl/grp_rgo/guan/meta_18s/00.data/database/SILVA_132_QIIME_release/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna"
silva_taxonomy="/biodata/dep_psl/grp_rgo/guan/meta_18s/00.data/database/SILVA_132_QIIME_release/taxonomy/16S_only/97/taxonomy_7_levels.txt"

### parameters
# number of threads for parallel steps
n_cores=30
# max. barcode errors for 454 data
bc_err=0
# max. barcode errors for miseq data
bc_err_miseq=2
# min. qual score for 454 data
qual=25
# phred quality threshold for miseq data
q_miseq=30
# min. number of reads for a given OTU
min_size=2
# identity threshold for OTU clustering
id_threshold=0.97
# min length of alingment required to sequences in the database
min_id_aln=0.75
# subsampling depth used for alpha-diversity
subsampling_depth=1000
# OTU clustering method
tax_method="uclust"
