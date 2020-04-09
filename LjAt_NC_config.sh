#!/bin/bash

# scripts to obtain an ASV table from 16S rRNA sequencing
# data from natural community experiments
# by Rui Guan <guan@mpipz.mpg.de>

working_dir="./"
data_dir="/biodata/dep_psl/grp_rgo/guan/ASV_At_Lj/00.data/"

## define the profiled kingdom
bacteria=true

# Get the list of miseq < run >
l_list_miseq=`ls -l $data_dir/*_forward_reads.fastq.gz | \
    awk '{print $9}' |sed 's/_forward_reads\.fastq\.gz//'| \
    sed 's/\//\t/g' |awk '{print $NF}' | xargs `

