#!/bin/bash

# scripts to obtain an ASV table from 16S rRNA sequencing
# data from natural community experiments
# by Rui Guan <guan@mpipz.mpg.de>
# modified from the DADA2 tutorial v.1.12

# exits whenever a function returns 1
set -e

# load functions
scripts_dir=`dirname $0`
source $scripts_dir/activate
source $scripts_dir/config.sh

log() {
    echo $(date -u)": "$1 >> $logfile
}

mkdir -p $working_dir

working_dir=$working_dir/01.split_fq

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
rm -f -r $working_dir/*

# demultiplexing reads1 and reads2 files
for l in $l_list_miseq
do
    # initialize lib. results directory
    log "["$l"] initializing the working_dir for miseq step 1..."
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

    # merge reads according to the overlap
    log "["$l"] decompressing forward/reverse/barcode files..."

    gzip -d -c $data_dir/"$l"_forward_reads.fastq.gz > \
        $working_dir/$l/forward_reads.fastq
    gzip -d -c $data_dir/"$l"_reverse_reads.fastq.gz > \
        $working_dir/$l/reverse_reads.fastq
    gzip -d -c $data_dir/"$l"_barcodes.fastq.gz  > \
        $working_dir/$l/barcodes.fastq

    ## for Bacteria
    log "["$l"] demultiplexing Bacteria reads..."
    bc_len=`less $data_dir/$l\_mapping.txt |tail -n1 |
            awk '{print $2}' |wc -c`
    let "bc_len=$bc_len - 1"
    split_libraries_fastq.py -i $working_dir/$l/forward_reads.fastq \
        -b $working_dir/$l/barcodes.fastq \
        -m $data_dir/"$l"_mapping.txt \
        --rev_comp_mapping_barcodes \
        --barcode_type $bc_len \
        --max_barcode_errors 0 \
        --store_demultiplexed_fastq \
        -q 0 -r 300 -p 0.01 -n 300 \
        --phred_offset 33 \
        -o $working_dir/$l/Bac_forward \
        &>> $output

    split_libraries_fastq.py -i $working_dir/$l/reverse_reads.fastq \
        -b $working_dir/$l/barcodes.fastq \
        -m $data_dir/"$l"_mapping.txt \
        --rev_comp_mapping_barcodes \
        --barcode_type $bc_len \
        --max_barcode_errors 0 \
        --store_demultiplexed_fastq \
        -q 0 -r 300 -p 0.01 -n 300 \
        --phred_offset 33 \
        -o $working_dir/$l/Bac_reverse \
        &>> $output

    split_sequence_file_on_sample_ids.py --file_type fastq \
        -i $working_dir/$l/Bac_forward/seqs.fastq \
        -o $working_dir/$l/Bac_forward/out \
        &>> $output

    split_sequence_file_on_sample_ids.py --file_type fastq \
        -i $working_dir/$l/Bac_reverse/seqs.fastq \
        -o $working_dir/$l/Bac_reverse/out \
        &>> $output

    log "["$l"] step 0 finished."
    `rm $working_dir/$l/forward_reads.fastq
        $working_dir/$l/reverse_reads.fastq \
        $working_dir/$l/barcodes.fastq`
done

working_dir_s1=$working_dir
working_dir=$working_dir/02.dada2

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
rm -f -r $working_dir/*

for l in $l_list_miseq
do
    # initialize lib. results directory
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

    mkdir -p $working_dir/$l/Bacteria
    $scripts_dir/dada2_bacteria.R $working_dir_s1/$l/ $working_dir/$l/
done

working_dir_s2=$working_dir
working_dir=$working_dir/03.ASV_tab

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
rm -f -r $working_dir/*

mkdir -p $working_dir/Bacteria
$scripts_dir/get_ASV_tab_bacteria.R $l_list_miseq \
        $working_dir_s2/ $working_dir/Bacteria
`less $working_dir/Bacteria/ASV_map.txt |
    sed 's/ASV/>ASV/; s/\t/\n/' >ASV.fasta`

