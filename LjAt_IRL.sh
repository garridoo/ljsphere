#!/bin/bash
# scripts for 16S data filter and OTU clustering for 454 data,
# truncate reads by the alignment of primer sequences
# Author guan@mpipz.mpg.de, modified from garridoo@mpipz.mpg.de
# exits whenever a function returns 1
set -e

# get path to scripts
scripts_dir=$(dirname $0)

# activate and load config file.
source $scripts_dir/activate.sh
source $scripts_dir/config.sh

log(){
    echo $(date -u)": "$1 >> $logfile
}

# prepare results director
prefix=`date +%d%m%Y`
working_dir=$working_dir/$prefix\_s1

output=$working_dir/"output.txt"
logfile=$working_dir/"log.txt"

mkdir -p $working_dir
rm -f -r $working_dir/*

# Get the list of 454 < run > used in step 1 and step 2
l_list=`ls -l $data_dir/s1/*.fasta |awk '{print $9}' |sed 's/\.fasta//'
    |sed 's/\//\t/g' |awk '{print $NF}' | xargs `

# Get reads with both forward/reverse primer sequences and truncate reads
# according to reverse primer position
for l in $l_list
do
    # initialize lib. results directory
    log "["$l"] initializing the workdir for 454 step 1..."
    rm -f -r $working_dir/"$l"
    mkdir $working_dir/"$l"

    # truncating reads to equal length
    log "["$l"] truncating reads..."
    perl $scripts_dir/truncate_by_length.pl -nr $scripts_dir/forward_primer.fasta \
                                                $scripts_dir/reverse_primer.fasta \
                                                $data_dir/s1/"$l".fasta \
                                                $data_dir/s1/"$l".qual \
                                                $working_dir/"$l"/filtered.fasta \
                                                $working_dir/"$l"/filtered.qual \
                                                &>> $output
done

log "Step1 454 data truncate DONE!"

## Step 2, demultiplexing reads according to barcodes
working_dir_s1=$working_dir/$prefix\_s1
working_dir_s2=$working_dir/$prefix\_s2

mkdir -p $working_dir_s2
rm -f -r $working_dir_s2/*

output=$working_dir_s2/"output.txt"
logfile=$working_dir_s2/"log.txt"

# Get the list of 454 < run > information.
l_list=`ls -l $data_dir/s1/*.fasta |awk '{print $9}'|sed 's/\.fasta//' |
    sed 's/\//\t/g' |awk '{print $NF}' |xargs`

# Demultiplexing the data according to barcodes
for l in $l_list
do
    # initialize results directory
    log "["$l"] initializing the 454 workdir for step 2..."
    rm -f -r $working_dir_s2/$l
    mkdir $working_dir_s2/$l

    # get the barcode length
    bc_length=`less $data_dir/s1/$l\_mapping.txt |tail -n1 |awk '{print $2}'|
               wc -c`
    let "bc_length=$bc_length-1"

    # quality filtering and demultiplexing
    log "["$l"] demultiplexing..."
    split_libraries.py -f $working_dir_s1/"$l"/filtered.fasta \
        -q $working_dir_s1/"$l"/filtered.qual \
        -m $data_dir/s1/"$l"_mapping.txt \
        -s $qual \
        -e $bc_err \
        -b $bc_length \
        -d \
        -o $working_dir_s2/"$l"/demux \
        &>>$output

    mv $working_dir_s2/"$l"/demux/* $working_dir_s2/"$l"
    rm -f -r $working_dir_s2/"$l"/demux

    # edit barcode label identifier for usearch compatibility
    cat $working_dir_s2/"$l"/seqs.fna | \
        sed 's/.*/;/g;s/>.*/&&/g;s/;>/;barcodelabel=/g;s/_[0-9]*;$/;/g' \
        >> $working_dir_s2/seqs.fasta

    log "["$l"] finish demultiplexing this  run"
done

log "Step2 454 data demultiplex DONE!"

# Step 3 OTU clustering
working_dir_s3="../03.OTU_cluster/"

output=$working_dir_s3/"output.txt"
logfile=$working_dir_s3/"log.txt"

# dereplication
log "dereplicating..."
usearch -derep_fulllength $working_dir_s3/seqs.fasta \
        -fastaout $working_dir_s3/seqs_unique.fasta \
        -sizeout \
        &>> $output

# abundance sort and discard singletons
log "sorting by abundance and discarding singletons..."
usearch -sortbysize $working_dir_s3/seqs_unique.fasta \
        -fastaout $working_dir_s3/seqs_unique_sorted.fasta \
        -minsize $min_size \
        &>> $output

# OTU clustering
log "OTU clustering using UPARSE..."
usearch -cluster_otus $working_dir_s3/seqs_unique_sorted.fasta \
        -otus $working_dir_s3/otus.fasta \
        &>> $output

# chimera detection
log "removing chimeras..."
usearch -uchime_ref $working_dir_s3/otus.fasta \
        -db $gold_db \
        -strand plus \
        -nonchimeras $working_dir_s3/otus_nc.fasta \
        -threads $n_cores \
        &>> $output

# align sequences to database using PyNAST and remove remaining
log "aligning OTU representative sequences to database..."
align_seqs.py -i $working_dir_s3/otus_nc.fasta \
              -t $gg_core_aligned_db \
              -p $min_id_aln \
              -o $working_dir_s3 \
              &>> $output

# rename OTUs and remove alignment gaps
log "renaming OTUs..."

sed -i 's/-//g' $working_dir_s3/otus_nc_aligned.fasta &>> $output

cat $working_dir_s3/otus_nc_aligned.fasta | \
    awk 'BEGIN {n=1}; />/ {print ">OTU_" n; n++} !/>/ {print}' \
    >> $working_dir_s3/rep_seqs.fasta

# generate OTU table
log "generating OTU table..."
usearch -usearch_global $working_dir_s3/seqs.fasta \
        -db $working_dir_s3/rep_seqs.fasta \
        -strand plus \
        -id $id_threshold \
        -uc $working_dir_s3/read_mapping.uc \
        &>> $output

# convert uc file to txt
log "converting uc OTU table file into text format..."
python $usearch_dir/uc2otutab.py $working_dir_s3/read_mapping.uc \
    1> $working_dir_s3/otu_table.txt \
    2>> $output

# taxonomy assignment
log "taxonomy assignment..."
assign_taxonomy.py -i $working_dir_s3/rep_seqs.fasta \
                   -r $silva_core_db \
                   -t $silva_taxonomy \
                   -m $tax_method \
                   -o $working_dir_s3/tax \
                   &>> $output

# cleanup
mv $working_dir_s3/tax/rep_seqs_tax_assignments.txt $working_dir_s3/taxonomy.txt
rm -f -r $working_dir_s3/tax

sed -i 's/; /\t/g' $working_dir_s3/taxonomy.txt

# convert OTU table to biom
log "converting OTU table to QIIME compatible biom format..."
biom convert -i $working_dir_s3/otu_table.txt \
             -o $working_dir_s3/otu_table.biom \
             --table-type="OTU table" \
             --to-json \
             &>> $output

log "Step3 OTU cluster DONE!"
