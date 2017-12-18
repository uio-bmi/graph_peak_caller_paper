#!/usr/bin/env bash

experiment_id=$1
replicate_number=$2
tf_name=$3
fragment_length=$4
read_length=$5
graph_dir=$6
graph_peak_caller=$7
base_dir=$(pwd)
chromosomes=$8
# Example
# ./pipeline.sh ENCSR000DUB 1 CTCF 135 36 /home/ivargry/dev/graph_peak_caller/tests/lrc_kir/ /home/ivargry/dev/graph_peak_caller/graph_peak_caller.py 1,2


vg_xg_index="$graph_dir/graph.xg"
vg_gcsa_index="$graph_dir/graph.gcsa"

work_dir="data/${tf_name}_${experiment_id}/$replicate_number"
mkdir -p $work_dir
cd $work_dir


# Step 1: Download data

if [ ! -f raw.fastq ]; then
    echo "Download fastq"
    encode_url=$(python3 $base_dir/download_encode_fastq.py $experiment_id $replicate_number)
    echo "Encode url: $encode_url"
    wget -qO- $encode_url > raw.fastq.gz
    gunzip -c raw.fastq.gz > raw.fastq
else
    echo "Raw fastq already exists. Not dowloading"
fi

# Step 2: Filter reads
# fastqc, trim_galore
if [ ! -f filtered.fastq ]; then
    head -n 1000 raw.fastq > filtered.fastq
else
    echo "Fastq already filtered"
fi


# Step 3: Map reads
if [ ! -f mapped.gam ]; then
    vg map -f filtered.fastq -g $vg_gcsa_index -x $vg_xg_index -M 2 > mapped.gam
else
    echo "Mapped reads exist. Not mapping"
fi

# Step 4: Filter mapped reads
vg filter -r 1.0 -s 2.0 -fu mapped.gam > filtered.gam
vg view -aj filtered.gam > filtered.json

# Step 5: Split filtered into chromosomes
python3 $graph_peak_caller split_vg_json_reads_into_chromosomes filtered.json $graph_dir

# Step 6 run peak caller on chromosomes
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    python3 $graph_peak_caller callpeaks \
        $graph_dir/$chromosome.json \
        $graph_dir/$chromosome.vg \
        $graph_dir/linear_map_$chromosome \
        filtered.json \
        filtered.json \
        False \
        "" \
        $fragment_length \
        $read_length > log_chr$chromosome.txt 2>&1
done