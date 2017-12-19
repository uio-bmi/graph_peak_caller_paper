#!/usr/bin/env bash
set -e

config_file=$1
experiment_id=$2
replicate_number=$3
tf_name=$4
base_dir=$(pwd)
chromosomes=$5
bam_alignments_url=$6

# Example
# ./pipeline.sh ENCSR000DUB 1 CTCF 135 36 /home/ivargry/dev/graph_peak_caller/tests/lrc_kir/ /home/ivargry/dev/graph_peak_caller/graph_peak_caller.py 1,2
# On server:
# ./pipeline.sh config_server.sh ENCSR000DUB 1 CTCF 16,17 https://www.encodeproject.org/files/ENCFF639IFG/@@download/ENCFF639IFG.bam

source $config_file
echo "Config file: $config_file"
echo "Experiment id: $experiment_id"

#vg_xg_index="$graph_dir/graph.xg"
#vg_gcsa_index="$graph_dir/graph.gcsa"
vg_xg_index=/home/ivar/data/whole_genome/graph.xg
vg_gcsa_index=/home/ivar/data/whole_genome/graph.gcsa

work_dir="data/${tf_name}_${experiment_id}/$replicate_number"
mkdir -p $work_dir
cd $work_dir
echo "Made working directory $work_dir"


# Step 1: Download data

if [ ! -f raw.fastq ]; then
    echo "Download fastq"
    encode_url=$(python3 $base_dir/download_encode_fastq.py $experiment_id $replicate_number)
    echo "Encode url: $encode_url"
    wget -O raw.fastq.gz --show-progress $encode_url 
    echo "Unzipping"
    gunzip -c raw.fastq.gz > raw.fastq
else
    echo "Raw fastq already exists. Not dowloading"
fi

# Downlaod bam alignments
if [ ! -f linear_alignments.bam ]; then
    echo "Downloading linear alignments (for macs2)"
    wget -O linear_alignments.bam $bam_alignments_url 
fi


# Run macs2 to get fragment length/read length
if [ ! -f macs_output.txt ]; then
	echo "Running macs2"
	macs2 predictd -g hs -i linear_alignments.bam > macs_output.txt 2>&1
else
	echo "Not running max. Already done"
fi

read_length=$(cat macs_output.txt | gawk 'match($0,  /tag size = ([0-9]+)/, ary) {print ary[1]}' )
echo "Found read length: $read_length"
fragment_length=$(cat macs_output.txt | gawk 'match($0,  /fragment length is ([0-9]+)/, ary) {print ary[1]}' )
echo "Found fragment length: $fragment_length"

# Step 2: Filter reads
# fastqc, trim_galore
if [ ! -f filtered.fastq ]; then
   head -n 10000 raw.fastq > filtered.fastq
    echo "Creating filtered.fastq from raw.fastq" 
    #mv raw.fastq filtered.fastq
else
    echo "Fastq already filtered"
fi


# Step 3: Map reads
echo "Mapping reads"
if [ ! -f mapped.gam ]; then
    echo "Using indices: $vg_gcsa_index and $vg_xg_index"
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
        filtered_$chromosome.json \
        filtered_$chromosome.json \
        False \
        "chr${chromosome}_" \
        $fragment_length \
        $read_length > log_chr$chromosome.txt 2>&1 &

    echo "Log output for chr $chromosome will be written to $work_dir/log_chr$chromosome.txt"

done

echo "Peak calling now running in background."
