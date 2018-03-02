#!/usr/bin/env bash
set -e

config_file=$1
experiment_id=$2
replicate_number=$3
tf_name=$4
base_dir=$(pwd)
chromosomes=$5
bam_alignments_url=$6
motif_url=$7
genome_size=$8


grch38_fasta_file=~/data/hg19/chromosomes/hg19_chr1-Y.fa
n_threads=$(grep -c ^processor /proc/cpuinfo)

source $config_file
echo "Config file: $config_file"
echo "Experiment id: $experiment_id"
echo "Hg 19 fasta file: $grch38_fasta_file"

vg_xg_index=/home/ivar/data/whole_genome/wg1.6.xg
vg_gcsa_index=/home/ivar/data/whole_genome/wg1.6.gcsa

work_dir="data/${tf_name}_${experiment_id}/$replicate_number"
mkdir -p $work_dir
cd $work_dir
echo "Made working directory $work_dir"


# Step 1: Download data

if [ ! -f raw.fastq.gz ]; then
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
if [ ! -f linear_alignments_from_encode.bam ]; then
    echo "Downloading linear alignments (for macs2)"
    wget -O linear_alignments_from_encode.bam $bam_alignments_url
fi

# Run macs2 to get fragment length/read length
if [ ! -f macs_output_predictd.txt ]; then
	echo "Running macs2"
	macs2 predictd -g hs -i linear_alignments_from_encode.bam > macs_output_predictd.txt 2>&1
else
	echo "Not running macs. Already done"
fi

read_length=$(cat macs_output_predictd.txt | gawk 'match($0,  /tag size = ([0-9]+)/, ary) {print ary[1]}' )
echo "Found read length: $read_length"
fragment_length=$(cat macs_output_predictd.txt | gawk 'match($0,  /fragment length is ([0-9]+)/, ary) {print ary[1]}' )
echo "Found fragment length: $fragment_length"



# Step 2: Filter reads
# fastqc, trim_galore
if [ ! -f raw_trimmed.fq ]; then
    echo "Trimming fastq using trim galore"
   #mv raw.fastq filtered.fastq
   trim_galore raw.fastq
else
    echo "Fastq already trimmed"
fi


# Step 3: Map reads
echo "Mapping reads"
if [ ! -f mapped.gam ]; then
    echo "Using indices: $vg_gcsa_index and $vg_xg_index"
    vg_mappering map -c 5000 -f raw_trimmed.fq -g $vg_gcsa_index -x $vg_xg_index > mapped.gam
    #vg mpmap --mq-method 2 -S -x $vg_xg_index -g $vg_gcsa_index -f raw_trimmed.fq > mapped.gam
else
    echo "Mapped reads exist. Not mapping"
fi

# Step 4: Filter mapped reads
echo "Filtering"
if [ ! -f filtered_low_qual_reads_removed.json ]; then
	vg_mappering filter -q 60 -r 0.95 -s 2.0 -fu -t 20 mapped.gam > filtered.gam
	vg_mappering view -aj filtered.gam > filtered_low_qual_reads_removed.json
else
	echo "Filtered exists. Not filtering"
fi

# Map linear reads (in order to remove bad linear mappings)
# Prepare linear reads for linear peak calling
if [ ! -f linear_alignments.bam ]; then

    echo "Mapping reads to linear genome"
    bwa aln -t $n_threads $grch38_fasta_file raw_trimmed.fq > reads.sai
    bwa samse $grch38_fasta_file reads.sai raw_trimmed.fq > alignments.sam

    # Convert to bam and sort
    echo "Converting to bam and filtering"
    samtools view -Su alignments.sam | samtools sort - alignments_sorted

    # Filter (removed duplicates and reads having low score)
    samtools view -F 1804 -q 37 -b alignments_sorted.bam > linear_alignments.bam
fi

# Get sequence id of reads that mapped with low mapq to linear genome
#echo "Finding reads that mapped bad to linear"
#awk '$5 < 37 { print $1  }' alignments.sam > low_qual.txt
#awk '$5 < 37 && $6 == "36M" { print $1  }' alignments.sam > low_qual.txt
#if [ ! -s "low_qual.txt" ]
#then 
   #echo "Something is probaly wrong. Found no low qual reads."
   #exit 0
#fi

echo "Removing low quality reads."
#python3 $base_dir/filter_json_alignments.py low_qual.txt filtered.json filtered_low_qual_reads_removed.json
#cp filtered.json filtered_low_qual_reads_removed.json



# Step 5: Split filtered into chromosomes
#if [ ! -f filtered_1.json ]; then
	graph_peak_caller split_vg_json_reads_into_chromosomes $chromosomes filtered_low_qual_reads_removed.json $graph_dir
#else
	#echo "Not splitting into chromosomes."
#fi

# Count unique reads in filtered files
#if [ ! -f count_unique_reads_output.txt ]; then
#    graph_peak_caller count_unique_reads $chromosomes $graph_dir/ filtered_low_qual_reads_removed_ > count_unique_reads_output.txt 2>&1
#else
#    echo "Unique reads already counted. Not counting"
#fi

#unique_reads=$(tail -n 1 count_unique_reads_output.txt)
echo "Counting unique reads"
unique_reads=$(pcregrep -o1 '"sequence": "([ACGTNacgtn]{20,})"' filtered_low_qual_reads_removed.json | sort | uniq | wc -l)
echo "$unique_reads unique reads in total"


# Step 6 run peak caller to get p-values for every chromosome
pids=""
RESULT=0
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    #if [ ! -f ${chromosome}_pvalues_values.npy ]; then
	graph_peak_caller callpeaks_whole_genome $chromosome \
		$graph_dir/ \
		$graph_dir/ \
		$graph_dir/linear_map_ \
		filtered_low_qual_reads_removed_ filtered_low_qual_reads_removed_ "" False $fragment_length $read_length \
		True $unique_reads $genome_size \
		> log_before_p_values_$chromosome.txt 2>&1 &
        pids="$pids $!"
	    echo "Peak calling for chr $chromosome started as process. Log will be written to $work_dir/log_before_p_values_$chromosome.txt"
    #else
        #echo "P values already computed for chromosome $chromosome."
    #fi
done

# Wait for all to finish between continuing
for pid in $pids; do
    wait $pid || let "RESULT=1"
done

if [ "$RESULT" == "1" ];
    then
       exit 1
fi


# Step 7 run from p values
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    #if [ -f ${chromosome}_max_paths.intervalcollection ]; then
    #	echo "Peaks already called for $chromosome. Not calling"
    #elif [ -f ${chromosome}_pvalues_values.npy ]; then
	graph_peak_caller callpeaks_whole_genome_from_p_values $chromosome $chromosome \
		$graph_dir/ \
		$graph_dir/ \
		$graph_dir/linear_map_ \
		filtered_low_qual_reads_removed filtered_low_qual_reads_removed_ "" False $fragment_length $read_length \
		 > log_after_p_values_$chromosome.txt 2>&1 &
	echo "Peak calling for chr $chromosome started as process. Log will be written to $work_dir/log_after_p_values_$chromosome.txt"
    #else
    #    echo "P values not computed for $chromosome. Will not call peaks now."
    #fi
done
