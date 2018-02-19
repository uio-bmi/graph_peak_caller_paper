#!/usr/bin/env bash

encode_id=$1
replicate=$2
tf=$3
chromosomes=$4
motif_url=$5

base_dir=$(pwd)
grch38_fasta_file=~/data/hg19/hg19.fasta
n_threads=$(grep -c ^processor /proc/cpuinfo)

echo "Will use $n_threads threads"

echo "RUNNING"
work_dir=data/${tf}_${encode_id}/${replicate}
cd $work_dir

echo "Changed dir to $work_dir"



# Prepare linear reads for linear peak calling
if [ ! -f linear_alignments.bam ]; then

    # Get raw fastq for this encode experiment and replicate number
    if [ ! -f raw.fastq.gz ]; then
        echo "Download fastq"
        encode_url=$(python3 $base_dir/download_encode_fastq.py $encode_id $replicate)
        echo "Encode url: $encode_url"
        wget -O raw.fastq.gz --show-progress $encode_url
    else
        echo "Raw fastq already exists. Not dowloading"
    fi

    if [ ! -f raw.fastq ]; then
        echo "Unzipping"
        gunzip -c raw.fastq.gz > raw.fastq
    else
        echo "Unzipped fastq alread exists. Not unzipping"
    fi

    if [ ! -f raw_trimmed.fq ]; then
        echo "Trimming fastq using trim galore"
       #mv raw.fastq filtered.fastq
       trim_galore raw.fastq
    else
        echo "Fastq already trimmed"
    fi


    echo "Mapping reads to linear genome"
    bwa aln -t $n_threads $grch38_fasta_file raw.fastq > reads.sai
    bwa samse $grch38_fasta_file reads.sai raw_trimmed.fastq > alignments.sam

    # Convert to bam and sort
    echo "Converting to bam and filtering"
    samtools view -Su alignments.sam | samtools sort - alignments_sorted


    # Filter (removed duplicates and reads having low score)
    # bwa aln has max mapping quality 37
    samtools view -F 0x904 -q 37 -b alignments_sorted.bam > linear_alignments.bam
fi




# Run macs with encode linear mapped reads
if [ ! -f macs_output_whole_run.txt ]; then
	echo "Running macs2"
	macs2 callpeak -g hs -t linear_alignments.bam -n macs > macs_output_whole_run.txt 2>&1
else
	echo "Not running max. Already done"
fi

# Extract macs2 sequences, write to fasta (for later comparison)
echo "Extracting macs sequences"
> macs_selected_chromosomes.bed  # Create empty file
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo "Getting macs peaks for chromosome $chromosome"
    # Add to single file for all chromosomems
    grep "chr${chromosome}\s" macs_peaks.narrowPeak >> macs_selected_chromosomes.bed

    # Also sort out specific chromosome, making later analysis faster
    grep "chr${chromosome}\s" macs_peaks.narrowPeak > macs_peaks_chr${chromosome}.bed
    graph_peak_caller linear_peaks_to_fasta macs_peaks_chr${chromosome}.bed $grch38_fasta_file macs_sequences_chr${chromosome}.fasta


done


# Fetch macs sequences for these peaks
echo "Fetch macs sequences for selected chromosomes"
graph_peak_caller linear_peaks_to_fasta macs_selected_chromosomes.bed $grch38_fasta_file macs_sequences.fasta


# Merge all graph peak caller result files into one single sorted sequence file
graph_peak_caller concatenate_sequence_files $chromosomes sequence_all_chromosomes.fasta


# Run motif enrichment analysis
$base_dir/plot_motif_enrichments.sh sequence_all_chromosomes.fasta macs_sequences.fasta $motif_url motif_enrichment.png $tf
cp motif_enrichment.png ../../../figures_tables/$tf.png

# Also run fimo for each chromosome
for chromosome in $(echo $chromosomes | tr "," "\n")
do
    echo ""
    echo "----- Running fimo separately for chr $chromosome --- "
    fimo -oc fimo_macs_chr$chromosome motif.meme macs_sequences_chr${chromosome}.fasta
    fimo -oc fimo_graph_chr$chromosome motif.meme ${chromosome}_sequences.fasta
done

# Analyse peak results
graph_peak_caller analyse_peaks_whole_genome $chromosomes ./ ~/data/whole_genome/ ../../../figures_tables/$tf
