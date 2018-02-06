#!/usr/bin/env bash

encode_id=$1
replicate=$2
tf=$3
chromosomes=$4
bam_alignments_url=$5
motif_url=$6
echo $bam_alignments_url
base_dir=$(pwd)
grch38_fasta_file=~/dev/hg38.fasta


echo "RUNNING"
work_dir=data/${tf}_${encode_id}/${replicate}
cd $work_dir

echo "Changed dir to $work_dir"


# Downlaod bam alignments
if [ ! -f linear_alignments.bam ]; then
    echo "Downloading linear alignments (for macs2)"
    wget -O linear_alignments.bam $bam_alignments_url
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
    grep "chr${chromosome}\s" macs_peaks.narrowPeak > macs_peaks_chr$chr.bed
    graph_peak_caller linear_peaks_to_fasta macs_peaks_chr$chr.bed $grch38_fasta_file macs_sequences_chr$chr.fasta

done


# Fetch macs sequences for these peaks
echo "Fetch macs sequences for selected chromosomes"
graph_peak_caller linear_peaks_to_fasta macs_selected_chromosomes.bed $grch38_fasta_file macs_sequences.fasta


# Merge all graph peak caller result files into one single sorted sequence file
graph_peak_caller concatenate_sequence_files $chromosomes sequence_all_chromosomes.fasta


# Run motif enrichment analysis
$base_dir/plot_motif_enrichments.sh sequence_all_chromosomes.fasta macs_sequences.fasta $motif_url motif_enrichment.png
cp motif_enrichment.png ../../../figures_tables/$tf.png
