#!/usr/bin/env nextflow

params.encode_experiment_id = "ENCSR000DUB"
params.replicate_number = "1"
params.tf_name = "CTCF"
params.work_dir = "/home/ivargry/dev/graph_chip_seq_pipeline/data"

params.graph_dir = "/home/ivargry/dev/graph_peak_caller/tests/lrc_kir/"
params.graph_peak_caller = "/home/ivargry/dev/graph_peak_caller/graph_peak_caller.py"

peak_caller = params.graph_peak_caller
encode_experiment_id = params.encode_experiment_id
tf_name = params.tf_name
replicate_number = params.replicate_number
vg_xg_index = params.graph_dir + "graph.xg"
vg_gcsa_index = params.graph_dir + "graph.gcsa"

process downloadData{
	
	input:
	val tf_name
	val replicate_number
	
	output:
	file raw_reads
	
    """
	encode_url=\$(python3 $baseDir/download_encode_fastq.py $encode_experiment_id $replicate_number)
    echo \$encode_url
    wget -qO- \$encode_url > raw.fastq.gz
    gunzip -c raw.fastq.gz > raw_reads
    """
}

process filterRawReads{
	input:
	file raw_reads
	
	output:
	file filtered_reads
	
	"""
	head -n 1000 raw_reads > filtered_reads
	"""
}

process mapReads{
	input:
	file filtered_reads
	
	output:
	file mapped_to_graph
		
	"""
	vg map -f filtered_reads -g $vg_gcsa_index -x $vg_xg_index -M 2 > mapped_to_graph
	"""
}

process filterMappedReads{
	input:
	file mapped_to_graph
	
	output:
	file filtered_json
	
	"""
	vg filter -r 1.0 -s 2.0 -fu mapped_to_graph > filtered.gam
	vg view -aj filtered.gam > filtered_json
	"""
}


process splitMappedReadsBychromosome{
	input:
	file filtered_json

	output:
	file "filtered_*.json" into result
	
	"""
	python3 $peak_caller split_vg_json_reads_into_chromosomes filtered_json $params.graph_dir 
	"""
}


result.flatMap().subscribe {
    println it.trim()
}


