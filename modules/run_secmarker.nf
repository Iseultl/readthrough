#!/usr/bin/env nextflow

process RUN_SECMARKER {
	tag "${genome_fasta.baseName}"
	label 'secmarker'
	publishDir "${params.output_dir}/secmarker", mode: 'copy'
	errorStrategy 'retry'
	maxRetries 2
	cpus 1
	memory '8GB'
	time '1h'

	input:
	path genome_fasta

	output:
	path "secmarker_output/*"

	script:
	"""
	mkdir -p secmarker_output

	# Equivalent to:
	# docker run -v \$(pwd):/data ileahy/secmarker:v1.0 \
	#   python2 /opt/secmarker-0.4/Secmarker.py -t /data/<input.fa> -o /data/<output_dir>
	python2 /opt/secmarker-0.4/Secmarker.py \
		-t ${genome_fasta} \
		-o secmarker_output
	"""
}