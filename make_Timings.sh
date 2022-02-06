#!/bin/bash

tmpfile=$(mktemp)

#Loop over all subjects and format timing files into FSL format
#for subj in `cat participants.tsv | tr -d '\r'` ; do
for subj in sub-FG25 sub-RG61R; do
	echo "generating for $subj"

	input="$subj/func/${subj}_task-Extinction_events.tsv"
	output_prefix="derivatives/afni/$subj"

	if [ -e $input ]; then
		cat $input | awk -F'\t' '{if ($3=="relief csm") {print $1, $2, "1"}}' | sort -n -k1 -t $'\t' > $tmpfile
		timing_tool.py -fsl_timing_files $tmpfile -write_timing $output_prefix/relief_csm.1D
		
		cat $input | awk -F'\t' '{if ($3=="relief csav") {print $1, $2, "1"}}' | sort -n -k1 -t $'\t' > $tmpfile
		timing_tool.py -fsl_timing_files $tmpfile -write_timing $output_prefix/relief_csav.1D
		
		cat $input | awk -F'\t' '{if ($3=="relief rating") {print $1, $2, "1"}}' | sort -n -k1 -t $'\t' > $tmpfile
		timing_tool.py -fsl_timing_files $tmpfile -write_timing $output_prefix/relief_rating.1D
	fi
done
