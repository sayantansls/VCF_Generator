#!/bin/bash
# author - sayantan(sayantan.ghosh@strandls.com)
# This script consumes multiple input TSV file, asks user to select genomic build, 
# checks file status and runs python script to create VCF file

HG_19="../../../data/human_genome_data/GrCh37/hg19.2bit"
HG_38="../../../data/human_genome_data/GrCh38/hg38.2bit"

GENES_37="../../../data/strandomics_input_data/GrCh37/genes.tsv"
GENES_38="../../../data/strandomics_input_data/GrCh38/genes_38.tsv"

FILES="$HG_19 $GENES_37 $HG_38 $GENES_38"

for argument in $@
do
	read -p "INPUT : Enter Genomic Build (37 or 38): "  genomic_build
	echo "INFO : You have chosen GrCh$genomic_build"

	
	files_absent="0"
	for file in $FILES
	do
		if [[ ! -f $file ]]; then
			echo "ERROR : $file is not present in location (run setup.sh)"
			files_absent="1"
		fi
	done

	if [[ $files_absent == "1" ]]; then
		echo "ERROR : Exiting script as one/more files absent"
		exit 0
	fi

	if [[ $genomic_build == '37' ]]; then
		GENOME=$HG_19
		GENES=$GENES_37
	elif [[ $genomic_build == '38' ]]; then
		GENOME=$HG_38
		GENES=$GENES_38
	else
		echo "ERROR : Incorrect genomic build provided"
	fi

	python vcf_generator_ver2.py $argument $GENOME $GENES
	echo -e "\n"
done
