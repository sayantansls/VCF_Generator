# VCF Generation Script (Version 2)

## Definition
The script 'vcf_generator_ver2.py' consumes a TSV file of variants and creates a VCF file (acc. to StrandOmics format)

## Input and Output Format

### Input
	Gene name	genomicHGVS
	PLP1	g.103041656G>A
	MFSD8	g.128859994C>T

### Output

	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	STRAN-404_S4_DMiSeq02-Run0059
	chr7	55249071	.	C	T	.	PASS	.	GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ	1/1:1000.00:608:54.28:54.28:0:0.41:36.11:254.01

## Execution

1. Run the 'setup.sh' script present in the home directory (required for setting up input files)

		./setup.sh

2. To generate VCF files, run the 'run.sh' script with input files as arguments

		./run.sh [input-file-path1] [input-file-path2] ....

## Limitations

1. This script requires a file and does not support single gene and genomic HGVS entry (implemented in version 3)

## Contact

Author - Sayantan Ghosh (sayantan.ghosh@strandls.com)