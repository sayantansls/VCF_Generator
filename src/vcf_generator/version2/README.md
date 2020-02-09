### The script 'vcf_generator_ver2.py' consumes a TSV file of variants anc creates a VCF (acc. to StrandOmics format)

#### The format of the input TSV file is as follows:

	Gene name	genomicHGVS
	PLP1	g.103041656G>A
	MFSD8	g.128859994C>T
	CEP290	g.88512260C>T
	PNKP	g.50365794C>A 

#### The format of the VCF file consumed by StrandOmics is as follows:

	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	STRAN-404_S4_DMiSeq02-Run0059
	chr7	55249071	.	C	T	.	PASS	.	GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ	1/1:1000.00:608:54.28:54.28:0:0.41:36.11:254.01

#### The following command runs the script and creates the VCF file for the variants in the input TSV file:

	DATA=../../../data
	PYTHONPATH=../../common python vcf_generator_ver2.py \
	$DATA/examples/input_file.tsv \
	$DATA/examples/output_file.vcf \
	$DATA/human_genome_data/hg19.2bit \
	$DATA/strandomics_input_data/genes.tsv

[Change the input and output file names and directories according to your convenience]