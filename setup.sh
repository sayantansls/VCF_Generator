echo "INFO : Starting the setup for vcf-generator"
echo "INFO : Checking scripts and files status"

cd data/

HG_19="human_genome_data/GrCh37/hg19.2bit"
HG_38="human_genome_data/GrCh38/hg38.2bit"

GENES_37="strandomics_input_data/GrCh37/genes.tsv"
GENES_38="strandomics_input_data/GrCh38/genes_38.tsv"

if [[ -f $GENES_37 ]]; then
	echo "SUCCESS : GrCh37 genes.tsv file present"
else
	echo "FAIL : GrCh37 genes.tsv file absent"
fi

if [[ -f $GENES_38 ]]; then
	echo "SUCCESS : GrCh38 genes_38.tsv file present"
else
	echo "FAIL : GrCh38 genes_38.tsv file absent"
fi

if [[ -f vcf_template/vcf_template.vcf ]]; then
	echo "SUCCESS : VCF template file present"
else
	echo "FAIL : VCF template file absent"
fi

echo "INFO : The binary human genome reference data to be downloaded"

if [[ ! -f $HG_19 ]]; then
	echo "INFO : Begining the download of GrCh37 data - hg19.2bit....."
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
	echo "INFO : Download of GrCh37 data - hg19.2bit ends"
	mv hg19.2bit human_genome_data/GrCh37/
fi

if [[ -f $HG_19 ]]; then
	echo "SUCCESS : hg19.2bit file downloaded and placed in appropriate location"
else
	echo "FAIL : hg19.2bit file absent"
fi

if [[ ! -f $HG_38 ]]; then
	echo "INFO : Begining the download of GrCh38 data - hg38.2bit....."
	wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
	echo "INFO : Download of GrCh38 data - hg38.2bit ends"
	mv hg38.2bit human_genome_data/GrCh38/
fi

if [[ -f $HG_38 ]]; then
	echo "SUCCESS : hg38.2bit file downloaded and placed in appropriate location"
else
	echo "FAIL : hg38.2bit file absent"
fi

echo "INFO : All input files successfully present"

