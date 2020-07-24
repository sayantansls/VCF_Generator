# VCF generation 

## Definition

It is a python-based script that can take in genomic HGVS nomenclature and creates a VCF in a format that can be parsed by StrandOmics. 

## Execution

1. Clone the repository with the command - 

		git clone https://github.com/sayantansls/vcf-generator.git

2. Once the repository is cloned, enter the directory vcf-generator and run the script *setup.sh*

		./setup.sh

This will ensure all the dependecy files (human genome reference data) are downloaded and placed in their appropriate location.
You can download the files separately as well
	* **GrCh37** - http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.2bit
	* **GrCh38** - http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit

3. To run and use the vcf-generator, run the script *run.sh*

		./run.sh

Once the script is run, it will ask for input to select certain parameters:
	* **Genomic Build** - Enter 'A' for GrCh37 OR 'B' for GrCh38 (this genomic build will be used for mapping in the VCF)
	* **Use Case** - The script either consumes a TSV file with headers *Gene name* and *genomic HGVS* with the data or it can take a single gene and genomic HGVS and create a VCF file. (Enter '1' for file option OR '2' for single gene and genomic HGVS option).

## Limitations

1. Currently, the script does not support duplications - any HGVS nomenclature like g.125456dup will not be handled.
2. The script only consumes genomic HGVS nomenclature, not cDNA or protein nomenclature.

*Both the limitations are currently being worked on and will be resolved in future versions*

## Contact

Author - Sayantan Ghosh (sayantan.ghosh@strandls.com)

