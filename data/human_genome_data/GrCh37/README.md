## Reference Genome Data (HUMAN) - GrCh37

The script consumes the hg19.2bit (human genome binary data) to set up the reference genome for mapping the position to base.
The binary file is of considerable size (~700MB) and is too large in size to be pushed into GitHub repository.
Download the required file from the link given below:

		http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.2bit
		
After downloading the file, move it to the location - '/vcf-generator/data/human_genome_data/GrCh37/'

[The script searches for the file in this given location - it will fail if hg19.2bit is not in the specified location]
