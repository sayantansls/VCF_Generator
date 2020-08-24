"""
author: sayantan (sayantan.ghosh@strandls.com)
This is an interactive python scripts which asks the user to select genomic build - GrCh37 or GrCh38
The script also asks the user to either provide a file as input or a single gene and genomic HGVS
"""
import time as tm
import csv, os
import twobitreader
import re, copy, sys
import vcf_utils

HEADERS = ['#CHROM',
           'POS',
           'ID',
           'REF',
           'ALT',
           'QUAL',
           'FILTER',
           'INFO',
           'FORMAT',
           'STRAN-404_S4_DMiSeq02-Run0059']

ENTRY_T = {'#CHROM': '',
           'POS': '',
           'ID': '.',
           'REF': '',
           'ALT': '',
           'QUAL': '.',
           'FILTER': 'PASS',
           'INFO': '.',
           'FORMAT': 'GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ',
           'STRAN-404_S4_DMiSeq02-Run0059': '1/1:1000.00:608:54.28:54.28:0:0.41:36.11:254.01'}

# A -- GrCh37 and B -- GrCh38
genomic_build_dict = {'A': {'genome_file': '../../../data/human_genome_data/GrCh37/hg19.2bit', 
                            'genesfile': '../../../data/strandomics_input_data/GrCh37/genes.tsv'}, 
                      'B': {'genome_file': '../../../data/human_genome_data/GrCh38/hg38.2bit', 
                            'genesfile': '../../../data/strandomics_input_data/GrCh38/genes_38.tsv'}}

sep = '\t'
OUTPUT_DIR = '../../../data/output/'

# Gets the pos, ref, alt for a Deletion Variant (g.1124566delG)
def del_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0]) - 1
        ref = genome[chrom][pos-1:int(positions[0])].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start) - 1
        ref = genome[chrom][pos-1:int(end)].upper()
        
    alt = genome[chrom][pos-1].upper()
    return pos, ref, alt

# Gets the pos, ref, alt for an Insertion Variant (g.14587156insAA)
def other_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    ins_bases = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    pos = positions[0]
    ref = genome[chrom][int(pos)-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

# Gets the pos, ref, alt for a Duplication Variant (g.98745415dupA)
def dup_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    refFromHgvs = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    if len(positions) == 1:
        pos = int(positions[0])
        refFromGenome = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        refFromGenome = genome[chrom][pos-1:int(end)].upper()
        
    ref = vcf_utils.check_ref_hgvs_genome(refFromHgvs, refFromGenome, genomicHGVS)
    alt = ref * 2
    return pos, ref, alt

# Gets the pos, ref, alt for a Delins Variant (g.54669745_54669748delinsCTGG)
def delins_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0])
        ref = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start)
        ref = genome[chrom][pos-1:int(end)].upper()

    alt = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    return pos, ref, alt

# Writes out the VCF entries for the Substitution Variants
def create_substitution_entries(output, subs):
    for variant in subs:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print('WARN: Gene name "{}" given is invalid'.format(gene))
        	print('WARN: This variant "{} {}" will be skipped'.format(gene, gHGVS))
        	continue

        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]
        chromosome = vcf_utils.get_chromosome(gene, gene_chrom_dict)
        (isPositionInvalid, actual_ref) = vcf_utils.validate_position(chromosome, position, ref, genome)

        if isPositionInvalid:
        	print('WARN: The genomic HGVS "{}" given is invalid'.format(gHGVS))
        	print('WARN: For position "{}", the ref provided is {} and actual ref is {}'.format(position, ref, actual_ref))
        	print('WARN: This variant "{} {}" will be skipped'.format(gene, gHGVS))
        	continue

        ENTRY['#CHROM'] = chromosome
        ENTRY['POS'] = position
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [ENTRY[i] for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

# Writes out the VCF entries for Non-Substitution Variants
def create_non_substitution_entries(output, others):
    for variant in others:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print('WARN: Gene name "{}" given is invalid'.format(gene))
        	print('WARN: This variant "{} {}" will be skipped'.format(gene, gHGVS))
        	continue
        	
        chrom = vcf_utils.get_chromosome(gene, gene_chrom_dict)

        if 'delins' in gHGVS:
            (pos, ref, alt) = delins_handling(gHGVS, chrom)
        elif 'dup' in gHGVS:
            (pos, ref, alt) = dup_handling(gHGVS, chrom)
        elif 'del' in gHGVS and 'ins' not in gHGVS:
            (pos, ref, alt) = del_handling(gHGVS, chrom)
        else:
            (pos, ref, alt) = other_handling(gHGVS, chrom)

        ENTRY['#CHROM'] = chrom
        ENTRY['POS'] = pos
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [str(ENTRY[i]) for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

# Checks the status of all the Input files required for running the script
def check_file_status(genomic_build):
    print('INFO : Checking all essential Files status....')
    genome_file = genomic_build_dict[genomic_build]['genome_file']
    genesfile = genomic_build_dict[genomic_build]['genesfile']

    vcf_template = '../../../data/vcf_template/vcf_template.vcf'
    global genome, gene_chrom_dict, genes_set

    if os.path.exists(genome_file):
        genome = vcf_utils.load_human_genome_sequence(genome_file)
    else:
        raise Exception('Human Genome 2bit file {} not present in location'.format(genome_file))

    if os.path.exists(vcf_template):
        metadata = vcf_utils.load_metadata(vcf_template)
    else:
        raise Exception('VCF template file {} not present in location'.format(vcf_template))
    
    if os.path.exists(genesfile):
        (gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(genesfile)
    else:
        raise Exception('Strand genes file {} not present in location'.format(genesfile))

    print('INFO : All input files present')
    return metadata
    
# Regulates workflow for single gene and genomic HGVS input option
def input_file_handling(input_file, genomic_build):
    metadata = check_file_status(genomic_build)
    variants_list = vcf_utils.process_input_file(input_file)
    output_file = os.path.join(OUTPUT_DIR, os.path.basename(input_file).replace('.tsv', '.vcf'))

    print("INFO : Total variants entered : {}".format(len(variants_list)))
    process_variants(variants_list, metadata, output_file)

# Regulates workflow for file input option
def single_entry_handling(gene, genomicHGVS, genomic_build):
    metadata = check_file_status(genomic_build)
    vcf_utils.validate_gene(gene, genes_set)
    variants_list = [[gene, genomicHGVS]]
    output_file = os.path.join(OUTPUT_DIR,''.join(['{}-{}'.format(gene, genomicHGVS), '.vcf']))
    process_variants(variants_list, metadata, output_file)

# Calls all functions required for generating VCFs 
def process_variants(variants_list, metadata, output_file):
    if len(variants_list) > 0:
        output = open(output_file, 'w+')
        for line in metadata:
            output.write(line)

        output.write(sep.join(HEADERS))
        output.write('\n')

        (subs, others) = vcf_utils.segregate_variants(variants_list)
        print("INFO : Substitution variants: {}".format(len(subs)))
        print("INFO : Non-Substitution variants: {}".format(len(others)))
        if others: vcf_utils.categorize_non_substitution_variants(others)
        create_substitution_entries(output, subs)
        create_non_substitution_entries(output, others)

if __name__ == '__main__':

    print('Start of code : {}'.format(tm.ctime(tm.time())))
    try:
        input = raw_input
    except NameError:
        pass

    genomic_build = input('INFO: Choose Genomic Build -- Enter A for GrCh37 and Enter B for GrCh38 --> ')

    print('INFO : Enter a filename containing gene and genomic HGVS or enter gene name and genomic HGVS here')
    option = input('INFO : Enter 1 for file option ; Enter 2 for gene and genomic HGVS --> ')
    if option == '1':
        print('INFO : You have chosen the file option')
        filename = input('INPUT : Enter the filename (file should be TSV format containing headers "Gene name" and "genomic HGVS"): ')
        input_file_handling(filename, genomic_build)
    elif option == '2':
        print('INFO : You have chosen to enter a single gene and genomic HGVS option')
        (gene, genomicHGVS) = input('INPUT : Enter gene and genomic HGVS (space separated) : ').split(' ')
        single_entry_handling(gene, genomicHGVS, genomic_build)
    else:
        print('ERROR : This is not a valid option')

    print('End of code : {}'.format(tm.ctime(tm.time())))