"""
author: sayantan (sayantan.ghosh@strandls.com)
This script consumes 'n' number of input files (TSV format with headers 'Gene name' and 'genomicHGVS') 
and creates separate VCF files for each input file
"""
import time as tm
import csv, os
import twobitreader
import re, copy, sys

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

sep = '\t'
OUTPUT_DIR = '../../../data/output/'

# Loads the genome information from the binary file of human genome reference
def load_human_genome_sequence(genome_file):
    global genome
    genome = twobitreader.TwoBitFile(genome_file)
    
# Loads the metadata to be written on the VCF file
def load_metadata(vcf_template):
    metadata = list()
    f = open(vcf_template, 'r')
    for line in f.readlines():
        if line.startswith('##'):
            metadata.append(line)
    return metadata

# Validates a gene provided by the user
def validate_gene(gene):
    isGeneInvalid = 0
    if gene not in genes_set:
        isGeneInvalid = 1
    return isGeneInvalid

# Validates a position entered by the user
def validate_position(chromosome, position, ref):
    isPositionInvalid = 0
    ref_from_genome = genome[chromosome][int(position) - 1].upper()
    if ref != ref_from_genome:
        isPositionInvalid = 1
    return isPositionInvalid, ref_from_genome

# Creates a gene to chromosome dictionary {'gene': 'chromosome'}                         
def create_gene_chromosome_map(genesfile):
    global gene_chrom_dict
    gene_chrom_dict = dict()

    global genes_set 
    genes_set = set()

    f = open(genesfile, 'r')
    gene_data = csv.DictReader(f, delimiter='\t')    
    for gene in gene_data:
        if gene['ChrName'] not in gene_chrom_dict:
            gene_chrom_dict[gene['ChrName']] = list()
        gene_chrom_dict[gene['ChrName']].append(gene['Symbol'])

        genes_set.add(gene['Symbol'])

# Process the user provided input file to create a list of variants
def process_input_file(input_file):
    variants_list = list()

    f = open(input_file)
    file_data = csv.DictReader(f, delimiter='\t')
    for variant in file_data:
        lst = [variant['Gene name'], variant['genomicHGVS']]
        variants_list.append(lst)
    return variants_list

# Segregates a list of variants into - Substitution and Non-Substitution Variants
def segregate_variants(variants_list):
    subs, others = [list(), list()]
    for variant in variants_list:
        if '>' in variant[1]:
            subs.append(variant)
        else:
            others.append(variant)
    return subs,others

# Returns the chromosome for a given gene
def get_chromosome(gene):
    for k,v in gene_chrom_dict.items():
        if gene in v:
            return k

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

# Checks whether the ref from HGVS and ref from Genome is same or not
def check_ref_hgvs_genome(refFromHgvs, refFromGenome, genomicHGVS):
    ref = ''
    if refFromHgvs:
        if refFromHgvs == refFromGenome:
            ref = refFromHgvs
        else:
            print('WARN : For Variant {}, ref in HGVS - {}, ref from genome - {}'.format(genomicHGVS, refFromHgvs, refFromGenome))
            print('WARN : Variant {} is faulty and will be skipped'.format(genomicHGVS))
    else:
        ref = refFromGenome
    return ref

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
        
    ref = check_ref_hgvs_genome(refFromHgvs, refFromGenome, genomicHGVS)
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
        isGeneInvalid = validate_gene(gene)

        if isGeneInvalid:
            print('WARN: Gene name "{}" given is invalid'.format(gene))
            print('WARN: This variant "{} {}" will be skipped'.format(gene, gHGVS))
            continue

        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]
        chromosome = get_chromosome(gene)
        (isPositionInvalid, actual_ref) = validate_position(chromosome, position, ref)

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
        isGeneInvalid = validate_gene(gene)

        if isGeneInvalid:
            print('WARN: Gene name "{}" given is invalid'.format(gene))
            print('WARN: This variant "{} {}" will be skipped'.format(gene, gHGVS))
            continue
            
        chrom = get_chromosome(gene)

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
def check_file_status(genome_file, genes_file):
    print('INFO : Checking all essential Files status....')

    vcf_template = '../../../data/vcf_template/vcf_template.vcf'

    if os.path.exists(vcf_template):
        metadata = load_metadata(vcf_template)
    else:
        raise Exception('VCF template file {} not present in location'.format(vcf_template))

    print('INFO : All input files present')
    return metadata

# Categorize the Non-Substitution Variants
def categorize_non_substitution_variants(others):
    dups, dels, indels, ins = [list(), list(), list(), list()]
    for non_sub in others:
        gene, genomicHGVS = [non_sub[0], non_sub[1]]
        if 'dup' in genomicHGVS:
            dups.append(genomicHGVS)
        elif 'del' in genomicHGVS and 'ins' not in genomicHGVS:
            dels.append(genomicHGVS)
        elif 'ins' in genomicHGVS and 'del' not in genomicHGVS:
            ins.append(genomicHGVS)
        elif 'delins' in genomicHGVS:
            indels.append(genomicHGVS)
        else:
            print('File object {} {} is of unknown category'.format(gene, genomicHGVS))

    print('--- Duplication Variants : {}'.format(len(dups)))
    print('--- Deletion Variants : {}'.format(len(dels)))
    print('--- Insertion Variants : {}'.format(len(ins)))
    print('--- Indels/Delins : {}'.format(len(indels)))

# Calls all functions required for generating VCFs 
def process_variants(variants_list, metadata, output_file):
    if len(variants_list) > 0:
        output = open(output_file, 'w+')
        for line in metadata:
            output.write(line)

        output.write(sep.join(HEADERS))
        output.write('\n')

        (subs, others) = segregate_variants(variants_list)
        print('INFO : Substitution variants: {}'.format(len(subs)))
        print('INFO : Non-Substitution variants: {}'.format(len(others)))
        if others:
            categorize_non_substitution_variants(others)

        create_substitution_entries(output, subs)
        create_non_substitution_entries(output, others)

# Regulates workflow for single gene and genomic HGVS input option
def input_file_handling(input_file, genome_file, genes_file):
    print('\nINFO : File provided : {}'.format(input_file))
    metadata = check_file_status(genome_file, genes_file)
    variants_list = process_input_file(input_file)
    output_file = os.path.join(OUTPUT_DIR, os.path.basename(input_file).replace('.tsv', '.vcf'))

    print('INFO : Total variants entered : {}'.format(len(variants_list)))
    process_variants(variants_list, metadata, output_file)

if __name__ == '__main__':
    input_file_handling(sys.argv[1], sys.argv[2], sys.argv[3])