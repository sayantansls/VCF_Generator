"""
author : sayantan (sayantan.ghosh@strandls.com)
This script contains all the utility functions used in vcf_generator_verX.py scripts
"""

import twobitreader
import csv, json

# Loads the genome information from the binary file of human genome reference
def load_human_genome_sequence(genome_file):
    genome = twobitreader.TwoBitFile(genome_file)
    return genome

# Loads the metadata to be written on the VCF file
def load_metadata(vcf_template):
    metadata = list()
    f = open(vcf_template, 'r')
    for line in f.readlines():
        if line.startswith('##'):
            metadata.append(line)
    return metadata

# Checks whether the ref from HGVS and ref from Genome is same or not
def check_ref_hgvs_genome(refFromHgvs, refFromGenome, genomicHGVS):
    ref = ''
    if refFromHgvs:
        if refFromHgvs == refFromGenome:
            ref = refFromHgvs
        else:
            print(messages['error_messages']['REF_MISMATCH'].format(genomicHGVS, refFromHgvs, refFromGenome))
            print(messages['error_messages']['VARIANT_SKIPPED'].format(genomicHGVS))
    else:
        ref = refFromGenome
    return ref

"""
Creates a gene to chromosome dictionary {'gene': 'chromosome'}                
The headers in the genes.tsv file are as follows:
0 -- Strand_gene_id
1 -- ChrName
2 -- Strand
3 -- Gene_start
4 -- Gene_end
5 -- Symbol
6 -- Entrez_id
"""                         
def create_gene_chromosome_map(genesfile):
    gene_chrom_dict, genes_set = [dict(), set()]

    f = open(genesfile, 'r')
    gene_data = csv.DictReader(f, delimiter='\t')    
    for gene in gene_data:
        if gene['ChrName'] not in gene_chrom_dict:
            gene_chrom_dict[gene['ChrName']] = list()
        gene_chrom_dict[gene['ChrName']].append(gene['Symbol'])

        genes_set.add(gene['Symbol'])
    return gene_chrom_dict, genes_set

# Returns the chromosome for a given gene
def get_chromosome(gene, gene_chrom_dict):
    for k,v in gene_chrom_dict.items():
        if gene in v:
            return k

"""
Process the user provided input file to create a list of variants        
The headers in the input file are as follows:
0 -- Gene name
1 -- genomicHGVS
"""
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

# Validates a gene provided by the user
def validate_gene(gene, genes_set):
    isGeneInvalid = 0
    if gene not in genes_set:
        isGeneInvalid = 1
    return isGeneInvalid

# Validates a position entered by the user
def validate_position(chromosome, position, ref, genome):
    isPositionInvalid = 0
    ref_from_genome = genome[chromosome][int(position) - 1].upper()
    if ref != ref_from_genome:
        isPositionInvalid = 1
    return isPositionInvalid, ref_from_genome

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
            print(messages['error_messages']['UNKNOWN_OBJECT'].format(gene, genomicHGVS))

    print(messages['success_messages']['DUPLICATION_VARIANTS'].format(len(dups)))
    print(messages['success_messages']['DELETION_VARIANTS'].format(len(dels)))
    print(messages['success_messages']['INSERTION_VARIANTS'].format(len(ins)))
    print(messages['success_messages']['DELINS_VARIANTS'].format(len(indels)))

# Reads json file which stores relevant information
def read_json(config_json):
	json_file = open(config_json, 'r')
	json_data = json.load(json_file)
	return json_data

messages = read_json('../../common/messages.json')