"""
author : sayantan (sayantan.ghosh@strandls.com)
This script contains all the utility functions used in vcf_generator_verX.py scripts
"""

import twobitreader
import csv

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
            print('WARN : For Variant {}, ref in HGVS - {}, ref from genome - {}'.format(genomicHGVS, refFromHgvs, refFromGenome))
            print('WARN : Variant {} is faulty and will be skipped'.format(genomicHGVS))
    else:
        ref = refFromGenome
    return ref

# Creates a gene to chromosome dictionary {'gene': 'chromosome'}                         
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
            print('File object {} {} is of unknown category'.format(gene, genomicHGVS))

    print('--- Duplication Variants : {}'.format(len(dups)))
    print('--- Deletion Variants : {}'.format(len(dels)))
    print('--- Insertion Variants : {}'.format(len(ins)))
    print('--- Indels/Delins : {}'.format(len(indels)))