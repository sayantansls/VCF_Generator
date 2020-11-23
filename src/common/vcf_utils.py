"""
author : sayantan (sayantan.ghosh@strandls.com)
This script contains all the utility functions used in vcf_generator_verX.py scripts
"""

import twobitreader
import csv, json, os

# Reads json file which stores relevant information
def read_json(config_json):
    json_file = open(config_json, 'r')
    json_data = json.load(json_file)
    return json_data

json_data = read_json('../../common/config.json')
messages = read_json(json_data["messages"])

# Loads the genome information from the binary file of human genome reference
def load_human_genome_sequence(genomic_build):
    genome_file = json_data[genomic_build]["genome"]
    if not os.path.exists(genome_file):
        raise Exception(messages['error_messages']['GENOME_NOT_PRESENT'].format(genome_file))

    return twobitreader.TwoBitFile(genome_file)

# Loads the metadata to be written on the VCF file
def load_metadata():
    vcf_template = json_data["vcf_template"]
    if not os.path.exists(vcf_template):
        raise Exception(messages['error_messages']['VCF_TEMPLATE_NOT_PRESENT'].format(vcf_template))

    metadata = list()
    f = open(vcf_template, 'r')
    for line in f.readlines():
        if line.startswith('##'):
            metadata.append(line)
    f.close()

    return metadata

# Checks whether the ref from HGVS and ref from Genome is same or not
def validate_reference(options):
    (refFromHgvs, refFromGenome, genomicHgvs) = options.values()
    ref = ''
    if refFromHgvs:
        if refFromHgvs == refFromGenome:
            ref = refFromHgvs
        else:
            print(messages['error_messages']['REF_MISMATCH'].format(genomicHgvs, refFromHgvs, refFromGenome))
            print(messages['error_messages']['VARIANT_SKIPPED'].format(genomicHgvs))
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
def create_gene_chromosome_map(genomic_build):
    genes_file = json_data[genomic_build]["genes"]
    if not os.path.exists(genes_file):
        raise Exception(messages['error_messages']['GENES_NOT_PRESENT'].format(genes_file))

    gene_chrom_dict, genes_set = [dict(), set()]
    f = open(genes_file, 'r')
    gene_data = csv.DictReader(f, delimiter='\t')    
    for gene in gene_data:
        if gene['ChrName'] not in gene_chrom_dict:
            gene_chrom_dict[gene['ChrName']] = list()
        gene_chrom_dict[gene['ChrName']].append(gene['Symbol'])

        genes_set.add(gene['Symbol'])
    f.close()

    return gene_chrom_dict, genes_set

# Returns the chromosome for a given gene
def get_chromosome(gene, gene_chrom_dict = {}):
    if not gene_chrom_dict:
        gene_chrom_dict = create_gene_chromosome_map()[0]
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
        variant_obj = {"gene": variant['Gene name'], "genomicHgvs": variant['genomicHGVS']}
        variants_list.append(variant_obj)
    f.close()

    return variants_list

# Segregates a list of variants into - Substitution and Non-Substitution Variants
def segregate_variants(variants_list):
    subs, others = [list(), list()]
    for variant in variants_list:
        if '>' in variant["genomicHgvs"]:
            subs.append(variant)
        else:
            others.append(variant)
    return subs,others

# Validates a gene provided by the user
def validate_gene(gene, genes_set = {}):
    isGeneInvalid = False
    if gene not in genes_set:
        isGeneInvalid = True
    return isGeneInvalid

# Validates a position entered by the user
def validate_position(options, genome):
    (chrom, pos, ref) = options.values()
    isPositionInvalid = False
    ref_from_genome = genome[chrom][int(pos) - 1].upper()
    if ref != ref_from_genome:
        isPositionInvalid = True
    return isPositionInvalid, ref_from_genome

# Categorize the Non-Substitution Variants
def categorize_non_substitution_variants(others):
    dups, dels, indels, ins = [list(), list(), list(), list()]
    for non_sub in others:
        (gene, genomicHGVS) = [non_sub["gene"], non_sub["genomicHgvs"]]
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