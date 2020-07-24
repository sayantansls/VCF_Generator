#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 21:59:56 2020

@author: sayantan
"""
import time as tm
import csv, errno, os
import twobitreader
import re, copy, sys
import datetime

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

def load_human_genome_sequence(genome_file):
    global genome
    genome = twobitreader.TwoBitFile(genome_file)
    
def load_metadata(vcf_template):
    metadata = list()
    f = open(vcf_template, 'r')
    for line in f.readlines():
        if line.startswith('##'):
            metadata.append(line)
    return metadata

def validate_gene(gene):
    if gene not in genes_set:
        raise Exception('ERROR : Given gene {} does not exist'.format(gene))
    else:
        print('INFO : Gene name {} is valid'.format(gene))
                
#The headers in the genes.tsv file are as follows:
#0 -- Strand_gene_id
#1 -- ChrName
#2 -- Strand
#3 -- Gene_start
#4 -- Gene_end
#5 -- Symbol
#6 -- Entrez_id
                           
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
        
#The headers in the input file are as follows:
#0 -- Gene name
#1 -- genomicHGVS

def process_input_file(input_file):
    variants_list = list()

    f = open(input_file)
    file_data = csv.DictReader(f, delimiter='\t')
    for variant in file_data:
        lst = [variant['Gene name'], variant['genomicHGVS']]
        variants_list.append(lst)
    return variants_list

def segregate_variants(variants_list):
    subs, others = [list(), list()]
    for variant in variants_list:
        if '>' in variant[1]:
            subs.append(variant)
        else:
            others.append(variant)
    return subs,others

def get_chromosome(gene):
    for k,v in gene_chrom_dict.items():
        if gene in v:
            return k

def del_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    if len(positions) == 1:
        pos = int(positions[0]) - 1
        ref = genome[chrom][pos-1:int(positions[0])].upper()
        alt = genome[chrom][pos-1].upper()
    else:
        start, end = [positions[0], positions[1]]
        pos = int(start) - 1
        ref = genome[chrom][pos-1:int(end)].upper()
        alt = genome[chrom][pos-1].upper()
    return pos, ref, alt

def other_handling(genomicHGVS, chrom):
    positions = re.findall(r'[0-9]+', genomicHGVS)
    ins_bases = ''.join(re.findall(r'[A-Z]+', genomicHGVS))
    pos = positions[0]
    ref = genome[chrom][int(pos)-1].upper()
    alt = ref + ins_bases
    return pos, ref, alt

def create_substitution_entries(output, subs):
    for variant in subs:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]

        ENTRY['#CHROM'] = get_chromosome(gene)
        ENTRY['POS'] = position
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [ENTRY[i] for i in HEADERS]
        output.write(sep.join(field_values))
        output.write('\n')

def create_non_substitution_entries(output, others):
    for variant in others:
        ENTRY = copy.deepcopy(ENTRY_T)
        gene, gHGVS = [variant[0], variant[1]]
        chrom = get_chromosome(gene)

        if "del" in gHGVS:
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

def check_file_status(genomic_build):
    print('INFO : Checking all essential Files status....')
    genome_file = genomic_build_dict[genomic_build]['genome_file']
    genesfile = genomic_build_dict[genomic_build]['genesfile']

    vcf_template = '../../../data/vcf_template/vcf_template.vcf'

    if os.path.exists(genome_file):
        load_human_genome_sequence(genome_file)
    else:
        raise Exception('Human Genome 2bit file {} not present in location'.format(genome_file))

    if os.path.exists(vcf_template):
        metadata = load_metadata(vcf_template)
    else:
        raise Exception('VCF template file {} not present in location'.format(vcf_template))
    
    if os.path.exists(genesfile):
        create_gene_chromosome_map(genesfile)
    else:
        raise Exception('Strand genes file {} not present in location'.format(genesfile))

    print('INFO : All input files present')
    return metadata
    

def input_file_handling(input_file, genomic_build):
    metadata = check_file_status(genomic_build)
    variants_list = process_input_file(input_file)

    print("INFO : Total variants entered : {}".format(len(variants_list)))
    process_variants(variants_list, metadata)

def single_entry_handling(gene, genomicHGVS, genomic_build):
    metadata = check_file_status(genomic_build)
    validate_gene(gene)
    variants_list = [[gene, genomicHGVS]]
    process_variants(variants_list, metadata)

def process_variants(variants_list, metadata):
    if len(variants_list) > 0:
        x = datetime.datetime.now()
        output_file = 'output_{}{}'.format(x, '.vcf')
        output = open(output_file, 'w+')
        for line in metadata:
            output.write(line)

        output.write(sep.join(HEADERS))
        output.write('\n')

        (subs, others) = segregate_variants(variants_list)
        print("INFO : Substitution variants: {}".format(len(subs)))
        print("INFO : INDEL variants: {}".format(len(others)))
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