"""
author: sayantan (sayantan.ghosh@strandls.com)
This is an interactive python scripts which asks the user to select genomic build - GrCh37 or GrCh38
The script also asks the user to either provide a file as input or a single gene and genomic HGVS
"""
import time as tm
import csv, os
import re, copy, sys
import vcf_utils, variant_functions

sep = '\t'

json_data = vcf_utils.read_json('../../common/config.json')
messages = vcf_utils.read_json(json_data['messages'])

# Writes out the VCF entries for the Substitution Variants
def create_substitution_entries(output, subs):
    for variant in subs:
        ENTRY = copy.deepcopy(json_data['vcf_entry'])
        gene, gHGVS = [variant[0], variant[1]]
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print(messages['error_messages']['INVALID_GENE'].format(gene))
        	print(messages['error_messages']['VARIANT_SKIPPED'].format(gene, gHGVS))
        	continue

        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]
        chromosome = vcf_utils.get_chromosome(gene, gene_chrom_dict)
        (isPositionInvalid, actual_ref) = vcf_utils.validate_position(chromosome, position, ref, genome)

        if isPositionInvalid:
        	print(messages['error_messages']['INVALID_GENOMIC_HGVS'].format(gHGVS))
        	print(messages['error_messages']['REF_MISMATCH'].format(position, ref, actual_ref))
        	print(messages['error_messages']['VARIANT_SKIPPED'].format(gene, gHGVS))
        	continue

        ENTRY['#CHROM'] = chromosome
        ENTRY['POS'] = position
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [ENTRY[i] for i in json_data['vcf_header']]
        output.write(sep.join(field_values))
        output.write('\n')

# Writes out the VCF entries for Non-Substitution Variants
def create_non_substitution_entries(output, others):
    for variant in others:
        ENTRY = copy.deepcopy(json_data['vcf_entry'])
        gene, gHGVS = [variant[0], variant[1]]
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print(messages['error_messages']['INVALID_GENE'].format(gene))
        	print(messages['error_messages']['VARIANT_SKIPPED'].format(gene, gHGVS))
        	continue
        	
        chrom = vcf_utils.get_chromosome(gene, gene_chrom_dict)

        if 'delins' in gHGVS:
            (pos, ref, alt) = variant_functions.delins_handling(gHGVS, chrom, genome)
        elif 'dup' in gHGVS:
            (pos, ref, alt) = variant_functions.duplication_handling(gHGVS, chrom, genome)
        elif 'del' in gHGVS and 'ins' not in gHGVS:
            (pos, ref, alt) = variant_functions.deletion_handling(gHGVS, chrom, genome)
        else:
            (pos, ref, alt) = variant_functions.insertion_handling(gHGVS, chrom, genome)

        ENTRY['#CHROM'] = chrom
        ENTRY['POS'] = pos
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [str(ENTRY[i]) for i in json_data['vcf_header']]
        output.write(sep.join(field_values))
        output.write('\n')

# Checks the status of all the Input files required for running the script
def check_file_status(genomic_build):
    print(messages['success_messages']['CHECKING_FILE_STATUS'])
    genome_file, genesfile, vcf_template = [json_data[genomic_build]['genome'], json_data[genomic_build]['genes'], json_data['vcf_template']]
    global genome, gene_chrom_dict, genes_set

    if os.path.exists(genome_file):
        genome = vcf_utils.load_human_genome_sequence(genome_file)
    else:
        raise Exception(messages['error_messages']['GENOME_NOT_PRESENT'].format(genome_file))

    if os.path.exists(vcf_template):
        metadata = vcf_utils.load_metadata(vcf_template)
    else:
        raise Exception(messages['error_messages']['VCF_TEMPLATE_NOT_PRESENT'].format(vcf_template))
    
    if os.path.exists(genesfile):
        (gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(genesfile)
    else:
        raise Exception(messages['error_messages']['GENES_NOT_PRESENT'].format(genesfile))

    print(messages['success_messages']['FILES_PRESENT'])
    return metadata
    
# Regulates workflow for single gene and genomic HGVS input option
def input_file_handling(input_file, genomic_build):
    metadata = check_file_status(genomic_build)
    variants_list = vcf_utils.process_input_file(input_file)
    output_file = os.path.join(json_data['output'], os.path.basename(input_file).replace('.tsv', '.vcf'))

    print(messages['success_messages']['TOTAL_VARIANTS'].format(len(variants_list)))
    process_variants(variants_list, metadata, output_file)

# Regulates workflow for file input option
def single_entry_handling(gene, genomicHGVS, genomic_build):
    metadata = check_file_status(genomic_build)
    vcf_utils.validate_gene(gene, genes_set)
    variants_list = [[gene, genomicHGVS]]
    output_file = os.path.join(json_data['output'],''.join(['{}-{}'.format(gene, genomicHGVS), '.vcf']))
    process_variants(variants_list, metadata, output_file)

# Calls all functions required for generating VCFs 
def process_variants(variants_list, metadata, output_file):
    if len(variants_list) > 0:
        output = open(output_file, 'w+')
        for line in metadata:
            output.write(line)

        output.write(sep.join(json_data['vcf_header']))
        output.write('\n')

        (subs, others) = vcf_utils.segregate_variants(variants_list)
        print(messages['success_messages']['SUB_VARIANTS'].format(len(subs)))
        print(messages['success_messages']['NON_SUB_VARIANTS'].format(len(others)))
        if others: vcf_utils.categorize_non_substitution_variants(others)
        create_substitution_entries(output, subs)
        create_non_substitution_entries(output, others)

if __name__ == '__main__':

    print('Start of code : {}'.format(tm.ctime(tm.time())))
    try:
        input = raw_input
    except NameError:
        pass

    genomic_build = input('INFO: Choose Genomic Build -- Enter 37 for GrCh37 and Enter 38 for GrCh38 --> ')

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