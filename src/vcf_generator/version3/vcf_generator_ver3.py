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
        gene, gHGVS = variant.values()
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print(messages['error_messages']['INVALID_GENE'].format(gene))
        	print(messages['error_messages']['VARIANT_SKIPPED'].format(gene, gHGVS))
        	continue

        ref_alt = re.findall(r'[A-Z]', gHGVS)
        position = ''.join(re.findall(r'[0-9]+', gHGVS))
        ref, alt = [ref_alt[0], ref_alt[1]]
        chromosome = vcf_utils.get_chromosome(gene, gene_chrom_dict)
        (isPositionInvalid, actual_ref) = vcf_utils.validate_position({"chrom": chromosome, "pos": position, "ref": ref}, genome)

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
        (gene, gHGVS) = variant.values()
        isGeneInvalid = vcf_utils.validate_gene(gene, genes_set)

        if isGeneInvalid:
        	print(messages['error_messages']['INVALID_GENE'].format(gene))
        	print(messages['error_messages']['VARIANT_SKIPPED'].format(gene, gHGVS))
        	continue
        	
        chrom = vcf_utils.get_chromosome(gene, gene_chrom_dict)
        variant_functions.set_human_genome(genomic_build)

        if 'delins' in gHGVS:
            (pos, ref, alt) = variant_functions.process_delins_variant(gHGVS, chrom)
        elif 'dup' in gHGVS:
            (pos, ref, alt) = variant_functions.process_duplication_variant(gHGVS, chrom)
        elif 'del' in gHGVS and 'ins' not in gHGVS:
            (pos, ref, alt) = variant_functions.process_deletion_variant(gHGVS, chrom)
        else:
            (pos, ref, alt) = variant_functions.process_insertion_variant(gHGVS, chrom)

        ENTRY['#CHROM'] = chrom
        ENTRY['POS'] = pos
        ENTRY['REF'] = ref
        ENTRY['ALT'] = alt

        field_values = [str(ENTRY[i]) for i in json_data['vcf_header']]
        output.write(sep.join(field_values))
        output.write('\n')

# Checks the status of all the Input files required for running the script
def check_file_status():
    print(messages['success_messages']['CHECKING_FILE_STATUS'])
    global genome, gene_chrom_dict, genes_set, metadata
    genome = vcf_utils.load_human_genome_sequence(genomic_build)
    metadata = vcf_utils.load_metadata()
    (gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(genomic_build)
    print(messages['success_messages']['FILES_PRESENT'])
    
# Regulates workflow for single gene and genomic HGVS input option
def input_file_handling(input_file):
    check_file_status()
    variants_list = vcf_utils.process_input_file(input_file)
    output_file = os.path.join(json_data['output'], os.path.basename(input_file).replace('.tsv', '.vcf'))

    print(messages['success_messages']['TOTAL_VARIANTS'].format(len(variants_list)))
    process_variants(variants_list, output_file)

# Regulates workflow for file input option
def single_entry_handling(gene, genomicHGVS):
    check_file_status()
    vcf_utils.validate_gene(gene, genes_set)
    variants_list = [[gene, genomicHGVS]]
    output_file = os.path.join(json_data['output'],''.join(['{}-{}'.format(gene, genomicHGVS), '.vcf']))
    process_variants(variants_list, output_file)

# Calls all functions required for generating VCFs 
def process_variants(variants_list, output_file):
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
    try:
        input = raw_input
    except NameError:
        pass

    global genomic_build
    genomic_build = input('INFO : Choose Genomic Build -- Enter 37 for GrCh37 and Enter 38 for GrCh38 --> ')

    print('INFO : Enter a filename containing gene and genomic HGVS or enter gene name and genomic HGVS here')
    option = input('INFO : Enter 1 for file option ; Enter 2 for gene and genomic HGVS --> ')
    if option == '1':
        print('INFO : You have chosen the file option')
        filename = input('INFO : Enter the filename (file should be TSV format containing headers "Gene name" and "genomic HGVS"): ')
        input_file_handling(filename)
    elif option == '2':
        print('INFO : You have chosen to enter a single gene and genomic HGVS option')
        (gene, genomicHGVS) = input('INFO : Enter gene and genomic HGVS (space separated) : ').split(' ')
        single_entry_handling(gene, genomicHGVS)
    else:
        print('ERROR : This is not a valid option')