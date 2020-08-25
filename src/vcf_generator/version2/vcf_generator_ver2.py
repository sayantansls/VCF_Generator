"""
author: sayantan (sayantan.ghosh@strandls.com)
This script consumes 'n' number of input files (TSV format with headers 'Gene name' and 'genomicHGVS') 
and creates separate VCF files for each input file
"""
import csv, os
import re, copy, sys
import vcf_utils, variant_functions

sep = '\t'

json_data = vcf_utils.read_json('../../common/config.json')
messages = vcf_utils.read_json(json_data['messages'])

def setting_up_data(genome_file, genes_file):
    global genome, gene_chrom_dict, genes_set
    genome = vcf_utils.load_human_genome_sequence(genome_file)
    (gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(genes_file) 

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
def check_file_status():
    print(messages['success_messages']['CHECKING_FILE_STATUS'])

    vcf_template = json_data['vcf_template']

    if os.path.exists(vcf_template):
        metadata = vcf_utils.load_metadata(vcf_template)
    else:
        raise Exception(messages['error_messages']['VCF_TEMPLATE_NOT_PRESENT'].format(vcf_template))

    print(messages['success_messages']['FILES_PRESENT'])
    return metadata

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
        if others:
            vcf_utils.categorize_non_substitution_variants(others)

        create_substitution_entries(output, subs)
        create_non_substitution_entries(output, others)

# Regulates workflow for single gene and genomic HGVS input option
def main(input_file, genome_file, genes_file):
    print('\nINFO : File provided : {}'.format(input_file))
    metadata = check_file_status()
    setting_up_data(genome_file, genes_file)
    variants_list = vcf_utils.process_input_file(input_file)
    output_file = os.path.join(json_data['output'], os.path.basename(input_file).replace('.tsv', '.vcf'))

    print(messages['success_messages']['TOTAL_VARIANTS'].format(len(variants_list)))
    process_variants(variants_list, metadata, output_file)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])