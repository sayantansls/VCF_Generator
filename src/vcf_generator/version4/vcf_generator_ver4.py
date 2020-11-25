"""
author: sayantan (sayantan.ghosh@strandls.com)
An object oriented approach for vcf-generation
"""

import sys, csv, re, os
import vcf_utils

GENOMIC_BUILD = '37'
json_data = vcf_utils.read_json('../../common/config.json')
messages = vcf_utils.read_json(json_data['messages'])
genome = vcf_utils.load_human_genome_sequence(GENOMIC_BUILD)
(gene_chrom_dict, genes_set) = vcf_utils.create_gene_chromosome_map(GENOMIC_BUILD)

# Datamodel(object) that defines a Variant
class Variant():
	def __init__(self, file_entry):
		self.file_entry = file_entry

	def getChromosome(self):
		return vcf_utils.get_chromosome(self.file_entry.gene, gene_chrom_dict)

# Datamodel(object) that defines a Substitution Variant
class SubstitutionVariant(Variant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

	def getPosition(self):
		return int(''.join(re.findall(r'[0-9]+', self.file_entry.genomicHgvs)))

	def getReference(self, chromosome, position):
		refFromHgvs = re.findall(r'[A-Z]', self.file_entry.genomicHgvs)[0]
		refFromGenome = genome[chromosome][int(position)-1].upper()
		return vcf_utils.validate_reference({"refFromHgvs": refFromHgvs, 
			"refFromGenome": refFromGenome, "genomicHgvs": self.file_entry.genomicHgvs})

	def getAlteration(self):
		return re.findall(r'[A-Z]', self.file_entry.genomicHgvs)[1]

# Datamodel(object) that defines a Non-Substitution Variant
class NonSubstitutionVariant(Variant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

# Datamodel(object) that defines a Deletion Variant (child of Non-Substitution Variant)
class DeletionVariant(NonSubstitutionVariant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

	def getPosition(self):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		return int(positions[0]) - 1

	def getReference(self, chromosome):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		if len(positions) == 1:
			pos = int(positions[0]) - 1
			reference = genome[chromosome][pos-1:int(positions[0])].upper()
		else:
			start, end = [positions[0], positions[1]]
			pos = int(start) - 1
			reference = genome[chromosome][pos-1:int(end)].upper()

		return reference

	def getAlteration(self, chromosome, position, reference):
		return genome[chromosome][position-1].upper()

# Datamodel(object) that defines a Duplication Variant (child of Non-Substitution Variant)
class DuplicationVariant(NonSubstitutionVariant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

	def getPosition(self):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		return int(positions[0])

	def getReference(self, chromosome):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		refFromHgvs = ''.join(re.findall(r'[A-Z]+', self.file_entry.genomicHgvs))
		if len(positions) == 1:
			pos = int(positions[0])
			refFromGenome = genome[chromosome][pos-1].upper()
		else:
			start, end = [positions[0], positions[1]]
			pos = int(start)
			refFromGenome = genome[chromosome][pos-1:int(end)].upper()
		return vcf_utils.validate_reference({"refFromHgvs": refFromHgvs, 
			"refFromGenome": refFromGenome, "genomicHgvs": self.file_entry.genomicHgvs})

	def getAlteration(self, chromosome, position, reference):
		return reference * 2

# Datamodel(object) that defines a Delins Variant (child of Non-Substitution Variant)
class DelinsVariant(NonSubstitutionVariant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

	def getPosition(self):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		return int(positions[0])

	def getReference(self, chromosome):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		if len(positions) == 1:
			pos = int(positions[0])
			reference = genome[chromosome][pos-1].upper()
		else:
			start, end = [positions[0], positions[1]]
			pos = int(start)
			reference = genome[chromosome][pos-1:int(end)].upper()
		return reference

	def getAlteration(self, chromosome, position, reference):
		return ''.join(re.findall(r'[A-Z]+', self.file_entry.genomicHgvs))

# Datamodel(object) that defines a Insertion Variant (child of Non-Substitution Variant)
class InsertionVariant(NonSubstitutionVariant):
	def __init__(self, file_entry):
		super().__init__(file_entry)

	def getPosition(self):
		return int(re.findall(r'[0-9]+', self.file_entry.genomicHgvs)[0])

	def getReference(self, chromosome):
		positions = re.findall(r'[0-9]+', self.file_entry.genomicHgvs)
		pos = int(positions[0])
		reference = genome[chromosome][pos-1].upper()
		return reference

	def getAlteration(self, chromosome, position, reference):
		ins_bases = ''.join(re.findall(r'[A-Z]+', self.file_entry.genomicHgvs))
		return reference + ins_bases

# Datamodel(object) that defines a VCF entry to be written to a file
class VCFEntry():
	def __init__(self, chrom, pos, ref, alt):
		self.chrom = chrom
		self.pos = pos
		self.ref = ref
		self.alt = alt
		self.vcf_id = '.'
		self.qual = '.'
		self.vcf_filter = 'PASS'
		self.info = '.'
		self.vcf_format = 'GT:GQ:DP:SR:VR:VA:SB:ABQ:AMQ'
		self.sample = '1/1:1000.00:608:54.28:54.28:0:0.41:36.11:254.01'

	def toString(self):
		print('CHROM : {} ; POS : {}, ID : {}, REF : {}, ALT : {}, QUAL : {}, FILTER : {}, INFO : {}, FORMAT : {}, \
			SAMPLE : {}'.format(self.chrom, self.pos, self.vcf_id, self.ref, self.alt, self.qual, self.vcf_filter, self.info, 
			self.vcf_format, self.sample));

	def writeEntry(self):
		write_output = '{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{format}\t{sample}'.format(chrom=self.chrom, \
			pos=self.pos, id=self.vcf_id, ref=self.ref, alt=self.alt, qual=self.qual, filter=self.vcf_filter, info=self.info, \
			format=self.vcf_format, sample=self.sample)
		return write_output

# Datamodel(object) that defines the File Entry (an entry from the input file)
class FileEntry():
	def __init__(self, gene, genomicHgvs):
		self.gene = gene
		self.genomicHgvs = genomicHgvs

	def toString(self):
		print('Gene : {} ; Genomic HGVS : {}'.format(self.gene, self.genomicHgvs))

# Reads a file and creates each entry from the file into a FileEntry object
def read_file_and_get_file_objects(input_file):
	file = open(input_file, 'r')
	data = csv.DictReader(file, delimiter='\t')
	file_objects = list()
	for row in data:
		if row['Gene name'] and row['genomicHGVS']:
			file_object = FileEntry(row['Gene name'], row['genomicHGVS'])
			file_objects.append(file_object)
		else:
			print(messages['error_messages']['GENE_HGVS_NOT_PRESENT'].format(row['Gene name'], row['genomicHGVS']))
	return file_objects

# Segregates FileEntry objects into Substitution and Non-Substitution variants
def segregate_file_objects(file_objects):
	sub_objects, non_sub_objects = [list(), list()]
	for file_object in file_objects:
		if ">" in file_object.genomicHgvs:
			sub_objects.append(file_object)
		else:
			non_sub_objects.append(file_object)
	return sub_objects, non_sub_objects

# Categorize Non-Substitution variant FileEntry Objects into types - Deletion, Duplication, Insertion and Delins
def categorize_non_sub_objects(non_sub_objects):
	dups, dels, indels, ins = [list(), list(), list(), list()]
	for file_object in non_sub_objects:
		if 'dup' in file_object.genomicHgvs:
			dups.append(file_object)
		elif 'del' in file_object.genomicHgvs and 'ins' not in file_object.genomicHgvs:
			dels.append(file_object)
		elif 'ins' in file_object.genomicHgvs and 'del' not in file_object.genomicHgvs:
			ins.append(file_object)
		elif 'del' and 'ins' in file_object.genomicHgvs:
			indels.append(file_object)
		else:
			print(messages['error_messages']['UNKNOWN_OBJECT'].format(file_object.gene. file_object.genomicHgvs))

	print(messages['success_messages']['DUPLICATION_VARIANTS'].format(len(dups)))
	print(messages['success_messages']['DELETION_VARIANTS'].format(len(dels)))
	print(messages['success_messages']['INSERTION_VARIANTS'].format(len(ins)))
	print(messages['success_messages']['DELINS_VARIANTS'].format(len(indels)))

# Creates a list of SubstitutionVariant objects from FileEntry objects
def create_substitution_variants(sub_objects):
	sub_variants = list()
	for substitution in sub_objects:
		sub_variant = SubstitutionVariant(substitution)
		sub_variants.append(sub_variant)
	return sub_variants

# Creates a list of NonSubstitutionVariant objects from FileEntry Objects
def create_non_substitution_variants(non_sub_objects):
	non_sub_variants = list()
	for non_substitution in non_sub_objects:
		genomic_hgvs = non_substitution.genomicHgvs
		if "delins" in genomic_hgvs:
			non_sub_variant = DelinsVariant(non_substitution)
		elif "del" in genomic_hgvs and "ins" not in genomic_hgvs:
			non_sub_variant = DeletionVariant(non_substitution)
		elif "ins" in genomic_hgvs and "del" not in genomic_hgvs:
			non_sub_variant = InsertionVariant(non_substitution)
		elif "dup" in genomic_hgvs:
			non_sub_variant = DuplicationVariant(non_substitution)
		else:
			print(messages['error_messages']['UNKNOWN_OBJECT'].format(non_substitution.gene, non_substitution.genomicHgvs))
		non_sub_variants.append(non_sub_variant)
	return non_sub_variants

# Converts SubstitutionVariant objects into VCFEntry objects
def create_vcf_entries_from_sub_variants(sub_variants):
	sub_vcf_entries = list()
	for sub_variant in sub_variants:
		chromosome, position = [sub_variant.getChromosome(), sub_variant.getPosition()]
		sub_vcf_entry = VCFEntry(chromosome, position, sub_variant.getReference(chromosome, position), sub_variant.getAlteration())
		sub_vcf_entries.append(sub_vcf_entry)
	return sub_vcf_entries

# Converts NonSubstitutionVariant objects into VCFEntry objects
def create_vcf_entries_from_non_sub_variants(non_sub_variants):
	non_sub_vcf_entries = list()
	for non_sub_variant in non_sub_variants:
		chromosome = non_sub_variant.getChromosome()
		position = non_sub_variant.getPosition()
		reference = non_sub_variant.getReference(chromosome)
		alteration = non_sub_variant.getAlteration(chromosome, position, reference)
		non_sub_vcf_entry = VCFEntry(chromosome, position, reference, alteration)
		non_sub_vcf_entries.append(non_sub_vcf_entry)
	return non_sub_vcf_entries

# Writes out the VCFEntry objects into an output file
def write_vcf_entry_objects(vcf_entry_objects, output):
	output.write('\t'.join(json_data['vcf_header']))
	output.write('\n')
	for vcf_entry_object in vcf_entry_objects:
		output.write(vcf_entry_object.writeEntry())
		output.write('\n')

def main(input_file):
	print('\nINFO : File - {}'.format(input_file))

	file_objects = read_file_and_get_file_objects(input_file)
	(sub_objects, non_sub_objects) = segregate_file_objects(file_objects)
	print(messages['success_messages']['SUB_VARIANTS'].format(len(sub_objects)))
	print(messages['success_messages']['NON_SUB_VARIANTS'].format(len(non_sub_objects)))
	
	if non_sub_objects: categorize_non_sub_objects(non_sub_objects)
	if sub_objects: sub_variants = create_substitution_variants(sub_objects)
	if non_sub_objects: non_sub_variants = create_non_substitution_variants(non_sub_objects)

	if sub_variants: sub_vcf_entries = create_vcf_entries_from_sub_variants(sub_variants)
	if non_sub_variants : non_sub_vcf_entries = create_vcf_entries_from_non_sub_variants(non_sub_variants)

	output_file = os.path.join(json_data['output'], os.path.basename(input_file).replace('.tsv', '.vcf'))
	output = open(output_file, 'w+')
	if sub_vcf_entries: write_vcf_entry_objects(sub_vcf_entries, output)
	if non_sub_vcf_entries: write_vcf_entry_objects(non_sub_vcf_entries, output)
	print(messages['success_messages']['OUTPUT'].format(output_file))

if __name__ == '__main__':
	for entry in sys.argv[1:]:
		main(entry)