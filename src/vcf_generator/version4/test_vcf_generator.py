"""
author : sayantan (sayantan.ghosh@strandls.com)
Unit testing for classes and functions in the script - vcf_generator_ver4.py
"""

import unittest
import vcf_generator_ver4
import vcf_utils

test_genomic_build = '37'
test_entry = {
	"gene": "BRCA1", 
	"genomicHgvs": "g.41215907C>T", 
	"chrom": "chr17",
	"position": 41215907,
	"ref": "C",
	"alt": "T"
}

def test_non_sub_variant_object(self, test_variant, test_entry):
	self.assertIsNotNone(test_variant)
	self.assertIsInstance(test_variant, vcf_generator_ver4.NonSubstitutionVariant)
	reference = test_variant.getReference(test_variant.getChromosome())
	self.assertEqual(test_variant.getPosition(), test_entry["position"])
	self.assertEqual(reference, test_entry["ref"])
	self.assertEqual(test_variant.getAlteration(test_variant.getChromosome(), 
		test_variant.getPosition(), reference), test_entry["alt"])

def create_file_objects(test_entries):
	test_file_objects = list()
	for test_entry in test_entries:
		test_file_objects.append(vcf_generator_ver4.FileEntry(test_entry["gene"], test_entry["genomicHgvs"]))
	return test_file_objects

class TestVcfGenerator(unittest.TestCase):

	def test_file_entry_object(self):
		test_file_entry = vcf_generator_ver4.FileEntry(test_entry["gene"], test_entry["genomicHgvs"])
		self.assertIsNotNone(test_file_entry)
		self.assertIsInstance(test_file_entry, vcf_generator_ver4.FileEntry)
		self.assertEqual(test_file_entry.gene, test_entry["gene"])
		self.assertEqual(test_file_entry.genomicHgvs, test_entry["genomicHgvs"])

	def test_variant_object(self):
		test_variant = vcf_generator_ver4.Variant(vcf_generator_ver4.FileEntry(test_entry["gene"], test_entry["genomicHgvs"]))
		self.assertIsNotNone(test_variant)
		self.assertIsInstance(test_variant, vcf_generator_ver4.Variant)
		self.assertEqual(test_variant.getChromosome(), test_entry["chrom"])

	def test_substitution_variant_object(self):
		test_file_entry = vcf_generator_ver4.FileEntry(test_entry["gene"], test_entry["genomicHgvs"])
		test_sub_variant = vcf_generator_ver4.SubstitutionVariant(test_file_entry)
		self.assertIsNotNone(test_sub_variant)
		self.assertIsInstance(test_sub_variant, vcf_generator_ver4.SubstitutionVariant)
		self.assertEqual(test_sub_variant.getPosition(), test_entry["position"])
		self.assertEqual(test_sub_variant.getReference(test_sub_variant.getChromosome(), test_entry["position"]), test_entry["ref"])
		self.assertEqual(test_sub_variant.getAlteration(), test_entry["alt"])

	def test_deletion_variant_object(self):
		test_del_entry = {
			"gene": "BRCA1",
			"genomicHgvs": "g.41209107delG",
			"position": 41209106,
			"ref": "TG",
			"alt": "T"
		}
		test_file_entry = vcf_generator_ver4.FileEntry(test_del_entry["gene"], test_del_entry["genomicHgvs"])
		test_del_variant = vcf_generator_ver4.DeletionVariant(test_file_entry)
		self.assertIsInstance(test_del_variant, vcf_generator_ver4.DeletionVariant)
		test_non_sub_variant_object(self, test_del_variant, test_del_entry)

	def test_insertion_variant_object(self):
		test_ins_entry = {
			"gene": "BRCA1",
			"genomicHgvs": "g.41244808_41244809insA",
			"position": 41244808,
			"ref": "C",
			"alt": "CA"
		}
		test_file_entry = vcf_generator_ver4.FileEntry(test_ins_entry["gene"], test_ins_entry["genomicHgvs"])
		test_ins_variant = vcf_generator_ver4.InsertionVariant(test_file_entry)
		self.assertIsInstance(test_ins_variant, vcf_generator_ver4.InsertionVariant)
		test_non_sub_variant_object(self, test_ins_variant, test_ins_entry)

	def test_duplication_variant_object(self):
		test_dup_entry = {
			"gene": "CFHR5",
			"genomicHgvs": "g.196963265dupA",
			"position": 196963265,
			"ref": "A",
			"alt": "AA"
		}
		test_file_entry = vcf_generator_ver4.FileEntry(test_dup_entry["gene"], test_dup_entry["genomicHgvs"])
		test_dup_variant = vcf_generator_ver4.DuplicationVariant(test_file_entry)
		self.assertIsInstance(test_dup_variant, vcf_generator_ver4.DuplicationVariant)
		test_non_sub_variant_object(self, test_dup_variant, test_dup_entry)

	def test_delins_variant_object(self):
		test_delins_entry = {
			"gene": "APC",
			"genomicHgvs": "g.112174576_112174577delinsAT",
			"position": 112174576,
			"ref": "GC",
			"alt": "AT"
		}
		test_file_entry = vcf_generator_ver4.FileEntry(test_delins_entry["gene"], test_delins_entry["genomicHgvs"])
		test_delins_variant = vcf_generator_ver4.DelinsVariant(test_file_entry)
		self.assertIsInstance(test_delins_variant, vcf_generator_ver4.DelinsVariant)
		test_non_sub_variant_object(self, test_delins_variant, test_delins_entry)		

	def test_create_substitution_variants(self):
		test_entries = [
			{"gene": "BRCA1", "genomicHgvs": "g.41215907C>T"}, 
			{"gene": "BRCA2", "genomicHgvs": "g.32929233C>G"}]
		test_file_objects = create_file_objects(test_entries)
		sub_variant_objects = vcf_generator_ver4.create_substitution_variants(test_file_objects)
		self.assertIsNotNone(sub_variant_objects)
		self.assertEqual(len(sub_variant_objects), len(test_entries))
		for index in range(len(sub_variant_objects)):
			self.assertIsInstance(sub_variant_objects[index], vcf_generator_ver4.SubstitutionVariant)
			self.assertIs(sub_variant_objects[index].file_entry, test_file_objects[index])

	def test_create_non_substitution_variants(self):
		test_entries = [ 
				{"gene": "BRCA2", "genomicHgvs": "g.32968884_32968885delTT"},
				{"gene": "BRCA1", "genomicHgvs": "g.41244808_41244809insA"}, 
				{"gene": "APC", "genomicHgvs": "g.g.112174576_112174577delinsAT"}]
		test_file_objects = create_file_objects(test_entries)
		non_sub_variant_objects = vcf_generator_ver4.create_non_substitution_variants(test_file_objects)
		self.assertIsNotNone(non_sub_variant_objects)
		self.assertEqual(len(non_sub_variant_objects), len(test_entries))
		for index in range(len(non_sub_variant_objects)):
			self.assertIsInstance(non_sub_variant_objects[index], vcf_generator_ver4.NonSubstitutionVariant)
			self.assertIs(non_sub_variant_objects[index].file_entry, test_file_objects[index])
		self.assertIsInstance(non_sub_variant_objects[0], vcf_generator_ver4.DeletionVariant)
		self.assertIsInstance(non_sub_variant_objects[1], vcf_generator_ver4.InsertionVariant)
		self.assertIsInstance(non_sub_variant_objects[2], vcf_generator_ver4.DelinsVariant)

if __name__ == '__main__':
	unittest.main()