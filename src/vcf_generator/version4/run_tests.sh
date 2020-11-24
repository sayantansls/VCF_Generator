echo "INFO : Running tests for VCF generator version 4"
PYTHONPATH=../../common/ python3 test_vcf_generator.py

if [[ -f vcf_generator_ver4.pyc ]]; then
	rm vcf_generator_ver4.pyc
fi

if [[ -d __pycache__ ]]; then
	rm -rf __pycache__
fi