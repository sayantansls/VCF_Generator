cd src/vcf_generator/version3/
PYTHONPATH=../../common/ python3 vcf_generator_ver3.py

echo 'INFO: output file present at location - data/output/'

cd ../../common/

if [[ -f vcf_utils.pyc ]]; then
	rm vcf_utils.pyc
fi

if [[ -f variant_functions.pyc ]]; then
	rm variant_functions.pyc
fi

if [[ -d __pycache__ ]]; then
	rm -rf __pycache__
fi