echo "INFO : Running Python Tests for common functions"
echo "INFO : Running tests for variant_functions.py script"
PYTHONPATH=../ python3 test_variant_functions.py

echo "INFO : Running tests for vcf_utils.py script"
PYTHONPATH=../ python3 test_vcf_utils.py

cd ..

if [[ -d __pycache__ ]]; then
	rm -rf __pycache__
fi