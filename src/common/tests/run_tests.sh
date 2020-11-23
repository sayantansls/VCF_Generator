echo "Running Python Tests for common functions"
echo "Running tests for variant_functions.py script"
PYTHONPATH=../ python3 test_variant_functions.py

echo "Running tests for vcf_utils.py script"
PYTHONPATH=../ python3 test_vcf_utils.py