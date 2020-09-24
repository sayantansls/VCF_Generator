if [[ ! -f vcf_generator_ver4.py ]]; then
	echo "ERROR : vcf_generator_ver4.py script is not present"
	exit 0
fi

PYTHONPATH=../../common python3 vcf_generator_ver4.py $@