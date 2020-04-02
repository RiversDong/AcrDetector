# Notice
    Do not change and delete any files in AcrFinder

# Requirement
    HMMER 3.3
    python 3.7
    pandas
    numpy
    joblib
    biopython
    scikit-learn
    
    Please install HMMER and add it to your PATH
    Please install packages listed above before using this script

# Running
    if you have add the script to your PATH, you can run it like
    python AcrFinder.py -i [or --infile] <infile> -o [or --output] <output>
    if the script don't in PATH, run it using the path of AcrFinder.py plus
    AcrFinder.py -i [or --infile] <infile> -o [or --output] <output>
    
    -i or --infile: specify the input file
    <infile>: file storing all CDS of whole genome in a fasta format. See the format in ./examples path
    -o or --outfile: specify the output file
    <output>: output filename
    
    For example
    python AcrFinder.py -i ./examples/GCF_001188875.1.fna -o re
    This command means seaching Acrs in GCF_001188875.1.fna and 
    The results will be stored in re*

# output
    After running python AcrFinder.py -i ./examples/GCF_001188875.1.fna -o re
    Four files prefixed with "re" are formed
        re, re_domain, re_protein and re_tbl
    re: potential Acrs
    re_domain and re_tbl: hth domain. hmmscan results in different formats
    re_protein: translated protein from CDS region

	There are 7 colums
	ID	length	function("1" means unknow function)	codon distance	dev	hth("1" hth in the downstream)

