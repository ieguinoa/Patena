Patena


--------------------------------------------------------------------------

Insallation guide:

clone this repo: git clone 
source install.sh
DONE!

REQUIREMENTS: python, perl, biopython, BLAST(blastp in path and BLASTDB defined with the path to DB)(or run web BLAST which considerably increases running time) 


NOT TESTED ON WINDOWS!


usage:   python bleach.py  [options]


GENERAL OPTIONS:
	
	--seq [sequence]	Define initial sequence
	--length [seq-length]   Define initial sequence length (and generate a random sequence)


 
OUTPUT FORMATING:
	--verbose		Write detailed output in stdout
	--minoutput		Write mutations history log in file (patena.log)
	--stepped 		Write detailed output and wait for user input at each step



