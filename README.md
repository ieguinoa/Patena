PATENA


--------------------------------------------------------------------------

Installation guide:

clone this repo: git clone https://github.com/ieguinoa/patena \\
install: source install.sh \\
DONE!

REQUIREMENTS: 
	-python
	-perl
	-biopython
	-BLAST: PATENA searches for possible functional elements in the sequences based on sequence homology evidence, for this it uses a BLAST search against a DB of known protein sequences. blastp should be in PATH and the env variable BLASTDB must be defined as the path to blast DBs dir. The DB of proteins to search against should preferentialy be UniprotKB. The name of the DB is (by default) uniprot_sprot.fasta, although this can be changed using parameter --db [DB_NAME]. An alternative option to the local blast is to run web BLAST (--blastweb), although this increases considerably the running time.


NOT TESTED ON WINDOWS!


USAGE:   
```
	python patena.py  [options]
```

OPTIONS:
	
	--seq [sequence]	Define initial sequence
	--length [seq-length]   Define initial sequence length (and generate a random sequence)


 
OUTPUT FORMATING:
```
	--verbose		Write detailed output in stdout
	--logoutput		Write mutations history log in file (patena.log)
	--stepped 		Write detailed output and wait for user input at each step
```



LOG FILE FORMAT:
	Each line has 4 columns. mutation, new sequence, score, loop id[1 | 2]

     

