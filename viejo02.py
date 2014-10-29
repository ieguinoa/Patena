import urllib2
import StringIO
import sys
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML
#from Bio import SearchIO
from array import array

import random




#****************************************
# CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))






#*********************
#******GLOBALS*******
#******************

endl = "\n"
tab = "\t"
cutoff=0.01
match=False
rand=True
change=True
verbose=False


#AA FREQUENCIES TO SELECT NEW RESIDUES FOR MUTATIONS (from http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html) 
weights= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]



#***********************************************************************



## PRINT PROGRAM HELP 
def printHelp():
  print endl + "Usage: "
  print "  python linkeado.py [options]"
  print endl + "   where options are:"
  print tab + "--length :" + tab + "Sequence lenght"
  print tab + "--db :" + tab + "swissprot | nr"
  print tab + "--composition :" + tab + "average | user_specified"
  print tab + "--seq :" + tab + "Initial sequence (not random)"
  
  
  
  
#***************************************************************************  
  
  
  
## AA MUTATION OF SEQUENCE - PROBABILITY OF MUTATION BASED ON WEIGHTS IN PONDERADA 
## PROBABILITY OF MUTATION BASED ON WEIGHTS ARRAY
## RETURNS MUTATED SEQUENCE
def mutar(sequence, mutationFreq):
    
    #weighted is a pairlist  (position,weight)
    weighted=[]
    
    #SELECT A POSITION 
    for x in range(len(mutationFreq)):
	weighted.append((x, mutationFreq[x]))
    print "Choose a position based on weights"
    mutatePosition= weighted_choice(weighted)  
    print "Mutate position: " + str(mutatePosition)
    
    #SELECT THE NEW AA FOR THAT POSITION (BASED ON LIST OF FREQUENCIES)
    
    previousResidue=sequence[mutatePosition]
    print "Mutate residue: " + previousResidue
    seleccionado= previousResidue
    
    #SELECT A NEWONE UNTIL THE RESIDUE IS DIFFERENT FROM PREVIOUS
    while previousResidue == seleccionado:
    	seleccionado = weighted_choice(weights)	
        print "Selected residue : " + seleccionado
        
    ##CREATE MUTATED SEQUENCE WITH NEW RESIDUE    
    mutatedSequence = sequence[0:mutatePosition]
    mutatedSequence += seleccionado
    mutatedSequence += sequence[mutatePosition+1:] 
    return mutatedSequence





##**********************************************

def getGlobalScore(mutationFreq):
  globalScore=0 
  for listIndex in range(len(mutationFreq)):
    globalScore=globalScore + mutationFreq[listIndex]
  return globalScore


##***************************************


def blastIt(sequence, mutationFreq, database):
        global match
        ##BLAST SEARCH
        result = NCBIWWW.qblast("blastp", database , sequence)
	records = NCBIXML.parse(result)
	first = records.next()


	#first.alignments contains all de alignments found
	if len(first.alignments) > 0:
	#get first alignment
	  firstAlign=first.alignments[0]
	  #print endl
	  
	  #print alignment stats 
	  print "Cutoff:" + str(cutoff)
	  for hsp in firstAlign.hsps:
	    if hsp.expect < cutoff:
	      match=True   #we have a matchi
	      print "****Alignment****"  
	      print "sequence:", firstAlign.title
	      
	      #length of the alignment (could be shorter than full sequence)
	      length=firstAlign.length
	      
	      #starting position of alignment in the sequence
	      start=hsp.query_start	
	      
	      #ending position of the alignment in the sequence
	      end=hsp.query_end
	      
	      #length = (end-start) ???
	     
	      print "E-Value:" + str(hsp.expect)
	      print "Query:   " + hsp.query 
	      print "Match:   " + hsp.match
	      print "Subject: " + hsp.sbjct 
	      print "Query Length:", len(sequence)
	      print "Query Start:", hsp.query_start
	      print "Query end:", hsp.query_end
	    else:
	      print "No hits found"
	      match=False
	else:
	  print "No hits found"
	  match=False


	if match:
	#TAKE THE mutationFreq LIST AND PUT "1" WHERE THERE WAS A MATCH AND "0" WHERE THERE WAS A GAP OR MISMATCH
		for j in range(len(sequence)):
			if j< (start-1) or j > (end-1):    
				#print sequence[j]
				mutationFreq[j]=1
			else:
				if hsp.match[j-start+1] <> "+" and hsp.match[j-start+1] <> " ":
					mutationFreq[j]=1
				else:
					mutationFreq[j]=0
					
					
					
					
					
					
					
					
					
  


#**************************
#******* MAIN *************
#**************************


#********DEFAULTS********
size=10
composition="average"
a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=0
database="swissprot"



#*******PROCESO LOS ARGUMENTOS******
if len(sys.argv) < 2:
  print "************************************************"
  print "************************************************"
  print "USING DEFAULT VALUES: length=10  - composition=average - sequence=RANDOM - db=swissprot"
  print "If you want a specific configuration, see the help using  --help"
  print "************************************************"
  print "************************************************"

else:
  for index in range(1,len(sys.argv)):
    arg = sys.argv[index]
    
    if (arg=='-H' or arg== '-help' or arg== '--help'):
      printHelp()
      exit()
    elif (arg=='--length') and (index < len(sys.argv)):
      size = int(sys.argv[index+1])
    elif (arg=='--seq') and (index < len(sys.argv)):
      sequence = sys.argv[index+1]
      rand=False
    elif (arg=='--db') and (index < len(sys.argv)):
      database = sys.argv[index+1]
    elif (arg=='--verbose'):
      verbose=True
    elif (arg== '--composition') and (index < len(sys.argv)):
      composition = sys.argv[index+1]
      if (composition=="user_specified"):  #DEBERIA RECIBIR LAS FRECUENCIAS DE TODOS LOS AMINOACIDOS
	for j in range(index+2,len(sys.argv),2):
	  if(sys.argv[j]=='-a'):
	    a=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-r'):
	    r=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-n'):
	    n=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-d'):
	    d=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-c'):
	    c=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-q'):
	    q=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-e'):
	    e=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-g'):
	    g=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-h'):
	    h=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-i'):
	    i=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-l'):
	    l=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-k'):
	    k=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-m'):
	    m=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-f'):
	    f=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-p'):
	    p=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-s'):
	    s=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-t'):
	    t=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-w'):
	    w=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-y'):
	    y=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-v'):
	    v=int(sys.argv[j+1])
	    


# FORMAT TO REQUEST RANDOM SEQUENCE:   
#  http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

if rand==True:
	#****************GET RANDOM SEQUENCE*************
        #print endl
	print "Generating random sequence..."    
	#print endl
	if not (composition=="user_specified"):
	  url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(size) + "&comp=" + composition + "&output=fasta"
	else:
	  print "fix this"
	response = urllib2.urlopen(url)
	html = response.read()
	i = html.index('\n')
	sequence = html[i+1:].replace('\n', '')
	print "*******************************" 
	

   
#print endl
#PRINT STARTING SEQUENCE (BEFORE ANY MUTATION)   
print "INITIAL SEQUENCE:" + sequence


#CREATE ARRAY TO SAVE MUTATION FREQUENCE
#mutationFreq=[]






################################
#### ITERATE OVER SEQUENCE ####
#############################


for iteration in range(2):
	#print endl
        print "*****************************"
	print "ITERATION NUMBER: " + str(iteration+1)
	print "*****************************"
	#BEFORE EACH ITERATION STEP, CLEAN M	UTATION FREQUENCE
	#THIS LIST SUMS UP THE PROBABILITY TO MUTATE EACH POSITION BASED ON ALL STEPS PERFORMED (SUMMING IT UP GIVES THE GLOBAL SCORE OF THE SEQUENCE)
	mutationFreq=[]
	for p in range(len(sequence)):
           mutationFreq.append(0)
	
	
	
	
	#########   FIRST STEP:  BLAST SEARCH  ######
	
	#print endl
	print "*************************************"
	print "STARTING BLAST SEARCH"
#	print endl
	blastIt(sequence,mutationFreq,database)
	print endl
	## NOW I HAVE mutationFreq ARRAY WITH 1s AND 0s
	#
#	print endl
	print "BLAST RESULTS:"
	print sequence
	print ''.join(map(str, mutationFreq))
 #       print endl
#        print "*************************************"



	##SECOND STEP: ????????
	#EXECUTE SECOND STEP AND INCREASE mutationFreq VALUES
	


 
	#AFTER ALL STEPS, MUTATE SEQUENCE
	
	###ONLY MUTATE IF SCORE GREATER THAN 0  (AT LEAST 1 POSITION HAS A MUTATION FREQUENCY GREATER THAN 0 )
	if getGlobalScore(mutationFreq) > 0:
	  print endl
	  print "*************************************"
	  print "MUTATION STEP"
	  print "Sequence before:            " + sequence 
	  mutatedSequence=mutar(sequence, mutationFreq)
	  print "Sequence after mutation:    " + mutatedSequence
	  sequence=mutatedSequence
	  print endl
	  print "*******************************************"




