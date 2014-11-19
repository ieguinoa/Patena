import urllib2
import StringIO
import sys
import os
import re
import time
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML
#from Bio import SearchIO
from array import array
import math
import random


#****************************************
# CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))






#*********************
#******GLOBALS*******
#******************

endl = "\n"
tab = "\t"
indent=""
cutoff=0.01
beta=0.1
match=False
rand=True
change=True
verbose=False
blastWeb=False  #BLAST SEARCH LOCAL OR WEB
blastIt=True
targetScore=0.0
output=True  ##print info to file
outputPath = "Output/"
mutationsIter=0
temperature=100
scoresFile="scores"   #save Scores Vs iteration number 
mutAttemptsFile="mutationsAttempt"  # save number of mutation attempts  Vs iteration number

#AA FREQUENCIES TO SELECT NEW RESIDUES FOR MUTATIONS (from http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html) 
weights= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]


#CLEAR ALL OUTPUT FILES
outputFile=outputPath + scoresFile
outFile=open(outputFile, "w")
outFile.close()
outputFile=outputPath + mutAttemptsFile
outFile=open(outputFile, "w")
outFile.close()
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
    global indent, mutationsIter
    #weighted is a pairlist  (position,weight)
    weighted=[]
    previousScore=getGlobalScore(mutationFreq)
    if previousScore > 0:
      mutationsIter=0
      for x in range(len(mutationFreq)):
	weighted.append((x, mutationFreq[x]+1))
      while 10000 > mutationsIter: 
	mutationsIter+=1
	#SELECT A POSITION 
	print endl
	print "*************************************"
	print "MUTATION ATTEMPT"
	print "*************************************"
	print "Score before:    " + str(previousScore)
	print "Choose a position based on weights"
	mutatePosition= weighted_choice(weighted)  
	print "Position chosen: " + str(mutatePosition)
	
	#SELECT THE NEW AA FOR THAT POSITION (BASED ON LIST OF FREQUENCIES)
	previousResidue=sequence[mutatePosition]
	print "Residue to mutate: " + previousResidue
	seleccionado= previousResidue
	
	#SELECT A NEWONE UNTIL THE RESIDUE IS DIFFERENT FROM PREVIOUS
	while previousResidue == seleccionado:
	    seleccionado = weighted_choice(weights)	
	    print "Selected residue : " + seleccionado
	    
	##CREATE MUTATED SEQUENCE WITH NEW RESIDUE    
	mutatedSequence = sequence[0:mutatePosition]
	mutatedSequence += seleccionado
	mutatedSequence += sequence[mutatePosition+1:]
	print "Secuencia original:" + sequence
	print "Secuencia mutada  :" + mutatedSequence
	#CREATE A LIST OF SCORES FOR THE MUTATED SEQUENCE
	mutatedFreq=[]
	for p in range(len(sequence)):
	    #mutatedFreq.append(0)
	    mutationFreq[p]=0
	#ACCEPT OR DENY THE MUTATION BASED ON A NEW EVALUATION OF THE SCORE
	indent=tab + tab
	sequenceEvaluation(mutatedSequence, mutationFreq, True)	 
	#IF THE GLOBAL SCORE DECREASED
	print "*************************************"
	print "MUTATION RESULTS"
	if previousScore > getGlobalScore(mutationFreq):
	    print "Previous score (" + str(previousScore) + ") > Mutation score (" + str(getGlobalScore(mutationFreq)) + ")" 
	    print "...ACCEPT MUTATION"
	    #mutationFreq=
	    return mutatedSequence
	else:
	    #DECISION BASED ON MONTE CARLO
	    diff=previousScore-getGlobalScore(mutationFreq)
	    print "LA DIFERENCIA ENTRE VALORES DE SCORE ES" + str(diff)
	    exponent=diff/beta
	    MCvalue=math.exp(exponent)
	    print "EL EXPONENTE ES  :" + str(exponent)
	    print "EL VALOR DE MC ES:" + str(MCvalue)
	    
	    #GENERATE RANDOM NUMBER BETWEEN 0 AND 1
	    randy=random.random()
	    print "EL VALOR RANDOM ELEGIDO ES:" + str(randy)
	    print "MONTE CARLO DECISION:"
	    if MCvalue > randy:
	      #ACCEPT MUTATION
	      print "...ACCEPT MUTATION"
	      return mutatedSequence
	    else:	    
	      #print "El score original " + str(getGlobalScore(mutatedFreq)) + " es mayor que "+ str(getGlobalScore(mutationFreq))
	      print " Mutation score (" + str(getGlobalScore(mutationFreq)) + ") >= Previous score (" + str(previousScore) + ")" 
	      print "...DENY MUTATION"
	      #return sequence
	print "*************************************"
      return sequence
    else:
      print "SCORE = 0 ******** QUIT MUTATION STEP"
    return sequence



###********************************************


  ######################################################################################
  ##########################       SEARCH ELMS     #####################################
  ######################################################################################




def makeElmSearch(sequence):
  ##MAKE THE SEARCH AND SAVE RESULTS IN outputELM
  elm_pattern_dict = {}
  with open("elm_patterns_20140701.txt", 'rU') as file_open :
	  patterns = file_open.readlines()
	  for line in patterns :
		  line = line.split("\t")
		  elm_id = line[0]
		  elm_pattern = line[1].replace("\n","")
		  elm_pattern_dict[elm_id] = elm_pattern

  #print "SEARCHING ELMS IN %s\n ..." % filename

  #output_file_name = "elms_search_in_%s.txt" % filename[:-4]
  output_file_name="outputELM"
  uniprot_list = []
  sequence_dict = {}

  #with open(filename, 'rU') as file_open :
	  #my_seq = file_open.read()
  
  with open(output_file_name, 'w') as file_write :
	  for elm_id in elm_pattern_dict :
		  where_to_start = []
		  elm_pos_dict = {}
		  pattern = re.compile('(?=%s)' % elm_pattern_dict[elm_id])
		  for matched_string in pattern.finditer('%s' % sequence) :
			  where_to_start.append(matched_string.start())
		  pattern = re.compile(elm_pattern_dict[elm_id])
		  for index in where_to_start :
			  match = re.search(pattern, '%s' % sequence[index:])
			  if match != None :
				  file_write.write("%s\t%s\n" % (index+1, index+len(match.group())))



##*************************************


def elmSearch(sequence, mutationFreq,verbose):
  makeElmSearch(sequence)	
  
  #FIRST SAVE HITS IN A LIST 
  elmSearchFreq=[]
  for p in range(len(sequence)):
	      elmSearchFreq.append(0)
  inputFile="outputELM"
  with open(inputFile, "r") as input_file:
    lines=input_file.readlines()
  for line in lines:
    line=line.split("\t")
    pattern_start=int(line[0])
    pattern_end=int(line[1])
    for x in range(pattern_start-1,pattern_end):
      elmSearchFreq[x]+=1

  if verbose:
    #print endl
    print indent + "ELM Search RESULTS:"
    print indent + sequence
    print indent + ''.join(map(str, elmSearchFreq))	  
  
  ##ADD hits to global score
  for i in range(len(sequence)):
    mutationFreq[i]+=elmSearchFreq[i]
  
  
                                         


##**********************************************

##getGlobalScore=List sum
def getGlobalScore(mutatFreqList):
  score=0.0
  for listIndex in range(len(mutatFreqList)):
    score= score + mutatFreqList[listIndex]
  return score


##***************************************


def blastIt(sequence, mutationFreq, database, verbose):
        global match
        ##BLAST SEARCH
        
        if blastWeb:       # WEB BLAST SEARCH
	  if verbose:
	    print indent + "WEB BLAST SEARCH IN PROGRESS..." 
	  result = NCBIWWW.qblast("blastp", database , sequence)
	  records = NCBIXML.parse(result)
	  first = records.next()
	else:     # LOCAL BLAST SEARCH
	  if verbose:
	    print indent + "LOCAL BLAST SEARCH IN PROGRESS..."
	  input=open('input', 'w')
	  #input.write("kqltqdddtdeveiaidntafmdeffseie\n")
	  input.write(sequence)
	  input.close()
	  commandLine=NcbiblastpCommandline(query="input", db=database, evalue=0.001, outfmt=5, out="output.xml")
	  #print commandLine
	  stdout, stderr = commandLine()
	  result_handle = open("output.xml")
	  blast_records = NCBIXML.parse(result_handle)
	  first = blast_records.next() 
	  
	
	
	#first.alignments contains all de alignments found
	if len(first.alignments) > 0:
	#get first alignment
	  firstAlign=first.alignments[0]
	  #print endl
	  
	  #print alignment stats 
	  if verbose:
	    print indent +"Cutoff:" + str(cutoff)
	  for hsp in firstAlign.hsps:
	    if hsp.expect < cutoff:
	      match=True   #we have a match
	      if verbose:
		print indent + "****Alignment****"  
		print indent + "sequence:", firstAlign.title
	      
	      #length of the alignment (could be shorter than full sequence)
	      length=firstAlign.length
	      
	      #starting position of alignment in the sequence
	      start=hsp.query_start	
	      
	      #ending position of the alignment in the sequence
	      end=hsp.query_end
	      
	      #length = (end-start) ???
	      if verbose:
		print indent + "E-Value:     " + str(hsp.expect)
		print indent + "Query:       " + hsp.query 
		print indent + "Match:       " + hsp.match
		print indent + "Subject:     " + hsp.sbjct 
		print indent + "Query Length:", len(sequence)
		print indent + "Query Start: ", hsp.query_start
		print indent + "Query end:   ", hsp.query_end
	    else:
	      if verbose:
		print indent + "No hits found"
	      match=False
	else:
	  if verbose:
	    print indent + "No hits found"
	  match=False

	    
	if match:
	#TAKE THE mutationFreq LIST AND PUT "1" WHERE THERE WAS A MATCH AND "0" WHERE THERE WAS A GAP OR MISMATCH
		for j in range(len(sequence)):
			if j< (start-1) or j > (end-1):    
				#print sequence[j]
				mutationFreq[j]=1
			else:
				if hsp.match[j-start+1] <> "+" and hsp.match[j-start+1] <> " ":
					mutationFreq[j] = 1
				else:
					mutationFreq[j]=0
					
	if verbose:	    
	    #print endl
	    print indent + "BLAST RESULTS:"
	    print indent + sequence
	    print indent + ''.join(map(str, mutationFreq))				
					
					
					
					

def iupred(sequence, mutationFreq, verbose):
	runCommand="iupred/iupredExe iupred/input long"
	input=open('iupred/input', 'w')
	input.write("Name" + endl)
	input.write(sequence)
	input.close()
	os.system(runCommand)	
	outputIUPred=open("salida", "r")
	
	#PRINT THE RESULTS OF IUPred
	if verbose:
	  iupredFreq=[]
	  iterOutputIUPred=iter(outputIUPred)
	  for p in range(len(sequence)):
	      iupredFreq.append(0)
	  for x in range(len(sequence)):
		  resultX=float(iterOutputIUPred.next())
		  if resultX < 0.5 :
			  iupredFreq[x] = 1
	  #print endl
	  print indent + "IUPred RESULTS:"
	  print indent + sequence
	  print indent + ''.join(map(str, iupredFreq))	  			  			
			
	outputIUPred.seek(0)
	rstFile_iter = iter(outputIUPred)
	#ADD 1 TO THE POSITION IN mutationFreq IF THE RESULT IS LESS THAN 0.5 (PREDICTING A GLOBULAR TENDENCY)
	for j in range(len(sequence)):
		resultJ=float(rstFile_iter.next())
		if resultJ > 0.5 :
			mutationFreq[j] += 0
		else:
			mutationFreq[j] += 1				
  


def anchor(sequence, mutationFreq, verbose):
	runCommand="ANCHOR/anchorExe ANCHOR/input"
	input=open('ANCHOR/input', 'w')
	input.write("Name" + endl)
	input.write(sequence)
	input.close()
	os.system(runCommand)	
	outputAnchor=open("outAnchor", "r")
	#print "aca evaluo anchor"
	#PRINT THE RESULTS OF ANCHOR
	if verbose:
	  anchorFreq=[]
	  iterOutputAnchor=iter(outputAnchor)
	  for p in range(len(sequence)):
	      anchorFreq.append(0)
	  for x in range(len(sequence)):
		  resultX=float(iterOutputAnchor.next())
		  if resultX >  0.5 :
			  anchorFreq[x] = 1
	  print indent + "ANCHOR RESULTS:"
	  print indent + sequence
	  print indent + ''.join(map(str, anchorFreq))	  			  			
			
	outputAnchor.seek(0)
	rstFile_iter = iter(outputAnchor)

	for j in range(len(sequence)):
		resultJ=float(rstFile_iter.next())
		if resultJ > 0.5 :
			mutationFreq[j] += 1
		else:
			mutationFreq[j] += 0				
  










def sequenceEvaluation(sequence, mutationFreq, verbose):
	#global indent
	
	##FIRST STEP: BLAST SEARCH
	if verbose:
	   print endl
	   print indent + "*************************************"
	   print indent + "STARTING BLAST SEARCH"
	blastIt(sequence,mutationFreq,database, verbose)
        
	    
        ##SECOND STEP: IUPred evaluation
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "STARTING IUPred"
	iupred(sequence, mutationFreq, verbose)
        
        ##THIRD STEP: ANCHOR evaluation
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "STARTING ANCHOR"
	anchor(sequence, mutationFreq, verbose)
	
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "STARTING ELM Search"
        elmSearch(sequence,mutationFreq, verbose)
	
        ##PRINT SCORE
        if verbose:
	  print indent + "*************************************"
	  print endl
	  print indent + "FINAL RESULTS:"
	  print indent + sequence
	  print indent + ''.join(map(str, mutationFreq))
	  print indent + "*************************************"
	



#**************************
#******* MAIN *************
#**************************


#********DEFAULTS********
size=10
composition="average"
a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=0
database="uniprot_sprot.fasta"



#*******PROCESO LOS ARGUMENTOS******
if len(sys.argv) < 2:
  print "************************************************"
  print "************************************************"
  print "USING DEFAULT VALUES: lenght=10  - composition=average - sequence=RANDOM - db=swissprot"
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
    elif (arg=='--blastweb'):
      blastWeb=True
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
#http://web.expasy.org/cgibin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

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
mutationFreq=[]
for p in range(len(sequence)):
  mutationFreq.append(0)
  #mutationFreq[p]=0



t0 = time.time()
################################
#### ITERATE OVER SEQUENCE ####
#############################
globalScore=10
iteration=1
#for iteration in range(20):
while globalScore > 0:
	#print endl
        print "*****************************"
	print "ITERATION NUMBER: " + str(iteration)
	print "*****************************"
	#BEFORE EACH ITERATION STEP, CLEAN M	UTATION FREQUENCE
	#THIS LIST SUMS UP THE PROBABILITY TO MUTATE EACH POSITION BASED ON ALL STEPS PERFORMED (SUMMING IT UP GIVES THE GLOBAL SCORE OF THE SEQUENCE)
	#mutationFreq=[]
	for p in range(len(sequence)):
           #mutationFreq.append(0)
	   mutationFreq[p]=0
	
	indent=""
	sequenceEvaluation(sequence,mutationFreq, True)	
	


 
	#AFTER ALL STEPS, MUTATE SEQUENCE
	
	###ONLY MUTATE IF SCORE GREATER THAN 0  (AT LEAST 1 POSITION HAS A MUTATION FREQUENCY GREATER THAN 0 )
	print endl
	print "*************************************"
	print "*************************************"
	print "MUTATION STEP"
	print "*************************************"
	print "*************************************"
	print "Sequence before: " + sequence 
	mutatedSequence=mutar(sequence, mutationFreq)
	print "Sequence after mutation:    " + mutatedSequence
	sequence=mutatedSequence
	print "Number of times before mutation accept:" + str(mutationsIter)
	print endl
	print "*******************************************"
      
        
        #else:     ##GLOBAL SCORE IS 0 - REACHED END OF LOOP
	#    break
	globalScore= getGlobalScore(mutationFreq)
	print "End of iteration: " + str(iteration)
	print "Global score :    " + str(globalScore)
	print endl
	#print "Score function: " +  str(getGlobalScore(mutationFreq)/len(sequence)) 
	#if ((getGlobalScore(mutationFreq)/len(sequence)) <= targetScore):
	#  break
	#Write score Vs iteration 
	outputFile=outputPath + scoresFile
	outFile=open(outputFile, "a")
	outFile.write(str(iteration)+ tab + str(globalScore) + endl)
	#Write mutations attempts
	outFile.close()
	outputFile=outputPath + mutAttemptsFile
	outFile=open(outputFile, "a")
	outFile.write(str(iteration)+ tab + str(mutationsIter) + endl)
        outFile.close()
        iteration=iteration+1
elapsedTime=time.time() - t0
print "Elapsed time:",elapsedTime, "Seconds"
if output:
  outputFile= outputPath + str(len(sequence)) 
  outFile=open(outputFile, "a")
  outFile.write(str(iteration-1)+ tab + str(elapsedTime) + endl)
   

