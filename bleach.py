import urllib2
import StringIO
import sys
import os
import re
import time
import errno
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML
#from Bio import SearchIO
from array import array
import math
import random
import subprocess




#****************************************
# CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))






#*********************
#******GLOBALS*******
#******************
maxIterations=10000
exeId=os.getpid()     #USE THIS TO IDENTIFY 
endl = "\n"
tab = "\t"
space=" "
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
mutAttempts=0
stepByStep=False
toolsPath=""   #SET THE PATH TO THE TOOL SET 
inputsPath="Input/"+ str(exeId) + "/" #SET PATH TO SAVE INPUTS FILES
outputsPath="Output/" + str(exeId) + "/"
try:
    os.makedirs("Input")
except OSError as exc: 
    if exc.errno == errno.EEXIST and os.path.isdir("Input"):
        pass




#AA FREQUENCIES TO SELECT NEW RESIDUES FOR MUTATIONS (from http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html) 
#aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]



  
#***************************************************************************  
  
  
def getGlobalScore(scoresList):
  score=0.0
  for listIndex in range(len(scoresList)):
    score= score + scoresList[listIndex]
  return score



#***************************************  



  




  ######################################################################################
  ##########################       SEARCH ELMS     #####################################
  ######################################################################################




def makeElmSearch(sequence,verbose):
  ##MAKE THE SEARCH AND SAVE RESULTS IN outputELM
  elm_pattern_dict = {}
  with open(toolsPath + "ELM/elm_patterns_20150301.txt", 'rU') as file_open :
	  patterns = file_open.readlines()
	  for line in patterns :
		  line = line.split("\t")
		  elm_id = line[0]
		  elm_pattern = line[1].replace("\n","")
		  elm_pattern_dict[elm_id] = elm_pattern

  #print "SEARCHING ELMS IN %s\n ..." % filename
  elm_pattern_desc_dict = {}
  if verbose:     #ONLY IF ITS REQUIRED TO PRINT A DETAILED A OUTPUT, READ THE DESCRIPTION OF EACH ELM FROM A DIFFERENT FILE
    with open(toolsPath + "ELM/elm_patterns_desc_20150301.txt", 'rU') as file_desc_open :
	    patterns_desc = file_desc_open.readlines()
	    for line in patterns_desc :
		    line = line.split("\t")
		    elm_id = line[0]
		    #print line[1]
		    elm_pattern_desc = line[1].replace("\n","")
		    elm_pattern_desc_dict[elm_id] = elm_pattern_desc
  #output_file_name = "elms_search_in_%s.txt" % filename[:-4]
  
  output_file_name=outputsPath + "outputELM"
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
				  if verbose:     #write description next to the indexes
				    #print elm_id
				    file_write.write(str(index+1) + tab + str(index+len(match.group())) + tab+ elm_id +tab+ elm_pattern_desc_dict[elm_id] + '\n')
				    #file_write.write("%s\t%s\t%s\t%s\n" % (index+1, index+len(match.group()) , elm_id , elm_pattern_desc_dict[elm_id] ))
				  else:
				    file_write.write("%s\t%s\n" % (index+1, index+len(match.group())))



##*************************************


def elmSearch(sequence, positionScores,verbose):
  
  makeElmSearch(sequence,verbose)   #SEARCH FOR MOTIFS IN MY SEQUENCE AND SAVE RESULTS IN A FILE	
  elmScores=[]
  for p in range(len(sequence)):
	      elmScores.append(0)
  #print "largoo:" + str(len(sequence))
  hitsFile=outputsPath + "outputELM"    #THE OUTPUT OF THE SEARCH CONTAINS THE LIST OF ELMs FOUND, NOW I HAVE TO PROCESS IT
  with open(hitsFile, "r") as input_file:
    lines=input_file.readlines()
  if verbose:  
    print indent + "ELM Search:"
  for line in lines:
    line=line.split("\t")
    pattern_start=int(line[0])
    pattern_end=int(line[1])
    if verbose:
      print indent + "Pattern found: " + line[2]
    #print "start:" + str(pattern_start)
    #print "end:" + str(pattern_end)
    for x in range(pattern_start-1,pattern_end):
      elmScores[x]+=1
      #print str(x)
	
  if verbose:
    print ""
    print indent + "ELM Search RESULTS:"
    data = [sequence,elmScores]
    col_width = max(len(str(word)) for row in data for word in row)  # padding  ***ACA IBA +1 AL FINAL PERO LO SAQUE
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
    #print indent + '\t'.join(map(str, sequence))
    #for u in range(len(sequence)):
      #if sequence[u]>=10:
	#sys.stdout.write(sequence[u])
	#sys.stdout.write(" ")
      #else:
	#sys.stdout.write(sequence[u])
	#sys.stdout.write("  ")
    #print endl	
    #print indent + '\t'.join(map(str, elmScores))	  
  
  ##ADD hits to global score
  for i in range(len(sequence)):
    positionScores[i]+=elmScores[i]
  
  
                                         


##*************************************



  ######################################################################################
  ##########################       PROSITE SEARCH      #####################################
  ######################################################################################



def prositeSearch(sequence, positionScores,verbose):
  
  
  #NEW LIST TO SAVE HITS 
  prositeScores=[]
  for p in range(len(sequence)):
	      prositeScores.append(0)
	      
  #SEARCH PROSITE USING PS_SCAN
  inputProsite=inputsPath + "sequenceFASTA"
  #input=open(inputsPath + "sequence" , "w")
  #input.write(">gi" + endl)
  #input.write(sequence)
  #input.close()
  proc = subprocess.Popen(['perl', toolsPath + 'ps_scan/ps_scan.pl','-r','-o', 'scan', inputProsite],stdout=subprocess.PIPE)
  while True:
    line = proc.stdout.readline()
    if line != '':
      pattern_start=int(line.split()[0])  
      pattern_end=int(line.split()[2])
      for x in range(pattern_start-1,pattern_end):
	prositeScores[x]+=1
      if verbose:
	  print indent + "Hit: " +line
      #print "Hit:" + line.split()[0] + "-" +line.split()[1] + space +  line.split()[2] + space + line.split()[3] + space + line.split()[4]
    else:
      break	      


  if verbose:
    #print endl
    print indent + "Prosite Search RESULTS:"
    print indent + sequence
    print indent + ''.join(map(str, prositeScores))	  
  
  ##ADD hits to global score
  for i in range(len(sequence)):
    positionScores[i]+=prositeScores[i]
  
  

                                         






  ######################################################################################
  ##########################      BLAST SEARCH    #####################################
  ######################################################################################



def blastIt(sequence, positionScores, database, verbose):
        global match
        ##BLAST SEARCH
        inputBlast=inputsPath+"inputBlast"
	outputBlast=outputsPath+"outputBlast"
        if blastWeb:       # WEB BLAST SEARCH
	  if verbose:
	    print indent + "WEB BLAST SEARCH IN PROGRESS..." 
	  result = NCBIWWW.qblast("blastp", database , sequence)
	  records = NCBIXML.parse(result)
	  first = records.next()
	else:     # LOCAL BLAST SEARCH
	  if verbose:
	    print indent + "LOCAL BLAST SEARCH IN PROGRESS..."
	  input=open(inputBlast, 'w')
	  input.write(sequence)
	  input.close()
	  commandLine=NcbiblastpCommandline(query=inputBlast, db=database, evalue=0.001, outfmt=5, out=outputBlast)
	  #print commandLine
	  stdout, stderr = commandLine()
	  result_handle = open(outputBlast)
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
		print indent + "Sequence name:", firstAlign.title
	      
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
	#TAKE THE positionScores LIST AND PUT "1" WHERE THERE WAS A MATCH AND "0" WHERE THERE WAS A GAP OR MISMATCH
		for j in range(len(sequence)):
			if j< (start-1) or j > (end-1):    
				#print sequence[j]
				positionScores[j]=1
			else:
				if hsp.match[j-start+1] <> "+" and hsp.match[j-start+1] <> " ":
					positionScores[j] = 1
				else:
					positionScores[j]=0
					
	if verbose:	    
	    #print endl
	    print indent + "BLAST RESULTS:"
	    print indent + sequence
	    print indent + ''.join(map(str, positionScores))				
					
					
					
					








  ######################################################################################
  ##########################      IUPRED SEARCH     #####################################
  ######################################################################################


def iupred(sequence, positionScores, verbose):
	runCommand=toolsPath + "iupred/iupredExe"+ space + inputsPath + "sequenceFASTA" +space+ "long" + space + outputsPath + "outIUPred"
	#input=open(inputsPath+"iupred/inputIupred"+exeid, 'w')
	#input.write("Name" + endl)
	#input.write(sequence)
	#input.close()
	os.system(runCommand)	
	outputIUPred=open(outputsPath + "outIUPred", "r")
	
	#PRINT THE RESULTS OF IUPred
	if verbose:
	  iupredScores=[]
	  iterOutputIUPred=iter(outputIUPred)
	  for p in range(len(sequence)):
	      iupredScores.append(0)
	  for x in range(len(sequence)):
		  resultX=float(iterOutputIUPred.next())
		  if resultX < 0.5 :
			  iupredScores[x] = 1
	  #print endl
	  print indent + "IUPred RESULTS:"
	  print indent + sequence
	  print indent + ''.join(map(str, iupredScores))	  			  			
			
	outputIUPred.seek(0)
	rstFile_iter = iter(outputIUPred)
	#ADD 1 TO THE POSITION IN positionScores IF THE RESULT IS LESS THAN 0.5 (PREDICTING A GLOBULAR TENDENCY)
	for j in range(len(sequence)):
		resultJ=float(rstFile_iter.next())
		if resultJ > 0.5 :
			positionScores[j] += 0
		else:
			positionScores[j] += 1				
  








  ######################################################################################
  ##########################    ANCHOR EVALUATION     #####################################
  ######################################################################################



def anchor(sequence, positionScores, verbose):
  inputAnchor=inputsPath + "sequenceFASTA"
  runCommand=toolsPath + "ANCHOR/anchor" + space + inputAnchor + space + outputsPath + "outAnchor"
  #input=open('ANCHOR/input', 'w')
  #input.write("Name" + endl)
  #input.write(sequence)
  #input.close()
  os.system(runCommand)	
  outputAnchor=open(outputsPath + "outAnchor", "r")
  #print "aca evaluo anchor"
  #PRINT THE RESULTS OF ANCHOR
  anchorScores=[]
  iterOutputAnchor=iter(outputAnchor)
  for p in range(len(sequence)):
      anchorScores.append(0)
  for x in range(len(sequence)):
	  resultX=float(iterOutputAnchor.next())
	  if resultX >  0.5 :
	    anchorScores[x] = 1
  if verbose:
    print indent + "ANCHOR RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, anchorScores))	  			  			
    data = [sequence,anchorScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
  outputAnchor.seek(0)
  rstFile_iter = iter(outputAnchor)
  #ADD 1 TO POSITIONS IN positionScores  IF THE RESULT IS GREATER THAN 0.5 (PREDICTING A DISORDERED BINDING REGION)
  for j in range(len(sequence)):
	  resultJ=float(rstFile_iter.next())
	  if resultJ > 0.5 :
		  positionScores[j] += 1
	  else:
		  positionScores[j] += 0				







  ######################################################################################
  ##########################       AMYLOID SEQUENCE DETERMINANTS  ######################
  ######################################################################################
  


def amyloidPatternSearch(sequence, positionScores,verbose):
    amyloidScore=[]
    for p in range(len(sequence)):
	amyloidScore.append(0)
    
    #with open(output_file_name, 'w') as file_write :
    #for elm_id in elm_pattern_dict :
    where_to_start = []
    elm_pos_dict = {}
    pattern = re.compile("[^P][PKRHW][VLSCWFNQE][ILTYWFNE][FIY][^PKRH]")
    for matched_string in pattern.finditer('%s' % sequence) :
	    where_to_start.append(matched_string.start())
    pattern = re.compile("[^P][PKRHW][VLSCWFNQE][ILTYWFNE][FIY][^PKRH]")
    for index in where_to_start :
	    match = re.search(pattern, '%s' % sequence[index:])
	    if match != None :
		    print "NEUTRAL pH SEQUENCE DETERMINANT FOUND "
		    for x in range(index,index+len(match.group())):
		      amyloidScore[x]+=1
		    #if verbose:     #write description next to the indexes
		      #print elm_id
		      #file_write.write(str(index+1) + tab + str(index+len(match.group())) + '\n')
		      
		      #file_write.write("%s\t%s\t%s\t%s\n" % (index+1, index+len(match.group()) , elm_id , elm_pattern_desc_dict[elm_id] ))
		    #else:
		      #file_write.write("%s\t%s\n" % (index+1, index+len(match.group())))
	
	
    if verbose:
      print indent + "AMYLOID SEQUENCE DETERMINANTS RESULTS:"
      data = [sequence,amyloidScore]
      col_width = max(len(str(word)) for row in data for word in row)   # padding
      for row in data:
	print indent + "|".join(str(word).ljust(col_width) for word in row)
      #print indent + sequence
      #print indent + ''.join(map(str, amyloidScore))







  ######################################################################################
  ##########################       TANGO EVALUATION     #####################################
  ######################################################################################


def tangoSearch(sequence, positionScores,verbose):
  outputTango= outputsPath+"outputTango"
  
  runCommand=toolsPath + 'tango/tango_x86_64_release tangoResults nt="N" ct="N" ph="7" te="298" io="0.05" seq="' + sequence + '" > ' + outputTango
  print runCommand 
  os.system(runCommand)
  outputTango=open(outputTango,'r')
  
  tangoScores=[]
  for p in range(len(sequence)):
    tangoScores.append(0)
  
  position=0
  for line in outputTango.readlines()[1:len(sequence)+1]:
    #print line.split()[1] 
    beta=float(line.split()[2])
    turn=float(line.split()[3])
    helix=float(line.split()[4])
    aggregation=float(line.split()[5])
    tangoCutff=1
    if beta > tangoCutff or turn > tangoCutff or helix > tangoCutff or aggregation > tangoCutff:
      tangoScores[position]=1
    position+=1
  #print str(beta) + tab + str(turn) + tab + str(helix) + tab + str(aggregation) 
  if verbose:
    print indent + "TANGO RESULTS:"
    print indent + sequence
    print indent + ''.join(map(str, tangoScores))	  
    
    
  for x in range(0,len(sequence)):
    if tangoScores[x] > 0:
      positionScores[x] += 1





  ######################################################################################
  ##########################       LIMBO EVALUATION     #####################################
  ######################################################################################
 
def limboEval(sequence, positionScores,verbose):
  #input=open(inputsPath + "sequenceLimbo" + exeId, "w")
  #input.write(">gi" + endl)
  #input.write(sequence)
  #input.close()
  outputLimbo= outputsPath + "outLimbo"
  
  #CALL LIMBO :   score.py + matrix + input + outpout
  runCommand="python" + space + toolsPath + "Limbo/score.py" + space + toolsPath+"Limbo/mergedmatrix.mat" + space + inputsPath + "sequenceFASTA" + space + outputLimbo
  os.system(runCommand)
  outputLimbo=open(outputLimbo,'r')
  limboScores=[]
  for p in range(len(sequence)):
    limboScores.append(0)
  for line in outputLimbo.readlines():
    #print line.split()[1] 
    hitStart=int(line.split()[0])  #first column is the start of the heptapeptide hit
    for y in range(hitStart-1,hitStart+7):
      limboScores[y] += 1

  for x in range(0,len(sequence)):
    positionScores[x] += limboScores[x]
  if verbose:
    print indent + "LIMBO RESULTS:"
    print indent + sequence
    print indent + ''.join(map(str, limboScores))	  






  ######################################################################################
  ##########################       TMHMM EVALUATION     #####################################
  ######################################################################################
 
def tmhmmEval(sequence, positionScores,verbose):
  #input=open("sequenceTmhmm", "w")
  #input.write(">gi" + endl)
  #input.write(sequence)
  #input.close()
  outputTmhmm=outputsPath + "outTmhmm"
  runCommand= toolsPath + "tmhmm/bin/tmhmm" + space + inputsPath + "sequenceFASTA" + space + ">" + outputTmhmm
  os.system(runCommand)
  outputTmhmm=open(outputTmhmm,'r')
  tmhmmScores=[]
  for p in range(len(sequence)):
    tmhmmScores.append(0)
  for line in outputTmhmm.readlines():
    if line.split()[0] == "TMhelix":
      #print "hit de estee"
      hitStart=int(line.split()[1])-1 
      hitEnd=int(line.split()[2]) 
      for y in range(hitStart,hitEnd):
	tmhmmScores[y] += 1
  if verbose:	
    print indent + "TMHMM RESULTS:"
    print indent + sequence
    print indent + ''.join(map(str, tmhmmScores))		
  for x in range(0,len(sequence)):
    positionScores[x] += tmhmmScores[x]










  ######################################################################################
  ##########################       GENERAL SEQUENCE EVALUATION     #####################################
  ######################################################################################



def sequenceEvaluation(sequence, positionScores, verbose):
	
	#SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
	input=open(inputsPath + "sequenceFASTA"  , "w")
	input.write(">gi" + endl)
	input.write(sequence)
	input.close()
	
	##: BLAST SEARCH
	if verbose:
	   print endl
	   print indent + "*************************************"
	   #print indent + "STARTING BLAST SEARCH"
	blastIt(sequence,positionScores,database, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	    
        ## IUPred evaluation
	if verbose:
	  print endl
	  print indent + "*************************************"
	  #print indent + "STARTING IUPred"
	iupred(sequence, positionScores, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	  
        ## ANCHOR evaluation
	if verbose:
	  print endl
	  print indent + "*************************************"
	  #print indent + "STARTING ANCHOR"
	anchor(sequence, positionScores, verbose)
	if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	  
	if verbose:
	  print endl
	  print indent + "*************************************"
	  #print indent + "STARTING ELM Search"
        elmSearch(sequence,positionScores, verbose)
	if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	  
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "Prosite Search in progress..."
        prositeSearch(sequence,positionScores, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "STARTING LIMBO Search"
        if (len(sequence) > 11):
	  limboEval(sequence,positionScores, verbose)
        else:
	  if verbose:
	    print indent + "LIMBO Search is only available for sequences of length >= 12 "
	print indent + "*************************************"
	  
	if verbose:
	  print endl
	  print indent + "*************************************"
	  print indent + "Search for transmembrane sections"
        tmhmmEval(sequence,positionScores, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	  
	if verbose:
	  print endl
	  print indent + "EVALUATE AMYLOID FIBRIL FORMATION "
	  print indent + "*************************************"
	  print indent + "Search for sequence determinants"
        amyloidPatternSearch(sequence,positionScores, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	  
        if verbose:
	  print endl
	  print indent + "*************************************"
	  #print indent + "STARTING TANGO Search"
        tangoSearch(sequence,positionScores, verbose)
        if stepByStep:
	  raw_input(indent + "Hit enter to continue with next evaluation")
	
	
	if stepByStep:
	  raw_input(indent + "Press enter to see final results...")
	
	
	
        ##PRINT SCORE
        if verbose:
	  print indent + "*************************************"
	  print endl
	  print indent + "FINAL RESULTS:"
	  print indent + sequence
	  print indent + ''.join(map(str, positionScores))
	  print indent + "GLOBAL SCORE:" + str(getGlobalScore(positionScores))
	  print indent + "*************************************"
	if stepByStep:
	  raw_input(indent + "....hit enter to continue")














	



#***********************************************************************



## PRINT PROGRAM HELP 
def printHelp():
  print endl + "Usage: "
  print "  python bleach.py [options]" + endl
  print "Options are:"
  print tab + "--length " + tab + "Sequence lenght"
  print tab + "--db " + tab + "swissprot | nr"        #blast database
  print tab + "--composition " + tab + "[average | user_specified]"
  print tab + "--seq " + tab + "initial-sequence"
  print tab + "--maxmutations " + tab + "Max amount of mutations"     # 
  #print tab + "--maxiterations"   ?????    SUM OF MUTATION ATTEMPTS ?????     ******MAINLY FOR PERFORMANCE REASONS
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 
#************************************************************************* 
 
#**************************
#******* MAIN *************
#**************************





#********DEFAULTS********
composition="average"
a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=0
database="uniprot_sprot.fasta"
sequence="RANDOM"
uvsilent=False




#*******PARAMETERS******
if len(sys.argv) < 2:
  print "************************************************"
  print "************************************************"
  print "USING DEFAULT VALUES - If you want a specific configuration, see the options using  --help"
  print "************************************************"
  print "************************************************"


else:
  for index in range(1,len(sys.argv)):
    arg = sys.argv[index]
    
    if (arg=='-H' or arg== '-help' or arg== '--help'):
      printHelp()
      exit()
    elif (arg=='--length') and (index < len(sys.argv)):
      length = int(sys.argv[index+1])
    elif (arg=='--beta') and (index < len(sys.argv)):
      beta = float(sys.argv[index+1])  
    elif (arg=='--seq') and (index < len(sys.argv)):
      sequence = sys.argv[index+1]
      length=len(sequence)
      rand=False
    elif (arg=='--db') and (index < len(sys.argv)):
      database = sys.argv[index+1]
    elif (arg=='--blastweb'):
      blastWeb=True
    elif (arg=='--uvsilent'):
      uvsilent=True  
    elif (arg=='--verbose'):
      verbose=True
    elif (arg=='--stepped'):
      stepByStep=True
      verbose=True   #MAKES NO SENSE TO GO STEP BY STEP IF CANT SEE A DETAILED OUTPUT
    elif (arg=='--maxiterations') and (index < len(sys.argv)):
      maxIterations=int(sys.argv[index+1])   
    elif (arg== '--composition') and (index < len(sys.argv)):
      composition = sys.argv[index+1]
      if (composition=="user_specified"):  #frequencies specified by parameter
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
	
  print "************************************************"
  print "************************************************"
  print "VALUES:"
  print "length=" + str(length) 
  print "composition=" + composition
  print "sequence=" + sequence +endl
  print "************************************************"
  print "************************************************"



#MAKE INPUT FOLDER
os.mkdir(inputsPath)
os.mkdir(outputsPath)

#OUTPUT
outputPath = "Output/"   
scoresFile="scores" + str(length)   #save Scores Vs iteration number 
mutAttemptsFile="mutationsAttempt" + str(length)  # save number of mutation attempts  Vs iteration number
timesFile='times' + str(length)   #save times Vs iteration number
totalTimesFile='totalTimes'  #save total time elapsed Vs sequence length


if output:
  totalTimesOutputFile=open( outputPath + totalTimesFile , "a")
  timesOutputFile=open( outputPath + timesFile , "a")
  scoresOutputFile=open( outputPath + scoresFile , "a")
  mutationsFile=open( outputPath + mutAttemptsFile , "a")
  


if uvsilent==True:   
  aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]
else:
  aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",0),("P",470),("S",656),("T",0),("W",108),("Y",0),("V",687) ]





# FORMAT TO REQUEST RANDOM SEQUENCE:   
#http://web.expasy.org/cgibin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

if rand==True:
	#****************GET RANDOM SEQUENCE*************
        #print endl
	print "Generating random sequence..."    
	#print endl
	if not (composition=="user_specified"):
	  url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(length) + "&comp=" + composition + "&output=fasta"
	else:
	  print "fix this"
	response = urllib2.urlopen(url)
	html = response.read()
	i = html.index('\n')
	sequence = html[i+1:].replace('\n', '')
	print "*******************************" 
	

   
#print endl
#PRINT STARTING SEQUENCE (BEFORE ANY MUTATION)   
if rand==True:
  print "INITIAL SEQUENCE:" + sequence


#CREATE ARRAY TO SAVE MUTATION FREQUENCY
positionScores=[]
mutatedScores=[]
for p in range(len(sequence)):
  positionScores.append(0)
  mutatedScores.append(0)



time0 = time.time()   #start time
timePrev=time0          #Previous iteration start time
################################
#### ITERATE OVER SEQUENCE ####
#############################

iteration=1
timePrev=time.time()

indent=""

for p in range(len(sequence)):
  #positionScores.append(0)
  positionScores[p]=0



if verbose:
  print "*****************************"
  print " INITIAL EVALUATION "
  print "*****************************"    
#EVALUATE INITIAL SEQUENCE
sequenceEvaluation(sequence,positionScores, verbose)	
globalScore= getGlobalScore(positionScores)

if verbose:
  print "*****************************"
  print "*****************************"
  print endl
  print endl

while globalScore > 0 and iteration <= maxIterations:
  if verbose:
    print "*****************************"
    print "STARTING ITERATION " + str(iteration)
    print "*****************************"
  
    print "Current sequence: " + sequence
    print "Current score:    " + ''.join(map(str, positionScores))
    print "Current global score:    " + str(globalScore)
  
  raw_input("")
  #BEFORE EACH ITERATION STEP, CLEAN MUTATION FREQUENCE
  #THIS LIST SUMS UP THE PROBABILITY TO MUTATE EACH POSITION BASED ON ALL STEPS PERFORMED (SUMMING IT UP GIVES THE GLOBAL SCORE OF THE SEQUENCE)
  #positionScores=[]
 
 
  
  ###ONLY MUTATE IF SCORE GREATER THAN 0  (AT LEAST 1 POSITION HAS A MUTATION FREQUENCY GREATER THAN 0 )
  #print endl
  #print "*************************************"
  #print "*************************************"
  #print "MUTATION STEP"   #START TRYING TO MUTATE
  #print "*************************************"
  #print "*************************************"
  #print "Sequence before: " + sequence 
 
    
  ## AA MUTATION OF SEQUENCE - PROBABILITY OF MUTATION BASED ON WEIGHTS IN PONDERADA 
  ## PROBABILITY OF MUTATION BASED ON WEIGHTS ARRAY
  ## RETURNS MUTATED SEQUENCE



#def mutar(sequence, positionScores):
    #global indent, mutAttempts
  
  #weighted IS A PAIRLIST(position,weight)
  weighted=[]
  
  previousScore=getGlobalScore(positionScores)    #CHANGE THIS!!!!
  #if previousScore > 0:
  mutAttempts=0       #COUNT MUTATIONS ATTEMPTS
  
  #WEIGHT LIST: CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. THE WEIGHT IS = SCORE + A BASE WEIGHT
  for x in range(len(positionScores)):
    weighted.append((x, positionScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
  while 10000 > mutAttempts:    ##JUST A SYMBOLIC MAX. AMOUNT OF MUTATIONS ATTEMPTS
      mutAttempts+=1
      
      indent = tab + tab    #output formatting 
      #SELECT A POSITION 
      if verbose:
	print endl
	print indent + "*************************************"
	print indent + "MUTATION ATTEMPT"
	print indent + "*************************************"
      #print indent + "Score before:    " + str(previousScore)
      #print indent + "Choose a position based on sequence weights"
      
      #CHOOSE A POSITION BASED ON WEIGHTS
      mutatePosition= weighted_choice(weighted) 
      if verbose:
	print indent + "Position chosen: " + str(mutatePosition)
      
      #SELECT THE NEW AA FOR THAT POSITION (BASED ON LIST OF FREQUENCIES)
      previousResidue=sequence[mutatePosition]
      if verbose:
	print indent + "Residue to mutate: " + previousResidue
      seleccionado= previousResidue
      
      #SELECT A NEWONE UNTIL THE RESIDUE IS DIFFERENT FROM PREVIOUS
      while previousResidue == seleccionado:
	  seleccionado = weighted_choice(aaFrequencies)	
      if verbose:	  
	print indent + "New residue : " + seleccionado

      ##BUILD MUTATED SEQUENCE WITH NEW RESIDUE    
      mutatedSequence = sequence[0:mutatePosition]
      mutatedSequence += seleccionado
      mutatedSequence += sequence[mutatePosition+1:]
      if verbose:
	print indent + "Original sequence: " + sequence
	print indent + "Mutated sequence : " + mutatedSequence

      ##RESET LIST OF SCORES FOR THE MUTATED SEQUENCE
      for p in range(len(sequence)):
	  mutatedScores[p]=0
      
      if stepByStep:
	raw_input(indent + "...Hit enter to start evaluation")
      if verbose:
	print ""
	indent=tab + tab + tab #output formatting stuff
	print indent + "STARTING PROPOSED MUTATION EVALUATION"
      sequenceEvaluation(mutatedSequence, mutatedScores, verbose)	 
      mutatedScore=getGlobalScore(mutatedScores)
      #if stepByStep:
	#raw_input(indent + "Hit enter to continue with mutation acceptance")
      #IF THE GLOBAL SCORE DECREASED
      if verbose:
	indent=tab + tab + tab + tab + tab
	print ""
	print indent + "*************************************"
	print indent + "DECISION"
	print indent + "Previous sequence"
	print indent + sequence
	print indent + ''.join(map(str, positionScores))
	print indent + "Global score: " + str(globalScore)
	print ""
	print indent + "Mutated sequence"
	print indent + mutatedSequence
	print indent + ''.join(map(str, mutatedScores))
	print indent + "Global score: " + str(mutatedScore)
	print ""
      if previousScore >= getGlobalScore(positionScores):
	  print indent + "Previous score (" + str(previousScore) + ") >= Mutated score (" + str(mutatedScore) + ")" 
	  print indent + "...ACCEPT MUTATION"
	  if stepByStep:
	    raw_input("")
	  break
	    #raw_input(indent + "Hit enter to continue with next iteration")
	  #return mutatedSequence
      else:
	  #DECISION BASED ON MONTE CARLO
	  print indent + "Previous score (" + str(previousScore) + ") < Mutated score (" + str(mutatedScore) + ")" 
	  diff=previousScore-getGlobalScore(mutatedScores)
	  print indent + "SCORE DIFFERENCE:" + str(diff)
	  exponent=diff/beta
	  MCvalue=math.exp(exponent)
	  #print "  :" + str(exponent)
	  print indent + "MC VALUE:" + str(MCvalue)     #e^(dif/beta)
	  
	  #GENERATE RANDOM NUMBER BETWEEN 0 AND 1
	  randy=random.random()
	  print indent + "RANDOM VALUE [0,1]:" + str(randy)
	  #print indent + "MONTE CARLO DECISION:"
	  if MCvalue > randy:
	    #ACCEPT MUTATION
	    print indent + "...ACCEPT MUTATION"
	    if stepByStep:
	      raw_input("")
	      #raw_input(indent + "Hit enter to continue with next iteration")
	    #return mutatedSequence
	    break
	  else:	    
	    #print "El score original " + str(getGlobalScore(mutatedScores)) + " es mayor que "+ str(getGlobalScore(positionScores))
	    print indent + " Mutation score (" + str(mutatedScore) + ") >= Previous score (" + str(previousScore) + ")" 
	    print indent + "...DENY MUTATION"
	    if stepByStep:
	      raw_input(indent + "Hit enter to continue with next attempt")
	    #return sequence
	    #break
      print "*************************************"
    #return sequence
  #else:
  #return sequence




###********************************************
 
  #mutatedSequence=mutar(sequence, positionScores)   *********esta funcion ya no se llama mas********
  
#/**********CUANDO SALE DE TODOS LOS INTENTOS DE MUTACION, mutatedSequence TIENE LA NUEVA SECUENCIA   
  
  
  #print endl
  print endl
  #print "Sequence after mutation:    " + mutatedSequence
  sequence=mutatedSequence
  positionScores=mutatedScores
  print "Attempts before mutation accept:" + str(mutAttempts)
  #print endl
  print "*******************************************"
  #sys.stdin.read(1) 
  
  #else:     ##GLOBAL SCORE IS 0 - REACHED END OF LOOP
  #    break
  globalScore= getGlobalScore(positionScores)
  print "End of iteration " + str(iteration)
  print "Global score :    " + str(globalScore)
  print "*******************************************"
  print endl
  if stepByStep:
    raw_input("Hit enter to continue with next iteration")
  if output:
    # MUTATION ATTEMPTS Vs. ITERATION NUMBER
    mutationsFile.write(str(iteration)+ tab + str(mutAttempts) + tab +str(beta) + endl)
    
    #SCORES Vs. ITERATION NUMBER
    scoresOutputFile.write(str(iteration)+ tab + str(globalScore) + tab + str(beta) + endl)
    
    #ITERATION ELAPSED TIME
    timeX=time.time()-timePrev   #ITERATION TIME
    timesOutputFile.write(str(iteration)+ tab + str(timeX) + tab + str(beta) + endl)
    
	    
  iteration=iteration+1   #NUMBER OF MUTATIONS ACCEPTED
  #indent=""


print "**END OF SEARCH**"
if globalScore==0:
  print "REACHED SCORE = 0"
  print sequence
else:
  print "REACHED LIMIT OF ITERATIONS"
  print "Final sequence: " + sequence
  print "Final score:    " + ''.join(map(str, positionScores))
  print "Global score: " + str(globalScore)
  



#PERFORMANCE TESTS OUTPUT
if output:        
  totalElapsedTime=time.time() - time0   ##
  print "Elapsed time:",totalElapsedTime, "Seconds"
  totalTimesOutputFile.write(str(length) + tab + str(iteration-1)+ tab + str(totalElapsedTime) + tab + str(beta) +  endl)

  ##CLOSE ALL OUTPUT FILES
  totalTimesOutputFile.close()
  timesOutputFile.close()
  scoresOutputFile.close()
  mutationsFile.close()
   

