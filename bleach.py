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



def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

#****************************************
# CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))






#*********************
#******GLOBALS*******
#******************

#  OUTPUT FORMATTING
endl = "\n"
tab = "\t"
space=" "
indent=""

#  EXECUTION Id . Identify runs
exeId=os.getpid() 

#  CUTOFFS, THRESHOLDS, LIMITS ....
maxIterations=10000
cutoff=0.01  #BLAST cutoff
waltzThreshold=79.0
beta=0.5   #MC 
targetScore=0.0
pastaThreshold=-5.5  #ENERGY Threshold
pastaProbabilityThreshold=0.05 #Aggregation probability threshold
iupredThreshold=0.5
anchorThreshold=0.5
tangoThreshold=1.0

#  EXECUTION PARAMETERS
minimalOutput=False  #True = only print global scores at end of iteration
verbose=False        #True = print detailed information of execution
stepByStep=False
match=False
rand=True
change=True
evaluateNetCharge=False
#targetNetCharge=0
blastWeb=False  #BLAST SEARCH LOCAL OR WEB
output=True  ##print info to file
mutAttempts=0
length=12   #defaul sequence length


#EVALUATION PARAMETERs (TEST MODE)



#All tools enabled by default 
runBlast=runTango=runPasta=runWaltz=runElm=runProsite=runLimbo=runTmhmm=runIupred=runAnchor=runAmyloidPattern=True

#FILES
logFileName='mutations' + str(exeId) + '.log'

#PATHS
basePath=getScriptPath() + '/'
toolsPath=basePath + 'Tools/'    #**************************TODO SET THE PATH TO THE TOOL SET 
inputsPath=basePath + "/Input/"+ str(exeId) + "/" #SET PATH TO SAVE INPUTS FILES
baseOutputPath=basePath + "/Output/" 
outputsPath=baseOutputPath + str(exeId) + "/"
testOutputPath=outputsPath   # DEFAULT OUTPUT FOR TESTs 


#####CREATE INPUT AND OUTPUT DIRS
try:
    os.makedirs("Input")
except OSError as exc: 
    if exc.errno == errno.EEXIST and os.path.isdir("Input"):
        pass

try:
    os.makedirs("Output")
except OSError as exc: 
    if exc.errno == errno.EEXIST and os.path.isdir("Output"):
        pass


#AA FREQUENCIES TO SELECT NEW RESIDUES FOR MUTATIONS (from http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html) 
#aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]



  
#***************************************************************************  
  

##    JUST SUM UP THE INDIVIDUAL SCORES  
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
      print indent + "start:" + str(pattern_start)
      print indent + "end:" + str(pattern_end)
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
  proc = subprocess.Popen(['perl', toolsPath + 'Prosite/ps_scan/ps_scan.pl','-r','-o', 'scan', inputProsite],stdout=subprocess.PIPE)
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
    print indent + "RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, prositeScores))	  
    data = [sequence,prositeScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
	print indent + "|".join(str(word).ljust(col_width) for word in row)
  ##ADD hits to global score
  for i in range(len(sequence)):
    positionScores[i]+=prositeScores[i]
  
  

                                         


##*************************************



  ######################################################################################
  ##########################       CHARGE SEARCH      #####################################
  ######################################################################################



def chargedSearch(sequence, positionScores,verbose):
   
  chargeScores=[]
  
  #FIRST EVALUATE THE NET CHARGE OF SEQUENCE
  netCharge=0
  for p in range(len(sequence)):
      	if (sequence[p]=='K') or (sequence[p]=='R') or (sequence[p]=='H'):
	  #IF POSITIVELY CHARGED
	  netCharge+=1
	else:
	  if (sequence[p]=='E') or (sequence[p]=='D'): 
	    #IF NEGATIVELY CHARGED 
	    netCharge-=1


  #SET SCORE VALUES BASED ON NET-CHARGE IMPULSE
  if netCharge > targetNetCharge:
    #AIMING FOR A MORE NEGATIVE NET CHARGE
    for p in range(len(sequence)):
	if (sequence[p]=='K') or (sequence[p]=='R') or (sequence[p]=='H'):
	  #POSITIVELY CHARGED AA
	  chargeScores.append(2)
	else:
	  if (sequence[p]=='E') or (sequence[p]=='D'): 
	    #IF NEGATIVELY CHARGED 
	    chargeScores.append(0)
	  else:
	    #NEUTRAL AA
	    chargeScores.append(1)
	    
  else:
    if netCharge < targetNetCharge:
      #AIMING FOR A MORE POSITIVE NET CHARGE
      for p in range(len(sequence)):
		if (sequence[p]=='K') or (sequence[p]=='R') or (sequence[p]=='H'):
		  chargeScores.append(0)
		else:
		  if (sequence[p]=='E') or (sequence[p]=='D'): 
		    #IF NEGATIVELY CHARGED 
		    chargeScores.append(2)
		  else:
		    #NEUTRAL AA
		    chargeScores.append(1)
    else:
      #THE NET CHARGE IS CORRECT
      for p in range(len(sequence)):
	chargeScores.append(0)
      
  
  if verbose:
    #print endl
    print indent + "Charge search RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, chargeScores))	  
    data = [sequence,chargeScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
	print indent + "|".join(str(word).ljust(col_width) for word in row)
  
  ##ADD hits to global score
  for i in range(len(sequence)):
    positionScores[i]+=chargeScores[i]
  
  

                                 




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

	blastScores=[]    
	for p in range(len(sequence)):
	      blastScores.append(0)
	if match:
		for j in range(len(sequence)):
			if j< (start-1) or j > (end-1):    
				blastScores[j]+=1
				#print sequence[j]
				positionScores[j]+=1
			else:
				if hsp.match[j-start+1] <> "+" and hsp.match[j-start+1] <> " ":
					positionScores[j] += 1
					blastScores[j]+=1
				#else:
					#positionScores[j]+=0
					
	if verbose:	    
	    #print endl
	    print indent + "BLAST RESULTS:"
	    #print indent + sequence
	    #print indent + ''.join(map(str, positionScores))				
	    data = [sequence,blastScores]
	    col_width = max(len(str(word)) for row in data for word in row)   # padding
	    for row in data:
	      print indent + "|".join(str(word).ljust(col_width) for word in row)
					
					
					








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
		  if resultX < iupredThreshold :
			  iupredScores[x] = 1
	  #print endl
	  print indent + "RESULTS:"
	  print indent + 'Threshold: '+ str(iupredThreshold)
	  #print indent + sequence
	  #print indent + ''.join(map(str, iupredScores))	  			  			
	  data = [sequence,iupredScores]
	  col_width = max(len(str(word)) for row in data for word in row)   # padding
	  for row in data:
	    print indent + "|".join(str(word).ljust(col_width) for word in row)
	outputIUPred.seek(0)
	rstFile_iter = iter(outputIUPred)
	#ADD 1 TO THE POSITION IN positionScores IF THE RESULT IS LESS THAN 0.5 (PREDICTING A GLOBULAR TENDENCY)
	for j in range(len(sequence)):
		resultJ=float(rstFile_iter.next())
		if resultJ < iupredThreshold :
			positionScores[j] += 1
				
  








  ######################################################################################
  ##########################    ANCHOR EVALUATION     #####################################
  ######################################################################################



def anchor(sequence, positionScores, verbose):
  inputAnchor=inputsPath + "sequenceFASTA"
  runCommand=toolsPath + "ANCHOR/anchor" + space + inputAnchor + space + outputsPath + "outAnchor"
 
  os.system(runCommand)	
  outputAnchor=open(outputsPath + "outAnchor", "r")
  
  if verbose:
    print indent + 'Threshold: ' + str(anchorThreshold)
  anchorScores=[]
  iterOutputAnchor=iter(outputAnchor)
  for p in range(len(sequence)):
      anchorScores.append(0)
  for x in range(len(sequence)):
	  resultX=float(iterOutputAnchor.next())
	  if resultX >  anchorThreshold :
	    anchorScores[x] = 1
  #PRINT THE RESULTS OF ANCHOR
  if verbose:
    print indent + "RESULTS:"
    data = [sequence,anchorScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
  outputAnchor.seek(0)
  rstFile_iter = iter(outputAnchor)
  #ADD 1 TO POSITIONS IN positionScores  IF THE RESULT IS GREATER THAN Threshold (PREDICTING A DISORDERED BINDING REGION)
  for j in range(len(sequence)):
	  resultJ=float(rstFile_iter.next())
	  if resultJ > anchorThreshold :
		  positionScores[j] += 1
	  #else:
		  #positionScores[j] += 0				







  #####################################################################################
  #########################  WALTZ EVALUATION #########################################
  #####################################################################################


def waltzSearch(sequence, positionScores,verbose):
  
  if verbose:
     #print indent + "WALTZ SEARCH"
     print indent + "Threshold: " + str(waltzThreshold)
  #EXECUTE WALTZ EVALUATION
  inputWaltz=inputsPath + "sequenceFASTA"
  proc = subprocess.Popen(['perl', toolsPath + 'waltz/scoreMatrixGT.pl', inputWaltz, toolsPath + 'waltz/616.mat', 'full'],stdout=subprocess.PIPE)

  #PROCESS OUTPUT  
  waltzScores=[]
  for p in range(len(sequence)):
    waltzScores.append(0)
 
  for q in range(0,len(sequence)-5):
    line = proc.stdout.readline()
    if float(line.split()[2])> waltzThreshold:
	  if verbose:
	     print indent + 'Subsequence above threshold: ' +  line.split()[0]
	  hit_start=int(line.split()[1]) - 1
          hit_end=hit_start + 6
          for x in range(hit_start-1,hit_end):
            waltzScores[x]+=1
    #else:
      #break

  #position=0
  if verbose:
    print indent + "RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, waltzScores))
    data = [sequence,waltzScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)

  for x in range(0,len(sequence)):
    if waltzScores[x] > 0:
      positionScores[x] += 1








  #####################################################################################
  #########################  PASTA EVALUATION #########################################
  #####################################################################################


def pastaSearch(sequence, positionScores,verbose):
  input=open(inputsPath + "seq.fasta" , "w")
  input.write(">gi" + endl)
  input.write(sequence)
  input.close()
  pastaPath=toolsPath + 'PASTA/pasta_exe/'
  runCommand = "perl " + pastaPath+'PastaPairs.pl' + space + pastaPath +'pot_pasta.dat '+ inputsPath + " 1 0 self " + str(pastaThreshold) +  " > /dev/null"
  #print 'comando:' + runCommand
  os.system(runCommand)
  ### CHECK THIS!!! THE OUTPUT IS IN THE SAME DIR AS INPUT (CHECK PASTA PERL SCRIPT)
  outputPasta=open(inputsPath + "seq-seq.best_pairings_list.pair.dat")
  pastaScores=[]
  for p in range(len(sequence)):
    pastaScores.append(0)
  if verbose:
    print indent + 'Threshold: ' + str(pastaThreshold)
  for line in outputPasta.readlines():
     fromP, dash,toP = (line.split()[9]).partition('-')
     for hits in range(int(fromP)-1,int(toP)):
       pastaScores[hits]=+1 
     if verbose: 
	print indent + 'Hit:'
	print indent + 'Energy value: ' + line.split()[4]	
	print indent + 'Positions:    ' + line.split()[9]
  #position=0
  if verbose:
    print indent + "RESULTS:"
    data = [sequence,pastaScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
  
  for x in range(0,len(sequence)):
    if pastaScores[x] > 0:
      positionScores[x] += 1






  ######################################################################################
  ##########################       AMYLOID SEQUENCE DETERMINANTS  ######################
  ######################################################################################
  


def amyloidPatternSearch(sequence, positionScores,verbose):
    amyloidScore=[]
    for p in range(len(sequence)):
	amyloidScore.append(0)
   
   #TODO READ PATTERN FROM FILE IN Tools DIR
    pattern = re.compile("[^P][PKRHW][VLSCWFNQE][ILTYWFNE][FIY][^PKRH]")
    where_to_start = []
    #elm_pos_dict = {}
    for matched_string in pattern.finditer('%s' % sequence) :
	    where_to_start.append(matched_string.start())
    #pattern = re.compile("[^P][PKRHW][VLSCWFNQE][ILTYWFNE][FIY][^PKRH]")
    hits=False
    for index in where_to_start :
	    match = re.search(pattern, '%s' % sequence[index:])
	    if match != None :
	            hits=True
		    if verbose:
		    	print indent + "NEUTRAL pH SEQUENCE DETERMINANT FOUND "
		    for x in range(index,index+len(match.group())):
		      amyloidScore[x]+=1
	    #else:
	      #if verbose:
		#print indent + 'NO MATCHES'
	
    if verbose:
      if hits:
	print indent + "RESULTS:"
	data = [sequence,amyloidScore]
	col_width = max(len(str(word)) for row in data for word in row)   # padding
	for row in data:
	  print indent + "|".join(str(word).ljust(col_width) for word in row)
      else:
	print indent + 'NO HITS'










  ######################################################################################
  ##########################       TANGO EVALUATION     #####################################
  ######################################################################################


def tangoSearch(sequence, positionScores,verbose):
  outputTango= outputsPath+"tangoResults.txt"
  #  32bits bin
  #runCommand=toolsPath + 'tango/tango_i386_release tangoResults nt="N" ct="N" ph="7" te="298" io="0.05" seq="' + sequence + '" > /dev/null' 
  
  #COULD NOT CHANGE THE OUTPUT PATH OF TANGO SO I HAVE TO CHANGE DIR MOMENTLY TO GET THE OUTPUT 
  os.chdir(outputsPath)
  runCommand=toolsPath + 'tango/tango_x86_64_release tangoResults nt="N" ct="N" ph="7" te="298" io="0.05" seq="' + sequence + '" > /dev/null'
  #print runCommand 
  os.system(runCommand)
  outputTango=open(outputTango,'r')
  os.chdir(basePath)
  tangoScores=[]
  for p in range(len(sequence)):
    tangoScores.append(0)
  if verbose:
    #print indent + 'TANGO SEARCH'
    print indent + 'Threshold: ' + str(tangoThreshold)
  position=0
  for line in outputTango.readlines()[1:len(sequence)+1]:
    #print line.split()[1] 
    beta=float(line.split()[2])
    turn=float(line.split()[3])
    helix=float(line.split()[4])
    aggregation=float(line.split()[5])
    if beta > tangoThreshold or turn > tangoThreshold or helix > tangoThreshold or aggregation > tangoThreshold:
      tangoScores[position]=1
    position+=1
  #print str(beta) + tab + str(turn) + tab + str(helix) + tab + str(aggregation) 
  if verbose:
    print indent + "RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, tangoScores))	  
    data = [sequence,tangoScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
    
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
  #if verbose:
    #print indent + "LIMBO SEARCH:"
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
    for y in range(hitStart-1,hitStart+6):
      limboScores[y] += 1

  for x in range(0,len(sequence)):
    positionScores[x] += limboScores[x]
  if verbose:
    print indent + "RESULTS:"
    #print indent + sequence
    #print indent + ''.join(map(str, limboScores))	  
    data = [sequence,limboScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)









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
    #print indent + sequence
    #print indent + ''.join(map(str, tmhmmScores))	
    data = [sequence,tmhmmScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
      print indent + "|".join(str(word).ljust(col_width) for word in row)
  for x in range(0,len(sequence)):
    positionScores[x] += tmhmmScores[x]










  ######################################################################################
  ##########################       GENERAL SEQUENCE EVALUATION     #####################################
  ######################################################################################



def firstPartialEvaluation(sequence, positionScores, verbose):
	
	#SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
	input=open(inputsPath + "sequenceFASTA"  , "w")
	input.write(">gi" + endl)
	input.write(sequence)
	input.close()
	if verbose:
	   print endl
	   print indent + "*************************************"
	   print indent + "FIRST PARTIAL EVALUATION"
	###: BLAST SEARCH
	#if verbose:
	   #print endl
	   #print indent + "*************************************"
	   ##print indent + "STARTING BLAST SEARCH"
	#blastIt(sequence,positionScores,database, verbose)
        #if stepByStep:
	  #raw_input(indent + "Hit enter to continue with next evaluation")
	    
        ## IUPred evaluation
        if runIupred:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    #print endl
	    print indent + "*************************************"
	    print indent + "STARTING IUPred SEARCH"
	  iupred(sequence, positionScores, verbose)
	    
	    
        ## ANCHOR evaluation
	if runAnchor:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "STARTING ANCHOR SEARCH"
	  anchor(sequence, positionScores, verbose)
	
	#ELM search
	if runElm:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "STARTING ELM Search"
	  elmSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	
	#Net charge evaluation
	if evaluateNetCharge:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "Starting sequence net charge evaluation"
	  chargedSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")


        #PASTA evaluation (self aggregation)
        if runPasta:
          if stepByStep and verbose:
            raw_input(indent + "Hit enter to continue with next evaluation")
          if verbose:
            print endl
            print indent + "*************************************"
            print indent + "Starting PASTA evaluation"
          pastaSearch(sequence,positionScores, verbose)
          #if stepByStep:
            #raw_input(indent + "Hit enter to continue with next evaluation")
	
	#Prosite search
	if runProsite:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation") 
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "Prosite Search in progress..."
	  prositeSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	
	
	#LIMBO evaluation
	if runLimbo:
	  if stepByStep and verbose:
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
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	  #print indent + "*************************************"
	
	
	#TMHMM evaluation
	if runTmhmm:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "Search for transmembrane sections"
	  tmhmmEval(sequence,positionScores, verbose)
	  
	
	
	#Search amyloid pattern
	if runAmyloidPattern:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    #print indent + "EVALUATE AMYLOID FIBRIL FORMATION "
	    print indent + "*************************************"
	    #print endl
	    print indent + "Search for Amyloid sequence determinants"
	  amyloidPatternSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	
	#Waltz evaluation
	if runWaltz:
	  if stepByStep and verbose:
             raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "Starting Waltz evaluation"
	  waltzSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	
	#Tango evaluation
	if runTango:
	  if stepByStep and verbose:
	    raw_input(indent + "Hit enter to continue with next evaluation")
	  if verbose:
	    print endl
	    print indent + "*************************************"
	    print indent + "STARTING TANGO Search"
	  tangoSearch(sequence,positionScores, verbose)
	  #if stepByStep:
	    #raw_input(indent + "Hit enter to continue with next evaluation")
	
	
	if stepByStep and verbose:
	  print endl
	  raw_input(indent + "Press enter to see final results...")
	
	
	
        ##PRINT SCORE
        if verbose:
	  print indent + "*************************************"
	  print endl
	  print indent + "RESULTS OF FIRST PARTIAL EVALUATION:"
	  #print indent + sequence
	  data = [sequence,positionScores]
	  col_width = max(len(str(word)) for row in data for word in row)   # padding
	  for row in data:
	      print indent + "|".join(str(word).ljust(col_width) for word in row)
	  #print indent + ''.join(map(str, positionScores))
	  print indent + "SCORE:" + str(getGlobalScore(positionScores))
	  print indent + "*************************************"
	if stepByStep:
	  raw_input(indent + "....hit enter to continue")













def secondPartialEvaluation(sequence, positionScores, verbose):
	
	#SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
	input=open(inputsPath + "sequenceFASTA"  , "w")
	input.write(">gi" + endl)
	input.write(sequence)
	input.close()
	
	##: BLAST SEARCH
	if verbose:
	   print indent + "*************************************"
	   print indent + "SECOND PARTIAL EVALUATION"
	   print endl
	   print indent + "*************************************"
	   #print indent + "STARTING BLAST SEARCH"
	blastIt(sequence,positionScores,database, verbose)
        #if stepByStep:
	  #raw_input(indent + "Hit enter to continue with next evaluation")
	
	
	if stepByStep:
	  raw_input(indent + "Press enter to see final results...")
	
	
	
        ##PRINT SCORE
        if verbose:
	  print indent + "*************************************"
	  print endl
	  print indent + "RESULTS OF SECOND PARTIAL EVALUATION:"
	  data = [sequence,positionScores]
	  col_width = max(len(str(word)) for row in data for word in row)   # padding
    	  for row in data:
      	     print indent + "|".join(str(word).ljust(col_width) for word in row)

 	  #print indent + sequence
	  #print indent + ''.join(map(str, positionScores))
	  print indent + "SCORE:" + str(getGlobalScore(positionScores))
	  print indent + "*************************************"
	if stepByStep:
	  raw_input(indent + "....hit enter to continue")




















#def mutation(evaluationFunction):






















	



#***********************************************************************



## PRINT PROGRAM HELP 
def printHelp():
  print endl + "Usage: "
  print "  python bleach.py [options]" + endl
  print "Options are:"
  print tab + "--length  sequence-length                " + tab + "Define random sequence length"
  print tab + "--db   [swissprot | nr]			"   + tab + ""    #blast database
  print tab + "--composition  [average | user_specified]"  + tab + "Composition used to select AA (if user_specified you must define AA frequencies)"
  print tab + "--seq   predefined-sequence		" + tab + "Define sequence to start with"
  print tab + "--maxmutations  max-number		"   + tab + "Limit in number of ACCEPTED mutations(NOT MUTTATION ATTEMPTS)"   # 
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
    elif (arg=='--netcharge') and (index < len(sys.argv)):
      targetNetCharge = int(sys.argv[index+1])
      evaluateNetCharge=True
    
  ##   OUTPUT DETAIL  
    elif (arg=='--verbose'):
      verbose=True
    elif (arg=='--minoutput'):
      minimalOutput=True
    elif (arg=='--testoutput') and (index < len(sys.argv)):
      testing=True
      testOutputPath = sys.argv[index+1]
      
      
      
 ##   SELECT WHICH TOOLS WONT ME EVALUATED
    elif (arg=='--noblast'):
      runBlast=False
    elif (arg=='--notango'):
      runTango=False
    elif (arg=='--noelm'):
      runElm=False
    elif (arg=='--noiupred'):
      runIupred=False
    elif (arg=='--noanchor'):
      runAnchor=False
    elif (arg=='--noprosite'):
      runProsite=False
    elif (arg=='--nolimbo'):
      runLimbo=False
    elif (arg=='--notmhmm'):
      runTmhmm=False
    elif (arg=='--nopasta'):
      runPasta=False
    elif (arg=='--nowaltz'):
      runWaltz=False  
    elif (arg=='--noamyloidpattern'):
      runAmyloidPattern=False
   
   
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


#CHECK IF TARGET NET CHARGE IS POSSIBLE BASED ON SEQUENCE LENGTH (AND PH??)
  if evaluateNetCharge:
    if abs(targetNetCharge) > length:
      print 'Net charge is impossible to reach with the specified sequence length'
      exit()
  
  
  
  print "************************************************"
  print "************************************************"
  print "EXECUTION PARAMETERS:"
  print 'Id=' + str(exeId) 
  print "Length=" + str(length) 
  print "Composition=" + composition
  print "Sequence=" + sequence
  if evaluateNetCharge:
    print "Target net charge=" + str(targetNetCharge)
  print "************************************************"
  print "************************************************"




#MAKE INPUT FOLDER
os.mkdir(inputsPath)
os.mkdir(outputsPath)

#OUTPUT
#outputPath = "Output/"   


###########THESE FILES SHOULD GO TO test/results DIRECTORY ***************************
#scoresFile="scores" + str(length)   #save Scores Vs iteration number 
#mutAttemptsFile="mutationsAttempt" + str(length)  # save number of mutation attempts  Vs iteration number
#timesFile='times' + str(length)   #save times Vs iteration number
#totalTimesFile='totalTimes'  #save total time elapsed Vs sequence length
#if output:
  #totalTimesOutputFile=open( testOutputPath + totalTimesFile , "a")
  #timesOutputFile=open( testOutputPath + timesFile , "a")
  #scoresOutputFile=open( testOutputPath + scoresFile , "a")
  #mutationsFile=open( testOutputPath  + mutAttemptsFile , "a")
#*************************************************

#TESTING MODE ON: WRITE EXECUTION PARAMETERS (TOTAL TIME, TIME PER BLOCK, MUT-ATTEMPTS, BETA VALUES, SCORES)
if testing: 
  testOutputFile=open(testOutputPath + '/'+ str(exeId), 'w')
  if rand:
    testOutputFile.write('RAND'+ tab +str(beta) + tab + str(length) + endl)
  else:
    testOutputFile.write('SEQ'+ tab + str(beta) + tab + str(length) + endl)

if minimalOutput:
  #CREATE .log FILE
  logFileStream=open( outputsPath+logFileName, 'w')


if uvsilent==True:   
  aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]
else:
  aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",0),("P",470),("S",656),("T",0),("W",108),("Y",0),("V",687) ]





# FORMAT TO REQUEST RANDOM SEQUENCE:   
#http://web.expasy.org/cgibin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

if rand==True:
	#****************GET RANDOM SEQUENCE*************
        #print endl
        if verbose:
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
	#if verbose:
	   #print "*******************************" 
	

   
#print endl
#PRINT STARTING SEQUENCE (BEFORE ANY MUTATION)   
if rand==True:
  if verbose:
    print "INITIAL SEQUENCE:    " + sequence
  if stepByStep:
    raw_input("Hit enter to start initial evaluation")
    

#CREATE ARRAY TO SAVE MUTATION FREQUENCY
positionScores=[]
partialScores=[]
mutatedScores=[]
for p in range(len(sequence)):
  positionScores.append(0)
  mutatedScores.append(0)
  partialScores.append(0)


time0 = time.time()   #start time



################################
#### ITERATE OVER SEQUENCE ####
#############################

iteration=1
globalIteration=0


indent=""

for p in range(len(sequence)):
  #positionScores.append(0)
  positionScores[p]=0
  partialScores[p]=0



if verbose:
  print endl
  print "*****************************"
  print " INITIAL EVALUATION "
  print "*****************************"    



##################################
#########EVALUATE INITIAL SEQUENCE
##################################

#MAKE BOTH PARTIAL EVALUATIONS TO GET A GLOBAL SCORE

#FIRST SET OF EVALUATIONS
firstPartialEvaluation(sequence, partialScores, verbose)

#SAVE RESULTS
firstPartialScore=getGlobalScore(partialScores)
firstPartialScores=partialScores

#ADD THE SCORE TO THE GLOBAL SCORE AND RESET PARTIAL LIST
for p in range(len(sequence)):
  positionScores[p]=positionScores[p]+partialScores[p]
  partialScores[p]=0
  
secondPartialScore = 0  
if runBlast:
  #SECOND PART OF EVALUATION
  secondPartialEvaluation(sequence, partialScores,verbose)
  secondPartialScore=getGlobalScore(partialScores)

  #ADD THE SCORE TO THE GLOBAL SCORE AND RESET PARTIAL LIST
  for p in range(len(sequence)):
    positionScores[p]=positionScores[p]+partialScores[p]
    partialScores[p]=0

##SUM OF SCORES LIST
globalScore=getGlobalScore(positionScores)

if verbose:
  print "*******************************************"
  print "INITIAL EVALUATION RESULTS"
  print "First partial score   : "  + str(firstPartialScore)
  print "Second partial score  : " + str(secondPartialScore) 
  print "Global score          : " + str(globalScore)
  print endl
  print "*******************************************"
  print "*****************************"
  print "*****************************"
  print endl
  print endl
  if stepByStep:
    raw_input("Hit enter to start mutations")

if minimalOutput:
  logFileStream.write('ISEQ' + tab + sequence + tab + str(globalScore) + endl)
    #print 'INITIAL SEQ:   ' + sequence + tab + str(globalScore)

if testing:
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'FIRST ' + str(firstPartialScore) + endl)
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'SECOND ' + str(secondPartialScore) + endl)
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'GLOBAL ' + str(globalScore) + endl)
  







######################################################
######################################################
#########  GLOBAL LOOP  ##############################
######################################################
######################################################


while globalScore > 0 and iteration <= maxIterations:
  if verbose:
    print "*****************************"
    print "STARTING GLOBAL ITERATION " + str(globalIteration)
    print "*****************************"
    data = [sequence,positionScores]
    col_width = max(len(str(word)) for row in data for word in row)   # padding
    for row in data:
    	print indent + "|".join(str(word).ljust(col_width) for word in row)
    #print "Current sequence:         " + sequence
    #print "Current total scores:     " + ''.join(map(str, positionScores))
    print "Current global score:     " + str(globalScore)
  if stepByStep: 
    if verbose:
	raw_input("Hit enter to start first round of evaluations")
    else:
	raw_input('')
 
 
    
  
 
 
 
 
 
  
  #################################################
  ######  FIRST ROUND OF EVALUATIONS - MUTATIONS
  #################################################
  
  if verbose:
    print "FIRST ROUND OF MUTATIONS: DECISION IS BASED ON (RESULTS OF)FIRST SET OF TOOLS" 
  
  partialScore=firstPartialScore
  while partialScore > 0 and iteration <= maxIterations:
      timePrev=time.time()
      weighted=[]
      #weighted IS A PAIRLIST(position,weight)
      #CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
      for x in range(len(positionScores)):
	weighted.append((x, positionScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
      
      
      
	
      mutAttempts=0       #COUNT MUTATIONS ATTEMPTS
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
	  firstPartialEvaluation(mutatedSequence, mutatedScores, verbose)	 
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
	    data = [sequence,positionScores]
	    col_width = max(len(str(word)) for row in data for word in row)   # padding
    	    for row in data:
      		print indent + "|".join(str(word).ljust(col_width) for word in row)
	    #print indent + sequence
	    #print indent + ''.join(map(str, positionScores))
	    print indent + "Global score: " + str(partialScore)
	    print ""
	    print indent + "Mutated sequence"
	    data = [sequence,mutatedScores]
    	    col_width = max(len(str(word)) for row in data for word in row)   # padding
    	    for row in data:
      		print indent + "|".join(str(word).ljust(col_width) for word in row)

	    #print indent + mutatedSequence
	    #print indent + ''.join(map(str, mutatedScores))
	    print indent + "Global score: " + str(mutatedScore)
	    print ""
	  if partialScore >= getGlobalScore(positionScores):
	      if verbose:
		    print indent + "Previous score (" + str(partialScore) + ") >= Mutated score (" + str(mutatedScore) + ")" 
		    print indent + "...ACCEPT MUTATION"
	      if stepByStep:
		raw_input("")
	      break
		#raw_input(indent + "Hit enter to continue with next iteration")
	      #return mutatedSequence
	  else:
	      #DECISION BASED ON MONTE CARLO
	      if verbose:	
		    print indent + "Previous score (" + str(partialScore) + ") < Mutated score (" + str(mutatedScore) + ")" 
	      diff=partialScore-getGlobalScore(mutatedScores)
	      if verbose:
		    print indent + "SCORE DIFFERENCE:" + str(diff)
	      exponent=diff/beta
	      MCvalue=math.exp(exponent)
	      if verbose:
		    #print "  :" + str(exponent)
		    print indent + "MC VALUE:" + str(MCvalue)     #e^(dif/beta)
	      
	      #GENERATE RANDOM NUMBER BETWEEN 0 AND 1
	      randy=random.random()
	      if verbose:
		    print indent + "RANDOM VALUE [0,1]:" + str(randy)
	      #print indent + "MONTE CARLO DECISION:"
	      if MCvalue > randy:
		#ACCEPT MUTATION
		if verbose:
		    print indent + "...ACCEPT MUTATION"
		if stepByStep:
		  raw_input("")
		  #raw_input(indent + "Hit enter to continue with next iteration")
		#return mutatedSequence
		break
	      else:	    
		if verbose:
		    #print "El score original " + str(getGlobalScore(mutatedScores)) + " es mayor que "+ str(getGlobalScore(positionScores))
		    print indent + " Mutation score (" + str(mutatedScore) + ") >= Previous score (" + str(partialScore) + ")" 
		    print indent + "...DENY MUTATION"
		if stepByStep:
		  raw_input(indent + "Hit enter to continue with next attempt")
		#return sequence
		#break
	  if verbose:
	    print indent + "*************************************"
	


      
    
      #####END OF MUTATION ITERATION
      if mutAttempts < 10000:   #MAKE SURE LOOP ENDED BY MUTATION ACCEPT
	    #print "Sequence after mutation:    " + mutatedSequence
	    ###NOW THE SEQUENCE IS THE MUTATED SEQUENCE
	    sequence=mutatedSequence
	    #AND THE POSITION SCORES ARE THE ONES CORRESPONDING TO THE MUTATED SEQUENCE
	    positionScores=mutatedScores
	    #AND THE GLOBAL SEQUENCE SCORE IS THE ONE CORRESPONDING TO THIS NEW SEQUENCE
	    partialScore= getGlobalScore(positionScores)
	    if verbose:
	      print endl
	      print "Attempts before mutation accept:" + str(mutAttempts)
	      #print endl
	      print "*******************************************"
	    
	      
	      print "End of (PARTIAL) iteration " + str(iteration)
	      print "(PARTIAL) score :    " + str(partialScore)
	      print "*******************************************"
	      print endl
	    if minimalOutput:
		logFileStream.write(str(mutatePosition) + '(' + previousResidue + ')->' + seleccionado + '   '+ tab + sequence + tab +str(partialScore) + tab + '1' + endl)
		  #print  mutatePosition + '(' + previousResidue + ') -> ' + seleccionado + tab + sequence
      else:
	   if verbose:
	      ### EXCEEDED THE NUMBER OF ATTEMPTS, SEQUENCE NOT CHANGED
	      print  " Mutations attempts exceeded " 
      
      if testing:
	timeX=time.time()-timePrev   #ITERATION TIME
	testOutputFile.write('LOOP1' + tab + str(iteration)+ tab + str(mutAttempts) + tab + str(partialScore) + tab + str(timeX) + endl ) 
	
      if stepByStep:
	raw_input("Hit enter to continue with next iteration")
		
      iteration=iteration+1   #TOTAL NUMBER OF ITERATIONS
      
      
      
      
      


  #RESET LIST OF PARTIAL SCORES
  for p in range(len(sequence)):
    partialScores[p]=0


  if runBlast:
    ####################################################	
    ##### SECOND ROUND OF EVALUATION - MUTATION  #######
    ####################################################
    
    
    secondPartialEvaluation(sequence,partialScores, verbose)	
    partialScore=getGlobalScore(partialScores)  
    while partialScore > 0 and iteration <= maxIterations:
	timePrev=time.time()
	weighted=[]
	#weighted IS A PAIRLIST(position,weight)
	#CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
	for x in range(len(partialScores)):
	  weighted.append((x, partialScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
	
	
	#MUTATION ATTEMPTS
	#INSIDE THIS LOOP, THE VARIABLE partialScore IS NOT MODIFIED. ONLY COMPARED TO THE MUTATED SCORES
	#THE LOOP ENDS WHEN A MUTATION IS ACCEPTED OR WHEN A LIMIT OF ATTEMPTS IS REACHED
	mutAttempts=0       #COUNT MUTATIONS ATTEMPTS
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
	    secondPartialEvaluation(mutatedSequence, mutatedScores, verbose)	 
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
	      data = [sequence,partialScores]
	      col_width = max(len(str(word)) for row in data for word in row)   # padding
    	      for row in data:
      		print indent + "|".join(str(word).ljust(col_width) for word in row)
	      #print indent + sequence
	      #print indent + ''.join(map(str, partialScores))
	      
	      print indent + "Global score: " + str(partialScore)
	      print ""
	      print indent + "Mutated sequence"
	      print indent + mutatedSequence
	      data = [sequence,mutatedScores]
	      col_width = max(len(str(word)) for row in data for word in row)   # padding
    	      for row in data:
      		print indent + "|".join(str(word).ljust(col_width) for word in row)
              #print indent + ''.join(map(str, mutatedScores))
	      #print indent + "Global score: " + str(mutatedScore)
	      print ""
	    if partialScore >= mutatedScore:
		if verbose:
		      print indent + "Previous score (" + str(partialScore) + ") >= Mutated score (" + str(mutatedScore) + ")" 
		      print indent + "...ACCEPT MUTATION"
		if stepByStep:
		  raw_input("")
		break
		  #raw_input(indent + "Hit enter to continue with next iteration")
		#return mutatedSequence
	    else:
		#DECISION BASED ON MONTE CARLO
		if verbose:	
		      print indent + "Previous score (" + str(partialScore) + ") < Mutated score (" + str(mutatedScore) + ")" 
		diff=partialScore-getGlobalScore(mutatedScores)
		if verbose:
		      print indent + "SCORE DIFFERENCE:" + str(diff)
		exponent=diff/beta
		MCvalue=math.exp(exponent)
		if verbose:
		      #print "  :" + str(exponent)
		      print indent + "MC VALUE:" + str(MCvalue)     #e^(dif/beta)
		
		#GENERATE RANDOM NUMBER BETWEEN 0 AND 1
		randy=random.random()
		if verbose:
		      print indent + "RANDOM VALUE [0,1]:" + str(randy)
		#print indent + "MONTE CARLO DECISION:"
		if MCvalue > randy:
		  #ACCEPT MUTATION
		  if verbose:
		      print indent + "...ACCEPT MUTATION"
		  if stepByStep:
		    raw_input("")
		    #raw_input(indent + "Hit enter to continue with next iteration")
		  #return mutatedSequence
		  break
		else:	    
		  if verbose:
		      #print "El score original " + str(getGlobalScore(mutatedScores)) + " es mayor que "+ str(getGlobalScore(positionScores))
		      print indent + " Mutation score (" + str(mutatedScore) + ") >= Previous score (" + str(partialScore) + ")" 
		      print indent + "...DENY MUTATION"
		  if stepByStep:
		    raw_input(indent + "Hit enter to continue with next attempt")
		  #return sequence
		  #break
	    if verbose:
	      print "*************************************"
	  


	
      
	#####END OF MUTATION ITERATION

	#MAKE SURE LOOP ENDED BY MUTATION ACCEPT
	if mutAttempts < 10000:   
	      #print "Sequence after mutation:    " + mutatedSequence
	      ###NOW THE SEQUENCE IS THE MUTATED SEQUENCE
	      sequence=mutatedSequence
	      #AND THE SCORES ARE THE ONES CORRESPONDING TO THE MUTATED SEQUENCE
	      partialScores=mutatedScores
	      #AND THE GLOBAL SEQUENCE SCORE IS THE ONE CORRESPONDING TO THIS NEW SEQUENCE
	      partialScore= getGlobalScore(partialScores)
	      if verbose:
		print endl
		print "Attempts before mutation accept:" + str(mutAttempts)
		#print endl
		print "*******************************************"
	      
		
		print "End of (PARTIAL) iteration " + str(iteration)
		print "(PARTIAL) score :    " + str(partialScore)
		print "*******************************************"
		print endl
	      if minimalOutput:
		  logFileStream.write(str(mutatePosition) + '(' + previousResidue + ')->' + seleccionado + tab + sequence + endl)
		  #print  mutatePosition + '(' + previousResidue + ') -> ' + seleccionado + tab + sequence
	else:
	      ### EXCEEDED THE NUMBER OF ATTEMPTS, SEQUENCE NOT CHANGED
	      if verbose:
		print  " Mutations attempts exceeded " 
	if testing:
	  timeX=time.time()-timePrev   #ITERATION TIME
	  testOutputFile.write('LOOP2' + tab + str(iteration)+ tab + str(mutAttempts) + tab + str(partialScore) + tab + str(timeX) + endl ) 
	
	  # MUTATION ATTEMPTS Vs. ITERATION NUMBER
	  #mutationsFile.write(str(iteration)+ tab + str(mutAttempts) + tab +str(beta) + endl)
	  
	  #SCORES Vs. ITERATION NUMBER
	  #scoresOutputFile.write(str(iteration)+ tab + str(partialScore) + tab + str(beta) + endl)
	  
	  #ITERATION ELAPSED TIME
	  #timesOutputFile.write(str(iteration)+ tab + str(timeX) + tab + str(beta) + endl)
	if stepByStep:
	  raw_input("Hit enter to continue with next iteration")
		  
	iteration=iteration+1   #TOTATL NUMBER OF ITERATIONS
	
  
    if verbose:
	print "SECOND ROUND OF EVALUATIONS FINISHED"
	if stepByStep:
	  raw_input("HIT ENTER TO GET CURRENT GLOBAL SCORE")
   
   
  ##AT THE END OF THE SECOND ROUND, THE SCORE CORRESPONDING TO THE SECOND PARTIAL EVALUATION IS 0 (EXCEPT WHEN THE NUMBER OF ITERATIONS IS EXCEEDED)....
  #THEN, GLOBAL SCORE IS DEFINED BY THE FIRST PARTIAL EVALUATION
       
  for p in range(len(sequence)):
    positionScores[p]=0
  firstPartialEvaluation(sequence, positionScores , False )
  globalScore=getGlobalScore(positionScores)
  globalIteration=globalIteration+1;  
  
  #PRINT RESULTS OF GLOBAL ITERATION
  if verbose:
    print "*******************************************"	     
    print "End of global iteration " + str(globalIteration)
    print "Global score :    " + str(globalScore)
    print "*******************************************"
    print endl
  #if minimalOutput:
   # logFileStream.write()
  if stepByStep:
	raw_input("Hit enter to continue with next iteration")
  if testing:
    testOutputFile.write('GLOBAL ' + tab + str(globalIteration) + endl ) 
	



###########################################
#### GLOBAL LOOP END ####################
#########################################

print "**END OF SEARCH**"
if globalScore==0:
  print "REACHED SCORE = 0"
  print sequence
else:
  print "REACHED LIMIT OF ITERATIONS"
  data = [sequence,positionScores]
  col_width = max(len(str(word)) for row in data for word in row)  # padding
  print 'FINAL SCORES'
  for row in data:
      print  "|".join(str(word).ljust(col_width) for word in row)
  print "Global score: " + str(globalScore)

  #print "Final sequence: " + sequence
  #print "Final score:    " + ''.join(map(str, positionScores))
  #print "Global score: " + str(globalScore)
  



#PERFORMANCE TESTS OUTPUT
totalElapsedTime=time.time() - time0   ##
#if output:        
  #print "Elapsed time:",totalElapsedTime, "Seconds"
  #totalTimesOutputFile.write(str(length) + tab + str(iteration-1)+ tab + str(totalElapsedTime) + tab + str(beta) +  endl)
if testing:
  testOutputFile.write('END' + tab + str(totalElapsedTime) + tab + str(globalScore) + endl) 
  testOutputFile.close()
  ##CLOSE ALL OUTPUT FILES
  #totalTimesOutputFile.close()
  #timesOutputFile.close()
  #scoresOutputFile.close()
  #mutationsFile.close()
   

