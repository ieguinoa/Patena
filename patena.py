#!/usr/bin/python
import urllib2
import StringIO
import sys
import os
# import re
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
#import subprocess
import shutil
import tool_functions

##TODO: CHECK NECESSARY ENVIRONMENT VARS:
#ANCHOR_PATH=$DIR/Tools/ANCHOR
#IUPred_PATH=$DIR/Tools/iupred
#PASTA_PATH=$DIR/Tools/PASTA/pasta_exe
#PROSITE=$DIR/Tools/Prosite/ps_scan


#****************************************

#get base path
def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


#****************************************

# WEIGHTED SELECTION
# ANON. FUNCTION: USE weighted_choide(param) TO CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
# RETURNS THE ID OF THE SELECTED ELEMENT
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
maxIterations=4000
cutoff=0.01  #BLAST cutoff
waltzThreshold=79.0
beta=0.5   #MC
targetScore=0.0
pastaThreshold=-5.5  #ENERGY Threshold
pastaProbabilityThreshold=0.05 #Aggregation probability threshold
iupredThreshold=0.5
anchorThreshold=0.5
tangoThreshold=1.0
targetNetCharge=0




### TODO:
# Make use of this  dictionary with the parameters to the set of tools and then get rid of the infividual vars (e.g pastaThershold, etc)
# step 1: when parsing parameters then modify this config_parms dict and not the individual vars.

##Default values
config_params= { 'iupredThreshold': iupredThreshold,
                'anchorThreshold':anchorThreshold,
                'waltzThreshold':waltzThreshold,
                'pastaThreshold':pastaThreshold,
                'targetNetCharge':targetNetCharge,
                'tangoThreshold':tangoThreshold
               }



#  EXECUTION PARAMETERS
testTimes=False   #True=print times of each part
minimalOutput=False  #True = only print global scores at end of iteration
verbose=False        #True = print detailed information of execution
stepByStep=False
match=False
rand=True
testing=False
change=True
evaluateNetCharge=False
blastWeb=False  #BLAST SEARCH LOCAL OR WEB
output=True  ##print info to file
mutAttempts=0
length=12   #defaul sequence length
global_evaluation=False  #when True, just make a general evaluation of the sequence: run all tools(from each loop) and print the results of each 
detailed_output=False

# EXECUTION TIMES OF DIFFERENT PARTS
times_dict={
    'Tango':0.0,
    'Pasta':0.0,
    'Waltz':0.0,
    'ELM':0.0,
    'Prosite':0.0,
    'Limbo':0.0,
    'Tmhmm':0.0,
    'IUpred':0.0,
    'Anchor':0.0,
    'Amyloid Pattern':0.0,
    'Net charge':0.0
}




#EVALUATION PARAMETERs (TEST MODE)



#All tools enabled by default 
runBlast=runTango=runPasta=runWaltz=runElm=runProsite=runLimbo=runTmhmm=runIupred=runAnchor=runAmyloidPattern=True

#TODO CHANGE TO USAGE OF EXEUTION SET
execution_set={'Tango','Pasta','Waltz','ELM','Prosite','Limbo','Tmhmm','IUpred','Anchor','Amyloid Pattern','Net charge'}

#FILES
logFileName='mutations' + str(exeId) + '.log'


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



#**************************************************************************

##    JUST SUM UP THE INDIVIDUAL SCORES
def getGlobalScore(scoresList):
  score=0.0
  for listIndex in range(len(scoresList)):
    score= score + scoresList[listIndex]
  return score



#**************************************

def print_evaluation_time(total_elapsed_time,times_dict):
    # before i was doing this calculation:
    # evaluationTime=totalElapsedTime-(pastaTime + anchorTime + tangoTime + blastTime + iupredTime + waltzTime + elmTime + prositeTime + tmhmmTime + limboTime)
    # not sure what was that idea, it doesnt makes sense
    print 'Total Elapsed Time: ' + total_elapsed_timed
    evaluation_time=0
    for key,value in times_dict.items():
        evaluation_time+=value
        print 'Time spent in '+ key +': ' + value
        # print 'Time spent in ANCHOR: ' + anchorTime
        # print 'Time spent in TANGO: ' + tangoTime
        # print 'Time spent in BLAST: ' + blastTime
        # print 'Time spent in IUPred: ' + iupredTime
        # print 'Time spent in Waltz: ' + waltzTime
        # print 'Time spent in ELM: ' + elmTime
        # print 'Time spent in Prosite: ' + prositeTime
        # print 'Time spent in TMHMM: ' + thmmTime
        # print 'Time spent in Limbo: ' + limboTime
    print 'Time spent on evaluations: ' +  evaluation_time





def print_execution_params(exeId,beta,length,composition,sequence,evaluateNetCharge,targetNetCharge)
    #####   ALWAYS PRINT GENERAL PARAMETERS OF EXECUTION
    #print endl
    print "************************************************"
    print "************************************************"
    print "EXECUTION PARAMETERS:"
    print 'Id=' + str(exeId) 
    print 'Beta= '+ str(beta)
    print "Length=" + str(length) 
    print "Composition=" + composition
    print "Sequence=" + sequence
    if evaluateNetCharge:
      print "Target net charge=" + str(targetNetCharge)
    print "************************************************"
    print "************************************************"





  ######################################################################################
  ##########################       GENERAL SEQUENCE EVALUATION     #####################################
  ######################################################################################



def firstPartialEvaluation(sequence, config_params,positionScores, verbose):
	#SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
	input=open(inputsPath + "sequenceFASTA"  , "w")
	input.write(">gi" + endl)
	input.write(sequence)
	input.close()
	if verbose:
	   print endl
	   print indent + "*************************************"
	   print indent + "FIRST PARTIAL EVALUATION"
        for tool_name in tool_functions.tool_functions_dict.keys():
            if tool_name in execution_set:  ## this set has the set of enabled tools to run
                timePrev=time.time()
                if stepByStep and verbose:
                    raw_input(indent + "Hit enter to continue with next evaluation")
                print indent + "STARTING "+ tool_name +" execution"
                tool_functions.tool_functions_dict[tool_name](sequence, positionScores,config_params, inputsPath,outputsPath,verbose,detailed_output)
                if detailed_output:
                    detailedOutFile.write('\n\n***********\n\n' )
                times_dict[tool_name]+=(time.time() - timePrev)
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
	timePrev=time.time()
	tool_functions.blastIt(sequence,positionScores,database,inputsPath, verbose,detailed_output)
	#blastTime+=(time.time() - timePrev)
        #if stepByStep:
	  #raw_input(indent + "Hit enter to continue with next evaluation")
	if detailed_output:
                detailedOutFile.write('\n\n***********\n\n' )
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








#***********************************************************************



## PRINT PROGRAM HELP 
def printHelp():
  print endl + "Usage: "
  print "  python patena.py [options]" + endl
  print "Options are:"
  print tab + "--length  sequence-length                " + tab + "Define random sequence length"
  #print tab + "--db   [swissprot | nr]			"   + tab + ""    #blast database
  #print tab + "--composition  [average | user_specified]"  + tab + "Composition used to select AA (if user_specified you must define AA frequencies)"
  print tab + "--seq   predefined-sequence		" + tab + "Define sequence to start with"
  print tab + "--maxmutations  max-number		"   + tab + "Limit in number of ACCEPTED mutations(NOT MUTTATION ATTEMPTS)"   # 
  #print tab + "--maxiterations"   ?????    SUM OF MUTATION ATTEMPTS ?????     ******MAINLY FOR PERFORMANCE REASONS











#**************************
#**************************
#******* MAIN *************
#**************************





#********DEFAULTS********
composition="average"
#a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=-1
userComposition={"A":-999 , "R":-999  , "N":-999  , "D":-999  , "C":-999  , "E":-999 , "Q":-999  , "G":-999  , "H":-999  , "I":-999  , "L":-999  , "K":-999  ,"M":-999  , "F":-999  , "P":-999  , "S":-999  , "T":-999  , "W":-999  , "Y":-999  , "V":-999 }
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
    elif (arg=='--noMutations'):
      global_evaluation=True
    elif (arg=='--db') and (index < len(sys.argv)):
      database = sys.argv[index+1]
    elif (arg=='--blastweb'):
      blastWeb=True
    elif (arg=='--uvsilent'):
      uvsilent=True 
    elif (arg=='--netcharge') and (index < len(sys.argv)):
      targetNetCharge = int(sys.argv[index+1])
      evaluateNetCharge=True


### TODO: CREATE A DICT OR SET WITH THE TOOLS THAT SHOULD BE APPLIED
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



    #CHECK IF THE USER DEFINED ANY OF THE AAs FREQUENCIES
    #elif (arg== '--composition') and (index < len(sys.argv)):
      #composition = sys.argv[index+1]
      #if (composition=="user_specified"):  #frequencies specified by parameter
	#for j in range(index+2,len(sys.argv),2):
    elif(arg=='-a'):
      userComposition['A']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-r'):
      userComposition['R']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-n'):
      userComposition['N']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-d'):
      userComposition['D']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-c'):
      userComposition['C']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-q'):
      userComposition['Q']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-e'):
      userComposition['E']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-g'):
      userComposition['G']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-h'):
      userComposition['H']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-i'):
      userComposition['I']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-l'):
      userComposition['L']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-k'):
      userComposition['K']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-m'):
      userComposition['M']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-f'):
      userComposition['F']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-p'):
      userComposition['P']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-s'):
      userComposition['S']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-t'):
      userComposition['T']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-w'):
      userComposition['W']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-y'):
      userComposition['Y']=int(sys.argv[index+1])
      composition=="user_specified"
    elif(arg=='-v'):
      userComposition['V']=int(sys.argv[index+1])
      composition=="user_specified"



###   OUTPUT DETAILS
    elif (arg=='--verbose'):
      verbose=True 
    elif (arg=='--detailed') and (index < len(sys.argv)): # print detailed output for each (to output file) 
      detailed_output=True
      detailedOutFile=open(sys.argv[index+1],'w') 	
    elif (arg=='--minoutput') and (index < len(sys.argv)):
      minimalOutput=True
      logsPath= sys.argv[index+1]
    elif (arg=='--testoutput') and (index < len(sys.argv)):
      testing=True
      testOutputPath = sys.argv[index+1]
    elif (arg=='--gettime'):
      testTimes=True   
    elif (arg=='--jobid') and (index < len(sys.argv)):
      exeId=sys.argv[index+1]



#CHECK IF TARGET NET CHARGE IS POSSIBLE BASED ON SEQUENCE LENGTH (AND PH??)
if evaluateNetCharge:
  if abs(targetNetCharge) > length:
      print 'Net charge is impossible to reach with the specified sequence length'
      exit()
  
  



#PATHS
basePath=getScriptPath() + '/'
# toolsPath=basePath + 'Tools/'    #**************************TODO SET THE PATH TO THE TOOL SET 
inputsPath=basePath + "/Input/"+ str(exeId) + "/" #SET PATH TO SAVE INPUTS FILES
baseOutputPath=basePath + "/Output/" 
outputsPath=baseOutputPath + str(exeId) + "/"
testOutputPath=outputsPath   # DEFAULT OUTPUT FOR TESTs 
logsPath=outputsPath #default path for log files



##### TRY TO CREATE INPUT AND OUTPUT DIRS
try:
    os.makedirs(inputsPath)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(inputsPath):
        # if exists, just remove it and create it again
        shutil.rmtree(inputsPath)
        os.makedirs(inputsPath)
        pass


try:
    os.makedirs(outputsPath)
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(inputsPath):
        # if exists, just remove it and create it again
        shutil.rmtree(outputsPath)
        os.makedirs(outputsPath)
        pass








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
  logFileStream=open( logsPath+'/'+logFileName, 'w')




#AMINOACIDS FREQUENCIES....
# THESE FREQUENCIES ARE USED TO SELECT REPLACEMENTS DURING MUTATIONS
# THE SELECTION IS MADE USING A WEIGTHED SELECTION. THE ONLY REQUERIMENT IS THAT THE FREQUENCIES (WEIGHTS) ARE WHOLE NUMBERS

# STANDARD COMPOSITION 
#if uvsilent==False:   
standardComposition={"A":825 , "R":553 , "N":406 , "D":545 , "C":137 , "E":393 , "Q":675 , "G":707 , "H":227 , "I":596 , "L":966 , "K":548 ,"M":242 , "F":386 , "P":470 , "S":656 , "T":534 , "W":108 , "Y":292 , "V":687}
  #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]
#else:
  #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",0),("P",470),("S",656),("T",534),("W",0),("Y",0),("V",687) ]

#print str(standardComposition)

if uvsilent:
# RESET Y,W,F FREQ = 0 
    userComposition['Y']=0
    userComposition['W']=0
    userComposition['F']=0


if (composition=="user_specified"):   
# USER HAS DEFINED AT LEAST ONE OF THE FREQUENCIES, THE FREQUENCIES DEFINED ARE IN
#FIRST CHECK IF THE SUM OF FREQUENCIES DEFINED IS LESS THAN 100 percent
  freqSum=0
  for key in userComposition:
    if userComposition[key] != -999:    #THE USER HAS DEFINED THIS FREQUENCE
      freqSum+= userComposition[key]
  if freqSum > 100:
    print 'Total defined frequencies exceeded 100%'
    exit()


#REDEFINE TOTAL aaFrequencies USING THE USER DEFINED FREQUENCIES
#TODO ***********************************************************
# TEMP. SOLUTION: JUST COPY THE FREQUENCIES DEFINED BY USER AND REPLACE UNKNOWN BY STANDARD
for key in userComposition: 
    if userComposition[key] == -999:
      userComposition[key] = standardComposition[key]
    else: 
      userComposition[key] = userComposition[key]*100  #change the frequencies to 10000 base



#CONVERT FREQUENCIES TO LIST OF PAIRS (AA, FREQ)
aaFrequencies=userComposition.items()
#print str(aaFrequencies)
#aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]




# FORMAT TO REQUEST RANDOM SEQUENCE:   
#http://web.expasy.org/cgibin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

#if rand==True:
	##****************GET RANDOM SEQUENCE*************
        ##print endl
        #if verbose:
	  #print "Generating random sequence..."    
	##print endl
	#if not (composition=="user_specified"):
	  #url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(length) + "&comp=" + composition + "&output=fasta"
	#else:
	  #print "fix this"
	##print url
  	#response = urllib2.urlopen(url)
	#html = response.read()
	#i = html.index('\n')
	#sequence = html[i+1:].replace('\n', '')
	##if verbose:
	   ##print "*******************************" 
	

#GENERATE RANDOM SEQUENCE WITH THE DEFINED COMPOSITION
if rand==True:
	sequence=[]
	for x in range(0,length):
	  sequence.append(weighted_choice(aaFrequencies)) 
	sequence=''.join(sequence)

  



print_execution_params(exeId,beta,length,composition,sequence,evaluateNetCharge,targetNetCharge)


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

timePrev=time.time()  #used to measure execution times of different parts (evaluation, mutations, etc)



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
firstPartialEvaluation(sequence, config_params,partialScores, verbose)

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

if detailed_output:
	index=0
        detailedOutFile.write('\n')
        detailedOutFile.write('*************************\n')
        detailedOutFile.write('\n')
        #detailedOutFile.write("First partial score"+ tab  + str(firstPartialScore) + '\n')
	#detailedOutFile.write("Second partial score" +tab+ str(secondPartialScore) + '\n')
	detailedOutFile.write("Global score"+tab + str(globalScore) + '\n')
 	#detailedOutFile.write('\n')
        #detailedOutFile.write('*************************\n')
	detailedOutFile.write('\n')
        detailedOutFile.write('Scores per position:\n')
        detailedOutFile.write('Pos' +tab+'AA' +tab+ 'Score\n')
	for aa in sequence:
		detailedOutFile.write(str(index) + tab + aa + tab + str(positionScores[index]))
		index+=1
		detailedOutFile.write('\n')
	#for score in positionScores:
        #        detailedOutFile.write(str(score) + tab)
        #data = [sequence,positionScores]
        #col_width = max(len(str(word)) for row in data for word in row)  # padding
        #for row in data:
        #      detailedOutFile.write("|".join(str(word).ljust(col_width) for word in row))


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
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'FIRST'  + tab + str(firstPartialScore) + endl)
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'SECOND' + tab + str(secondPartialScore) + endl)
  testOutputFile.write('ISEQ' + tab + sequence + tab + 'GLOBAL' + tab + str(globalScore) + endl)
  

















######################################################
######################################################
#########  GLOBAL LOOP  ##############################
######################################################
######################################################


while globalScore > 0 and iteration <= maxIterations and (not global_evaluation):
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
      weights=[]
      #weights IS A PAIRLIST(position,weight)
      #CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
      for x in range(len(positionScores)):
	weights.append((x, positionScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
      
      
      
	
      mutAttempts=0       #COUNT MUTATIONS ATTEMPTS
      while 1000 > mutAttempts:    ##JUST A SYMBOLIC MAX. AMOUNT OF MUTATIONS ATTEMPTS
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
	  mutatePosition= weighted_choice(weights) 
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
	  firstPartialEvaluation(mutatedSequence, config_params, mutatedScores, verbose)	 
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
            
              ###  THIS SHOULD BE EXTRACTED to a function monte_carlo_eval() that returns true or false 
              # IT IS ALSO REPEATED IN THE 2ND ROUND OF EVALUATIONS
              
	      #DECISION BASED ON MONTE CARLO
	      if verbose:	
		    print indent + "Previous score (" + str(partialScore) + ") < Mutated score (" + str(mutatedScore) + ")" 
	      diff=partialScore-getGlobalScore(mutatedScores)
	      if verbose:
		    print indent + "SCORE DIFFERENCE:" + str(diff)
	      exponent=diff/beta
	      if exponent<-100:   #SATURATION
		MCvalue=-100
	      else:
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
        ##########  until here goes the new monte_carlo_eval() function
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
	weights=[]
	#weights IS A PAIRLIST(position,weight)
	#CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
	for x in range(len(partialScores)):
	  weights.append((x, partialScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
	
	
	#MUTATION ATTEMPTS
	#INSIDE THIS LOOP, THE VARIABLE partialScore IS NOT MODIFIED. ONLY COMPARED TO THE MUTATED SCORES
	#THE LOOP ENDS WHEN A MUTATION IS ACCEPTED OR WHEN A LIMIT OF ATTEMPTS IS REACHED
	mutAttempts=0       #COUNT MUTATIONS ATTEMPTS
	while 10000 > mutAttempts:    ##JUST A SYMBOLIC MAX. AMOUNT OF MUTATIONS ATTEMPTS
	    mutAttempts+=1



         ####################  from here it can be extracted to a function mutate_sequence() that returnes the mutated sequence
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
	    mutatePosition= weighted_choice(weights) 
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

            ########### UNTIL HERE THE new function mutate_sequence()
              



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


                #############THIS COULD BE AN EXTRACTED FUNCTION called monte_carlo() that returns true or false wheter it accepts or not
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
  firstPartialEvaluation(sequence, config_params,positionScores , False )
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

if not global_evaluation:
	print "**END OF SEARCH**"
	if globalScore==0:
	  print "REACHED SCORE = 0"
	  print 'FINAL SEQUENCE: ' + sequence
	else:
	  print "REACHED LIMIT OF ITERATIONS"
	  data = [sequence,positionScores]
	  col_width = max(len(str(word)) for row in data for word in row)  # padding
	  print 'LAST SEQUENCE: ' + sequence
	  print 'FINAL SCORES'
	  for row in data:
	      print  "|".join(str(word).ljust(col_width) for word in row)
	  print "Global score: " + str(globalScore)

	  #print "Final sequence: " + sequence
	  #print "Final score:    " + ''.join(map(str, positionScores))
	  #print "Global score: " + str(globalScore)

	if minimalOutput:
	    logFileStream.write('END' + tab + str(globalScore) + tab + sequence + endl)  



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
	if testTimes:
            print_evaluation_time(totalElapsedTime,times_dict)



# if __name__ == '__main__':
    # main()
