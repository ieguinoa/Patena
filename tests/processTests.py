import operator
import os
import errno
from math import *
import sys
from os.path import isfile
from itertools import *
import numpy as np
import string
import matplotlib.pyplot as plt
#import statistics
sys.path.insert(0, 'graphics')   #GRAPHICS FUNCTIONS
import glob
from graphs import *


#  	FROM BETA TEST, PLOT:
#		-SCATTER EACH TOTAL TIME IN TIME vs BETA PLOT (DIFFERENT COLOR FOR RAND/SEQ)
#		-SCATTER SCORE vs STEP# FOR A FEW DIFFERENT BETA VALUES(SHOULD I PUT A LIMIT IN STEP# AXIS?)
#		-SCATTER MUTAttempts vs 
#def processBetaTestResults2(basePath):
  #xSeqValues,ySeqValues=[],[]
  #xRandValues,yRandValues=[],[]
  #iterations,mutAttempts,stepScores, colorValues=[],[],[],[]
  #for files in os.listdir(basePath):
    #randomSeq=False
    #if not(files.endswith('.png')):
  ##CHANGE THIS LOOP FOR  
  ##for filename in glob.glob(os.path.join(basePath, '*.out')):
  ##AND CHANGE THE TESTOUTPUT IN BLEACH
      #with open(basePath+'/'+files,'r') as rawdata:
	#betaValue=0.0
	#time=0.0 #default value(in case the execution has not finished)
	#for lines in rawdata.readlines():
	  #cols=lines.split('\t')
	  ##print cols
	  #if cols[0]=='RAND' or cols[0]=='SEQ':
	    #if cols[0]=='RAND':
	      #randomSeq=True
	    #else:
	      #randomSeq=False
	    #betaValue=float(cols[1])
	  #if cols[0]=='END':
	    #time=float(cols[1])
	  #if cols[0]=='LOOP1' or cols[0]=='LOOP2':
	    #if betaValue==0.5:
	      #mutAttempts.append(int(cols[2]))
	      #stepScores.append(float(cols[3]))
	      #iterations.append(int(cols[1]))
	      ##colorValues.append(int(betaValue*35))
	      #colorValues.append('red')
	      ##the iteration time is stored in cols[4] 
	    ##if betaValue==1.0:
	      ##mutAttempts.append(int(cols[2]))
	      ##stepScores.append(float(cols[3]))
	      ##iterations.append(int(cols[1]))
	      ###colorValues.append(int(betaValue*35))
	      ##colorValues.append('blue')
	    #if betaValue==1.5:
	      #mutAttempts.append(int(cols[2]))
	      #stepScores.append(float(cols[3]))
	      #iterations.append(int(cols[1]))
	      ##colorValues.append(int(betaValue*35))
	      #colorValues.append('green')
	    #if cols[0]=='LOOP2':
		#print 'AT LEAST 1 REACHED SECOND LOOP'
	#if time>7000:
	  #time=7000
	##if betaValue==0.4 and time>3000:
	  ##print files
	#if randomSeq:
	  #xRandValues.append(betaValue)	  
	  #yRandValues.append(time) 
	#else:
	  #xSeqValues.append(betaValue)
	  #ySeqValues.append(time)
	##if betaValue==0.0:
	  ##print files
    #params = {
	  #'xlabel': u"Beta values",
	  #'ylabel': u"Total time [s]",
	  #'xRandValues': xRandValues,
	  #'xSeqValues': xSeqValues,
	  #'yRandValues':yRandValues,
	  #'ySeqValues':ySeqValues,
	  #'filename': basePath+'/beta-vs-time-length50.png',
	  #'title': 'Sequence length=50'
      #}	
    #scatterGraphBinaryColour(**params)
   



#  	FROM TIME TEST, PLOT:
#		-SCATTER EACH TOTAL TIME IN TIME vs LENGTH PLOT (DIFFERENT COLOR FOR RAND/SEQ)
#		-
#		-
def processTimeTestResults(basePath):
  xSeqValues,ySeqValues=[],[]   
  xRandValues,yRandValues=[],[]
  iterations,colorValues,mutAttempts,stepScores
  for files in os.listdir(basePath):
    randomSeq=False
    if not(files.endswith('.png')):
      with open(basePath+'/'+files,'r') as rawdata:
	time=0.0 #default value(in case the execution has not finished)
	length=0
	for lines in rawdata.readlines():
	  cols=lines.split('\t')
	  #print cols
	  if cols[0]=='RAND' or cols[0]=='SEQ':
	    if cols[0]=='RAND':
	      randomSeq=True
	    else:
	      randomSeq=False
	    betaValue=float(cols[1])  
	    length=float(cols[2])
	  if cols[0]=='END':
	    time=float(cols[1])
	  if cols[0]=='LOOP1' or cols[0]=='LOOP2':
	    colorValues.append(int(betaValue*35))
	    iterations.append(int(cols[1]))
	    mutAttempts.append(int(cols[2]))
	    stepScores.append(float(cols[3]))
	    #the iteration time is stored in cols[4] 
	    if cols[0]=='LOOP2':
	      print 'AT LEAST 1 REACHED SECOND LOOP'
	if randomSeq:
	  xRandValues.append(length)	  
	  yRandValues.append(time) 
	else:
	  xSeqValues.append(length)
	  ySeqValues.append(time)
    params = {
	  'xlabel': u"Sequence length",
	  'ylabel': u"Total time [s]",
	  'xRandValues': xRandValues,
	  'xSeqValues': xSeqValues,
	  'yRandValues':yRandValues,
	  'ySeqValues':ySeqValues,
	  'filename': basePath+'/seqLength-vs-time-beta1.png',
	  'title': 'Beta value = 1.0'
      }	
    #scatterGraphBinaryColour(**params)


###################################################
####    PROCESS BETA TESTS RESULTS AND PLOT: ######
####	BETA(x)   vs   EXEC.TIME(y)log?  ##########
####	MEAN TIMES, INCLUDES ERROR	###########
###################################################

def makeBetaVsTime(basePath):
  betaSeqDictTime,betaRandDictTime={},{}  #save pairs (beta-[list of times])
  iterations,mutAttempts,stepScores, colorValues=[],[],[],[]
  excludeList=[0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.6,1.7,1.8,1.9,2.1,2.2,2.4]
  for files in os.listdir(basePath):
    randomSeq=False
    if not(files.endswith('.png')):
      with open(basePath+'/'+files,'r') as rawdata:
	betaValue=0.0
	time=16000.0 #default value(in case the execution has not finished)
	
	#PROCESS FILE: SAVE BETA VALUE, RAND/NATURAL AND ELAPSED TIME
	for lines in rawdata.readlines():
	  cols=lines.split('\t')
	  #print cols
	  if cols[0]=='RAND' or cols[0]=='SEQ':
	    betaValue=float(cols[1])
	    if cols[0]=='RAND':
	      randomSeq=True
	    else:
	      randomSeq=False
	  if cols[0]=='END':
	    time=float(cols[1])
	    #if time>5000:
	      #time=5000
	
	#SAVE EXECUTION DATA IN DICTIONARY
	if betaValue not in excludeList:
	  if randomSeq:
	    if betaValue in betaRandDictTime:
	      betaRandDictTime[betaValue].append(time)
	    else:
	      betaRandDictTime[betaValue]=[time]
	  else:
	    if betaValue in betaSeqDictTime:
	      betaSeqDictTime[betaValue].append(time)
	    else:
	      betaSeqDictTime[betaValue]=[time]



  
  #BUILD X,Y LISTS OF EVERYTHING (VALUES, MEAN, STDEV )
  
  #RANDOM EXECUTIONS LISTS
  yMeanTimeValuesRand,yTimeValuesRand,yStdevTimeValuesRand=[],[],[]
  xMeanBetaValuesRand=[]  #mean and stdev x values are the same (length = different beta values)
  xBetaTimeValuesRand=[]  #one beta value for each execution (length = number of executions ) 
  
  #NATURAL SEQs EXECUTION LISTS	
  yMeanTimeValuesSeq,yTimeValuesSeq,yStdevTimeValuesSeq=[],[],[]
  xMeanBetaValuesSeq=[]  #mean and stdev x values are the same (one for each beta value)
  xBetaTimeValuesSeq=[]  #one beta value for each execution (length = number of executions ) 
  
  
  #ITERATE OVER THE DICTIONARIES AND FILL LISTS
  
  #RANDOM DICTTIONARIES
  keylist = betaRandDictTime.keys()
  keylist.sort()
  for key in keylist:
    xMeanBetaValuesRand.append(key)
    yMeanTimeValuesRand.append(np.mean(np.asarray(betaRandDictTime[key])))
    yStdevTimeValuesRand.append(np.std(np.asarray(betaRandDictTime[key])))
    for index in range(len(betaRandDictTime[key])):
      xBetaTimeValuesRand.append(key)	  				#append the key
      yTimeValuesRand.append(betaRandDictTime[key][index])		#append the value
 
  #NATURAL SEQs DICTIONARIES
  keylist = betaSeqDictTime.keys()
  keylist.sort()
  #for key in keylist:
  for key in keylist:
  #for key in betaSeqDictTime:
    xMeanBetaValuesSeq.append(key)
    yMeanTimeValuesSeq.append(np.mean(np.asarray(betaSeqDictTime[key])))
    yStdevTimeValuesSeq.append(np.std(np.asarray(betaSeqDictTime[key])))
    for index in range(len(betaSeqDictTime[key])):
      xBetaTimeValuesSeq.append(key)	  				#append the key
      yTimeValuesSeq.append(betaSeqDictTime[key][index])			#append the value

  params = {
	'xlabel': u"Beta=%Aceptacion score+1",
	'ylabel': u"Tiempo ejec. [s]",
	#'xRandValues': xRandValues,
	#'xSeqValues': xSeqValues,
	#'yRandValues':yRandValues,
	#'ySeqValues':ySeqValues,
	'xMeanValuesRand': xMeanBetaValuesRand,
	'yMeanValuesRand': yMeanTimeValuesRand,
	'xMeanValuesSeq': xMeanBetaValuesSeq,
	'yMeanValuesSeq': yMeanTimeValuesSeq,
	'yErrorValuesRand': yStdevTimeValuesRand,
	'yErrorValuesSeq': yStdevTimeValuesSeq,
	'filename': basePath+'/beta-vs-time-length50.png',
	#'ymax':5000,
	'title': 'Largo secuencia=50'
    }	
  meanErrorLines(**params)











###################################################
####  PROCESS BETA TESTS RESULTS AND PLOT:   ######
####       ITERATION NUMBER vs % ACEPTACION  ######
####   FOR EACH EXECUTION, MEAN  AND ERROR   ######
###################################################
def makeIterationVsAcceptRate(basePath):
  betaSeqDictTime,betaRandDictTime={},{}  #save pairs (beta-[list of times])
  
  betaList=[0.5,1.5,2.4]  #LIST OF BETAS TO PRINT
  #betaColourList=['red','green']
  betaSeqDictExecutions,betaRandDictExecutions={},{}  ##DICT OF EXECUTIONS BY BETA VALUE
  #betaValuesExecutionsRand=[]
  maxIterations=350
  step=10

  #iterationList=range(0,10,1)+range(10,100,20)+range(200,1000,50)
  betaValues=[]
  random=[]
  executions=[]  # this lists saves the executions I want to print
  iterations,mutAttempts,stepScores, colorValues=[],[],[],[]
  for files in os.listdir(basePath):
    randomSeq=False
    if not(files.endswith('.png')):
      with open(basePath+'/'+files,'r') as rawdata:
	execution=[]  #list of %accepted per step
	betaValue=0.0
	time=15000.0 #default value(in case the execution has not finished)
	for lines in rawdata.readlines():
	  cols=lines.split('\t')
	  #print cols
	  if cols[0]=='RAND' or cols[0]=='SEQ':
	    if cols[0]=='RAND':
	      randomSeq=True
	    else:
	      randomSeq=False
	    betaValue=float(cols[1])
	  if cols[0]=='END':
	    time=float(cols[1])
	  if cols[0]=='LOOP1' or cols[0]=='LOOP2':
	    if betaValue in betaList:
	      iteration=int(cols[1])
	      #if iteration in iterationList:
	      if ((iteration%step)==0):
	      #if int(iteration) in np.logspace(0,3):
		#print iteration
		mutAttempts.append(int(cols[2]))
		stepScores.append(float(cols[3]))
		iterations.append(int(cols[1]))
		#acceptRate=(1/int(cols[2]))
		#acceptRate=((1.0/int(cols[2]))*100.0)
		#execution.append(acceptRate)
		#SAVE MUTATTEMPTS
		execution.append(int(cols[2]))
		#SAVE SCORES
		#execution.append(float(cols[3]))
	      #colorValues.append(int(betaValue*35))
	      #colorValues.append('red')
	      #the iteration time is stored in cols[4] 
	    #if betaValue==1.0:
	      #mutAttempts.append(int(cols[2]))
	      #stepScores.append(float(cols[3]))
	      #iterations.append(int(cols[1]))
	      ##colorValues.append(int(betaValue*35))
	      #colorValues.append('blue')
	    #if betaValue==1.5:
	      #mutAttempts.append(int(cols[2]))
	      #stepScores.append(float(cols[3]))
	      #iterations.append(int(cols[1]))
	      ##colorValues.append(int(betaValue*35))
	      #colorValues.append('green')
	    #if cols[0]=='LOOP2':
		#print 'AT LEAST 1 REACHED SECOND LOOP'
	#if time>7000:
	  #time=7000
	#if betaValue==0.4 and time>3000:
	  #print files1
	
	#SAVE WHOLE EXECUTION
	if betaValue in betaList:
	  betaValues.append(betaValue)
	  executions.append(execution)
	  #print betaValue
	  random.append(randomSeq)
	  
	  
	    
	#if randomSeq:
	  ##xRandValues.append(betaValue)	  
	  ##yRandValues.append(time) 
	  #if betaValue in betaRandDictTime:
	    #betaRandDictTime[betaValue].append(time)
	    #betaRandDictExecutions[betaValue].append(execution)
	    ##cantBetaRand[xRandValues[x]] += 1
	  #else:
	    #betaRandDictTime[betaValue]=[time]
	    #betaRandDictExecutions[betaValue][execution]
	    ##cantBetaRand[xRandValues[x]] = 1	
	#else:
	  ##xSeqValues.append(betaValue)
	  ##ySeqValues.append(time)
	  #if betaValue in betaSeqDictTime:
	    #betaSeqDictTime[betaValue].append(time)
	    #betaSeqDictExecutions[betaValue].append(execution)
	    ##cantBetaRand[xRandValues[x]] += 1
	  #else:
	    #betaSeqDictTime[betaValue]=[time]
	    #betaSeqDictExecutions=[execution]
	##if betaValue==0.0:
	  ##print files

 
  params = { 
      'executionsList':executions,
      'random':random,
      'beta':betaValues,
      'logScale': True,
      'maxIterations': maxIterations,
      'step': step,
      'xlabel':'Iteration',
      'ylabel':'Intentos de mutac.',
	
   }
 
  iterationVsX(**params)	    
 
 
 
 
 
 
  #NOW PROCESS EXECUTIONS DATA TO GET MEAN AND STDEV 
  
  
  #FIRST MAKE DICT WITH LIST OF MUTATTEMPTS PER ITERATION #
  iterRandDict={}
  iterSeqDict={}
  for exeIndex in range(len(executions)):
    if random[exeIndex]:
      if betaValues[exeIndex] in iterRandDict:
	    iteration=0
	    while iteration<len(executions[exeIndex]) and iteration<maxIterations:
	      #AGREGO EL VALOR
	      iterRandDict[betaValues[exeIndex]][iteration].append(executions[exeIndex][iteration]) 
	      iteration+=1
	      
      else:   #IF IT IS THE FIRST EXECUTION FOR THIS BETA VALUE
	iterRandDict[betaValues[exeIndex]]= [[] for i in range(maxIterations)]
    
    else:
      if betaValues[exeIndex] in iterSeqDict:
	    iteration=0
	    while iteration<len(executions[exeIndex]) and iteration<maxIterations:
	      #AGREGO EL VALOR
	      iterSeqDict[betaValues[exeIndex]][iteration].append(executions[exeIndex][iteration])
	      iteration+=1
	      
      else:   #IF IT IS THE FIRST EXECUTION FOR THIS BETA VALUE
	#print 'entro 1 vez'
	iterSeqDict[betaValues[exeIndex]]= [[] for i in range(maxIterations)]
  
  
  #NOW PROCESS PREVIOUS LISTS TO GET MEANS AND STDEV
  meanRandDict,meanSeqDict,errorRandDict,errorSeqDict={},{},{},{}
  for betas in iterRandDict.keys():
    #if betas in betaList:
      meanRandDict[betas]=[]
      errorRandDict[betas]=[]
      iteration=0
      while iteration<len(iterRandDict[betas]) and iteration<maxIterations:
        #print iterRandDict[betas][iteration]
        if len(iterRandDict[betas][iteration]) > 0:
	  mean=np.mean(np.asarray(iterRandDict[betas][iteration]))
	  #print mean
	  std=np.std(np.asarray(iterRandDict[betas][iteration]))
	  #print std
	  meanRandDict[betas].append(mean)
	  errorRandDict[betas].append(std)
	iteration+=1
      
  for betas in iterSeqDict.keys():
    #if betas in betaList:
      meanSeqDict[betas]=[]
      errorSeqDict[betas]=[]
      iteration=0
      while iteration<len(iterSeqDict[betas]) and iteration<maxIterations:
	if len(iterSeqDict[betas][iteration]) > 0:
	  mean=np.mean(np.asarray(iterSeqDict[betas][iteration]))
	  #print mean
	  std=np.std(np.asarray(iterSeqDict[betas][iteration]))
	  #print std
	  meanSeqDict[betas].append(mean)
	  errorSeqDict[betas].append(std)
	iteration+=1
     
    
  betaValues, random, executionsMean, executionsError =[],[],[],[]
  for betas in meanRandDict.keys():
    betaValues.append(betas)
    random.append(True)
    executionsMean.append(meanRandDict[betas])
    executionsError.append(errorRandDict[betas])
    
  for betas in meanSeqDict.keys():
    betaValues.append(betas)
    random.append(False)
    executionsMean.append(meanSeqDict[betas])
    executionsError.append(errorSeqDict[betas])
  
  
  
  params = { 
      'executionsList':executionsMean,
      'executionsErrorList':executionsError,
      'random':random,
      'beta':betaValues,
      'logScale': True,
      'maxIterations': maxIterations,
      'step': step,
      'ylabel': 'Intentos de mutac.',
      'xlabel': 'Iteration'
    }
  
  iterationVsXError(**params)
  
  
  
   ##proceso los datos recolectados
    #betaMeanRand,cantBetaRand = {},{}
    #betaMeanNat,cantBetaNat = {},{}
    #for x in len(xRandValues)
      #if xRandValues[x] in :
	#betaMeanRand[xRandValues[x]].append(yRandValues[x])
	##cantBetaRand[xRandValues[x]] += 1
      #else:
	#betaMeanRand[xRandValues[x]]=[yRandValues[x]]
        ##cantBetaRand[xRandValues[x]] = 1	
    
    
    #for y in betaMeanRand:
      #betaMeanRand.append(mean(betaMeanRand[y]))
      #betaMeanStdv.append(pstdev(betaMeanRand[y]))
      ##betaMeanRand[y]=betaMeanRand[y]/cantBetaRand[y]
    
  #print 'paso todo'
  
  ##BUILD X,Y LISTS OF EVERYTHING (VALUES, MEAN, STDEV )
  #yMeanTimeValuesRand,yTimeValuesRand,yStdevTimeValuesRand=[],[],[]
  #xMeanBetaValuesRand,xBetaTimeValuesRand=[],[]  #mean and stdev x values are the same (one for each beta value)
  
  #yMeanTimeValuesSeq,yTimeValuesSeq,yStdevTimeValuesSeq=[],[],[]
  #xMeanBetaValuesSeq,xBetaTimeValuesSeq=[],[]  #mean and stdev x values are the same (one for each beta value)
  
  ##iterate over rand dict and get values together in a signle list
  
  #keylist = betaRandDictTime.keys()
  #keylist.sort()
  ##for key in keylist:
  #for key in keylist:
    #xMeanBetaValuesRand.append(key)
    #yMeanTimeValuesRand.append(np.mean(np.asarray(betaRandDictTime[key])))
    #yStdevTimeValuesRand.append(np.std(np.asarray(betaRandDictTime[key])))
    #for index in range(len(betaRandDictTime[key])):
      #xBetaTimeValuesRand.append(key)	  				#append the key
      #yTimeValuesRand.append(betaRandDictTime[key][index])		#append the value
  
  #keylist = betaSeqDictTime.keys()
  #keylist.sort()
  ##for key in keylist:
  #for key in keylist:
  ##for key in betaSeqDictTime:
    #xMeanBetaValuesSeq.append(key)
    #yMeanTimeValuesSeq.append(np.mean(np.asarray(betaSeqDictTime[key])))
    #yStdevTimeValuesSeq.append(np.std(np.asarray(betaSeqDictTime[key])))
    #for index in range(len(betaSeqDictTime[key])):
      #xBetaTimeValuesSeq.append(key)	  				#append the key
      #yTimeValuesSeq.append(betaSeqDictTime[key][index])			#append the value
  #print 'paso todo'

  #params = {
	#'xlabel': u"Beta value",
	#'ylabel': u"Total time [s]",
	##'xRandValues': xRandValues,
	##'xSeqValues': xSeqValues,
	##'yRandValues':yRandValues,
	##'ySeqValues':ySeqValues,
	#'xMeanValuesRand': xMeanBetaValuesRand,
	#'yMeanValuesRand': yMeanTimeValuesRand,
	#'xMeanValuesSeq': xMeanBetaValuesSeq,
	#'yMeanValuesSeq': yMeanTimeValuesSeq,
	#'yErrorValuesRand': yStdevTimeValuesRand,
	#'yErrorValuesSeq': yStdevTimeValuesSeq,
	#'filename': basePath+'/beta-vs-time-length50.png',
	##'ymax':5000,
	#'title': 'Sequence length=50'
    #}	
  #meanErrorLines(**params)
  #scatterGraphBinaryColour(**params)
  #params = {
	#'xlabel': u"Iteration number",
	#'ylabel': u"Score [patena Units]",
	#'xValues': iterations,
	#'yValues': stepScores,
	#'colorValues': colorValues ,
	#'filename': basePath+'/iterationNum-vs-score.png',
	#'title': ''
    #}	
  #scatterGraphNaryColour(**params)
  #params = {
	#'xlabel': u"Iteration number",
	#'ylabel': u"Mutation Attempts",
	#'xValues': iterations,
	#'yValues': mutAttempts,
	#'colorValues': colorValues ,
	#'filename': basePath+'/iterationNum-vs-mutAtt.png',
	#'title': ''
    #}	
  #scatterGraphNaryColour(**params)
  
  
 






####  OPEN ALL LOGS IN PATH AND ALIGN RESULT SEQUENCES 
####   SAVE RESULTS IN results FILE
def getDivergenceList(basePath):
   multiAlignment=open(basePath+'/resultsList','w')
   resultNumber=1
   for files in os.listdir(basePath):
    if (files.endswith('.log')):
       with open(basePath+'/'+files,'r') as rawdata:
	for line in rawdata:
	  pass
	last = line
	#print last
	resultSequence=last.split('\t')[1]
	#print resultSequence
	#multiAlignment.write('>RESULTSEQ_'+str(resultNumber) + '\n')
        multiAlignment.write(resultSequence + '\n')
        resultNumber+=1



def getIdentityPercent(seq1,seq2):
	#print seq1
	#print seq2
	hits=0
	for aa in seq1:
		if aa == seq2[seq1.index(aa)]:
			hits+=1
	#print (float(hits)/(len(seq1)))*100
	return (float(hits)/float((len(seq1))))*100


def processDivergenceTest(basePath):
	getDivergenceList(basePath)
	sequenceList=[]
	initialSeq='MALWMRLLPLLALLALWGPDPAAAFVNQHL'
	#FIRST, LOAD SEQUENCE LIST FROM FILE
	#sequenceList = (with open(basePath+'/'+'resultsList','r')).read().splitlines()
	with open(basePath+'/'+'resultsList','r') as rawdata:
		for line in rawdata.readlines():
			sequenceList.append(line.rstrip())
	
	idPercentInitial, idPercentAll=[],[]
	#FIRST, GET IDENTITY % AGAINST INITIAL SEQUENCE:
	for seq in sequenceList:
		idPercentInitial.append(getIdentityPercent(initialSeq,seq))
	
	
	#print len(idPercentInitial)
	plt.hist(np.asarray(idPercentInitial),10)
	plt.title('Identity against starting sequence')	
	plt.ylabel('Count')
	plt.xlabel('Sequence identity %')	
	#plt.ylim(0,15)
	#plt.savefig('againstStart.png', bbox_inches='tight', frameon=True)
	#plt.show()
	#GET IDENTITY % ALL AGAINST ALL
	for seq1 in sequenceList:
		for seq2 in sequenceList:
			if sequenceList.index(seq1) != sequenceList.index(seq2):
				idPercentAll.append(getIdentityPercent(seq1,seq2))
	plt.hist(np.asarray(idPercentAll),18)
	plt.title('Identity against equivalent results')
	plt.ylabel('Count')
	plt.xlabel('Sequence identity %')
	plt.savefig('againstAll.png', bbox_inches='tight', frameon=True)	
	#plt.show()





#print sys.argv[1]
#makeIterationVsAcceptRate(sys.argv[1])
#makeIterationVsAcceptRate('/home/ieguinoa/results/0.5-2.0')        
#makeBetaVsTime(sys.argv[1])        
#getDivergenceList(sys.argv[1])
processDivergenceTest(sys.argv[1])
#processBetaTestResults('/home/ieguinoa/results/todos')    
#processBetaTestResults('/home/ieguinoa/results/beta-0.1-1.9')
#processBetaTestResults('/home/ieguinoa/results/beta-0.5-1.5-step-0.1')
#processTimeTestResults('/home/ieguinoa/results/timeTestBeta1')
#processBetaTestResults('/home/ieguinoa/results/betaTests')
