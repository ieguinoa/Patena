import operator
import os
import errno
from math import *
import sys
from os.path import isfile
from itertools import *
import numpy as np
import string
sys.path.insert(0, 'graphics')   #GRAPHICS FUNCTIONS
import glob


def processBetaTestResults(basePath):
  xRandValues,yRandValues=[],[]
  xSeqValues,ySeqValues=[],[]
  for files in os.listdir(basePath):
  #CHANGE THIS LOOP FOR  
  #for filename in glob.glob(os.path.join(basePath, '*.out')):
  #AND CHANGE THE TESTOUTPUT IN BLEACH
    with open(basePath+'/'+files,'r') as rawdata:
      for lines in rawdata.readlines():
	#TODO MODIFICAR BLEACH PARA QUE IMPRIMA EN LA PRIMER LINEA SI LA SEQ INICIAL ES RANDOM O NO
	cols=lines.split('\t')
	if cols[0]=='RAND' or cols[0]=='SEQ':
	  if cols[0]=='RAND':
	    randomSeq=True
	  else:
	    randomSeq=False
	  betaValue=cols[1]
	#print cols[0]
	if cols[0]=='END':
	  time=cols[1]
      if randomSeq:
	xRandValues.append(betaValue)	  
	yRandValues.append(time) #APPEND 1 IF INITIAL SEQUENCE FOR THIS TEST WAS RANDOM(WILL PRINT IN DIFFERENT COLOR) 
      else:
	ySeqValues.append(time)
	xSeqValues.append(betaValue)
  scatterGraphBinaryColour(xlabel, ylabel, xRandvalues,xSeqValues, yRandValues,ySeqValues,filename,
                          ylegend=u'',fitlegend='u', ticks='', title=u""):
	
processBetaTestResults('/home/ieguinoa/results/betaTests')