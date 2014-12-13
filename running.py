import operator
import os
from math import *
import sys
from os.path import isfile
from itertools import *
import numpy as np
import string

#***********************************************************************
#************************** Global Variables ***************************
#***********************************************************************

endl = "\n"
tab = "\t"

#************************
#Folders and files names:
#************************






#***********************************************************************

def runLinkeado(length):
  print endl + "****************************************************"
  print "Running Linkeado (v1.2)"
  print "****************************************************"

  runCommand = "python link1.2.py --length " + str(length)
  print runCommand
  os.system(runCommand)


#***********************************************************************




#***********************************************************************

def printHelp():
  print endl + "The usage mode is: "


#***********************************************************************



#*********************
#***DEFAULT VALUES***
#*******************

beta=False
score=False
mutAttempt=False
time=False
length=30




#********************************
#***** CHECK INPUT PARAMETERS **
#******************************

for i in range(1,len(sys.argv)):
  arg = sys.argv[i]
  if (arg=='--beta' or arg== '-b'):
    #TEST BETA vs ITERATIONS-TIME
    beta=True
  elif (arg=='--score' or arg== '-s'):
    #TEST score vs iteration number
    score=True
  elif (arg=='--mutations' or arg== '-m'):
    #TEST mutation attempts  Vs. iteration number
    mutAttempt=True
  elif (arg=='--time' or arg== '-t'):
    #TEST total time vs sequence length
    time=True
  elif (arg=='--length' or arg== '-l'):
    #DEFINE SEQUENCE LENGTH
    length= sys.argv[i+1] 

sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"





#TEST AVARAGE SCORE-MUT.ATTEMPT  FOR EACH STEP IN DIFFERENT RUNS WITH SAME LENGTH
#THIS TEST RESULTS IN A FILE WITH SEVERAL  SCORE-STEP  or  MutATTEMPT-STEP  
#TO GET THE AVARAGE YOU NEED TO PROCESS THAT FILE AND AVARAGE VALUES FOR THE SAME STEP NUMBER
#THIS RUN ALSO GIVES A FILE WITH  TIMES FOR EACH STEP WICH CAN BE PROCESSED IN THE SAME WAY
if score or mutAttempt:
  runCommand = "python link1.2.py --beta 0.5 --length " + str(length)
  for x in range(0,10):  #hago 10 pruebas identicas
    	print runCommand
	os.system(runCommand)


#TEST AVARAGE RUN TIME (TOTAL)
#THIS TEST GIVES A FILE WITH TOTAL TIME vs LENGTH
#TO GET THE AVARAGE FOR EACH LENGHT YOU NEED TO PROCESS AND AVARAGE ALL TIMES ASSOCIATED WITH THE SAME LENGTH
if time:
  #para diferentes longitudes hago 10 pruebas identicas
  for largo in range(30,60,5):
    runCommand = "python link1.2.py --beta 2.0 --length " + str(largo)
    for x in range(0,10):  #hago 10 pruebas identicas
      os.system(runCommand)


#TEST TOTAL NUMBER OF STEPS vs BETA VALUE FOR THE SAME SEQUENCE LENGTH 
#THIS TEST RESULTS IN A FILE WITH BETA - TOTAL STEPS VALUES. THESE SHOULD BE PROCESSED TO GET AVARAGE FOR 
#IT ALSO PRINTS THE TIME, SCORE AND MUTATION ATTEMPS FOR EACH STEP WICH CAN BE PROCESSED IN A SIMILAR WAY AS BEFORE
if beta:
	for value in range (4):
		betaValue= 1.0 + (0.5*value)
		runCommand = "python link1.2.py --beta " + str(betaValue) + " --length " + str(length)
		for x in range (0,10):
			print runCommand
			os.system(runCommand)
