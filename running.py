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

  #makeCommand = "make"
  runCommand = "python link1.2.py --length " + str(length)
  #print makeCommand
  #os.system(makeCommand)
  print runCommand
  os.system(runCommand)


#***********************************************************************




#***********************************************************************

def printHelp():
  print endl + "The usage mode is: "
  print "  python amber.py [options]"
  print endl + "   where options are:"
  print tab + "-P or -p :" + tab + "TopFile"
  print tab + "-C or -c :" + tab + "RstFile"
  print tab + "-I or -i :" + tab + "inputFile"





#***********************************************************************
#**************************** Main Program *****************************
#***********************************************************************

#top = input_amber_path + defaultTopFile
#rst = input_amber_path + defaultRstFile
#mdin = input_amber_path + defaultMdinFile



#for i in range(10,40,1):
  #for j in range(0,5):
    #runLinkeadoSeq(sequence[:i])







 #len(sys.argv) < 2:
  #print "\n Using all the default input files:"
  #print tab + top
  #print tab + rst
  #print tab + mdin
  #print "\n If you want a specific input file, see the help running:"
  #print tab + "python amber.py -h"
#else:






beta=False
score=False
mutAttempt=False
time=False




#CLEAR ALL OUTPUT FILES

#outputFile=outputPath + scoresFile
#outFile=open(outputFile, "w")
#outFile.close()
#outputFile=outputPath + mutAttemptsFile
#outFile=open(outputFile, "w")
#outFile.close()


length=30

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




if score or mutAttempt:
  runCommand = "python link1.2.py --beta 2.0 --length " + str(length)
  for x in range(0,10):  #hago 10 pruebas identicas
    os.system(runCommand)


if time:
  #para diferentes longitudes hago 10 pruebas identicas
  for largo in range(30,60,5):
    runCommand = "python link1.2.py --beta 2.0 --length " + str(largo)
    for x in range(0,10):  #hago 10 pruebas identicas
      os.system(runCommand)

#openLennardTable()

#parseTop(top)
#parseRst(rst)
#parseMdin(mdin)


#makeParticlesInputFile()
#makeLennardTable()


#runMacheAmber()
#runSander(top, rst, mdin)
  
  
  
  







