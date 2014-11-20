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
defaultTopFile = "input.prmtop"
defaultRstFile = "input.rst7"
defaultMdinFile = "input.mdin"

input_amber_path = "Input_Amber/"
output_path = "Input_Mache/"
output_sander_path = "Output_Sander/"

lennardTableFilename = "TablaLennardOriginal"
lennardOutFilename = "TablaCoeficientesLennard"
particlesFilename = "particles.in"
outSanderFilename = "outSander"



inputPath="Input/"
inputFile="input"


#************************
#Top and Rst variables:
#************************
TYPE = {}

TIME_i = 0

NATOM = 0
NTYPES = 0

ATOM_NAME = []
CHARGE = []
ATOMIC_NUMBER = []
MASS = []

ACOEF = []
BCOEF = []

sigma = []
epsilon = []

POSITION = []
VELOCITY = []

#************************
#Default Input Values:
#************************
imin = 0
ntb = 0
ntp=0
ntt = 1
nstlim =  20
dt = 0.001
vlimit = 1
ntf = 1
ibelly = 0
ntx = 5
irest = 1
ntpr = 1
ntwx = 1
ntwe = 1
temp0 =100.0
tempi = 0
tautp= 2.0
cut=12.0






#***********************************************************************
#****************************** Functions ******************************
#***********************************************************************
def openLennardTable():
  global lennardTableFilename, TYPE

  print endl + "****************************************************"
  print "Parsing Lennard Table"
  print "****************************************************"

  if not isfile(output_path + lennardTableFilename):
    print "The file ", output_path + lennardTableFilename, " doesn't exist."
    exit()
  
  lennard = open(output_path + lennardTableFilename, 'r')
  for l in lennard:
    l = l.strip()
    if not l.split() or l.startswith('<end>'):
      print l
      break
    if l.startswith('#') or l.startswith('<start>'):
      print l
      continue
    f = l.split()
    TYPE[f[0]] = {'sigma': float(f[1]), 'epsilon': float(f[2]),\
     'charge': float(f[3])}
  
  lennard.close()

#***********************************************************************


def parseTop(_top):
  global NATOM, NTYPES, ATOM_NAME, CHARGE, \
  ATOMIC_NUMBER, MASS, ACOEF, BCOEF, particlesFilename

  print endl + "****************************************************"
  print "Parsing the top file"
  print "****************************************************"
  
  if not isfile(_top):
    print "The file ", _top, " doesn't exist."
    exit()
  
  topFile = open(_top, 'r')
  
  topFile_iter = iter(topFile)
  for t in topFile_iter:
    if t.startswith("%FLAG POINTERS"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(10I8)"):
	l = (topFile_iter.next()).split()
	NATOM = int(l[0])
	NTYPES = int(l[1])
	#Here we can get all necesary parameters
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "NATOM: ", NATOM
      print "NTYPES: ", NTYPES
    
    if t.startswith("%FLAG ATOM_NAME"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(20a4)"):
	l = (topFile_iter.next()).split()
	for name in l:
	  ATOM_NAME.append(name)
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "ATOM_NAME: ", ATOM_NAME
    
    if t.startswith("%FLAG CHARGE"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(5E16.8)"):
	l = (topFile_iter.next()).split()
	for ch in l:
	  CHARGE.append(float(ch))
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "CHARGE: ", CHARGE
    
    
    if t.startswith("%FLAG ATOMIC_NUMBER"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(10I8)"):
	l = (topFile_iter.next()).split()
	for atn in l:
	  ATOMIC_NUMBER.append(int(atn))
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "ATOMIC NUMBER: ", ATOMIC_NUMBER
    
    
    if t.startswith("%FLAG MASS"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(5E16.8)"):
	l = (topFile_iter.next()).split()
	for m in l:
	  MASS.append(float(m))
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "MASS: ", MASS
    
    
    
    
    # %FLAG LENNARD_JONES_ACOEF
    # %FORMAT(5E16.8)  (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
    # CN1  : Lennard Jones r**12 terms for all possible atom type interactions,
    #       indexed by ICO and IAC; for atom i and j where i < j, the index into
    #       this array is as follows (assuming the value of ICO(index) is positive):
    #       CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
    
    if t.startswith("%FLAG LENNARD_JONES_ACOEF"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(5E16.8)"):
	l = (topFile_iter.next()).split()
	for acoef in l:
	  ACOEF.append(float(acoef))
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "ACOEF: ", ACOEF
    
    
    if t.startswith("%FLAG LENNARD_JONES_BCOEF"):
      f = topFile_iter.next()
      if f.startswith("%FORMAT(5E16.8)"):
	l = (topFile_iter.next()).split()
	for bcoef in l:
	  BCOEF.append(float(bcoef))
      else:
	print "error: The file", _top, " isnt a prmtop file"
      print "BCOEF: ", BCOEF
    
  
  #for debug
  for i in range(0, NTYPES*(NTYPES+1)/2):
    A = ACOEF[i]
    B = BCOEF[i]
    sigma.append((A/B)**(1.0/6))
    epsilon.append((B**2)/(4*A))
  print "sigma = " + str(sigma)
  print "epsilon = " + str(epsilon)
  
  topFile.close()

#***********************************************************************


def parseRst(_rst):
  global NATOM, NTYPES, ATOM_NAME, CHARGE, POSITION, \
  ATOMIC_NUMBER, MASS, ACOEF, BCOEF, particlesFilename, \
  VELOCITY, TIME_i

  print endl + "****************************************************"
  print "Parsing the rst file"
  print "****************************************************"
  
  if not isfile(_rst):
    print "The file ", _rst, " doesn't exist."
    exit()
  
  rstFile = open(_rst, 'r')
  
  rstFile_iter = iter(rstFile)
  rstFile_iter.next()
  l = (rstFile_iter.next()).split()
  if int(l[0]) != NATOM:
    print "ERROR: the file ", _rst, " doesn\'t match with the top file"
    exit()
  TIME_i = float(l[1])
  
  #POSITIONS
  for i in range(0, NATOM/2):
    f = (rstFile_iter.next()).split()
    POSITION.append({'x': float(f[0]), 'y': float(f[1]), 'z': float(f[2])})
    POSITION.append({'x': float(f[3]), 'y': float(f[4]), 'z': float(f[5])})
  if NATOM/2 < NATOM/2.0:
    f = (rstFile_iter.next()).split()
    POSITION.append({'x': float(f[0]), 'y': float(f[1]), 'z': float(f[2])})
  
  #VELOCITIES
  if rstFile.tell == os.fstat(rstFile.fileno()).st_size:
    exit()
  for i in range(0, NATOM/2):
    f = (rstFile_iter.next()).split()
    VELOCITY.append({'x': float(f[0]), 'y': float(f[1]), 'z': float(f[2])})
    VELOCITY.append({'x': float(f[3]), 'y': float(f[4]), 'z': float(f[5])})
      
  if NATOM/2 < NATOM/2.0:
    f = f.split()
    VELOCITY.append({'x': float(f[0]), 'y': float(f[1]), 'z': float(f[2])})
  
  #BOX
  if rstFile.tell != os.fstat(rstFile.fileno()).st_size:
    print "There isn\'t Box"
    #TODO
  
  print POSITION
  print VELOCITY
  rstFile.close()

#***********************************************************************
def parseMdin(_mdin):
  global NATOM, NTYPES, ATOM_NAME, CHARGE, POSITION, \
  ATOMIC_NUMBER, MASS, ACOEF, BCOEF, particlesFilename, \
  VELOCITY, TIME_i

  print endl + "****************************************************"
  print "Parsing the mdin file"
  print "****************************************************"
  
  if not isfile(_mdin):
    print "The file ", _mdin, " doesn't exist."
    exit()
  
  mdinFile = open(_mdin, 'r')
  
  for l in mdinFile:
    if '&end' in l:
      break;
    noSpaces = l.replace(' ', '')
    parms = noSpaces.split(',')
    for parm in parms:
      setting = parm.split('=')
      for i in range(0,len(setting)):
	if setting[i] == 'temp0':
	  temp0 = setting[i+1]
	if setting[i] == 'tempi':
	  tempi = setting[i+1]
	if setting[i] == 'dt':
	  dt = setting[i+1]
	if setting[i] == 'nstlim':
	  nstlim = setting[i+1]
	
	
  print "dt = ", dt
  print "temp0 = ", temp0
  print "tautp = ", tautp
  print "nstlim = ", nstlim
  mdinFile.close()
    
    
#***********************************************************************

def makeParticlesInputFile():
  global NATOM, NTYPES, ATOM_NAME, CHARGE, \
  ATOMIC_NUMBER, MASS, ACOEF, BCOEF, particlesFilename

  print endl + "****************************************************"
  print "Making the inputFile"
  print "****************************************************"

  particlesFile = open(output_path + particlesFilename, 'w')
  
  #type, position, velocitie, charge
  particlesFile.write(str(NATOM)+endl)
  for i in range(0, NATOM):
    particlesFile.write(ATOM_NAME[i] + tab + 
    str(POSITION[i]['x']) + tab + str(POSITION[i]['y']) + tab + str(POSITION[i]['z']) + tab +
    str(VELOCITY[i]['x']) + tab + str(VELOCITY[i]['y']) + tab + str(VELOCITY[i]['z']) + tab +
    str(CHARGE[i]) + endl)
    
  #box
  
  #parameters for run
  particlesFile.write(str(nstlim) + endl +
		      str(dt) + endl +
		      str(temp0) + endl +
		      str(tempi) + endl +
		      str(tautp) + endl )
  
  particlesFile.write( endl + "#Format:" + endl +
			"NATOM" + endl +
			"TYPE" + tab + "POS(x)" + tab + "POS(y)" + tab + "POS(z)" +
			tab + "VEL(x)" + tab + "VEL(y)" + tab + "VEL(z)"+
			tab + "CHARGE" + endl +
			"dt" + endl +
			"temp0" + endl +
			"tempi" + endl +
			"tautp" + endl)
  particlesFile.close()



#***********************************************************************

def makeLennardTable():
  global NATOM, NTYPES, ATOM_NAME, CHARGE, TYPE, \
  ATOMIC_NUMBER, MASS, ACOEF, BCOEF, lennardOutFilename, \
  lennardInputFilename

  print endl + "****************************************************"
  print "Making the lennardTable"
  print "****************************************************"

  lennard = open(output_path + lennardOutFilename, 'w')
  
  used_types = set(ATOM_NAME)
  
  lennard.write(str(NTYPES) + endl)
  for t in used_types:
    lennard.write(t + tab + 
		  str(TYPE[t]['sigma']) + tab + 
		  str(TYPE[t]['epsilon']) + tab + 
		  str(TYPE[t]['charge']) + endl)
  lennard.close()

#***********************************************************************

def runMacheAmber():
  print endl + "****************************************************"
  print "Running MacheAmber"
  print "****************************************************"

  makeCommand = "make"
  runCommand = "./amberMache"
  print makeCommand
  os.system(makeCommand)
  print runCommand
  os.system(runCommand)


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

def runLinkeadoSeq(sequence):
  print endl + "****************************************************"
  print "Running Linkeado (v1.2)"
  print "****************************************************"

  #makeCommand = "make"
  runCommand = "python link1.2.py --seq " + sequence
  #print makeCommand
  #os.system(makeCommand)
  print runCommand
  os.system(runCommand)


#***********************************************************************

def runSander(_top, _rst, _mdin):
  print endl + "****************************************************"
  print "Running Sander"
  print "****************************************************"

  runCommand = "$P_PATH/sander -O"
  runCommand += " -i " + _mdin
  runCommand += " -o " + outSanderFilename
  runCommand += " -p " + _top 
  runCommand += " -c " + _rst
  
  print runCommand
  os.system(runCommand)
  
  moveCommand = "mv "
  for file in os.listdir("."):
    if file.endswith("mdcrd") or file.endswith("mdinfo") or file.endswith("mden") or file.endswith("restrt") or file.endswith(outSanderFilename):
      os.system(moveCommand + file + " " + output_sander_path)


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
    

sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"




if 




#openLennardTable()

#parseTop(top)
#parseRst(rst)
#parseMdin(mdin)


#makeParticlesInputFile()
#makeLennardTable()


#runMacheAmber()
#runSander(top, rst, mdin)
  
  
  
  







