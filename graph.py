#!/usr/bin/python

import operator
import os
from math import *
import sys
from os.path import isfile
from itertools import *
import numpy as np
import matplotlib.mlab as mlab
from matplotlib.pyplot import *
import matplotlib.pyplot as plt

import pylab as plot
params = {'legend.fontsize': 20,
          'legend.linewidth': 4}
plot.rcParams.update(params)

file=sys.argv[1]  #REcibo el nombre del archivo por parametro 

machePath = "Output/"
#sanderPath = "Output_Sander/"

#dirs = [x for x in os.listdir(machePath) if os.path.isdir(sanderPath+x) and os.path.isdir(sanderPath+x)]
#dirs=[machePath]
#print dirs
fig = plt.figure()
#length = []  ##file name = length of sequence      
#iterations = []
#for f in sorted(os.listdir(machePath)):
  #if f.endswith(".out") and f.startswith("time"):
    #print str(int(f))
    #inputFile = open(machePath + f)
    #matches=0.0  #number of matches (more than 1 iterations)
    #for l in inputFile:
      ##if int(l.split()[0]) > 1:
      #matches=matches+ int(l.split()[0])
      ##length.append(int(f))
      ##iterations.append(int(l.split()[0]))  ##append iterations in a list
    #length.append(int(f))   ## X index = length of sequence
    #iterations.append(matches/5)
    #inputFile.close()

value=[]  ##ACA GUARDO LOS VALORES QUE GRAFICO Vs ITERACION (PUEDEN SER SCORE, INTENTOS DE MUTACIONES, ETC)
iterations = []
inputFile = open(machePath + file)
#matches=0.0  #number of matches (more than 1 iterations)
for l in inputFile:
  iterations.append(int(l.split()[0]))  ##append iterations in a list
  value.append(float(l.split()[1]))  ##agrego el valor a graficar vs iteraciones
#length.append(int(f))   ## X index = length of sequence
#iterations.append(matches/5)
inputFile.close()


#print length
#print iterations
np_time = np.array(iterations)
np_step = np.array(value)
p1 = plt.plot(np_time, np_step, lw=2, label="hemoglobin")

#for f in os.listdir(sanderPath+d):
  #if f.endswith(".out") and f.startswith("time"):
    #print f
    #inputFile = open(sanderPath + d + "/" + f)
    #steps = []      
    #time = []
    
    #for l in inputFile:
      #stp = int(l.split()[0])
      #steps.append(stp)
      #time.append(float(l.split()[1]))		#Time is in S and is an average
    #inputFile.close()
    
    #np_time = np.array(time)
    #np_step = np.array(steps)
    #p1 = plt.plot(np_step, np_time, lw=2, label="snd_"+f)

legend=plt.legend(bbox_to_anchor=(0.38, 0.52), loc=4, borderaxespad=0.,
		    fancybox=True,shadow=True,title='Method')
plt.grid()
plt.show()
fig.savefig(machePath + "out.png", dpi = 100)