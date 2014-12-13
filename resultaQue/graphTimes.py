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

#file=sys.argv[1]  #REcibo el nombre del archivo por parametro 

outputPath = "Output/"

if sys.argv[1] == '--time':
  time=True
else:
  time=False
  

#lleno las listas con 0
#for x in range(0,3000):
  #valores.append(0.)
  #cant.append(0)
  
#time=False
  
  
fig = plt.figure()
#
for files in os.listdir(outputPath):
      valores=[]    #tiempo total / iteraciones totales 
      cant=[]   #muestras
      muestras=0
      inputFile = open(outputPath + files)
      for l in inputFile:
	muestras+=1
	cant.append(muestras)
	if time:
	  valores.append(float(l.split()[2]))
	else:
	    valores.append(int(l.split()[1]))
	#valorAsociado=float(l.split()[1])
	#if valorAsociado > 0:   #NO TENGO EN CUENTA EL ULTIMO  SCORE
	  #cant[iteracion]+=1
	  #valores[iteracion]+=valorAsociado
      inputFile.close()


      #iterations = []
      #value=[]  ##ACA GUARDO LOS VALORES QUE GRAFICO Vs ITERACION (PUEDEN SER SCORE, INTENTOS DE MUTACIONES, ETC)

      #for x in range(muestras):
	#if valores[x] > 0:
	    #value.append(valores[x] / cant[x])
	    #iterations.append(x)
      np_time = np.array(cant)
      np_step = np.array(valores)
      p1 = plt.plot(np_time, np_step, lw=2, label=files)




legend=plt.legend(bbox_to_anchor=(0.8, 0.8), loc=4, borderaxespad=0.,
		    fancybox=True,shadow=True,title='Method')
plt.grid()
plt.show()
fig.savefig(outputPath + "out.png", dpi = 100)
