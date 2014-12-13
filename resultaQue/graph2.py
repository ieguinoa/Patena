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

inputPath = "Output/"


##el mode puede ser:
#        --time = grafico Tiempo de corrida VS corrida (puntos)  a partir de los archivos totalTimesXXXX
#        --iterations = grafico iteraciones Vs corrida (puntos)   a partir de los archivos totalTimesXXXX
#	 --scores = grafico de promedios de scores vs numero de iteracion a partir de los archivos scoreXXX
# 	 --mutAttempts = grafico de intentos de mutaciones vs numero de iteracion a partir de los archivos mutAttemptsXXX


time=False
iterations=False
scores=False
mutAttempts=False


##procesar los parametros
for x in range(len(sys.argv)):
  mode=sys.argv[x]
  if mode == '--time':
    time=True
  if mode == '--iterations':
    iterations=True
  if mode == '--mutAttempts':
    mutAttempts=True
  if mode == '--scores':  
    scores=True




if time:
  fig = plt.figure()
  for files in os.listdir(inputPath): 
    if files.startswith("totalTimes"):
      cant=[]   #muestras
      muestras=0
      inputFile = open(inputPath + files)
      for l in inputFile:
	muestras+=1
	cant.append(muestras)
	valores.append(float(l.split()[2]))
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
  fig.savefig(inputPath + "graficoTimes.png", dpi = 100)



if iterations:
  fig = plt.figure()
  for files in os.listdir(inputPath): 
    if files.startswith("totalTimes"):
      cant=[]   #muestras
      muestras=0
      inputFile = open(inputPath + files)
      for l in inputFile:
	muestras+=1
	cant.append(muestras)
	valores.append(float(l.split()[1]))
      inputFile.close()
      np_time = np.array(cant)
      np_step = np.array(valores)
      p1 = plt.plot(np_time, np_step, lw=2, label=files)
  legend=plt.legend(bbox_to_anchor=(0.8, 0.8), loc=4, borderaxespad=0.,
		    fancybox=True,shadow=True,title='Method')
  plt.grid()
  plt.show()
  fig.savefig(inputPath + "graficoIterations.png", dpi = 100)




    


if scores:
  fig = plt.figure()
  for files in os.listdir(inputPath): 
    if files.startswith("scores"): 
	cant=[]   #cantidad de valores
	valores=[]    #acumula los valores 
	#lleno las listas con 0
	for x in range(0,3000):
	  valores.append(0.)
	  cant.append(0)
	#
	inputFile = open(inputPath + files)
	for l in inputFile:
	  iteracion=int(l.split()[0])
	  valorAsociado=float(l.split()[1])
	  if valorAsociado > 0:   #NO TENGO EN CUENTA EL ULTIMO  SCORE
	    cant[iteracion]+=1
	    valores[iteracion]+=valorAsociado
	#inputFile.close()

	iterations = []
	value=[]  ##ACA GUARDO LOS VALORES QUE GRAFICO Vs ITERACION (PUEDEN SER SCORE, INTENTOS DE MUTACIONES, ETC)

	for x in range(0,100):
	  if valores[x] > 0:
	      value.append(valores[x] / cant[x])
	      iterations.append(x)

	np_time = np.array(iterations)
	np_step = np.array(value)
	p1 = plt.plot(np_time, np_step, lw=2, label=files)
	inputFile.close()



  legend=plt.legend(bbox_to_anchor=(0.38, 0.52), loc=4, borderaxespad=0.,
		      fancybox=True,shadow=True,title='Method')
  plt.grid()
  plt.show()
  fig.savefig(inputPath + "graficoScores.png", dpi = 100)



if mutAttempts:
  fig = plt.figure()
  for files in os.listdir(inputPath): 
    if files.startswith("mutationsAttempt"): 
	cant=[]   #cantidad de valores
	valores=[]    #acumula los valores 
	#lleno las listas con 0
	for x in range(0,3000):
	  valores.append(0.)
	  cant.append(0)
	#
	inputFile = open(inputPath + files)
	for l in inputFile:
	  iteracion=int(l.split()[0])
	  valorAsociado=float(l.split()[1])
	  if valorAsociado > 0:   #NO TENGO EN CUENTA EL ULTIMO  SCORE
	    cant[iteracion]+=1
	    valores[iteracion]+=valorAsociado
	#inputFile.close()

	iterations = []
	value=[]  ##ACA GUARDO LOS VALORES QUE GRAFICO Vs ITERACION (PUEDEN SER SCORE, INTENTOS DE MUTACIONES, ETC)

	for x in range(0,100):
	  if valores[x] > 0:
	      value.append(valores[x] / cant[x])
	      iterations.append(x)

	np_time = np.array(iterations)
	np_step = np.array(value)
	p1 = plt.plot(np_time, np_step, lw=2, label=files)
	inputFile.close()



  legend=plt.legend(bbox_to_anchor=(0.38, 0.52), loc=4, borderaxespad=0.,
		      fancybox=True,shadow=True,title='Method')
  plt.grid()
  plt.show()
  fig.savefig(inputPath + "graficoMutAttempts.png", dpi = 100)
