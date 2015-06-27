#!/usr/bin/env python
# -*- coding: utf-8 -*-
##
## Ejemplo de uso
##
##from graphs import *
##
##params =  {'title': u"Speedup de uso de memorias de textura",
##    'xlabel':u"Arquitectura GPU",
##    'ylabel':u"Aceleración (en veces)",
##    'values':(1.0,1.21,1.38),
##    'ticks':(u'Referencia', u'Fermi', u'Kepler'),
##    'filename':"texture.png"}
##barGraph(**params)

import matplotlib.lines as mlines
#import matplotlib.cm as cm

hasScipy = False
try:
  from scipy.interpolate import interp1d
  hasScipy = True
except ImportError:
  print 'quilombo'
  pass

try:
  import seaborn as sns
  palette = sns.color_palette()
  sns.set(palette=palette)
  sns.set_context("paper",font_scale=1.7)
  sns.set_style("ticks", {"font.family": "serif"})
except ImportError:
  palette = ["b"] * 1000
  print 'quilombo 2'
  pass

import matplotlib as mpl
import matplotlib.pyplot as plt

from matplotlib.patches import Rectangle

import pylab
import numpy as np

class Palette(object):
  def __init__(self, colors):
    self.colors = colors
    self.index = 0

  def next_color(self):
    tmp = self.index
    self.index = self.index + 1
    if self.index == len(self.colors):
      self.index = 0;
    return self.colors[tmp]

def repeatPaletteColors(sample_size, total_size):
    pieces = total_size / sample_size
    result = []
    for sample in xrange(0,sample_size):
        for piece in xrange(0, pieces):
            result.append(palette[sample])
    return result

def barGraph(ylabel, yvalues, ticks, filename, xlabel=None, ylim=None, title=u"", rotation=0, colors=palette):
  assert(len(yvalues) == len(ticks))
  N = len(yvalues)
  barLocations = np.arange(N)    # the x locations for the groups
  barWidth = 0.5       # the width of the bars: can also be len(x) sequence
  pylab.title(title)
  pylab.ylabel(ylabel)
  if xlabel:
    pylab.xlabel(xlabel)
  plt.xticks(np.arange(0.25, N), ticks, rotation=rotation)
  p1 = plt.bar(barLocations, yvalues, barWidth, color=colors)
  if ylim:
    plt.ylim(ylim)
  pylab.legend()
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def lineGraph(yvalues, filename,
               xlabel=u"", ylabel=u"", scale=u'linear',xvalues=None,
               ylegend=u'',ticks=None, title=u"", ylim=None, xlim=None):
  if xvalues is None:
    xvalues = range(len(yvalues))
  fig, ax = plt.subplots()

  if type(yvalues) != type([]):
    yvalues = [yvalues]
    ylegend = [ylegend]

  for _yvalues,_ylegend in zip(yvalues,ylegend):
    if hasScipy:
      f = interp1d(xvalues, _yvalues)
      ax.plot(xvalues,_yvalues,'o')
      ax.plot(xvalues,f(xvalues),'-', label=_ylegend)
    else:
      ax.plot(xvalues,_yvalues,'o', label=_ylegend)

  if ticks:
    locs, labels = plt.xticks()
    plt.xticks(locs,ticks)

  if ylim:
      plt.ylim(ylim)
  if xlim:
      plt.xlim(xlim)
  pylab.title(title)
  ax.set_ylabel(ylabel)
  ax.set_xlabel(xlabel)
  pylab.legend(loc='best')
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def stackGraph(xlabel, ylabel, yvalues, filename,
               scale=u'linear',xvalues=None,ylegend=u'',ticks='', title=u""):
  pylab.title(title)
  fig, ax = plt.subplots()
  stack = ax.stackplot(xvalues, yvalues, label=ylegend)
  if(ylegend):
    proxy_rects = [Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0]) for pc in stack]
    label_list = [ylegend]
    # make the legend
    ax.legend(proxy_rects, label_list, loc=2)
    ax.set_xscale(scale)
  plt.xlim((0,xvalues[-1]))
  ax.set_ylabel(ylabel)
  ax.set_xlabel(xlabel)
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def normalized(elements):
  total = sum(elements)
  return [ v / total for v in elements ]


def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
    """
    Minimize chartjunk by stripping out unnecesasry plot borders and axis ticks

    The top/right/left/bottom keywords toggle whether the corresponding plot border is drawn
    """
    ax = axes or plt.gca()
    ax.spines['top'].set_visible(top)
    ax.spines['right'].set_visible(right)
    ax.spines['left'].set_visible(left)
    ax.spines['bottom'].set_visible(bottom)

    #turn off all ticks
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticks_position('none')

    #now re-enable visibles
    if top:
        ax.xaxis.tick_top()
    if bottom:
        ax.xaxis.tick_bottom()
    if left:
        ax.yaxis.tick_left()
    if right:
        ax.yaxis.tick_right()

def stackedBarChart( ylabel, data, labels, filename, ticks, xlabel=None):
  rows = len(data)
  data = np.array([ normalized(v) for v in data ]).T
  bars = len(data)
  colors = Palette(palette)
  for i in xrange(0, bars):
    bot =np.add.reduce(data[:i])
    plt.bar(range(0,rows), data[i], bottom=bot,
        label=labels[i], color=colors.next_color())

  plt.ylim(0, 1.0)
  plt.xticks([0.35 + v for v in xrange(0,rows+2)], ticks)
  plt.ylabel(ylabel)
  if xlabel:
    plt.xlabel(xlabel)
  legend = plt.legend(loc="upper right")
  remove_border()

  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()



def scatterGraphNaryColour(xlabel, ylabel, xValues,yValues,colorValues,filename,title=u""):
  #A = np.vstack([xvalues, np.ones(len(xvalues))]).T
  #m, c = np.linalg.lstsq(A, yvalues)[0]
  pylab.title(title)
  pylab.ylabel(ylabel)
  pylab.xlabel(xlabel)
  plt.scatter(xValues,yValues,color=colorValues, marker='.',linewidth='0.4',label='Random Seq.')
  #plt.scatter(xRandValues,yRandValues,color='red',marker='x',linewidth='0.6',label='Random Seq.')
  #plt.scatter(xSeqValues,ySeqValues,color='blue',marker='x',linewidth='0.6',label='Natural Seq.')
  plt.xlim(0, 300)
  plt.ylim(0, 300)
  #p1 = plt.scatter(xvalues, yvalues, label=ylegend)
  #plt.plot(xvalues, m*xvalues + c,  label=fitlegend)
  #plt.ylim((0, 10000)
  #plt.xlim((0, max(yvalues)*1.1))
  blue_dot = mlines.Line2D([], [], color='blue', marker='.', markersize=15, label='Beta 1.0')
  red_dot = mlines.Line2D([], [], color='red', marker='.', markersize=15, label='Beta 1.0')
  green_dot = mlines.Line2D([], [], color='green', marker='.', markersize=15, label='Beta 1.0')
  plt.legend([green_dot, red_dot],["Beta 1.5",'Beta 0.5'])
    #plt.legend(['r','b','g'],['sdfsdf''sdfsd''aaa'],loc='upper left')
  #plt.legend(['Beta{}'.format(i) for i in range(len(colorValues)+1)], loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0., fontsize=11)
  #pylab.legend(loc='upper right', frameon=True, fontsize=7)
  pylab.savefig(filename, bbox_inches='tight', frameon=True)
  pylab.close()


def scatterGraphBinaryColour(xlabel, ylabel, xRandValues,xSeqValues, yRandValues,ySeqValues,filename,
                          title=u""):
  #A = np.vstack([xvalues, np.ones(len(xvalues))]).T
  #m, c = np.linalg.lstsq(A, yvalues)[0]
  pylab.title(title)
  pylab.ylabel(ylabel)
  pylab.xlabel(xlabel)
  plt.scatter(xRandValues,yRandValues,color='red',marker='x',linewidth='0.6',label='Random Seq.')
  plt.scatter(xSeqValues,ySeqValues,color='blue',marker='x',linewidth='0.6',label='Natural Seq.')
  plt.xlim(0.0, 2.6)
  plt.ylim(0, 7100)
  #p1 = plt.scatter(xvalues, yvalues, label=ylegend)
  #plt.plot(xvalues, m*xvalues + c,  label=fitlegend)
  #plt.ylim((0, 10000)
  #plt.xlim((0, max(yvalues)*1.1))
  pylab.legend(loc='upper left', frameon=True, fontsize=7)
  pylab.savefig(filename, bbox_inches='tight', frameon=True)
  pylab.close()
  
def scatterGraphFitLineal(xlabel, ylabel, xvalues, yvalues, filename,
                          ylegend=u'',fitlegend='u', ticks='', title=u""):
  A = np.vstack([xvalues, np.ones(len(xvalues))]).T
  m, c = np.linalg.lstsq(A, yvalues)[0]
  pylab.title(title)
  pylab.ylabel(ylabel)
  pylab.xlabel(xlabel)
  p1 = plt.scatter(xvalues, yvalues, label=ylegend)
  plt.plot(xvalues, m*xvalues + c,  label=fitlegend)
  plt.xlim((0, max(xvalues)*1.1))
  plt.ylim((0, max(yvalues)*1.1))
  pylab.legend(loc='best')
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def piechart(labels, values, filename, title):
  def my_autopct(val):
    total=sum(values)
    pct = 100.0*(float(val)/float(total))
    val = int(val)
    return '{p:.2f}%  ({v:d}ms)'.format(p=pct,v=val)

  patches, texts = plt.pie(values, labels=labels, startangle=90, labeldistance=1.05, colors=palette)
  edited_labels = [labels[i] +" "+ my_autopct(values[i]) for i in range(len(values))]
  pylab.axis('equal')

  legend = plt.legend(patches, edited_labels,loc='upper center' , bbox_to_anchor=(0.5, -0.05),
                fancybox=True, shadow=True)
  for pie_wedge in patches:
    pie_wedge.set_edgecolor('white')

  pylab.savefig(filename, bbox_inches="tight")
  pylab.close()

def comparisonBarGraph(ylabel, xvalues, yvalues, filename, rotation=0, xlabel=None,  colors=palette):
  if xlabel:
    pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  ticks = range(0,len(xvalues))
  pylab.bar(ticks, yvalues, align="center", color=colors)
  pylab.xticks(ticks, xvalues, rotation=rotation)
  pylab.legend(loc="best")
  pylab.savefig(filename, bbox_inches="tight")
  pylab.close()

def multiComparativeBarChart(ticks, values, filename, ylabel, width=0.1, rotation=0, loc="best"):
  fig, ax = plt.subplots()
  indexes = np.arange(max(len(v) for u,v in values.items()))
  rects, count = [], 0
  for key, vals in values.items():
      bar = ax.bar(indexes + width * count, vals, width, color=palette[count])
      rects.append(bar)
      count += 1

  ax.set_ylabel(ylabel)
  ax.set_xticks(indexes)
  ax.set_xticklabels(ticks, rotation=rotation)
  ax.legend(rects, values.keys(), loc="best")

  plt.savefig(filename, bbox_inches="tight")
  plt.close()

def initialize():
  mpl.rcParams['savefig.dpi'] = 150

def histogram(xlabel, ylabel, values, nbins, title, filename):
  pylab.hist(values, nbins, histtype="bar",  alpha=0.5)
  pylab.title(title)
  pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def scatter(xlabel, ylabel, filename, xdata, ydata):
  pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  pylab.plot(xdata, ydata, "o-")
  pylab.savefig(filename, bbox_inches='tight')
  pylab.close()

def comparativeScatter(xlabel, ylabel, filename, xdata, ydata, label, ymax=None,xmin=None):
  pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  #frmt = ['-s','-^','-s']
  frmt = ['-o','-s','-^','-s']
  c=0
  for yd, l in zip(ydata, label):
      pylab.plot(xdata, yd, frmt[c],  label=l)
      c+=1
  pylab.legend(loc="best")
  x1,x2,y1,y2 = pylab.axis()
  if ymax:
	pylab.axis((x1,x2,0,ymax))
  if xmin:
	pylab.axis((xmin,x2,0,ymax))
  pylab.savefig(filename, bbox_inches="tight")
  pylab.close()
  


#plot mean, error bars and individual values(with points)
def meanErrorLines(xlabel, ylabel, filename, yErrorValuesSeq, yErrorValuesRand,xMeanValuesRand,xMeanValuesSeq, yMeanValuesRand, yMeanValuesSeq, title,ymax=None,xmin=None):
  pylab.xlabel(xlabel)
  pylab.ylabel(ylabel)
  #frmt = ['-s','-^','-s']
  frmt = ['-o','-s','-^','-s']
  #plt.plot(xMeanValuesSeq,yMeanValuesSeq,label='Natural seq.')
  
  #plt.scatter(xValuesRand,yValuesRand,label='Random seq.')
  plt.errorbar(xMeanValuesSeq,yMeanValuesSeq,yErrorValuesSeq, label='Natural seq.')
  plt.errorbar(xMeanValuesRand,yMeanValuesRand,yErrorValuesRand, label='Random seq.')
  #plt.plot(xMeanValuesRand,yMeanValuesRand,label='Random seq.')
  
  x1,x2,y1,y2 = pylab.axis()
  #pylab.yscale('log')
  if ymax:
	pylab.axis((x1,x2,0,ymax))
  if xmin:
	pylab.axis((xmin,x2,0,ymax))
  pylab.legend(loc="best")
  #pylab.savefig(filename, bbox_inches="tight")
  pylab.show()
  #pylab.close()



def iterationVsX(executionsList,beta,random,maxIterations,logScale,step):
  iterations=range(1,maxIterations,step)
  for x in range(len(executionsList)):
    betaValue=beta[x]
    if betaValue==0.5:
      color='red'
    if betaValue==2.0:
      color='green'  
    if random[x]==True:
      symbol='^'
    else:
      symbol='o'
    pairList=zip(iterations,executionsList[x])
    xvalues=[]
    yvalues=[]
    for p in pairList:
      xvalues.append(p[0])
      yvalues.append(p[1])
    plt.plot(xvalues,yvalues, marker=symbol,color=color,linestyle='-',linewidth=0.5, markersize=7)
    #plt.plot(xvalues,yvalues)
  
  #if logScale:
    #pylab.xscale('log')  
  #x1,x2,y1,y2 = pylab.axis()
  #pylab.axis((x1,maxI,y1,y2))
  pylab.show()
  



def iterationVsXError(executionsList,executionsErrorList,beta,random,maxIterations,logScale,step):
  iterations=range(1,maxIterations,step)
  for x in range(len(executionsList)):
    print 'iteration'
    betaValue=beta[x]
    if betaValue==0.5:
      color='red'
    if betaValue==2.0:
      color='green'  
    if random[x]==True:
      symbol='*'
    else:
      symbol='.'
    pairList=zip(iterations,executionsList[x])
    errorPairList=zip(iterations, executionsErrorList[x])
    xvalues=[]
    yvalues=[]
    yerrorValues=[]
    for p in pairList:
      xvalues.append(p[0])
      yvalues.append(p[1])
    for p in errorPairList:  
      yerrorValues.append(p[1])
    plt.errorbar(xvalues,yvalues,yerrorValues,color=color,marker=symbol,capsize=10,markeredgewidth=1)
    #plt.plot(xvalues,yvalues, marker=symbol,color=color,linestyle='-',linewidth=0.2)
    #plt.plot(xvalues,yvalues)
  
  #if logScale:
    #pylab.xscale('log')  
  #x1,x2,y1,y2 = pylab.axis()
  #pylab.axis((x1,maxI,y1,y2))
  pylab.show()
  
  
initialize()
