import urllib2
import StringIO
import sys
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML
#from Bio import SearchIO

#*********************
#******GLOBALS*******
#******************

endl = "\n"
tab = "\t"



#***********************************************************************

def printHelp():
  print endl + "Usage: "
  print "  python linkeado.py [options]"
  print endl + "   where options are:"
  print tab + "--length :" + tab + "Sequence lenght"
  print tab + "--composition :" + tab + "average | user_specified"
  
  
  
  
  
  
#**************************
#******* MAIN *************
#**************************


#********DEFAULTS********
size=10
composition="average"
a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=0


#*******PROCESO LOS ARGUMENTOS******
if len(sys.argv) < 2:
  print "\n Use default input: size=10 composition=average "
  print "\n If you want a specific configuration, see the help running:"
  print tab + "python linkeado.py -h"
else:
  
  for indice in range(1,len(sys.argv)):
    arg = sys.argv[indice]
    
    if (arg=='-H' or arg== '-help' or arg== '--help'):
      printHelp()
      exit()
    elif (arg=='--length') and (indice < len(sys.argv)):
      size = int(sys.argv[indice+1])
    elif (arg== '--composition') and (indice < len(sys.argv)):
      composition = sys.argv[indice+1]
      if (composition=="user_specified"):  #DEBERIA RECIBIR LAS FRECUENCIAS DE TODOS LOS AMINOACIDOS
	for j in range(indice+2,len(sys.argv),2):
	  if(sys.argv[j]=='-a'):
	    a=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-r'):
	    r=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-n'):
	    n=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-d'):
	    d=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-c'):
	    c=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-q'):
	    q=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-e'):
	    e=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-g'):
	    g=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-h'):
	    h=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-i'):
	    i=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-l'):
	    l=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-k'):
	    k=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-m'):
	    m=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-f'):
	    f=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-p'):
	    p=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-s'):
	    s=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-t'):
	    t=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-w'):
	    w=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-y'):
	    y=int(sys.argv[j+1])
	  elif(sys.argv[j]=='-v'):
	    v=int(sys.argv[j+1])
	    
	    
#http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   


#****************OBTENGO LA SECUENCIA RANDOM *************    
if not (composition=="user_specified"):
  url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(size) + "&comp=" + composition + "&output=fasta"
else:
  url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(size) + "&comp=user_specified" + "&A=" + str(a) + "&R=" + str(r) + "&N=" + str(n) + "&D=" + str(d) + "&C=" + str(c) + "&Q=" + str(q) + "&E="+ str(e) + "&G=" + str(g) +"&H="+ str(h) + "&I=" + str(i) + "&L="+ str(l) + "&K="+ str(k) + "&M="+ str(m) + "&F="+ str(f) + "&P="+ str(p) + "&S=" + str(s) +"&T=" + str(t) +"&W="+ str(w) + "&Y=" + str(y) +"&V=" +str(v) + "&output=fasta" 

#print url
response = urllib2.urlopen(url)
html = response.read()
i = html.index('\n')
sequence = html[i+1:].replace('\n', '')
print sequence


#handle = open("D.rerio_calcineurin.gb")
#records = SeqIO.parse(handle,"genbank")
#first_seq = records.next().seq
result = NCBIWWW.qblast("blastp", "nr", "AAMVLSEGEWQLVLHVWAKVEADVAGHGQDIIILLFKSHPETLEKFDRFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKK")
#print result.read()
#result = NCBIWWW.qblast("blastp", "nr", "AA")
records = NCBIXML.parse(result)
#print records.next()
first = records.next()
#if first is None:
##  print "No hits found"
#else:
#for alignment in blast_record.alignments:
if len(first.alignments) > 0:
 firstAlign=first.alignments[0]
 for hsp in firstAlign.hsps:
  if hsp.expect < 0.1:
   print "****Alignment****"
   print "sequence:", firstAlign.title
   print "length:", firstAlign.length
   length=firstAlign.length
   print "E-value:", hsp.expect
   print hsp.query[0:firstAlign.length] 
   print hsp.match[0:firstAlign.length] 
   print hsp.sbjct[0:firstAlign.length] 
else:
 print "No hits found"


#analisis de las secuencias
for j in range(len(hsp.query)):
  #if hsp.query[j] == hsp.sbjct[j]:
	print j + ": " + hsp.query[j]
        #print hsp.sbjct[j]
  
 
