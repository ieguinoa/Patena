#import os 

#os.system("perl ps_scan/ps_scan.pl -o pff sequence > salidita")

import subprocess
import os
endl = "\n"
#os.chdir('ps_scan')
sequence="MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLS"
input=open("sequence", "w")
input.write(">gi" + endl)
input.write(sequence)
input.close()
proc = subprocess.Popen(['perl', 'ps_scan/ps_scan.pl','-o', 'pff', 'sequence'],stdout=subprocess.PIPE)
while True:
  line = proc.stdout.readline()
  if line != '':
    #the real code does filtering here
    #print "test:", line.rstrip()
    desde=int(line.split()[1])  
    hasta=int(line.split()[2])
    print "DESDE=" + str(desde)
    print "HASTA=" + str(hasta)
  else:
    break
