#!/usr/bin/python
import argparse
import urllib2
import StringIO
import sys
import os
# import re
import time
import errno
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import Seq
from Bio.Blast import NCBIXML
#from Bio import SearchIO
from array import array
import math
import random
#import subprocess
import shutil
import tool_functions

##TODO: CHECK NECESSARY ENVIRONMENT VARS:
#ANCHOR_PATH=$DIR/Tools/ANCHOR
#IUPred_PATH=$DIR/Tools/iupred
#PASTA_PATH=$DIR/Tools/PASTA/pasta_exe
#PROSITE=$DIR/Tools/Prosite/ps_scan


#****************************************

#get base path
def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))


#****************************************

# WEIGHTED SELECTION
# ANON. FUNCTION: USE weighted_choide(param) TO CHOSE AN ITEM OF A PAIRLIST (ID, WEIGHT) , BASED ON WEIGHTS
# RETURNS THE ID OF THE SELECTED ELEMENT
weighted_choice = lambda s : random.choice(sum(([v]*wt for v,wt in s),[]))





### TODO: move all global values inside main function

#*********************
#******GLOBALS*******
#******************

#  OUTPUT FORMATTING
endl = "\n"
tab = "\t"
space=" "
indent=""

#  EXECUTION Id . Identify runs
exeId=os.getpid()





class Mutation: 
    def __init__(self,mutated_sequence,position,previous,replacement):
        self.sequence=mutated_sequence
        self.mutated_position=position
        self.prev_aa=previous
        self.replacement_aa=replacement


#**************************************************************************

##    JUST SUM UP THE INDIVIDUAL SCORES
def get_global_score(scoresList):
  score=0.0
  for listIndex in range(len(scoresList)):
    score= score + scoresList[listIndex]
  return score



#**************************************

def print_evaluation_time(total_elapsed_time,times_dict):
    # before i was doing this calculation:
    # evaluationTime=total_elapsed_time-(pastaTime + anchorTime + tangoTime + blastTime + iupredTime + waltzTime + elmTime + prositeTime + tmhmmTime + limboTime)
    # not sure what was that idea, it doesnt makes sense
    print 'Total Elapsed Time: ' + total_elapsed_timed
    evaluation_time=0
    for key,value in times_dict.items():
        evaluation_time+=value
        print 'Time spent in '+ key +': ' + value
    print 'Time spent on evaluations: ' +  evaluation_time





def print_execution_params(exe_id,beta,length,composition,sequence,evaluate_netcharge,config_params):
    #####   ALWAYS PRINT GENERAL PARAMETERS OF EXECUTION
    #print endl
    print "************************************************"
    print "************************************************"
    print "EXECUTION PARAMETERS:"
    print 'Id=' + str(exe_id) 
    print 'Beta= '+ str(beta)
    print "Length=" + str(length) 
    print "Composition=" + composition
    print "Sequence=" + sequence
    if evaluate_netcharge:
      print "Target net charge=" + str(config_params['targetNetCharge'])
    print "************************************************"
    print "************************************************"





  ######################################################################################
  ##########################       GENERAL SEQUENCE EVALUATION     #####################################
  ######################################################################################



def firstPartialEvaluation(sequence, config_params,position_scores,execution_set,times_dict,inputsPath,job_out_path,step_by_step,detailed_output,verbose):
	#SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
	input=open(inputsPath + "sequenceFASTA"  , "w")
	input.write(">gi" + endl)
	input.write(sequence)
	input.close()
	if verbose:
	   # print endl
	   print "*************************************"
	   print "FIRST PARTIAL EVALUATION"
        for tool_name in tool_functions.tool_functions_dict.keys():
            if tool_name in execution_set:  ## this set has the set of enabled tools to run
                timePrev=time.time()
                if step_by_step:
                    raw_input("Hit enter to continue with next evaluation")
                if verbose:
                    print indent + "*************************************"
                    print "STARTING "+ tool_name +" execution"
                tool_functions.tool_functions_dict[tool_name](sequence, position_scores,config_params, inputsPath,job_out_path,verbose,detailed_output)
                if detailed_output:
                    details_out.write('\n\n***********\n\n' )
                times_dict[tool_name]+=(time.time() - timePrev)
        score_result=get_global_score(position_scores)
        ##PRINT SCORE
        if verbose:
	  print indent + "*************************************"
	  # print endl
	  print indent + "RESULTS OF FIRST PARTIAL EVALUATION:"
	  #print indent + sequence
          print_formatted_scores(sequence,position_scores)
	  print indent + "SCORE:" + str(score_result)
	  print indent + "*************************************"
	if step_by_step:
	  raw_input(indent + "....hit enter to continue")
        return score_result













def secondPartialEvaluation(sequence, position_scores, verbose):
    #SAVE SEQUENCE TO EVALUATE(FASTA FORMAT) IN A FILE
    input=open(inputsPath + "sequenceFASTA"  , "w")
    input.write(">gi" + endl)
    input.write(sequence)
    input.close()
    ##: BLAST SEARCH
    if verbose:
        print indent + "*************************************"
        print indent + "SECOND PARTIAL EVALUATION"
        # print endl
        print indent + "*************************************"
        #print indent + "STARTING BLAST SEARCH"
    timePrev=time.time()
    tool_functions.blastIt(sequence,position_scores,database,inputsPath, verbose,detailed_output)
    #blastTime+=(time.time() - timePrev)
    #if step_by_step:
      #raw_input(indent + "Hit enter to continue with next evaluation")
    if detailed_output:
        details_out.write('\n\n***********\n\n' )
    if step_by_step:
        raw_input(indent + "Press enter to see final results...")
    ##PRINT SCORE
    if verbose:
        print indent + "*************************************"
        # print endl
        print indent + "RESULTS OF SECOND PARTIAL EVALUATION:"
        print_formatted_scores(sequence,position_scores)
        print indent + "SCORE:" + str(get_global_score(position_scores))
        print indent + "*************************************"
    if step_by_step:
        raw_input(indent + "....hit enter to continue")











def mutation_attempt(sequence,weights,aaFrequencies,eval_step,iteration,mutation_attempts,verbose):
    indent = tab + tab    #output formatting
    #SELECT A POSITION
    if verbose:
        # print endl
        print "*************************************"
        print "MUTATION ATTEMPT"
        print "*************************************"
    #CHOOSE A POSITION BASED ON WEIGHTS
    mutate_position= weighted_choice(weights)
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+ eval_step + tab + "SEL POSITION:"  + tab +  str(mutate_position)
    #SELECT THE NEW AA FOR THAT POSITION (BASED ON LIST OF FREQUENCIES)
    previous_residue=sequence[mutate_position]
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "AA MUTATE:"  + tab +  previous_residue
    selected= previous_residue
    #SELECT A NEWONE UNTIL THE RESIDUE IS DIFFERENT FROM PREVIOUS
    while previous_residue == selected:
        selected = weighted_choice(aaFrequencies)
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "NEW:       "  + tab +  selected
    ##BUILD MUTATED SEQUENCE WITH NEW RESIDUE
    mutated_sequence = sequence[0:mutate_position]
    mutated_sequence += selected
    mutated_sequence += sequence[mutate_position+1:]
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+ eval_step + tab +"PREV-SEQ:"  + tab +  sequence
        print str(iteration) + tab + str(mutation_attempts) +tab+ eval_step + tab +"NEW-SEQ:" + tab +  mutated_sequence
    #create Mutation instance with all the info
    return Mutation(mutated_sequence,mutate_position,previous_residue,selected)




def mc_evaluation(beta,partialScore,mutated_score,eval_step,iteration,mutation_attempts,verbose):
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "PREV-SCORE (" + str(partialScore) + ") < MUT-SCORE (" + str(mutated_score) + ")" 
    diff=partialScore-mutated_score
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "SCORE-DIFF:" + str(diff)
    exponent=diff/beta
    if exponent<-100:   #SATURATION
        MCvalue=-100
    else:
        MCvalue=math.exp(exponent)
    if verbose:
        #print "  :" + str(exponent)
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "MC VALUE:" + str(MCvalue)     #e^(dif/beta)
    #GENERATE RANDOM NUMBER BETWEEN 0 AND 1
    randy=random.random()
    if verbose:
        print str(iteration) + tab + str(mutation_attempts) +tab+  eval_step + tab + "RANDOM VALUE [0,1]:" + str(randy)
    #print indent + "MONTE CARLO DECISION:"
    if MCvalue > randy:
        #ACCEPT MUTATION
        return True
    else:
        return False



def print_formatted_scores(sequence, scores):
    data = [sequence,scores]
    col_width = max(len(str(word)) for row in data for word in row)  # padding
    for row in data:
        print  "|".join(str(word).ljust(col_width) for word in row)

#***********************************************************************













#***********************************************************************






#**************************
#**************************
#******* MAIN *******
#**************************





def main():
    #  CUTOFFS, THRESHOLDS, LIMITS , DEFAULTS
    targetScore=0.0 ### not working????
    composition="average"
    #a=r=n=d=c=q=e=g=h=i=l=k=m=f=p=s=t=w=y=v=-1
    userComposition={"A":-999 , "R":-999  , "N":-999  , "D":-999  , "C":-999  , "E":-999 , "Q":-999  , "G":-999  , "H":-999  , "I":-999  , "L":-999  , "K":-999  ,"M":-999  , "F":-999  , "P":-999  , "S":-999  , "T":-999  , "W":-999  , "Y":-999  , "V":-999 }
    database="uniprot_sprot.fasta"
    # sequence="RANDOM"




    #  EXECUTION PARAMETERS
    # testTimes=False   #True=print times of each part
    min_logging=False  #True = only print global scores at end of iteration
    # verbose=False        #True = print detailed information of execution
    # step_by_step=False
    # rand=True
    testing=False
    # change=True
    output=True  ##print info to file
    # global_evaluation=False  #when True, just make a general evaluation of the sequence: run all tools(from each loop) and print the results of each 

    #FILES
    logFileName='mutations' + str(exeId) + '.log'



    #AA FREQUENCIES TO SELECT NEW RESIDUES FOR MUTATIONS (from http://web.expasy.org/protscale/pscale/A.A.Swiss-Prot.html) 
    #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]


    # EXECUTION TIMES OF DIFFERENT PARTS
    times_dict={
        'Tango':0.0,
        'Pasta':0.0,
        'Waltz':0.0,
        'ELM':0.0,
        'Prosite':0.0,
        'Limbo':0.0,
        'Tmhmm':0.0,
        'IUpred':0.0,
        'Anchor':0.0,
        'Amyloid Pattern':0.0,
        'Net charge':0.0
    }

    execution_set={'Tango','Pasta','Waltz','ELM','Prosite','Limbo','Tmhmm','IUpred','Anchor','Amyloid Pattern'}

    evaluateNetCharge = False

    config_params = { 'iupredThreshold': 0.5,
                    'anchorThreshold': 0.5,
                    'waltzThreshold': 79.0,
                    'pastaThreshold': -5.5,
                    'tangoThreshold': 1.0
                   }


    #**************************
    #******* PARAMS *******

    parser = argparse.ArgumentParser(description='Evaluate/Generate linker sequences')
    parser.add_argument('--evaluation-only', dest="global_evaluation",action='store_true',help='Only perform evaluation steps on the sequence. Do NOT attempt mutations.')
    parser.add_argument('--blast-web', dest="blast_web",action='store_true',help='Use BLAST web service instead of local installation.')
    parser.add_argument('--uvsilent', dest="uvsilent",action='store_true',help='UV silent........')
    parser.add_argument('--gettime', dest="testTimes",action='store_true',help='Print execution times for each evaluation.')
    parser.add_argument('--verbose', dest="verbose",action='store_true',help='Verbose output.')
    parser.add_argument('--stepped', dest="step_by_step",action='store_true',help='Ask for user input after each step.')
    parser.add_argument("--detailed-output", action='store', type=argparse.FileType('w'), dest='details_out',help="Directs the output to a name of your choice")
    parser.add_argument('--max-iterations', nargs=1, type=int, default=4000, help='Max amount of iterations')
    parser.add_argument('--length',  type=int, default=12, help='Sequence length')
    parser.add_argument('--beta', nargs=1, type=float, default=0.5, help='Monte Carlo Beta value')
    parser.add_argument('--net-charge', nargs=1, type=int, help='Net charge of the final sequence')
    parser.add_argument('--seq', nargs=1, help='Starting sequence')



    args = parser.parse_args()
    global_evaluation = args.global_evaluation
    blastWeb = args.blast_web
    uvsilent = args.uvsilent
    max_iterations = args.max_iterations
    verbose = args.verbose
    testTimes = args.testTimes
    step_by_step = args.step_by_step
    length = args.length
    beta = args.beta

    ### TODO:check sequence with re
    ## (ARNDCQEGHILKMFPSTWYV)
    sequence=args.seq  #sequence could be None if it was not defined by user


    if step_by_step:
        verbose = True   #MAKES NO SENSE TO GO STEP BY STEP IF CANT SEE A DETAILED OUTPUT

    details_out = args.details_out
    detailed_output = False
    if details_out:
        detailed_output = True

    if args.net_charge:  # != None
        evaluateNetCharge=True
        config_params['targetNetCharge'] = args.net_charge
        execution_set.add('Net charge')  # 'Net charge' is not in the execution_set by default, only included by user request






    ## parameters that may be removed
    # elif (arg=='--minoutput') and (index < len(sys.argv)):
      # min_logging=True
      # logsPath= sys.argv[index+1]
    # elif (arg=='--testoutput') and (index < len(sys.argv)):
      # testing=True
      # testOutputPath = sys.argv[index+1]
        # elif (arg=='--jobid') and (index < len(sys.argv)):
          # exeId=sys.argv[index+1]
        # elif (arg=='--db') and (index < len(sys.argv)):
          # database = sys.argv[index+1]




    runBlast=False
    ### NEED TO ADD:
    # # execution_set={'Tango','Pasta','Waltz','ELM','Prosite','Limbo','Tmhmm','IUpred','Anchor','Amyloid Pattern'}
        # elif (arg=='--noblast'):
          # runBlast=False
        # elif (arg=='--notango'):
          # execution_set.remove('Tango')
        # elif (arg=='--noelm'):
          # execution_set.remove('ELM')
        # elif (arg=='--noiupred'):
          # execution_set.remove('IUpred')
        # elif (arg=='--noanchor'):
          # execution_set.remove('Anchor')
        # elif (arg=='--noprosite'):
          # execution_set.remove('Prosite')
        # elif (arg=='--nolimbo'):
          # execution_set.remove('Limbo')
        # elif (arg=='--notmhmm'):
          # execution_set.remove('Tmhmm')
        # elif (arg=='--nopasta'):
          # execution_set.remove('Pasta')
        # elif (arg=='--nowaltz'):
          # execution_set.remove('Waltz')
        # elif (arg=='--noamyloidpattern'):
          # execution_set.remove('Amyloid Pattern')




        # elif(arg=='-a'):
          # userComposition['A']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-r'):
          # userComposition['R']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-n'):
          # userComposition['N']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-d'):
          # userComposition['D']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-c'):
          # userComposition['C']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-q'):
          # userComposition['Q']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-e'):
          # userComposition['E']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-g'):
          # userComposition['G']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-h'):
          # userComposition['H']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-i'):
          # userComposition['I']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-l'):
          # userComposition['L']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-k'):
          # userComposition['K']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-m'):
          # userComposition['M']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-f'):
          # userComposition['F']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-p'):
          # userComposition['P']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-s'):
          # userComposition['S']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-t'):
          # userComposition['T']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-w'):
          # userComposition['W']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-y'):
          # userComposition['Y']=int(sys.argv[index+1])
          # composition=="user_specified"
        # elif(arg=='-v'):
          # userComposition['V']=int(sys.argv[index+1])
          # composition=="user_specified"















    #***********************************************************************
    ### PRE-RUN SETUP AND PARAMS CHECK
    # ****************************************



    #CHECK IF TARGET NET CHARGE IS POSSIBLE BASED ON SEQUENCE LENGTH (AND PH??)
    if evaluateNetCharge:
        if abs(config_params['targetNetCharge']) > length:
            print 'Net charge is impossible to reach with the specified sequence length'
            exit()



    #PATHS
    #####CREATE INPUT AND OUTPUT DIRS
    try:
        os.makedirs("Input")
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir("Input"):
            pass
    try:
        os.makedirs("Output")
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir("Output"):
            pass

    base_path=get_script_path() + '/'
    # toolsPath=base_path + 'Tools/'    #**************************TODO SET THE PATH TO THE TOOL SET 
    inputsPath=base_path + "/inputs/"+ str(exeId) + "/" #SET PATH TO SAVE INPUTS FILES
    job_out_path=base_path + "/outputs/" + str(exeId) + "/"
    testOutputPath=job_out_path   # DEFAULT OUTPUT FOR TESTs 
    logsPath=job_out_path #default path for log files


    ##### TRY TO CREATE INPUT AND OUTPUT DIRS
    try:
        os.makedirs(inputsPath)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(inputsPath):
            # if exists, just remove it and create it again
            shutil.rmtree(inputsPath)
            os.makedirs(inputsPath)
            pass
    try:
        os.makedirs(job_out_path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(inputsPath):
            # if exists, just remove it and create it again
            shutil.rmtree(job_out_path)
            os.makedirs(job_out_path)
            pass


    #OUTPUT
    #outputPath = "Output/"   

    ###########THESE FILES SHOULD GO TO test/results DIRECTORY ***************************
    #scoresFile="scores" + str(length)   #save Scores Vs iteration number 
    #mutation_attemptsFile="mutationsAttempt" + str(length)  # save number of mutation attempts  Vs iteration number
    #timesFile='times' + str(length)   #save times Vs iteration number
    #totalTimesFile='totalTimes'  #save total time elapsed Vs sequence length
    #if output:
      #totalTimesOutputFile=open( testOutputPath + totalTimesFile , "a")
      #timesOutputFile=open( testOutputPath + timesFile , "a")
      #scoresOutputFile=open( testOutputPath + scoresFile , "a")
      #mutationsFile=open( testOutputPath  + mutation_attemptsFile , "a")
    #*************************************************

    #TESTING MODE ON: WRITE EXECUTION PARAMETERS (TOTAL TIME, TIME PER BLOCK, MUT-ATTEMPTS, BETA VALUES, SCORES)
    if testing: 
        testOutputFile=open(testOutputPath + '/'+ str(exeId), 'w')
        if rand:
            testOutputFile.write('RAND'+ tab +str(beta) + tab + str(length) + endl)
        else:
            testOutputFile.write('SEQ'+ tab + str(beta) + tab + str(length) + endl)

    if min_logging:
        #CREATE .log FILE
        log_file_stream=open( logsPath+'/'+logFileName, 'w')

    #AMINOACIDS FREQUENCIES....
    # THESE FREQUENCIES ARE USED TO SELECT REPLACEMENTS DURING MUTATIONS
    # THE SELECTION IS MADE USING A WEIGTHED SELECTION. THE ONLY REQUERIMENT IS THAT THE FREQUENCIES (WEIGHTS) ARE WHOLE NUMBERS

    # STANDARD COMPOSITION 
    #if uvsilent==False:   
    standardComposition={"A":825 , "R":553 , "N":406 , "D":545 , "C":137 , "E":393 , "Q":675 , "G":707 , "H":227 , "I":596 , "L":966 , "K":548 ,"M":242 , "F":386 , "P":470 , "S":656 , "T":534 , "W":108 , "Y":292 , "V":687}
      #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]
    #else:
      #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",0),("P",470),("S",656),("T",534),("W",0),("Y",0),("V",687) ]

    #print str(standardComposition)

    if uvsilent:
    # RESET Y,W,F FREQ = 0 
        userComposition['Y']=0
        userComposition['W']=0
        userComposition['F']=0

    if (composition=="user_specified"):   
    # USER HAS DEFINED AT LEAST ONE OF THE FREQUENCIES, THE FREQUENCIES DEFINED ARE IN
    #FIRST CHECK IF THE SUM OF FREQUENCIES DEFINED IS LESS THAN 100 percent
        freqSum=0
        for key in userComposition:
            if userComposition[key] != -999:    #THE USER HAS DEFINED THIS FREQUENCE
                freqSum+= userComposition[key]
        if freqSum > 100:
            print 'Total defined frequencies exceeded 100%'
            exit()


    #REDEFINE TOTAL aaFrequencies USING THE USER DEFINED FREQUENCIES
    #TODO ***********************************************************
    # TEMP. SOLUTION: JUST COPY THE FREQUENCIES DEFINED BY USER AND REPLACE UNKNOWN BY STANDARD
    for key in userComposition: 
        if userComposition[key] == -999:
          userComposition[key] = standardComposition[key]
        else: 
          userComposition[key] = userComposition[key]*100  #change the frequencies to 10000 base



    #CONVERT FREQUENCIES TO LIST OF PAIRS (AA, FREQ)
    aaFrequencies=userComposition.items()
    #print str(aaFrequencies)
    #aaFrequencies= [("A",825), ("R",553),("N",406),("D",545),("C",137),("E",393),("Q",675),("G",707),("H",227),("I",596),("L",966),("K",548),("M",242),("F",386),("P",470),("S",656),("T",534),("W",108),("Y",292),("V",687) ]

    # FORMAT TO REQUEST RANDOM SEQUENCE:   
    #http://web.expasy.org/cgibin/randseq/randseq.pl?size=100&comp=user_specified&A=10&R=10&N=10&D=10&C=10&Q=10&E=10&G=10&H=0&I=0&L=0&K=0&M=0&F=0&P=0&S=0&T=0&W=0&Y=10&V=0&output=fasta   

    #if rand==True:
            ##****************GET RANDOM SEQUENCE*************
            ##print endl
            #if verbose:
              #print "Generating random sequence..."    
            ##print endl
            #if not (composition=="user_specified"):
              #url="http://web.expasy.org/cgi-bin/randseq/randseq.pl?size=" + str(length) + "&comp=" + composition + "&output=fasta"
            #else:
              #print "fix this"
            ##print url
            #response = urllib2.urlopen(url)
            #html = response.read()
            #i = html.index('\n')
            #sequence = html[i+1:].replace('\n', '')
            ##if verbose:
               ##print "*******************************" 

    #GENERATE RANDOM SEQUENCE WITH THE DEFINED COMPOSITION
    if not sequence:
        sequence=[]
        for x in range(0,length):
            sequence.append(weighted_choice(aaFrequencies)) 
        sequence=''.join(sequence)

    print_execution_params(exeId,beta,length,composition,sequence,evaluateNetCharge,config_params)
    if step_by_step:
        raw_input("Hit enter to start initial evaluation")
    #CREATE ARRAY TO SAVE MUTATION FREQUENCY
    position_scores=[]
    partialScores=[]
    mutatedScores=[]
    for p in range(len(sequence)):
        position_scores.append(0)
        mutatedScores.append(0)
        partialScores.append(0)

    # TIME
    time0 = time.time()   #start time
    timePrev=time.time()  #used to measure execution times of different parts (evaluation, mutations, etc)


    ##################################
    #########EVALUATE INITIAL SEQUENCE
    ##################################

    for p in range(len(sequence)):
        #position_scores.append(0)
        position_scores[p]=0
        partialScores[p]=0
    if verbose:
        print endl
        print "*****************************"
        print " INITIAL EVALUATION "
        print "*****************************"    
    #MAKE BOTH PARTIAL EVALUATIONS TO GET A GLOBAL SCORE

    #FIRST SET OF EVALUATIONS
    firstPartialScore=firstPartialEvaluation(sequence, config_params,partialScores,execution_set,times_dict,inputsPath,job_out_path,step_by_step,detailed_output,verbose)

    #SAVE RESULTS
    # firstPartialScore=get_global_score(partialScores)
    firstPartialScores=partialScores

    #ADD THE SCORE TO THE GLOBAL SCORE AND RESET PARTIAL LIST
    for p in range(len(sequence)):
        position_scores[p]=position_scores[p]+partialScores[p]
        partialScores[p]=0
      
    secondPartialScore = 0  
    if runBlast:
        #SECOND PART OF EVALUATION
        secondPartialEvaluation(sequence, partialScores,verbose)
        secondPartialScore=get_global_score(partialScores)

        #ADD THE SCORE TO THE GLOBAL SCORE AND RESET PARTIAL LIST
        for p in range(len(sequence)):
            position_scores[p]=position_scores[p]+partialScores[p]
            partialScores[p]=0

    ##SUM OF SCORES LIST
    global_score=get_global_score(position_scores)

    if detailed_output:
        index=0
        details_out.write('\n')
        details_out.write('*************************\n')
        details_out.write('\n')
        #details_out.write("First partial score"+ tab  + str(firstPartialScore) + '\n')
        #details_out.write("Second partial score" +tab+ str(secondPartialScore) + '\n')
        details_out.write("Global score"+tab + str(global_score) + '\n')
        #details_out.write('\n')
        #details_out.write('*************************\n')
        details_out.write('\n')
        details_out.write('Scores per position:\n')
        details_out.write('Pos' +tab+'AA' +tab+ 'Score\n')
        for aa in sequence:
            details_out.write(str(index) + tab + aa + tab + str(position_scores[index]))
            index+=1
            details_out.write('\n')
        #for score in position_scores:
        #        details_out.write(str(score) + tab)
        #data = [sequence,position_scores]
        #col_width = max(len(str(word)) for row in data for word in row)  # padding
        #for row in data:
        #      details_out.write("|".join(str(word).ljust(col_width) for word in row))


    if verbose:
        print "*******************************************"
        print "INITIAL EVALUATION RESULTS"
        print "First partial score   : "  + str(firstPartialScore)
        print "Second partial score  : " + str(secondPartialScore) 
        print "Global score          : " + str(global_score)
        print endl
        print "*******************************************"
        print "*****************************"
        print "*****************************"
        print endl
        print endl
        if step_by_step:
            raw_input("Hit enter to proceed with mutations")

    if min_logging:
        log_file_stream.write('ISEQ' + tab + sequence + tab + str(global_score) + endl)
        #print 'INITIAL SEQ:   ' + sequence + tab + str(global_score)

    if testing:
        testOutputFile.write('ISEQ' + tab + sequence + tab + 'FIRST'  + tab + str(firstPartialScore) + endl)
        testOutputFile.write('ISEQ' + tab + sequence + tab + 'SECOND' + tab + str(secondPartialScore) + endl)
        testOutputFile.write('ISEQ' + tab + sequence + tab + 'GLOBAL' + tab + str(global_score) + endl)








    ######################################################
    ######################################################
    #########  GLOBAL LOOP  ##############################
    ######################################################
    ######################################################

    iteration=1
    global_iteration=0

    while global_score > 0 and iteration <= max_iterations and (not global_evaluation):
        if verbose:
            print "*****************************"
            print "STARTING GLOBAL ITERATION " + str(global_iteration)
            print "*****************************"
            print_formatted_scores(sequence,position_scores)
            print "Current global score:     " + str(global_score)
        if step_by_step:
            raw_input("Hit enter to start first round of evaluations")
            raw_input('')


        #################################################
        ######  FIRST ROUND OF EVALUATIONS - MUTATIONS
        #################################################
        if verbose:
            print "FIRST ROUND OF MUTATIONS: DECISION IS BASED ON (RESULTS OF)FIRST SET OF TOOLS" 
        partialScore=firstPartialScore
        while partialScore > 0 and iteration <= max_iterations:
            timePrev=time.time()
            weights=[]
            #weights IS A PAIRLIST(position,weight)
            #CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
            for x in range(len(position_scores)):
                weights.append((x, position_scores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
            mutation_attempts=0       #COUNT MUTATIONS ATTEMPTS
            while 1000 > mutation_attempts:    ##JUST A SYMBOLIC MAX. AMOUNT OF MUTATIONS ATTEMPTS
                mutation_attempts+=1
                # get a possible mutation
                mutation=mutation_attempt(sequence,weights,aaFrequencies, "EVAL1",iteration,mutation_attempts,verbose)
                ##RESET LIST OF SCORES FOR THE MUTATED SEQUENCE
                for p in range(len(sequence)):
                    mutatedScores[p]=0
                if step_by_step:
                    raw_input("...Hit enter to start evaluation")
                if verbose:
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "STARTING PROPOSED MUTATION EVALUATION"
                mutatedScore=firstPartialEvaluation(mutation.sequence, config_params, mutatedScores,execution_set,times_dict,inputsPath,job_out_path,step_by_step,detailed_output,verbose)
                #if step_by_step:
                  #raw_input(indent + "Hit enter to continue with mutation acceptance")
                if verbose:
                    print "*************************************"
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "MUTATION ATTEMPT RESULTS"
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "PREV-SEQ:"
                    print_formatted_scores(sequence,position_scores)
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "PREV-SCORE:" + str(partialScore)
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "MUT-SEQ:"
                    print_formatted_scores(mutation.sequence,mutatedScores)
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "MUT-SCORE:" + str(mutatedScore)
                    # print "*************************************"
                    # print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "ATTEMPT RESULTS"
                if partialScore >= mutatedScore: #get_global_score(position_scores):
                    if verbose:
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "PREV-SCORE (" + str(partialScore) + ") >= MUT-SCORE (" + str(mutatedScore) + ")" 
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "ACCEPT MUTATION"
                    if step_by_step:
                        raw_input("")
                    break
                else:
                    #decide using MC whether the mutation is accepted
                    mc_decision =  mc_evaluation(beta,partialScore, mutatedScore,"EVAL1",iteration, mutation_attempts,verbose)
                    if mc_decision:
                        if verbose:
                            #if true then accept the mutation...no need to keep trying new mutations
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "PREV-SCORE (" + str(partialScore) + ") < MUT-SCORE (" + str(mutatedScore) + ")" 
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "ACCEPT MUTATION - MONTE CARLO DECISION"
                        break
                    else:
                        if verbose:
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "PREV-SCORE (" + str(partialScore) + ") < MUT-SCORE (" + str(mutatedScore) + ")" 
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "DENY MUTATION - MONTE CARLO DECISION"

                if verbose:
                    print "*************************************"


            #####END OF MUTATION ITERATION
            if mutation_attempts < 10000:   #MAKE SURE LOOP ENDED BY MUTATION ACCEPT
                ###NOW THE SEQUENCE IS THE MUTATED SEQUENCE
                sequence=mutation.sequence
                #AND THE POSITION SCORES ARE THE ONES CORRESPONDING TO THE MUTATED SEQUENCE
                position_scores=mutatedScores
                #AND THE GLOBAL SEQUENCE SCORE IS THE ONE CORRESPONDING TO THIS NEW SEQUENCE
                partialScore = mutatedScore # get_global_score(position_scores)
                if verbose:
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "ATTEMPTS MADE:" +tab+ str(mutation_attempts)
                    print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL1" + tab + "END"
                    # print "(PARTIAL) score :    " + str(partialScore)
                    print "*******************************************"
                    print endl
                if min_logging:
                    ## TODO : CHECK THIS FIX!!!! ***fix this, i should print the new selected AA but after refactoring to mutation_attempt() method I dont get that anymore
                    # change mutation_attempt() to retun an Object of a class that has the 
                    log_file_stream.write(str(mutation.mutated_position) + '(' + mutation.prev_aa + ')->' + mutation.replacement_aa + tab + sequence + endl)
                    # log_file_stream.write(str(mutatePosition) + '(' + previousResidue + ')->' + tab + sequence + tab +str(partialScore) + tab + '1' + endl)
                    #log_file_stream.write(str(mutatePosition) + '(' + previousResidue + ')->' + seleccionado + '   '+ tab + sequence + tab +str(partialScore) + tab + '1' + endl)
            else:
                if verbose:
                    ### EXCEEDED THE NUMBER OF ATTEMPTS, SEQUENCE NOT CHANGED
                    print  str(iteration) + tab + str(mutation_attempts) +tab+ "EVAL1" + tab + "ATTEMPTS EXCEEDED" 

            if testing:
                timeX=time.time()-timePrev   #ITERATION TIME
                testOutputFile.write('LOOP1' + tab + str(iteration)+ tab + str(mutation_attempts) + tab + str(partialScore) + tab + str(timeX) + endl ) 
            if step_by_step:
                raw_input("Hit enter to continue with next iteration")
            iteration=iteration+1   #TOTAL NUMBER OF ITERATIONS


        #RESET LIST OF PARTIAL SCORES
        for p in range(len(sequence)):
            partialScores[p]=0


        if runBlast:
            ####################################################
            ##### SECOND ROUND OF EVALUATION - MUTATION  #######
            ####################################################

            secondPartialEvaluation(sequence,partialScores, verbose)
            partialScore=get_global_score(partialScores)
            while partialScore > 0 and iteration <= max_iterations:
                timePrev=time.time()
                weights=[]
                #weights IS A PAIRLIST(position,weight)
                #CONTAINS THE WEIGHT USED TO SELECT THE MUTATION POSITION. EACH ELEMENT IS A PAIR (X, WEIGHT), WHERE X= POSITION AND EIGHT IS = (SCORE(X) + A BASE WEIGHT)
                for x in range(len(partialScores)):
                  weights.append((x, partialScores[x]+1))    #the weight is score+1 - this gives a slight chance to all the position to suffer mutation
                #MUTATION ATTEMPTS
                #INSIDE THIS LOOP, THE VARIABLE partialScore IS NOT MODIFIED. ONLY COMPARED TO THE MUTATED SCORES
                #THE LOOP ENDS WHEN A MUTATION IS ACCEPTED OR WHEN A LIMIT OF ATTEMPTS IS REACHED
                mutation_attempts=0       #COUNT MUTATIONS ATTEMPTS
                while 10000 > mutation_attempts:    ##JUST A SYMBOLIC MAX. AMOUNT OF MUTATIONS ATTEMPTS
                    mutation_attempts+=1
                    indent = tab + tab    #output formatting

                    ## Attempt mutation
                    mutation=mutation_attempt(sequence,weights,aaFrequencies,"EVAL2",iteration,mutation_attempts,verbose)
                    ##RESET LIST OF SCORES FOR THE MUTATED SEQUENCE
                    for p in range(len(sequence)):
                        mutatedScores[p]=0
                    if step_by_step:
                        raw_input("...Hit enter to start evaluation")
                    if verbose:
                        print str(iteration) + tab + str(mutation_attempts) +tab+ "EVAL2" + tab + "START MUT EVAL"
                    secondPartialEvaluation(mutation.sequence, mutatedScores, verbose)
                    mutatedScore=get_global_score(mutatedScores)
                    #if step_by_step:
                      #raw_input(indent + "Hit enter to continue with mutation acceptance")
                    #IF THE GLOBAL SCORE DECREASED
                    if verbose:
                        print "*************************************"
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "MUTATION ATTEMPT RESULTS"
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SEQ:"
                        print_formatted_scores(sequence,position_scores)
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SCORE:" + str(partialScore)
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "MUT-SEQ:"
                        print_formatted_scores(mutation.sequence,mutatedScores)
                        print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "MUT-SCORE:" + str(mutatedScore)
                        # print "*************************************"
                        # print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "DECISION"
                        # print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SEQ"
                        # print_formatted_scores(sequence,partialScores)
                        # print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab +"SCORE: " + str(partialScore)
                        # print str(iteration) + tab + str(mutation_attempts) +tab+ "EVAL2" + tab + "Mutated sequence"
                        # print indent + mutation.sequence
                        # print_formatted_scores(sequence,mutatedScores)
                    if partialScore >= mutatedScore:
                        if verbose:
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SCORE (" + str(partialScore) + ") >= MUT-SCORE (" + str(mutatedScore) + ")" 
                            print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "ACCEPT MUTATION"
                            # print indent + "Previous score (" + str(partialScore) + ") >= Mutated score (" + str(mutatedScore) + ")" 
                            # print indent + "...ACCEPT MUTATION"
                        if step_by_step:
                            raw_input("")
                        break
                          #raw_input(indent + "Hit enter to continue with next iteration")
                        #return mutatedSequence
                    else:
                        mc_decision =  mc_evaluation(beta,partialScore, mutatedScore, "EVAL2",iteration, mutation_attempts,verbose)
                        if mc_decision:
                            if verbose:
                                print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SCORE (" + str(partialScore) + ") < MUT-SCORE (" + str(mutatedScore) + ")" 
                                print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "ACCEPT MUTATION - MONTE CARLO DECISION"
                            #if true then accept the mutation...no need to keep trying new mutations
                            break
                        else:
                            if verbose:
                                print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "PREV-SCORE (" + str(partialScore) + ") < MUT-SCORE (" + str(mutatedScore) + ")" 
                                print str(iteration) + tab + str(mutation_attempts) +tab+  "EVAL2" + tab + "DENY MUTATION - MONTE CARLO DECISION"
                    if verbose:
                        print "*************************************"
                #####END OF MUTATION ITERATION

                #MAKE SURE LOOP ENDED BY MUTATION ACCEPT
                if mutation_attempts < 10000:
                    #print "Sequence after mutation:    " + mutatedSequence
                    ##NOW THE SEQUENCE IS THE MUTATED SEQUENCE
                    sequence=mutatedSequence
                    #AND THE SCORES ARE THE ONES CORRESPONDING TO THE MUTATED SEQUENCE
                    partialScores=mutatedScores
                    #AND THE GLOBAL SEQUENCE SCORE IS THE ONE CORRESPONDING TO THIS NEW SEQUENCE
                    partialScore= get_global_score(partialScores)
                    if verbose:
                        # print endl
                        print "Attempts before mutation accept:" + str(mutation_attempts)
                        #print endl
                        print "*******************************************"
                        print "End of (PARTIAL) iteration " + str(iteration)
                        print "(PARTIAL) score :    " + str(partialScore)
                        print "*******************************************"
                        # print endl
                    if min_logging:
                        log_file_stream.write(str(mutation.mutated_position) + '(' + mutation.prev_aa + ')->' + mutation.replacement_aa + tab + sequence + endl)
                        #log_file_stream.write(str(mutatePosition) + '(' + previousResidue + ')->' + seleccionado + tab + sequence + endl)
                        #print  mutatePosition + '(' + previousResidue + ') -> ' + seleccionado + tab + sequence
                else:
                    ### EXCEEDED THE NUMBER OF ATTEMPTS, SEQUENCE NOT CHANGED
                    if verbose:
                        print  " Mutations attempts exceeded " 
                if testing:
                    timeX=time.time()-timePrev   #ITERATION TIME
                    testOutputFile.write('LOOP2' + tab + str(iteration)+ tab + str(mutation_attempts) + tab + str(partialScore) + tab + str(timeX) + endl ) 
                    # MUTATION ATTEMPTS Vs. ITERATION NUMBER
                    #mutationsFile.write(str(iteration)+ tab + str(mutation_attempts) + tab +str(beta) + endl)
                    #SCORES Vs. ITERATION NUMBER
                    #scoresOutputFile.write(str(iteration)+ tab + str(partialScore) + tab + str(beta) + endl)
                    #ITERATION ELAPSED TIME
                    #timesOutputFile.write(str(iteration)+ tab + str(timeX) + tab + str(beta) + endl)
                if step_by_step:
                    raw_input("Hit enter to continue with next iteration")
                iteration=iteration+1   #TOTATL NUMBER OF ITERATIONS
            if verbose:
                print "SECOND ROUND OF EVALUATIONS FINISHED"
                if step_by_step:
                    raw_input("HIT ENTER TO GET CURRENT GLOBAL SCORE")


        ##AT THE END OF THE SECOND ROUND, THE SCORE CORRESPONDING TO THE SECOND PARTIAL EVALUATION IS 0 (EXCEPT WHEN THE NUMBER OF ITERATIONS IS EXCEEDED)....
        #THEN, GLOBAL SCORE IS DEFINED BY THE FIRST PARTIAL EVALUATION

        for p in range(len(sequence)):
            position_scores[p]=0
        global_score=firstPartialEvaluation(sequence, config_params,position_scores,execution_set,times_dict,times_dict,inputsPath,job_out_path,step_by_step,detailed_output,False )
        #global_score=get_global_score(position_scores)
        global_iteration=global_iteration+1;  
        #PRINT RESULTS OF GLOBAL ITERATION
        if verbose:
            print "*******************************************"	     
            print "End of global iteration " + str(global_iteration)
            print "Global score :    " + str(global_score)
            print "*******************************************"
            print endl
        if step_by_step:
            raw_input("Hit enter to continue with next iteration")
        if testing:
            testOutputFile.write('GLOBAL ' + tab + str(global_iteration) + endl )



    ###########################################
    #### GLOBAL LOOP END ####################
    #########################################

    if not global_evaluation:
        print "**END OF SEARCH**"
        if global_score==0:
            print "REACHED SCORE = 0"
            print 'FINAL SEQUENCE: ' + sequence
        else:
            print "REACHED LIMIT OF ITERATIONS"
            print 'LAST SEQUENCE: ' + sequence
            print 'FINAL SCORES'
            print_formatted_scores(sequence, position_scores)
            print "Global score: " + str(global_score)
        if min_logging:
            log_file_stream.write('END' + tab + str(global_score) + tab + sequence + endl)  

        #PERFORMANCE TESTS OUTPUT
        total_elapsed_time=time.time() - time0   ##
        #if output:        
          #print "Elapsed time:",total_elapsed_time, "Seconds"
          #totalTimesOutputFile.write(str(length) + tab + str(iteration-1)+ tab + str(total_elapsed_time) + tab + str(beta) +  endl)
        if testing:
            testOutputFile.write('END' + tab + str(total_elapsed_time) + tab + str(global_score) + endl) 
            testOutputFile.close()
            ##CLOSE ALL OUTPUT FILES
            #totalTimesOutputFile.close()
            #timesOutputFile.close()
            #scoresOutputFile.close()
            #mutationsFile.close()
        if testTimes:
            print_evaluation_time(total_elapsed_time,times_dict)



if __name__ == '__main__':
    main()
