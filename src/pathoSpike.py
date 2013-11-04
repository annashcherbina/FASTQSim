'''
This class spikes in silico reads from provided source files into a background FASTQ file.
 @author        Anna Shcherbina (mailto: anna.shcherbina@ll.mit.edu)
License:          GNU GPL license (http://www.gnu.org/licenses/gpl.html)  

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''

import sys 
import io 
import os 
from Params import *
from random import *
from helpers import *  
import numpy
from MultiprocessingHelpers import * 
from multiprocessing import Process, Queue, current_process,freeze_support,Manager

def main():    
    print "Parsing input arguments...\n"
    if len(sys.argv)<5:
        raise SyntaxError("usage: python pathoSpike.py <background fastq> <parameters.csv> [<coverage> <pathogen.fasta>] [-o output prefix] [-threads number of threads, default=1] [-plothistogram]")
    #determine whether to perform plotting 
    plothist=False
    plotCollection=None
    if sys.argv.__contains__('-plothistogram'): 
        plothist=True
        import matplotlib.pyplot as plt
        from plotSpiked import *
        plotCollection=plotSpiked();
        plotindex=sys.argv.index('-plothistogram') 
        sys.argv.pop(plotindex) 

    background_reads = sys.argv[1]
    parameter_file_set = sys.argv[2]
    print "background_reads: " + background_reads + "\n"
    print "parameter_file_set: " + str(parameter_file_set) + "\n"
    


    #check to see if user specified an output file. 
    if "-o" in sys.argv: 
        outindex=sys.argv.index('-o') 
        outputfilename=sys.argv[outindex+1]
        sys.argv.pop(outindex+1)
        sys.argv.pop(outindex)

    else: 
        outputfilename= background_reads.split('.')[0]


    #check to see the number of threads to use 
    if "-threads" in sys.argv: 
        threadindex=sys.argv.index('-threads') 
        numthreads=int(sys.argv[threadindex+1]) 
        sys.argv.pop(threadindex+1) 
        sys.argv.pop(threadindex) 
    else: 
        numthreads=1 

    
    #may have multiple sources of pathogen DNA (i.e. example uses E.coli and plasmid) 
    pathogen_coverage = [] 
    pathogen_source = [] 
    is_circular = [] 
   
    for i in range(0,len(sys.argv)-3,3): 
        pathogen_coverage.append(float(sys.argv[3+i]))
        pathogen_source.append(sys.argv[4+i])
        is_circular.append(sys.argv[5+i])
        
    params = Params(parameter_file_set,background_reads,plotCollection)
    
    readSource = []  
    revcompSource = [] #the reverse complement of the source files. 
    pathogenSourceLength = [] 
    for i in range(len(pathogen_source)):  
        readSrcF = open(pathogen_source[i],'r') 
        readSource.append(readSrcF.read()) 
        readSrcF.close() 
        firstnewline = readSource[i].find('\n') 
        readSource[i] = readSource[i][firstnewline+1:]
        readSource[i] = readSource[i].replace('\n','')
        revcompSource.append(getRevComp(readSource[i])) 
        pathogenSourceLength.append(len(readSource[i])) 
        
    
    spikedReads = []

#If -plothistogram is specified, keep track of actual position/size for insertions/deletions/mutations to ensure that the proper distribution/ parameters are being used to generate these values. This step should be ommitted if greater speed is desired. 
    if plothist:
        actualReadLengths=dict() 
        actualDelPos = dict() 
        actualInsertPos = dict() 
        actualMutPos = dict()         
        actualDelSize = dict() 
        actualInsertSize = dict()
        actualPosCount = dict() 

    task_queue=Queue() 
    done_queue=Queue() 



    #start worker processes 
    procs=[] 
    for i in range(numthreads): 
        p=Process(target=read_gen_worker,args=(task_queue,done_queue,params))
        p.daemon=True
        p.start() 
        procs.append(p) 

    
    print "generating reads\n"
    
    for i in range(len(readSource)):    
        src = readSource[i]
        revComp = revcompSource[i]
        circularity = bool(is_circular[i]) 
        coverage = pathogen_coverage[i] 
        basesToGenerate = round(coverage*len(src)) 
        task_queue.put((src,revComp,circularity,coverage,basesToGenerate,plothist))

    
    #get items from the output queue 
    returned = 0 
    while True: 
        if done_queue.empty()==False: 
            [posCountDict_instance,spikedReads_instance]=done_queue.get() 
            if plothist: 
                actualPosCount=mergeDict(actualPosCount,posCountDict_instance)
            spikedReads=spikedReads+spikedReads_instance
            returned+=1
        elif returned >=len(readSource):
            break; 
        

    #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
    for i in range(numthreads): 
        task_queue.put('STOP')     
    del readSource # we no longer need this variable. 
    del revcompSource
    

    #make sure the daemon processes are dead 
    for p in procs: 
        p.join() 


    print "processing reads (generating insertions, mutations, deletions)...\n" 
    processedReads = [] 


    #dispatch each read to be processed in parallel. 
    task_queue=Queue() 
    done_queue=Queue() 
    #start worker processes 
    procs=[] 
    for i in range(numthreads): 
        p=Process(target=read_mutate_worker,args=(task_queue,done_queue,params))
        p.daemon=True
        p.start() 
        procs.append(p) 
    putval=0;
    for entry in spikedReads:  
        task_queue.put((entry,plothist))
        putval+=1 
    returned = 0; 
    oldval=0; 
    while True: 
        if done_queue.empty()==False: 
            [processed_read,stat_dict]=done_queue.get() 
            returned+=1;
            processedReads.append(processed_read)
            if plothist:
                #update all relevant dictionaries 
                if 'readLength' in stat_dict: 
                    if stat_dict['readLength'] not in actualReadLengths: 
                        actualReadLengths[stat_dict['readLength']]=1 
                    else: 
                        actualReadLengths[stat_dict['readLength']]+=1 
                if 'posCount' in stat_dict: 
                    for posval in stat_dict['posCount']: 
                        if posval not in actualPosCount: 
                            actualPosCount[posval] = 1
                        else: 
                            actualPosCount[posval]+=1 

                if 'mutPos' in stat_dict: 
                    if stat_dict['mutPos'] not in actualMutPos: 
                        actualMutPos[stat_dict['mutPos']]=1 
                    else:
                        actualMutPos[stat_dict['mutPos']]+=1 
                if 'insertPos' in stat_dict: 
                    if stat_dict['insertPos'] not in actualInsertPos: 
                        actualInsertPos[stat_dict['insertPos']]=1 
                    else: 
                        actualInsertPos[stat_dict['insertPos']]+=1 
                if 'insertSize' in stat_dict: 
                    if stat_dict['insertSize'] not in actualInsertSize: 
                        actualInsertSize[stat_dict['insertSize']]=[1,0]  
                    else: 
                        actualInsertSize[stat_dict['insertSize']][0]+=1 
                if 'insertRepeat' in stat_dict: 
                    if stat_dict['insertRepeat'] not in actualInsertSize: 
                        actualInsertSize[stat_dict['insertRepeat']]=[0,1] 
                    else: 
                        actualInsertSize[stat_dict['insertRepeat']][1]+=1 
                if 'delPos' in stat_dict: 
                    if stat_dict['delPos'] not in actualDelPos: 
                        actualDelPos[stat_dict['delPos']]=1 
                    else: 
                        actualDelPos[stat_dict['delPos']]+=1 
                if 'delSize' in stat_dict: 
                    if stat_dict['delSize'] not in actualDelSize: 
                        actualDelSize[stat_dict['delSize']]=[1,0] 
                    else: 
                        actualDelSize[stat_dict['delSize']][0]+=1 
                if 'delRepeat' in stat_dict: 
                    if stat_dict['delRepeat'] not in actualDelSize: 
                        actualDelSize[stat_dict['delRepeat']]=[0,1] 
                    else: 
                        actualDelSize[stat_dict['delRepeat']][1]+=1
        elif returned >=len(spikedReads): 
            print "finished mutating reads!"
            break; 
        else:
            sleep(1)
            if returned==oldval:                 
                print "processsed:"+str(returned)+' spiked reads out of ' + str(len(spikedReads))+'\n' 
            oldval=returned; 



   #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
    for i in range(numthreads): 
        task_queue.put('STOP') 
    #make sure the daemon processes are dead 
    for p in procs: 
        p.join() 

    del spikedReads 

    #plot statistics for readLength, insertions, deletions, and mutations 

    if plotCollection:
        plotCollection.actualReadLengths=actualReadLengths;         
        plotCollection.actualPosCount=actualPosCount; 
        plotCollection.actualDelPos=actualDelPos; 
        plotCollection.actualMutPos=actualMutPos; 
        plotCollection.actualInsertPos=actualInsertPos; 
        plotCollection.actualInsertSize=actualInsertSize; 
        plotCollection.actualDelSize=actualDelSize; 


    #insert spiked reads into the background 
    print "reading in background read block\n"                   
    #read in the portion of the background reads to be spiked. 
    backgroundReadsf = open(background_reads,"r") 
    backgroundSize = os.stat(background_reads).st_size 
    maxToRead = 1048576/10 #.1 MB chunk 
    numIterations = math.ceil(float(backgroundSize)/maxToRead) 
    numToSpike = math.ceil(float(len(processedReads))/numIterations)
     
    allSpikedFastq = [] 
    allInsertedQual=[] 
    allInsertedFasta = [] 
    allInsertedNames=[] 
    
    #parallelization code for read spiking 
    task_queue=Queue() 
    done_queue=Queue() 
    
    procs=[] 
    for i in range(numthreads): 
        p=Process(target=insert_spiked_reads_worker,args=(task_queue,done_queue))
        p.daemon=True 
        p.start() 
        procs.append(p) 
        

    for iter in range(int(numIterations)): 
        backgroundReads= backgroundReadsf.read(maxToRead) 
        startIndex = backgroundReads.index('@') 
        if iter < (int(numIterations)-1):
            endIndex = -1*backgroundReads[::-1].index('@')-1  
        else: 
            endIndex = len(backgroundReads)
        backgroundReads = backgroundReads[startIndex:endIndex] 
        backgroundReads = backgroundReads.split("\n") 
        
        #get the fastq reads from this file 
        #print "starting fastq extract step\n"
        qualwithname,qualDict = fastqextract(backgroundReads)
        #print "completedfastq extract step\n"
                
        #pick a subset of the reads ( randomly) to spike into this portion of the file. 
        subset_processed=dict()
        longest_allowed = max(qualDict.keys()) 

        for i in range(int(numToSpike)): 
            if len(processedReads) > 0: 
                pos = random.randint(0,len(processedReads)-1) 
                new_read = processedReads.pop(pos) 
                if len(new_read)> longest_allowed: 
                    new_read=new_read[0:longest_allowed] 
                if subset_processed.__contains__(len(new_read)): 
                    subset_processed[len(new_read)].append(new_read) 
                else: 
                    subset_processed[len(new_read)]=[new_read]  
        
        #print "inserting reads into background slice " + str(iter)+"\n" 
        #print "starting insertSpikedreads\n"
        task_queue.put((backgroundReads,subset_processed,qualwithname,qualDict))

    #add a 'STOP' to each task queue to indicate that no more tasks will be passed in 
    for i in range(numthreads): 
        task_queue.put('STOP') 

    returned =0;
    while True: 
        if done_queue.empty()==False: 
            [spikedFastq,insertedQual,insertedFasta,insertedReadNames]= done_queue.get() 
            allSpikedFastq=allSpikedFastq + spikedFastq
            allInsertedQual = allInsertedQual+insertedQual 
            allInsertedFasta = allInsertedFasta+insertedFasta 
            allInsertedNames = allInsertedNames+insertedReadNames
            returned+=1; 
            print "spiked "+str(returned) + " file blocks, out of " + str(numIterations)+'\n' 
        elif (task_queue.empty()==False) or (returned < int(numIterations)):
            sleep(1) 
        else: 
            break;     
    backgroundReadsf.close()


    for p in procs: 
        p.join() 
    print "writing spiked reads to separate file\n"
    spikedOutputFastaf = open(outputfilename+"_Spiked.fq","w") 

    spikedreadlengthdict=dict() 
    for pindex in range(len(allInsertedFasta)):
        p = allInsertedFasta[pindex] 
        if spikedreadlengthdict.__contains__(len(p)):
            spikedreadlengthdict[len(p)]+=1 
        else: 
            spikedreadlengthdict[len(p)]=1 
        q = allInsertedQual[pindex]
        name = allInsertedNames[pindex] 
        spikedOutputFastaf.write(name+"\n")
        spikedOutputFastaf.write(p.upper()+"\n") 
        spikedOutputFastaf.write("+\n") 
        spikedOutputFastaf.write(q.upper()+'\n')
    spikedOutputFastaf.close() 

    
    print "writing spiked output file"
    fullOutputFastaf = open(outputfilename+"_Full.fq","w")
    for f in allSpikedFastq: 
        fullOutputFastaf.write(f.upper()+"\n") 
    fullOutputFastaf.close()    

    #generate summary graphs for background distributions, regression curves, and spiked reads 
    if plotCollection:
        plotCollection.plotAll(); 
if __name__ == '__main__':
    main() 
