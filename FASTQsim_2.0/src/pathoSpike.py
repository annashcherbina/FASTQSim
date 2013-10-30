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


def main():    
    print "Parsing input arguments...\n"
    if len(sys.argv)<5:
        raise SyntaxError("usage: python pathoSpike.py <background fastq> <numSlices> <targetSlice> <parameters.csv> [<coverage> <pathogen.fasta>]")
    #determine whether to perform plotting 
    plothist=False
    plotCollection=None
    if sys.argv.__contains__('-plothistogram'): 
        plothist=True
        import matplotlib.pyplot as plt
        from plotSpiked import *
        plotCollection=plotSpiked();
    toremove=None; 
    for i in sys.argv:
        if i.__contains__('plothistogram'):
            toremove=i
            break; 
    if toremove:
        sys.argv.remove(toremove)
    background_reads = sys.argv[1]
    numSlices = int(sys.argv[2]) 
    targetSlice = int(sys.argv[3])
    parameter_file_set = sys.argv[4]
    print "background_reads: " + background_reads + "\n"
    print "numSlices: " + str(numSlices) + "\n"
    print "targetSlice: " + str(targetSlice) + "\n"
    print "parameter_file_set: " + str(parameter_file_set) + "\n"
    

    #check to see if user specified an output file. 
    if "-o" in sys.argv: 
        outindex=sys.argv.index('-o') 
        outputfilename=sys.argv[outindex+1]
        sys.argv.remove(outputfilename)
        sys.argv.remove('-o') 
    else: 
        outputfilename= background_reads.split('.')[0]

    
    #may have multiple sources of pathogen DNA (i.e. example uses E.coli and plasmid) 
    pathogen_coverage = [] 
    pathogen_source = [] 
    is_circular = [] 
    for i in range(0,len(sys.argv)-5,3): 
        pathogen_coverage.append(float(sys.argv[5+i]))
        pathogen_source.append(sys.argv[6+i])
        is_circular.append(sys.argv[7+i])
        
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
    actualPosCount = dict() 
    print "generating reads\n"
    for i in range(len(readSource)): 
        
        src = readSource[i]
        revComp = revcompSource[i]
        circularity = bool(is_circular[i]) 
        coverage = pathogen_coverage[i] 
        basesToGenerate = round(coverage*len(src)) 
        basesGenerated =0 
        
        while basesGenerated < basesToGenerate: 
            #choose read length and position to insert the read within the background fastq file.             
            sourceToUse=src
            prob = 1 
            readLength = 0 
            prob = round(numpy.random.uniform(),3) 
            readLength = int(cdfToLength(prob, params.readLengthsCDF))

            #choose read position randomly
            startPos = random.randint(0,len(src)-Params.minLength)                       

            
            #determine whether we use the source sequence or its reverse complement. 
            if numpy.random.random() < 0.5:
                sourceToUse = revComp
            
            if circularity == False: 
                read = sourceToUse[startPos:min(startPos+readLength,len(sourceToUse))] 
                readLength = len(read) 
                #keep track of bases covered in the log file 
                s = startPos 
                e = min(startPos+readLength,len(sourceToUse)) 
                if sourceToUse==revComp: 
                    enew = len(sourceToUse) - s -1 
                    snew = len(sourceToUse)- e -1
                    e =enew 
                    s = snew  
            else:    
                if len(sourceToUse) > startPos+readLength: 
                    read = sourceToUse[startPos:startPos+readLength] 
                    s = startPos 
                    e = startPos+readLength 
                    if sourceToUse==revComp: 
                        snew = len(sourceToUse)-e-1 
                        enew = len(sourceToUse)-s -1
                        s = snew 
                        e = enew  
                else:
                    overflow = startPos+readLength - len(sourceToUse) 
                    read = sourceToUse[startPos:]+sourceToUse[0:overflow]
                    s = startPos 
                    e = len(sourceToUse) 
                    if sourceToUse==revComp: 
                        snew= len(sourceToUse)-e-1 
                        enew = len(sourceToUse)-s-1
                        s = snew 
                        e = enew  
                    s = 0 
                    e= overflow 
                    if sourceToUse==revComp: 
                        snew = len(sourceToUse)-e -1 
                        enew = len(sourceToUse)-s -1
                        e= enew 
                        s = snew  
                #increment number of bases generated so far 
            basesGenerated+=readLength 
            #create a buffer zone for the read (to be used in the case a deletion causes the read to become shorter than the specified length). 
            if circularity == False: 
                startPos+=readLength 
                numRemainingBases = max(0,len(sourceToUse) - startPos)
                bufferSize = min(numRemainingBases,readLength) 
                if bufferSize ==0:
                    buffer = ""
                else: 
                    buffer = sourceToUse[startPos:startPos+bufferSize] 
            else:
                startPos+=readLength
                if len(sourceToUse) > startPos+readLength: 
                    buffer = sourceToUse[startPos:startPos+readLength] 
                else:
                    overflow = startPos+readLength - len(sourceToUse)
                    buffer = sourceToUse[startPos:] + sourceToUse[0:overflow] 
            spikedReads.append(tuple([read+buffer,readLength]))
            
            #update actual position count statistics. 
            for i in range(readLength): 
                if actualPosCount.__contains__(i): 
                    actualPosCount[i]+=1 
                else: 
                    actualPosCount[i] = 1  
    del readSource # we no longer need this variable. 
    del revcompSource
    del sourceToUse
    processedReads = [] 
    c= 0 
    
    #keep track of actual position/size for insertions/deletions/mutations to ensure that the proper distribution/ parameters are being used to generate these values. 
    actualDelPos = dict() 
    actualInsertPos = dict() 
    actualMutPos = dict() 
    
    actualDelSize = dict() 
    actualInsertSize = dict()
    insertposcheck = [] 
    delposcheck = [] 
    print "Processing reads\n"
    for entry in spikedReads:    
        c+=1 
        if c%10000 ==0:
            print "processing read: " + str(c)+ " out of: " + str(len(spikedReads))+"\n"        
        read = entry[0] 
        if len(read)==0: 
            continue; 
        desiredLength = entry[1] 
        #generating deletions
        
        numdel = int(round(random.gauss(params.meanDeletionsPerRead,params.stdDevDeletionsPerRead))) 
        if numdel < 0: 
            numdel = 0 
        for d in range(numdel): 
            #size of deletion 
            isRepeat = False 
            prob=round(numpy.random.uniform(),5); 
            delsize=cdfToLength(prob,params.delSizeCDF)
            #is the deletion removing a repeat? 
            pdelRepeat= params.delSize[delsize][1]/float(params.delSize[delsize][0]) 
            if numpy.random.uniform()<pdelRepeat:
                isRepeat = True 
                #find a repeat to delete, if a repeat of the desired size is not available, delete a smaller repeat.  
                delPos = findRepeat(read,delsize,desiredLength)         
                while delPos == -1:
                    delsize-=1  
                    delPos = findRepeat(read,delsize,desiredLength)
            else: 
                #print "deleting a non-repeat\n"
                #calculate position of deletion based on probability of deletion at a given position 
                prob = 1 
                delKeys = params.delPosCDFreversed.keys() 
                delKeys.sort() 
                for key in delKeys: 
                    if key > desiredLength:
                        maxprob = params.delPosCDFreversed[key]
                        break 
                while prob > maxprob :
                    prob = round(numpy.random.uniform(),3) 
                delPos = cdfToLength(prob,params.delPosCDF)
            read = read[0:delPos]+read[delPos+delsize:]
            delposcheck.append(delPos)
            #update deletion statistics 
            if actualDelPos.__contains__(delPos): 
                actualDelPos[delPos]+=1 
            else: 
                actualDelPos[delPos]=1 
            if actualPosCount.__contains__(delPos)==False: 
                raise SyntaxError("actualPosCount doesn't have position: " + str(delPos)+'\n')
            if actualDelSize.__contains__(delsize): 
                actualDelSize[delsize][0]+=1 
            else: 
                actualDelSize[delsize] = [1,0]
            if isRepeat: 
                actualDelSize[delsize][1] +=1
            
        if len(read) ==0:
            print "deletion step generated empty read!\n"
        
        #generate insertions 
        numinsert = int(round(random.gauss(params.meanInsertionsPerRead,params.stdDevInsertionsPerRead))) 
        for i in range(numinsert): 
            #size of insertion
            prob=round(numpy.random.uniform(),5); 
            insertsize = cdfToLength(prob,params.insertSizeCDF) 
            prob=1
            
            insertKeys = params.insertPosCDFreversed.keys() 
            insertKeys.sort() 
            maxprob=1; 
            for key_index in range(len(insertKeys)): 
                if insertKeys[key_index] > desiredLength:
                    maxprob = params.insertPosCDFreversed[insertKeys[key_index-1]]
                    break 
            while prob > maxprob :
                prob = round(numpy.random.uniform(),3) 
            pos = cdfToLength(prob,params.insertPosCDF) 
            #is the insertion a repeat?
            isRepeat = False  
            pinsertRepeat = params.insertSize[insertsize][1]/float(params.insertSize[insertsize][0]) 
            if numpy.random.uniform() <pinsertRepeat:
                toInsert = read[pos-insertsize:pos] 
                isRepeat= True 
            else: 
                #insert a random sequence 
                toInsert = "" 
                for i in range(insertsize):
                    baseProb = random.random()  
                    if baseProb < .25: 
                        toInsert+="a" 
                    elif baseProb < 0.50: 
                        toInsert+="t" 
                    elif baseProb < 0.75: 
                        toInsert+="c"
                    else: 
                        toInsert+="g" 
            for i in range(len(read),len(read)+len(toInsert)): 
                if actualPosCount.__contains__(i)==False: 
                    actualPosCount[i]=1 
                else: 
                    actualPosCount[i]+=1 
            
            read = read[0:pos]+toInsert+read[pos:]
            insertposcheck.append(pos)
            #update insertion statistics 
            if actualInsertPos.__contains__(pos): 
                actualInsertPos[pos]+=1 
            else: 
                actualInsertPos[pos] =1 
            if actualInsertSize.__contains__(insertsize): 
                actualInsertSize[insertsize][0]+=1 
            else:
                actualInsertSize[insertsize] = [1,0] 
            if isRepeat: 
                actualInsertSize[insertsize][1]+=1 
           
        #adjust read length to the desired value. 
        if len(read) > desiredLength: 
            toRemove = len(read) - desiredLength 
            read = read[0:-1*toRemove] 
        
        #generate mutation: 
        for pos in range(0,len(read)):
            if pos not in params.pos: 
                probMut=0; 
            elif pos not in params.mutationCount: 
                probMut=0; 
            else:
                probMut=params.mutationCount[pos]/float(params.pos[pos]); 
            if params.mutationCount.keys().__contains__(pos) and numpy.random.uniform()<probMut:
                #mutate 
                start_base = read[pos] 
                if start_base.lower() not in ['a','t','c','g','n']: 
                    start_base = 'n' 
                mutOptions = params.mutationType[start_base.lower()].keys()
                mutProb = random.random() 
                probSum = 0 
                for m in mutOptions: 
                    probSum +=params.mutationType[start_base.lower()][m] 
                    if probSum >= mutProb: 
                        if pos==0 and len(read)>1:
                            read = m+read[1::] 
                        elif pos==0:
                            read = m 
                        else:
                            read = read[0:pos-1]+m+read[pos:]  
                        break
                
                #update mutation statistics 
                if actualMutPos.__contains__(pos): 
                    actualMutPos[pos]+=1 
                else: 
                    actualMutPos[pos]=1 
                
        if len(read) != desiredLength:
            print "desired length of read: " + str(desiredLength)+"\n"
            print "actual length of read: "+ str(len(read)) + "\n"      
            print "read: "+ str(read)+"\n"
            #assert len(read) == desiredLength 
        processedReads.append(read) 
        
    del spikedReads 

    #plot statistics for readLength, insertions, deletions, and mutations 
    #get actual read length distribution 
    actualReadLengths=dict()
    for r in processedReads: 
        readLength=len(r) 
        if readLength not in actualReadLengths: 
            actualReadLengths[readLength]=1; 
        else: 
            actualReadLengths[readLength]+=1; 
    if plotCollection:
        plotCollection.actualReadLengths=actualReadLengths;         
        plotCollection.actualPosCount=actualPosCount; 
        plotCollection.actualDelPos=actualDelPos; 
        plotCollection.actualMutPos=actualMutPos; 
        plotCollection.actualInsertPos=actualInsertPos; 
        plotCollection.actualInsertSize=actualInsertSize; 
        plotCollection.actualDelSize=actualDelSize; 

    print "reading in background read block\n"                   
    #read in the portion of the background reads to be spiked. 
    backgroundReadsf = open(background_reads,"r") 
    backgroundSize = os.stat(background_reads).st_size 
    maxToRead = 104857600 #100 MB chunk 
    numIterations = math.ceil(float(backgroundSize)/maxToRead) 
    numToSpike = math.ceil(float(len(processedReads))/numIterations)
     
    allSpikedFastq = [] 
    allInsertedQual=[] 
    allInsertedFasta = [] 
    allInsertedNames=[] 
    
    for iter in range(int(numIterations)): 
        print "block number: " + str(iter)+"\n"
        backgroundReads= backgroundReadsf.read(maxToRead) 
        startIndex = backgroundReads.index('@') 
        if iter < (int(numIterations)-1):
            endIndex = -1*backgroundReads[::-1].index('@')-1  
        else: 
            endIndex = len(backgroundReads)
        backgroundReads = backgroundReads[startIndex:endIndex] 
        backgroundReads = backgroundReads.split("\n") 
        
        #get the fastq reads from this file 
        print "starting fastq extract step\n"
        qualwithname,qualDict = fastqextract(backgroundReads)
        print "completedfastq extract step\n"
                
        #pick a subset of the reads ( randomly) to spike into this portion of the file. 
        subset_processed=dict()
        longest_allowed = max(qualDict.keys()) 
        print "made length figure\n"
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
        print "starting insertSpikedreads\n"
        spikedFastq,insertedQual,insertedFasta,insertedReadNames= insertSpikedReads(backgroundReads,subset_processed,qualwithname,qualDict)
        print "complete insertSpikedReads\n"
        allSpikedFastq=allSpikedFastq + spikedFastq 
        allInsertedQual = allInsertedQual+insertedQual 
        allInsertedFasta = allInsertedFasta+insertedFasta 
        allInsertedNames = allInsertedNames+insertedReadNames
        
    backgroundReadsf.close()      
    print "writing spiked reads to separate file\n"
    spikedOutputFastaf = open(outputfilename+"_Spiked.fq","w") 
    c=0
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
        c+=1  
        spikedOutputFastaf.write(p+"\n") 
        spikedOutputFastaf.write("+\n") 
        spikedOutputFastaf.write(q+'\n')
    spikedOutputFastaf.close() 

    
    print "writing spiked output file"
    fullOutputFastaf = open(outputfilename+"_Full.fq","w")
    for f in allSpikedFastq: 
        fullOutputFastaf.write(f+"\n") 
    fullOutputFastaf.close()    

    if plotCollection:
    #generate summary graphs for background distributions, regression curves, and spiked reads 
        plotCollection.plotAll(); 
if __name__ == '__main__':
    main() 
