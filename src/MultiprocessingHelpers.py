import numpy; 
from helpers import *; 
from random import *; 
from curve_fit_params import *; 



def insert_spiked_reads_worker(input,output): 
    for args in iter(input.get,'STOP'): 
        backgroundReads=args[0] 
        subset_processed=args[1] 
        qualwithname=args[2] 
        qualDict=args[3] 
        spikedFastq,insertedQual,insertedFasta,insertedReadNames= insertSpikedReads(backgroundReads,subset_processed,qualwithname,qualDict)
        output.put([spikedFastq,insertedQual,insertedFasta,insertedReadNames]) 



def generic_worker(input,output): 
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

def read_mutate_worker(input,output,params): 
    
    for entry in iter(input.get,'STOP'):
        read = entry[0][0]
        desiredLength = entry[0][1]
        plothist=entry[1]
        if plothist:
            statDict=dict()         
        else:
            statDict=None
        
        #generating deletions
        numdel = int(round(random.gauss(params.meanDeletionsPerRead,params.stdDevDeletionsPerRead))) 
        #numdel=max(numdel,len(read))
        for d in range(numdel):
            if len(read)<2:
                break; 
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
            
            #update deletion statistics
            if plothist:
                statDict["delPos"]=delPos; 
                statDict["delSize"]=delsize 
                if isRepeat: 
                    statDict["delRepeat"]=delsize; 
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
            if plothist:
                statDict['posCount']=[]; 
                for i in range(len(read),len(read)+len(toInsert)):
                    statDict['posCount'].append(i) 
            
            read = read[0:pos]+toInsert+read[pos:]
            
            #update insertion statistics
            if plothist:
                statDict['insertPos']=pos; 
                statDict['insertSize']=insertsize 
                if isRepeat: 
                    statDict['insertRepeat']=insertsize 
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
                if plothist:
                    statDict['mutPos']=pos; 
        if len(read) != desiredLength:
            print "desired length of read: " + str(desiredLength)+"\n"
            print "actual length of read: "+ str(len(read)) + "\n"      
            print "read: "+ str(read)+"\n"
        if plothist:
            statDict['readLength']=len(read) 
        output.put([read,statDict])
        
        
        
def read_gen_worker(input,output,params): 
    for args in iter(input.get,'STOP'):
        src=args[0] 
        revComp=args[1] 
        circularity=args[2] 
        coverage=args[3] 
        basesToGenerate=args[4]
        plothist=args[5]
        print "bases to generate:"+str(basesToGenerate)+"\n" 
        if plothist==False:
            actualPosCount=None
        else:
            actualPosCount=dict() 

        basesGenerated=0; 
        spikedReads=[]

        while basesGenerated < basesToGenerate: 
            #choose read length and position to insert the read within the background fastq file.        
            prob = 1 
            readLength = 0 
            prob = round(numpy.random.uniform(),3) 
            readLength = int(cdfToLength(prob, params.readLengthsCDF))
            #choose read position randomly
            startPos = random.randint(0,len(src)-params.minLength)                       
            #determine whether we use the source sequence or its reverse complement. 
            if numpy.random.random() < 0.5:
                sourceToUse = revComp
            else:
                sourceToUse=src
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
            if plothist:
                for i in range(readLength): 
                    if actualPosCount.__contains__(i): 
                        actualPosCount[i]+=1 
                    else: 
                        actualPosCount[i] = 1  
        output.put([actualPosCount,spikedReads]) 
