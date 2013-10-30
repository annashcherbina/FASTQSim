# helper functions for spiking 
import math 
import random 
import numpy
import string 

#parses contents of MD tag in a SAM file
#returns list of reference bases for SNP positions 
def parseMD(md): 
    indeletion=False
    numbers=['0','1','2','3','4','5','6','7','8','9']; 
    ref_snp_bases=[] 
    for i in range(len(md)): 
        if md[i] in numbers: 
            indeletion=False
        elif md[i]=='^':
            indeletion=True 
        elif indeletion==False:
            ref_snp_bases.append(md[i])
    return ref_snp_bases 

#finds the longest repeating substring in the passed argument. 
def getLongestRepeat(strand): 
    maxpossible=len(strand)/2
    if maxpossible==0: 
        return 0 
    for l in range(maxpossible,0,-1): 
        for startpos in range(len(strand)-2*l+1): 
            firstblock=strand[startpos:startpos+l]
            secondblock=strand[startpos+l:startpos+2*l]
            if firstblock==secondblock:
                return len(firstblock) 
    return 0 

#parses the CIGAR string in a SAM file 
def parseCIGAR(strand):
    operations=['m','i','d','n','s','h','p','=','x']
    strand=strand.lower() 
    cur=''
    blocks=[]
    for i in range(len(strand)): 
        if strand[i] in operations: 
            blocks.append(cur+strand[i])
            cur=''
        else:
            cur=cur+strand[i] 
    return blocks 

#identifies SNP positions and sequence bases, returns a list of tuples [(SNP position, sequence base)]
def identifySNPs(cigarlist,strand): 
    #we only care about operations "i","m","=","x" 
    filtered_cigar=[] 
    for i in cigarlist: 
        if i.__contains__('s') or i.__contains__('h') or i.__contains__('d') or i.__contains__('n') or i.__contains__('p'):
            continue
        else: 
            filtered_cigar.append(i) 
    snplist=[] 
    curpos=0
    for block in filtered_cigar: 
        if 'i' in block:
            op='i'
            length=int(block.replace('i',''))
        elif 'm' in block: 
            op='m'
            length=int(block.replace('m',''))
        elif 'x' in block:
            op='x'
            length=int(block.replace('x',''))
        elif '=' in block: 
            op='='
            length=int(block.replace('=',''))
        else:
            continue #unknown operation. 
        #we care about the M case 
        if (op=='m') or (op=='x'): 
            for basepos in range(curpos,curpos+length): 
                if strand[basepos]!="=":
                    mutationval=strand[basepos]
                    snplist.append(tuple([basepos,mutationval]))
        curpos+=length; 
    return snplist



#removes the soft/hard clipping to get the position counts for a sequence 
def getPosCount(cigarlist,strandlength):
    hstart=None
    hend=None
    sstart=None
    send=None
    doneWithStart=False
    for b in cigarlist:
        if (b.__contains__('h')) and (hstart==None) and (doneWithStart==False):
            hstart=int(b.replace('h',''))
        elif b.__contains__('h') :
            hend=int(b.replace('h',''))
        elif b.__contains__('s') and (sstart==None) and (doneWithStart==False): 
            sstart=int(b.replace('s',''))
        elif b.__contains__('s'):
            send=int(b.replace('s',''))
        else:
            doneWithStart=True; 
    subtract_from_front=0
    if hstart:
        subtract_from_front+=hstart
    if sstart:
        subtract_from_front+=sstart
    subtract_from_end=0
    if hend: 
        subtract_from_end+=hend 
    if send:
        subtract_from_end+=send 
    startpos=subtract_from_front
    endpos=strandlength-subtract_from_end-1
    return startpos,endpos 



#find the reverse complement of the input
def getRevComp(source):
    rev = source[::-1] 
    revComp = "" 
    for base in rev: 
        if base.lower() =="a":
            revComp+="t" 
        elif base.lower() =="t":
            revComp+="a" 
        elif base.lower() =="c":
            revComp+="g"
        elif base.lower() =="g":
            revComp+="c"
        else:
            revComp+=base 
    return revComp

#pick the first key that exceeds the probability input "prob" 
def cdfToLength(prob,valDict):
    keys = valDict.keys() 
    keyDiff = [abs(prob-key) for key in keys]
    posKeyDiff=[]; 
    for k in keyDiff:
        if k > 0: 
            posKeyDiff.append(k); 
    closestkey = posKeyDiff.index(min(posKeyDiff)) 
    return valDict[keys[closestkey]]  

#find a repeat unit of a specified size within a given read.
#return -1 if a repeat of the desired length is not present in the read. 
def findRepeat(read,delSize,desiredLength): 
    for startPos in range(desiredLength): 
        repeat =read[startPos:startPos+delSize] 
        tester = read[startPos+delSize:startPos+2*delSize] 
        if repeat.lower() == tester.lower(): 
            return startPos 
    return -1 

def mean(values):
    return sum(values)/len(values) 

#calculate standard deviation of a list of values with a known mean 
def SD(values, mean):
    size = len(values)
    sum = 0.0
    for n in range(0, size):
        sum+=((values[n])-mean)**2 
    return math.sqrt((1.0/(size-1))*(sum))

#split fasta file into slices, each having "readsPerSlice" number of reads  
def chopFasta(src,readsPerSlice,origin):
    slices= dict() #tag each slice with a number (these are the keys). The values are the lists of reads associated with each block.
    i =0 
    block = 0 
    numReads = len(src) 
    while i < len(src): 
        if origin == 'b': 
            slices[block]=src[i:min(numReads,i+2*readsPerSlice)] #each read has a header of the form "> ...", so we use 2*readsPerSlice instead of readsPerSlice  
            i+=2*readsPerSlice 
        else: 
            assert(origin == 'p') 
            slices[block] = src[i:min(numReads,i+readsPerSlice)] 
            i+=readsPerSlice 
        
        block+=1 
    return slices
    
def generateReadName(prior_read_name,following_read_name, readlength):
    
    letters = list(string.letters)
    numbers = list(string.digits)
    new_read_name = "" 
    prior_read_name= list(prior_read_name)
    following_read_name = list(following_read_name)
    length_to_use = min(len(prior_read_name),len(following_read_name)) 
    prior_read_name = prior_read_name[0:length_to_use] 
    following_read_name= prior_read_name[0:length_to_use] 
    for pos in range(length_to_use): 
        if prior_read_name[pos] == following_read_name[pos]: 
            new_read_name= new_read_name + prior_read_name[pos]
        elif prior_read_name[pos] in numbers and following_read_name[pos] in numbers: 
            new_read_name = new_read_name + str(int(.5*(int(prior_read_name[pos])+int(following_read_name[pos]))))
        else: 
            decider = random.random() 
            if decider < 0.5: 
                new_read_name = new_read_name + prior_read_name[pos]
            else: 
                new_read_name = new_read_name + following_read_name[pos] 
    #new_read_name = prior_read_name
    return new_read_name 
            
    
#accepts a fastq file separated into lines, these are passed in as a list 
def fastqextract(fastqreads):
    qualDict= dict() 
    qualname=[] 
    seqName = None
    qual=False
    curseq=""
    curqual=""
    first=True
    for line in fastqreads:
        if line.startswith('@'): #mke sure it is not a qual field that starts with @
            if curseq!="" and curqual!="" and seqName!=None:
                qualname.append(seqName) 
                if qualDict.__contains__(len(curqual)): 
                    qualDict[len(curqual)].append(curqual)
                else:
                    qualDict[len(curqual)]=[curqual]
                
                curseq=""
                curqual=""
                seqName = line.rstrip()
                qual=False 
            elif first==True: 
                seqName = line.rstrip() 
                first=False 
            else: 
                curqual=curqual+line.rstrip()             
        elif line.startswith('+') and len(line.rstrip())==1:
            qual=True 
        elif qual==False:
            curseq=curseq+line.rstrip() 
        else:
            curqual=curqual+line.rstrip() 
    qualname.append(seqName)
    if qualDict.__contains__(len(curqual)):
        qualDict[len(curqual)].append(curqual)
    else: 
        qualDict[len(curqual)]=[curqual]
    #interpolate to fill in any "missing" lengths
    maxlength = max(qualDict.keys())
    minlength = min(qualDict.keys())
    lastdefined = qualDict[maxlength] 
    for i in range(maxlength,minlength-1,-1): 
        if qualDict.__contains__(i) == False: 
            qualDict[i] = lastdefined[0:i] 
        else: 
            lastdefined=qualDict[i] 
    return qualname,qualDict
    
    
    
def insertSpikedReads(qualSource,insertedList,qualname,qualSourceDict):
    
    insertedQuality = [] 
    insertedFinalized = [] 
    readNames = [] 
    randomupper = len(qualname)-2 
    counter=0;
    #quality scores are generated for the spiked reads by choosing a set of quality scores for the base position in the background reads. Each quality score for the base position in the backround is used in a rotation format -- all scores are used once before re-using the same score again. 
    for readlength in insertedList.keys(): 
        if readlength==0:
            continue;
        counter+=1 
        if (counter %len(insertedList.keys()))==1000:
            print "inserted "+str(counter) + " reads out of "+str(len(insertList.keys()))+"; "+str(counter*100.0/len(insertList.keys())) +" % done\n" 
        if readlength not in qualSourceDict:
            continue; 
        numoptions = len(qualSourceDict[readlength])-1 
        if numoptions==0:
            continue;
        positions = [int(randomupper*random.random()) for i in xrange(len(insertedList[readlength]))] 
        for i in range(len(insertedList[readlength])): 
            read = insertedList[readlength][i]
            scores = qualSourceDict[readlength][i%numoptions] 
            insert_pos = positions[i] 
            
            #generate a read name for the new pathogen sequence.  
            prior_read_name = qualname[insert_pos]
            if (insert_pos+1) >= len(qualname): 
                following_read_name = qualname[insert_pos-1]   
            following_read_name =qualname[insert_pos+1]  
            readName = generateReadName(prior_read_name,following_read_name,len(read))
            insert_pos = 4*insert_pos 

            qualSource.insert(insert_pos,readName)
            qualSource.insert(insert_pos+1,read)  
            qualSource.insert(insert_pos+2,"+")
            qualSource.insert(insert_pos+3,scores)
            insertedQuality.append(scores)  
            insertedFinalized.append(read) 
            readNames.append(readName)
        #avoid re-using the same read for future quality scores/ titles 
    return qualSource,insertedQuality,insertedFinalized,readNames
                 
            
    
