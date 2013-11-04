'''
Created on Jul 28, 2012
%stores parameters used for spiking the dataset
@author: AN22471
'''
from helpers import * 
from helpers import *
from curve_fit_params import * 
from time import sleep 
import thread 
    

class Params(object):
    #parameter dictionaries
    readLengths = dict() 
    quality = dict()
    pos = dict() 
    insertCount = dict() 
    insertSize = dict() 
    insertsByRead = [] 
    delCount = dict() 
    delSize = dict() 
    delsByRead = [] 
    mutationCount = dict() 
    mutationType = dict() 
    minLength = 50 #we are not interested in short "noise" reads. No read shorter than 50 bases should be generated. 
    maxLength = None 
    
    #if a generated read falls into these ranges, the empirical distribution data should be used. Outside of this range, we interpolate for missing values. 
    readLengthsRange =[] 
    insertPosRange = [] 
    delPosRange = [] 
    mutationPosRange = [] 
    delSizeRange = [] 
    insertSizeRange = [] 
    delRepeatRange = [] 
    insertRepeatRange=[]
    
    #get mean and standard deviation of number of insertions/deletions for each read in the background data. 
    meanInsertionsPerRead = None 
    meanDeletionsPerRead = None 
    stdDevInsertionsPerRead = None 
    stdDevDeletionsPerRead = None 
    
    #cdf equivalents for parameter dictionaries 
    readLengthsCDF = dict() 
    insertSizeCDF = dict() 
    delSizeCDF = dict() 
    insertPosCDF = dict() 
    delPosCDF = dict() 
    
    plotCollection=None

    #populate the specified dictionary with information from the specified parameter file.
    def populate(self,srcFile,destDict):
        srcFile = srcFile.strip('\n') 
        f = open(srcFile,'r') 
        for line in f: 
            line = line.lstrip(",")
            line = line.replace("\n","") 
            line = line.split(",")
            if type(destDict).__name__=='list': 
                for l in line: 
                    destDict.append(int(l))          
            elif len(line) <2:
                raise SyntaxError("line must have at least two values!: " + line + " in file: " + srcFile +"\n")  
            elif len(line) ==2: 
                destDict[int(line[0])] = float(line[1]) 
            elif srcFile.__contains__("mutationType"): 
                key = line[0] 
                self.mutationType[key]= dict() 
                for i in range(1,len(line),2): 
                    self.mutationType[key][line[i]] = float(line[i+1]) 
            else:                 
                value = [float(i) for i in line] 
                key = int(value[0]) 
                value = value[1::] 
                destDict[key] = value 
        f.close() 
    
    #find the total number of fasta reads in the source file. 
    def srcSize(self):
        if self.readLengths == None: 
            raise SyntaxError("populate the parameters list first! \n") 
        else: 
            return sum(self.readLengths.values())
        
    def srcTotalInsertions(self):
        if len(self.insertSize) ==0: 
            raise SyntaxError("self.insertSize has not been populated yet\n") 
        return sum([i[0] for i in self.insertSize.values()])
    def srcTotalDeletions(self):
        if len(self.delSize) == 0: 
            raise SyntaxError("self.delSize has not been populated yet\n")
        return sum([i[0] for i in self.delSize.values()])  
            
        
    #make a CDF for values in the source dict store this cdf in the destination dict
    def makeCDF(self,srcDict,destDict):
        
        if type(srcDict.values()[0]).__name__ == 'int'or type(srcDict.values()[0]).__name__ == 'float64' or type(srcDict.values()[0]).__name__ == "str" or  type(srcDict.values()[0]).__name__ == "float": 
            total = float(sum(srcDict.values())) 
        elif type(srcDict.values()[0]).__name__ == 'list': 
            total = float(sum([i[0] for i in srcDict.values()]))
        else: 
            raise SyntaxError('cdf source must have list,string,float or int values! \n')
        keys = srcDict.keys() 
        keys.sort() 
        sumSoFar = 0 
        for key in keys: 
            if type(srcDict[key]).__name__ == "list": 
                sumSoFar +=srcDict[key][0]
            else: 
                sumSoFar+=srcDict[key]  
            destDict[round(sumSoFar/total,5)] = key 
    
    def getRangeForEmpiricalFit(self,inputfile):   
        fdata = inputfile.read().split('\n')
        empties = fdata.count('') 
        for i in range(empties): 
            fdata.remove('')
        values = [int(i.split(',')[0]) for i in fdata]
        lowval=values[0] 
        highval=math.floor(values[-1]*.95)
        return lowval,highval

    def __init__(self,param_file,backgroundReads,plotCollection):
        if plotCollection:
            self.plotCollection=plotCollection; 
            plotCollection.plotName=backgroundReads.split('.')[0]
            plotCollection.plotName=plotCollection.plotName.split('/')[-1]
        categoryMap = {"readHist": self.readLengths,"qualHist":self.quality,"posCount":self.pos,"insertCount":self.insertCount,"insertSize":self.insertSize,"insertsByRead":self.insertsByRead,"delCount":self.delCount,"delSize":self.delSize,"delsByRead":self.delsByRead,"mutationCount": self.mutationCount,"mutationType":self.mutationType}   
        limitMap = dict() 
        paramfilelist= open(param_file,'r') 
        for line in paramfilelist:   
            #find the range of read lengths that should be fit empirically (within 2 std dev of the mean). 
            #assumes input parameter .csv files are sorted
            if line.__contains__('readHist'): 
                fin=open(line.rstrip(),'r')  
                self.minLength,self.maxLength = self.getRangeForEmpiricalFit(fin)
                self.readLengthsRange=[self.minLength,self.maxLength]  
            elif line.__contains__('insertCount'): 
                fin=open(line.rstrip(),'r')
                insertPosMin,insertPosMax=self.getRangeForEmpiricalFit(fin)
                self.insertPosRange=[insertPosMin,insertPosMax]
            elif line.__contains__('delCount'):
                fin=open(line.rstrip(),'r')
                delPosMin,delPosMax=self.getRangeForEmpiricalFit(fin)
                self.delPosRange=[delPosMin,delPosMax]
            elif line.__contains__('insertSize'):
                fin=open(line.rstrip(),'r')
                insertSizeMin,insertSizeMax=self.getRangeForEmpiricalFit(fin)
                self.insertSizeRange=[insertSizeMin,insertSizeMax] 
                self.insertRepeatRange=self.insertSizeRange 
            elif line.__contains__('delSize'):
                fin=open(line.rstrip(),'r')
                delSizeMin,delSizeMax = self.getRangeForEmpiricalFit(fin)
                self.delSizeRange=[delSizeMin,delSizeMax]
                self.delRepeatRange=self.delSizeRange 
            elif line.__contains__('mutationCount'): 
                fin=open(line.rstrip(),'r')
                mutPosMin,mutPosMax=self.getRangeForEmpiricalFit(fin)
                self.mutationPosRange=[mutPosMin,mutPosMax] 
                
            for category in categoryMap.keys(): 
                if line.__contains__(category): 
                    self.populate(line,categoryMap[category])  
        paramfilelist.close() 
        
        #fit curves to empirical data to  account for noise effects and very long reads and the "cliffs" observed for PacBio indels/mutations 
        #read lengths
        if plotCollection:
            plotCollection.pos=self.pos
            plotCollection.readLengths=self.readLengths; 
        A= max(self.readLengths.values())
        total = sum(self.readLengths.values()) 
        meanlength = sum([float(self.readLengths.keys()[i]*self.readLengths.values()[i])/total for i in range(len(self.readLengths.keys()))]) 
        stdev = .2*(max(self.readLengths.values())-min(self.readLengths.values()))
        gauss_params=[A,meanlength,stdev]
        self.readLengths = fitToData(self.readLengths,self.readLengthsRange,gauss_params,3) 
        if plotCollection:
            plotCollection.curveFitReadLengths=self.readLengths; 
            plotCollection.insertProb=self.insertCount; 
        self.insertCount = fitToData(self.insertCount,self.insertPosRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitInsertProb=self.insertCount
         
        #insertion size/ insertion repeats
        insertsize_tmp = dict()
        insertrepeat_tmp = dict() 
        for key in self.insertSize:
            insertsize_tmp[key] = self.insertSize[key][0]
            insertrepeat_tmp[key] = self.insertSize[key][1] 
        if plotCollection:
            plotCollection.insertSize=insertsize_tmp; 
            plotCollection.insertRepeat=insertrepeat_tmp; 


        insertsize_tmp = fitToData(insertsize_tmp,self.insertSizeRange,gauss_params,3)
        insertrepeat_tmp=fitToData(insertrepeat_tmp,self.insertSizeRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitInsertSize=insertsize_tmp; 
            plotCollection.curveFitInsertRepeat=insertrepeat_tmp; 

        self.insertSize=dict() 
        for i in insertsize_tmp.keys(): 
            self.insertSize[i] = [insertsize_tmp[i],insertrepeat_tmp[i]] 

        if plotCollection:
            plotCollection.delProb=self.delCount; 
        #deletion probability   
        self.delCount = fitToData(self.delCount,self.delPosRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitDelProb=self.delCount; 
        
        #deletion size/ deletion repeats
        delsize_tmp = dict() 
        delrepeat_tmp=dict() 
        for key in self.delSize: 
            delsize_tmp[key]= self.delSize[key][0]
            delrepeat_tmp[key]=self.delSize[key][1]
        if plotCollection:
            plotCollection.delSize=delsize_tmp; 
            plotCollection.delRepeat=delrepeat_tmp; 

            
        
        delsize_tmp = fitToData(delsize_tmp,self.delSizeRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitDelSize=delsize_tmp; 
        
        delrepeat_tmp = fitToData(delrepeat_tmp,self.delSizeRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitDelRepeat=delrepeat_tmp; 
        self.delSize=dict() 
        for i in delsize_tmp.keys(): 
            self.delSize[i] = [delsize_tmp[i],delrepeat_tmp[i]] 
            
        #mutation probability
        if plotCollection:
            plotCollection.mutationProb=self.mutationCount; 
        
        self.mutationCount = fitToData(self.mutationCount,self.mutationPosRange,gauss_params,3)
        if plotCollection:
            plotCollection.curveFitMutationProb=self.mutationCount; 
                
        
        self.makeCDF(self.readLengths,self.readLengthsCDF) 
        self.makeCDF(self.delSize,self.delSizeCDF) 
        self.makeCDF(self.insertSize,self.insertSizeCDF)
        self.makeCDF(self.insertCount,self.insertPosCDF) 
        self.makeCDF(self.delCount,self.delPosCDF)
        
        #readlength to cdf keys 
        self.delPosCDFreversed =dict() 
        for key in self.delPosCDF:
            self.delPosCDFreversed[self.delPosCDF[key]] = key 
            
        self.insertPosCDFreversed = dict() 
        for key in self.insertPosCDF: 
            self.insertPosCDFreversed[self.insertPosCDF[key]]=key 
                        
        keys = self.insertPosCDF.keys() 
        keys.sort()
        iscdf = [] 
        for k in keys: 
            iscdf.append(self.insertPosCDF[k]) 
        
        keys = self.delPosCDF.keys() 
        keys.sort()
        iscdf = [] 
        for k in keys: 
            iscdf.append(self.delPosCDF[k]) 

        self.meanInsertionsPerRead = sum(self.insertsByRead)/float(len(self.insertsByRead))   
        self.meanDeletionsPerRead= sum(self.delsByRead)/float(len(self.delsByRead)) 
        self.stdDevDeletionsPerRead = SD(self.delsByRead,self.meanDeletionsPerRead) 
        self.stdDevInsertionsPerRead = SD(self.insertsByRead,self.meanInsertionsPerRead) 
