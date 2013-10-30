#This file does three things: 
#1. calculate the frequency of mutations, insertions, and deletions at each position in an
#NGS read.
#2. find the degree of coverage of the reference genome
#3. find the fraction of each subject read that was aligned to a corresponding sequence
#in the refenece genome. 
# @author        Anna Shcherbina (mailto: anna.shcherbina@ll.mit.edu)
#License:          GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
#
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.


#takes a sam file as input and produces characterization csv files as an output
from helpers import *; 
from pandas import DataFrame;  
import sys; 
inputf=open(sys.argv[1],'r') 
outputprefix=sys.argv[2]
plothistogram=False 
if len(sys.argv) > 3: 
    if sys.argv[3]=="-plothistogram": 
        plothistogram=True 
        import matplotlib.pyplot as plt 

posCount=dict()
delCount=dict()
insertCount=dict()
mutCount=dict()
mutType=dict() 
insertSize=dict()
delSize=dict()
Usage=dict()
insertsByRead=[] 
delsByRead=[] 

for line in inputf:
    if line.startswith('@'): 
        continue 
    line=line.strip(); 
    line=line.split('\t')
    #print str(line)+'\n' 
    seqname=line[0] 
    aligned=int(line[1]) 
   # print str(aligned)+'\n'
    if aligned !=0:
        continue 
    refname=line[2]
    refstartpos=int(line[3]) 
    cigar=line[5]
    strand=line[9] 
    strandlength=len(strand) 
    md=""
    #account for single base pair mutations 
    for i in range(11,len(line)): 
        if 'MD' in line[i]: 
            md=line[i].replace('MD','') 
            md=md.replace('Z','') 
            md=md.replace(':','') 
            break 
    ref_snp_bases=parseMD(md) 

    

    cigarlist=parseCIGAR(cigar)
    
    #posCount
    startpos,endpos=getPosCount(cigarlist,strandlength)
        
    Usage[seqname]=[strandlength,endpos-startpos+1,startpos,endpos]
    strandofinterest=strand[startpos:endpos+1]
    for p in range(startpos,endpos+1):
        if p not in posCount:
            posCount[p]=1
        else: 
            posCount[p]+=1
    
    alignmentpos=0; 
    curpos=startpos;
    insertcountforread=0; 
    delcountforread=0;
    #get Insertions/Deletions 
    for block in cigarlist:
        if 'h' in block:
            continue
        elif 'p' in block:
            continue 
        elif 's' in block:
            continue 
        elif ('m' in block) or ('=' in block) or ('x' in block):
            #match or single-base mutation 
            numbases=int(block.replace('m',''))
            alignmentpos+=numbases 
            curpos+=numbases 
        elif 'i' in block: 
            insertcountforread+=1
            insertionsize=int(block.replace('i',''))
            insertionbases=strand[curpos+1:curpos+2+insertionsize]
            for p in range(curpos+1,curpos+2+insertionsize): 
                if p not in insertCount: 
                    insertCount[p]=1
                else: 
                    insertCount[p]+=1
                    
                    
            longestRepeat=getLongestRepeat(insertionbases)
            
            if insertionsize not in insertSize: 
                insertSize[insertionsize]=[1,0] 
            else:
                insertSize[insertionsize][0]+=1
            if (longestRepeat > 0) and (longestRepeat not in insertSize):
                insertSize[longestRepeat]=[0,1]
            elif longestRepeat > 0:
                insertSize[longestRepeat][1]+=1 
            
            curpos+=insertionsize
            alignmentpos+=insertionsize 
        elif 'd' in block:
            delcountforread+=1
            deletionsize=int(block.replace('d',''))
            deletionbases=strand[curpos+1:curpos+2+deletionsize]
            for p in range(curpos+1,curpos+2+deletionsize): 
                if p not in delCount:
                    delCount[p]=1
                else:
                    delCount[p]+=1 
            longestRepeat=getLongestRepeat(deletionbases) 
            if deletionsize not in delSize: 
                delSize[deletionsize]=[1,0] 
            else:
                delSize[deletionsize][0]+=1 
            if (longestRepeat > 0) and (longestRepeat not in delSize): 
                delSize[longestRepeat]=[0,1]
            elif longestRepeat > 0:
                delSize[longestRepeat][1]+=1 
            
            curpos+=deletionsize 
            alignmentpos+=deletionsize
        else: 
            print "unknown Op in CIGAR:"+str(block)
    #print "strandofinterest:"+str(strandofinterest)+'\n'
    #Handle single point mutations from MD tag 
    insertsByRead.append(insertcountforread)
    delsByRead.append(delcountforread)

    mutation_pos=identifySNPs(cigarlist,strandofinterest); 
    #offset each mutation position by the start of the alignment 
    if len(ref_snp_bases)!=len(mutation_pos): 
        print "MD does not agree with MN for strand:\n"
        print str(strand)+'\n' 
        print "using conservative estimate for SNPs\n"
    for i in range(min(len(mutation_pos),len(ref_snp_bases))):  
        ref_base=(ref_snp_bases[i]).lower() 
        snp_pos=int(mutation_pos[i][0])+startpos
        seq_base=(mutation_pos[i][1]).lower()
        if snp_pos not in mutCount: 
            mutCount[snp_pos]=1 
        else:
            mutCount[snp_pos]+=1
        if ref_base not in mutType: 
            mutType[ref_base]=dict() 
        if seq_base not in mutType[ref_base]: 
            mutType[ref_base][seq_base]=1
        else: 
            mutType[ref_base][seq_base]+=1 

#generate the output files 
fout=open(outputprefix+"posCount.csv",'w') 
for entry in posCount: 
    fout.write(str(entry)+','+str(posCount[entry])+'\n') 

fout=open(outputprefix+"delCount.csv",'w') 
for entry in delCount: 
    fout.write(str(entry)+','+str(delCount[entry])+'\n') 

fout=open(outputprefix+"insertCount.csv",'w') 
for entry in insertCount: 
    fout.write(str(entry)+','+str(insertCount[entry])+'\n') 

fout=open(outputprefix+"mutationCount.csv","w") 
for entry in mutCount: 
    fout.write(str(entry)+','+str(mutCount[entry])+'\n') 

fout=open(outputprefix+"mutationType.csv","w") 
for entry in mutType:
    fout.write(entry) 
    for subentry in mutType[entry]: 
        fout.write(','+subentry+','+str(mutType[entry][subentry]))
    fout.write('\n') 

fout=open(outputprefix+"insertsByRead.csv","w") 
if len(insertsByRead)>0:
    fout.write(str(insertsByRead[0]))

for i in range(1,len(insertsByRead)): 
    fout.write(','+str(insertsByRead[i]))

fout=open(outputprefix+"delsByRead.csv","w") 
if len(delsByRead) > 0:
    fout.write(str(delsByRead[0]))
for i in range(1,len(delsByRead)):
    fout.write(','+str(delsByRead[i]))


fout=open(outputprefix+"insertSize.csv","w") 
for entry in insertSize: 
    fout.write(str(entry)+","+str(insertSize[entry][0])+','+str(insertSize[entry][1])+'\n') 

fout=open(outputprefix+"delSize.csv","w") 
for entry in delSize:
    fout.write(str(entry)+','+str(delSize[entry][0])+","+str(delSize[entry][1])+'\n') 

fout=open(outputprefix+"Usage.csv","w") 
for entry in Usage: 
    fout.write(str(entry)+","+str(Usage[entry][0])+','+str(Usage[entry][1])+','+str(Usage[entry][2])+','+str(Usage[entry][3])+'\n')


#write a summary file 
fout=open(outputprefix+"summary.csv",'w')
fout.write(outputprefix+"posCount.csv\n")
fout.write(outputprefix+"delCount.csv\n")
fout.write(outputprefix+"insertCount.csv\n")
fout.write(outputprefix+"mutationCount.csv\n") 
fout.write(outputprefix+"mutationType.csv\n")
fout.write(outputprefix+"insertSize.csv\n")
fout.write(outputprefix+"delSize.csv\n") 
fout.write(outputprefix+"Usage.csv\n") 
fout.write(outputprefix+"readHist.csv\n") 
fout.write(outputprefix+"qualHist.csv\n") 
fout.write(outputprefix+"insertsByRead.csv\n") 
fout.write(outputprefix+"delsByRead.csv\n") 

#if plotting is enabled, generate plots of the characterization statistics 
if plothistogram: 

    #insertion probability from insertion count. 
    insertprob=dict() 
    for entry in posCount: 
        if entry in insertCount: 
            insertprob[entry]=float(insertCount[entry])/float(posCount[entry])
    plotname=outputprefix+"CharacterizationInsertCount.png" 
    fig=plt.figure(); 
    ax=fig.add_subplot(111) 
    ax.bar(insertprob.keys(), insertprob.values(),width=1,color='r',log=True) 
    ax.set_ylim([10e-5,1])
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Base Position Along a Read',fontsize=20) 
    ax.set_ylabel('Probability of Insertion',fontsize=20) 
    ax.set_title('Probability of Insertion as a Function of Base Position\n Dataset '+ str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True) 
    plt.savefig(plotname,bbox_inches=0) 
        
    #plot the insertion size distribution 
    plotname=outputprefix+"CharacterizationInsertSize.png"
    fig=plt.figure() 
    ax=fig.add_subplot(111)
    barwidth=0.2
    overallSize=[i - barwidth for i in insertSize.keys()] 
    overallCount=[i[0] for i in insertSize.values()] 
    repeatCount=[i[1] for i in insertSize.values()] 
    bar1=ax.bar(overallSize,overallCount,width=barwidth,color='b',align='center',log=True,label='Total Insertions') 
    bar2=ax.bar(insertSize.keys(),repeatCount,width=barwidth,color='r',align='center',log=True,label='Repeat Insertions')        
    ax.legend(loc=1,borderaxespad=0)
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Insert Size',fontsize=20) 
    ax.set_ylabel('Insert Count',fontsize=20) 
    ax.set_title('Insertion Size \n Dataset '+sys.argv[1].split('/')[-1],fontsize=20)
    plt.grid(True) 
    plt.savefig(plotname,bbox_inches=0) 


    #deletion probability from deletion count. 
    delprob=dict() 
    for entry in posCount: 
        if entry in delCount: 
            delprob[entry]=float(delCount[entry])/float(posCount[entry])
    plotname=outputprefix+"CharacterizationDelCount.png"
    fig = plt.figure() 
    ax=fig.add_subplot(111) 
    ax.bar(delprob.keys(), delprob.values(),color='r',log=True) 
    ax.set_ylim([10e-5,1])
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Base Position Along a Read',fontsize=20) 
    ax.set_ylabel('Probability of Deletion',fontsize=20) 
    ax.set_title('Probability of Deletion as a Function of Base Position\n Dataset '+ str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True) 
    plt.savefig(plotname,bbox_inches=0) 

    #plot the deletion size distribution 
    plotname=outputprefix+"CharacterizationDelSize.png"
    fig=plt.figure() 
    ax=fig.add_subplot(111) 
    ax.set_yscale('log') 
    overallSize=[i - 0.1 for i in delSize.keys()] 
    overallCount=[i[0] for i in delSize.values()] 
    repeatCount=[i[1] for i in delSize.values()] 
    bar1=ax.bar(overallSize,overallCount,width=0.2,color='b',align='center',label='Total Deletions',log=True) 
    bar2=ax.bar(delSize.keys(),repeatCount,width=0.2,color='r',align='center',label='Repeat Deletions',log=True)      
    ax.legend(loc=1,borderaxespad=0)
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Deletion Size',fontsize=20) 
    ax.set_ylabel('Deletion Count',fontsize=20) 
    ax.set_title('Deletion Size \n Dataset'+str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True) 
    plt.savefig(plotname,bbox_inches=0) 


    #plot the probability from mutation count 
    mutprob=dict() 
    for entry in posCount: 
        if entry in mutCount: 
            mutprob[entry]=float(mutCount[entry])/posCount[entry]
        plotname=outputprefix+"CharacterizationMutationCount.png" 
    fig = plt.figure() 
    ax=fig.add_subplot(111) 
    ax.bar(mutprob.keys(), mutprob.values(),color='r',log=True) 
    ax.set_ylim([10e-5,1])
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Base Position Along a Read',fontsize=20) 
    ax.set_ylabel('Probability of Mutation',fontsize=20) 
    ax.set_title('Probability of Mutation as a Function of Base Position\n Dataset '+ str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True) 
    plt.savefig(plotname,bbox_inches=0) 

    #plot mutation type histogram 
    mtp=dict() 
    bases=['a','t','c','g','n'] 
    for b1 in bases:
        mtp[b1]=dict() 
        for b2 in bases: 
            if b1==b2:
                continue
            mtp[b1][b2]=0 
    
    for key in mutType:
        totalmuts=float(sum(mutType[key].values()))
        for subkey in mutType[key]: 
            mtp[key][subkey]=mutType[key][subkey]/max(0.01,totalmuts) 
    mutationDF=DataFrame([[0,mtp['a']['t'],mtp['a']['c'],mtp['a']['g'],mtp['a']['n']],[mtp['t']['a'],0,mtp['t']['c'],mtp['t']['g'],mtp['t']['n']],[mtp['c']['a'],mtp['c']['t'],0,mtp['c']['g'],mtp['c']['n']],[mtp['g']['a'],mtp['g']['t'],mtp['g']['c'],0,mtp['g']['n']],[mtp['n']['a'],mtp['n']['t'],mtp['n']['c'],mtp['n']['g'],0]],columns=['A','T','C','G','N'])
    plotname=outputprefix+'CharacterizationMutType.png'
    mutplot=mutationDF.plot(kind='bar',stacked=True)
    group_labels=['A','T','C','G','N']
    mutplot.set_xticklabels(group_labels) 
    mutplot.set_ylim([0,1])
    plt.setp(mutplot.get_xticklabels(),fontsize=18)
    plt.setp(mutplot.get_yticklabels(),fontsize=18)        
    mutplot.set_xlabel('Original Base',fontsize=20) 
    mutplot.set_ylabel('Mutated Base',fontsize=20) 
    mutplot.set_title('Probability of Mutation by Base\n Dataset '+str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True)        
    plt.savefig(plotname,bbox_inches=0) 

    #plot read usage 
    readUsageDist=[]
    for entry in Usage: 
        totalReadLength=Usage[entry][0] 
        usedReadLength=Usage[entry][1] 
        fractionUsed=float(usedReadLength)/float(totalReadLength) 
        readUsageDist.append(fractionUsed) 
        plotname=outputprefix+'CharacterizationUsage.png'
    fig=plt.figure() 
    ax = fig.add_subplot(111) 
    ax.set_yscale('log') 
    ax.set_xlim([0,1])
    ax.hist(readUsageDist,100,histtype='bar',color='r') 
    plt.setp(ax.get_xticklabels(),fontsize=18)
    plt.setp(ax.get_yticklabels(),fontsize=18)        
    ax.set_xlabel('Fraction of Aligned Bases in Read',fontsize=20)
    ax.set_ylabel('Number of Reads',fontsize=20)
    ax.set_title('Read Alignment Quality \n Dataset '+str(sys.argv[1].split('/')[-1]),fontsize=20)
    plt.grid(True)
    plt.savefig(plotname,bbox_inches=0) 
