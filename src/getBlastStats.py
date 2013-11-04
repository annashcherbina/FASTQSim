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


import sys
import io
import os 
from pandas import DataFrame 

#mark columns in blast CSV files
read_name_pos = 0  
query_start = 2
query_end = 3
ref_start = 13
ref_end = 14
query_seq = 18
ref_seq = 19 

#return the reverse complement of a given string 
def getrevcomp(startstring): 
    rev = startstring[::-1] 
    revc = "" 
    for base in rev: 
        if base.lower() == "a": 
            revc+="t" 
        elif  base.lower() == "t": 
            revc+="a" 
        elif base.lower() == "c": 
            revc+="g" 
        elif base.lower() == "g": 
            revc+="c"
        else: 
            revc+=base.lower() 
    return revc 


def main():
    if len(sys.argv)<3:
        raise SyntaxError("usage: python getBlastStats.py <fasta_base_name> <parsed blast base name> \n")
    base_fasta = sys.argv[1]
    base_blast = sys.argv[2]
    plothistogram=False 
    if len(sys.argv) > 3: 
        if sys.argv[3]=="-plothistogram": 
            plothistogram=True 
            import matplotlib.pyplot as plt 
    base_output = sys.argv[1].replace('.fasta','') 


    poscount = dict() # read position -> number of reads this position occurs in 


    insertcount = dict() # read position -> number of insertions relative to reference genome 
    insertprob=dict() 
    insertSize = dict() #bp -> [number of insertions, number of insertions that are repeats]
    insertsByRead = [] 
    
    delcount = dict() # read position -> number of delettions relative to reference genome 
    delprob=dict() 
    delSize = dict() #bp -> [number of deletions, number of deletions that are repeats]
    delsByRead = [] 

    mutationcount = dict() # read position -> number of mutations relative to reference genome
    mutprob=dict() 
    mutationtype= {'a':{'t':0,'c':0,'g':0,'n':0},
                   't':{'a':0,'c':0,'g':0,'n':0},
                   'c':{'a':0,'t':0,'g':0,'n':0},
                   'g':{'a':0,'t':0,'c':0,'n':0},
                   'n':{'a':0,'t':0,'c':0,'g':0}}


    coverage = dict() # reference genome pos -> number of query sequences that overlap this pos 
    readUsage = dict() # readname->[read length, length of aligned part, coverstart, coverend]
    readUsageDist=[]

    #check for start/end primers and their reverse complements. 
    primercheck = dict() 

    fasta = open(base_fasta,'r')
    currentReadLength = 0
    currentRead = None     
    for line in fasta:
        if line.startswith('>'): 
            line = line.replace('>','') 
            line = line.replace('\n','') 
            line = line.replace(' ','') 
            line = line.replace('"','') 
        
            seq_name = str(line)  
            
            if currentRead == None:
                currentRead = seq_name
            else:
                readUsage[str(currentRead)] = [currentReadLength]
                currentReadLength = 0
                currentRead = seq_name            
            #find the read length
        else:
            line = line.strip('\n')
            currentReadLength+=len(line)
    #write the final line 
    readUsage[str(currentRead)] = [currentReadLength] 
    
    fasta.close() 


    blast = open(base_blast,'r')
    for line in blast:
        line = line.split(",")
        readName = str(line[read_name_pos])
        readName = readName.replace('\n','') 
        readName = readName.replace(' ','')
        readName = readName.replace('\t','') 
        readName = readName.replace('"','') 


        query_start_pos = int(line[query_start])
        query_end_pos = int(line[query_end])
        ref_start_pos = int(line[ref_start])
        ref_end_pos = int(line[ref_end])

        query = line[query_seq].strip('0123456789')
        query = query.replace('"','') 
        ref = line[ref_seq].strip('0123456789')
        ref = ref.replace('"','') 

        for i in range(0,10): 
            query = query.replace(str(i),'') 
            ref = ref.replace(str(i),'') 
  
        #sometimes parser truncates read names. account for this possibility. 
        for k in readUsage.keys(): 
            if k.__contains__(readName): 
                readName = k 
                break 

        readUsage[readName].append(query_end_pos - query_start_pos+1) 
        readUsage[readName].append(query_start_pos) 
        readUsage[readName].append(query_end_pos) 
        
        if float(readUsage[readName][1])/float(readUsage[readName][0]) < .80: 
            #don't calculate statistics because the read is  most likely overlapping a transposon/ other junk sequence 
            continue 

        for i in range(query_start_pos,query_end_pos+1):
            if poscount.__contains__(i) == False:
                poscount[i]=1
            else:
                poscount[i]+=1
        #we are not aiming to fully cover the reference genome, so this is commented out. Uncomment it to calculate the number 
        #of times that the reads hit each position in the reference genome. 

        for i in range(ref_start_pos,ref_end_pos+1):
            if coverage.__contains__(i) == False:
                coverage[i] = 1
            else:
                coverage[i]+=1 
        if len(query)!= len(ref): 
            print "query: " + str(query)+"\n" 
            print "ref: " + str(ref) + "\n" 

        shorterLen = min(len(query),len(ref))
        delLen = 0 
        insertLen = 0 
        numDels = 0 
        numInserts = 0
        if query !="":
            q_start = query[0:4].lower()
            q_end = query[-4:].lower()
            q_start_revc = getrevcomp(q_start) 
            q_end_revc = getrevcomp(q_end)
            primercandidates = [q_start,q_end,q_start_revc, q_end_revc] 
            for candidate in primercandidates: 
                if primercheck.__contains__(candidate): 
                    primercheck[candidate]+=1 
                else: 
                    primercheck[candidate] = 1 
        for i in range(shorterLen):
            assert((query[i] =="-" and ref[i] == "-")==False) #there should never be a case with a gap in both sequences at the same base

            if query[i] =="-": #deletion
                delLen+=1 
                if delcount.__contains__(i) == False:
                    delcount[i] =1
                else:
                    delcount[i]+=1
            elif ref[i] =='-': #insertion
                insertLen+=1
                if insertcount.__contains__(i) == False:
                    insertcount[i] = 1
                else:
                    insertcount[i]+=1

            if ref[i].lower()!=query[i].lower():
                if ref[i].lower() not in ['a','t','c','g','n']: 
                    continue 
                if query[i].lower() not in ['a','t','c','g','n']:
                    continue

                if mutationcount.__contains__(i) == False:
                    mutationcount[i] = 1
                else:
                    mutationcount[i]+=1
                if delLen!=0: 
                    if delSize.__contains__(delLen): 
                        delSize[delLen][0]+=1 
                    else: 
                        delSize[delLen] = [1,0]
                    numDels+=1 
                    #is this a repeat deletion?
                    bases_deleted = ref[i-delLen:i] 
                    flank_left = ref[i-2*delLen:i-delLen] 
                    flank_right = ref[i:min(shorterLen,i+delLen)]
                    if flank_left == bases_deleted or flank_right == bases_deleted:  
                        delSize[delLen][1]+=1 
                    
                if insertLen!=0: 
                    if insertSize.__contains__(insertLen): 
                        insertSize[insertLen][0]+=1 
                    else: 
                        insertSize[insertLen] = [1,0] 
                    numInserts+=1 
                    #is this a repeat insertion? 
                    bases_inserted = query[i-insertLen:i] 
                    flank_left = query[i-2*insertLen:i-insertLen] 
                    flank_right = query[i:min(shorterLen,i+insertLen)]
                    if flank_left == bases_inserted or flank_right == bases_inserted: 
                        insertSize[insertLen][1]+=1 
                delLen = 0 
                insertLen = 0 
                                
                mutationtype[ref[i].lower()][query[i].lower()]+=1

            else: 
                if delLen!=0: 
                    if delSize.__contains__(delLen): 
                        delSize[delLen][0]+=1 
                    else: 
                        delSize[delLen] = [1,0]
                    numDels+=1
                    #is this a repeat deletion?
                    bases_deleted = ref[i-delLen:i] 
                    flank_left = ref[i-2*delLen:i-delLen] 
                    flank_right = ref[i:min(shorterLen,i+delLen)]
                    if flank_left == bases_deleted or flank_right == bases_deleted:  
                        delSize[delLen][1]+=1 
                    
                if insertLen!=0: 
                    if insertSize.__contains__(insertLen): 
                        insertSize[insertLen][0]+=1 
                    else: 
                        insertSize[insertLen] = [1,0]
                    numInserts+=1 
                    #is this a repeat insertion? 
                    bases_inserted = query[i-insertLen:i] 
                    flank_left = query[i-2*insertLen:i-insertLen] 
                    flank_right = query[i:min(shorterLen,i+insertLen)]
                    if flank_left == bases_inserted or flank_right == bases_inserted: 
                        insertSize[insertLen][1]+=1 
                delLen = 0 
                insertLen = 0
        insertsByRead.append(numInserts)
        delsByRead.append(numDels) 

    #write all outputs to file.
    blast.close() 

    fprimerCheck = open(base_output+"primerCheck.csv","w") 
    keys = primercheck.keys() 
    keys.sort() 
    for c in keys: 
        fprimerCheck.write(str(c)+","+str(primercheck[c])+"\n") 
    fprimerCheck.close() 

 
    fposcount = open(base_output+"posCount.csv",'w') 
    for p in poscount:
        fposcount.write(str(p)+","+str(poscount[p])+'\n')
    fposcount.close() 
    

    finsertcount = open(base_output + "insertCount.csv",'w')
    
    for i in insertcount:
        finsertcount.write(str(i) + "," + str(insertcount[i]) + '\n')
        if i in poscount: 
            insertprob[i]= min(1,insertcount[i]/float(poscount[i]))
    finsertcount.close()
    #plot the insertion histogram, if specified in input arguments. 
    if plothistogram: 
        plotname=base_output+"CharacterizationInsertCount.png" 
        fig = plt.figure() 
        ax=fig.add_subplot(111) 
        ax.bar(insertprob.keys(), insertprob.values(),width=1,color='r',log=True) 
        #ax.set_yscale('log')
        ax.set_ylim([10e-5,1])
        plt.setp(ax.get_xticklabels(),fontsize=18)
        plt.setp(ax.get_yticklabels(),fontsize=18)        
        ax.set_xlabel('Base Position Along a Read',fontsize=20) 
        ax.set_ylabel('Probability of Insertion',fontsize=20) 
        ax.set_title('Probability of Insertion as a Function of Base Position\n Dataset '+ str(sys.argv[1].split('/')[-1]),fontsize=20)
        plt.grid(True) 
        plt.savefig(plotname,bbox_inches=0) 
        
    
   
    finsertSize = open(base_output +"insertSize.csv","w") 
    for it in insertSize: 
        finsertSize.write(str(it)+","+str(insertSize[it][0])+ ","+str(insertSize[it][1])+"\n") 
    finsertSize.close() 
    #plot insertion size histogram, if specified in input arguments 
    if plothistogram: 
        plotname=base_output+"CharacterizationInsertSize.png"
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

    finsertsByRead = open(base_output.rstrip() +"insertsByRead.csv","w") 
    if len(insertsByRead) > 0:
        finsertsByRead.write(str(insertsByRead[0]))
        if len(insertsByRead) > 1: 
            for ir in range(1,len(insertsByRead)): 
                finsertsByRead.write(","+str(insertsByRead[ir])) 
        finsertsByRead.write("\n") 
    finsertsByRead.close() 


    fdelcount = open(base_output  + "delCount.csv",'w')
    
    for d in delcount:
        fdelcount.write(str(d) + ","+ str(delcount[d])+'\n')
        if d in poscount:
            delprob[d]=min(1,delcount[d]/float(poscount[d])) 
    fdelcount.close() 

    #plot deletion probability as a function of position, if specified in input arguments. 
    if plothistogram: 
        plotname=base_output+"CharacterizationDelCount.png"
        fig = plt.figure() 
        ax=fig.add_subplot(111) 
        ax.bar(delprob.keys(), delprob.values(),color='r',log=True) 
        #ax.set_yscale('log')
        ax.set_ylim([10e-5,1])
        plt.setp(ax.get_xticklabels(),fontsize=18)
        plt.setp(ax.get_yticklabels(),fontsize=18)        
        ax.set_xlabel('Base Position Along a Read',fontsize=20) 
        ax.set_ylabel('Probability of Deletion',fontsize=20) 
        ax.set_title('Probability of Deletion as a Function of Base Position\n Dataset '+ str(sys.argv[1].split('/')[-1]),fontsize=20)
        plt.grid(True) 
        plt.savefig(plotname,bbox_inches=0) 
               
 
    fdelSize = open(base_output +"delSize.csv","w") 
    for dt in delSize: 
        fdelSize.write(str(dt) + "," + str(delSize[dt][0])+"," + str(delSize[dt][1])+"\n") 
    fdelSize.close() 

    #plot deletion size histogram, if specified in input arguments 
    if plothistogram: 
        plotname=base_output+"CharacterizationDelSize.png"
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

    fdelsByRead = open(base_output +"delsByRead.csv","w") 
    if len(delsByRead) > 0: 
        fdelsByRead.write(str(delsByRead[0]))
        if len(delsByRead) > 1: 
            for dr in range(1,len(delsByRead)): 
                fdelsByRead.write(","+ str(delsByRead[dr])) 
        fdelsByRead.write("\n") 
    fdelsByRead.close() 

    fmutationcount = open(base_output+ "mutationCount.csv",'w')
    for m in mutationcount:
        fmutationcount.write(str(m) + "," + str(mutationcount[m])+"\n") 
        if m in poscount: 
            mutprob[m]=min(1,mutationcount[m]/float(poscount[m]))
    fmutationcount.close()
    #plot the mutation histogram histogram, if specified in input arguments. 
    if plothistogram: 
        plotname=base_output+"CharacterizationMutationCount.png" 
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


    fmutationType = open(base_output + "mutationType.csv","w") 
    for m in mutationtype: 
        fmutationType.write(m) 
        for mprime in mutationtype[m].keys(): 
            fmutationType.write(","+mprime+"," + str(mutationtype[m][mprime])) 
        fmutationType.write("\n") 
    fmutationType.close() 

    mtp=dict() 
    for key in mutationtype:
        mtp[key]=dict() 
        totalmuts=float(sum(mutationtype[key].values()))
        for subkey in mutationtype[key]: 
            mtp[key][subkey]=mutationtype[key][subkey]/max(0.01,totalmuts) 
    mutationDF=DataFrame([[0,mtp['a']['t'],mtp['a']['c'],mtp['a']['g'],mtp['a']['n']],[mtp['t']['a'],0,mtp['t']['c'],mtp['t']['g'],mtp['t']['n']],[mtp['c']['a'],mtp['c']['t'],0,mtp['c']['g'],mtp['c']['n']],[mtp['g']['a'],mtp['g']['t'],mtp['g']['c'],0,mtp['g']['n']],[mtp['n']['a'],mtp['n']['t'],mtp['n']['c'],mtp['n']['g'],0]],columns=['A','T','C','G','N'])
    #plot mutation type histogram, if specified in input arguments 
    if plothistogram: 
        plotname=base_output+'CharacterizationMutationType.png'
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

    fusage = open( base_output + "Usage.csv",'w')
    for u in readUsage:
        fusage.write(str(u)) 
        if len(readUsage[u])>1: 
            readUsageDist.append(min(1,float(readUsage[u][1])/float(readUsage[u][0])))
        for i in readUsage[u]: 
            fusage.write(","+str(i)) 
        fusage.write('\n') 
    fusage.close()
    #plot read usage distribution
    if plothistogram: 
        plotname=base_output+'CharacterizationUsage.png'
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


if __name__ == '__main__':
    main() 
