import os 
import sys 


def main():
    if len(sys.argv) < 2: 
        raise SyntaxError("specify the .qual file for which quality score stats are to be calculated") 
    meanVals = dict()
    maxVals = dict()
    minVals = dict() 
    modeVals = dict()
    positions  = dict() 
    
    f= open(sys.argv[1],'r') 
    plothist=False 
    if len(sys.argv) > 2: 
        if sys.argv[2] == '-plothistogram': 
            plothist=True 
            import matplotlib.pyplot as plt 
    
    pos = -1 
    first = True
    #iterate through the quality scores at each position, keeping track of the 
    #min, max, mean, and mode values. 
    counter =-1 
    for line in f:
         counter+=1
         #print "line " + str(counter) + "\n"
         if line.startswith('>'): 
             if first == True: 
                 first = False 
             else: 
                 pos = -1 
         else: 
             line = [int(i) for i in line.split(" ")]
             for base in range(len(line)): 
                 pos+=1 
                 if positions.__contains__(pos): 
                     positions[pos]+=1 
                 else: 
                     positions[pos] = 1 
                 if meanVals.__contains__(pos): 
                     meanVals[pos]+=line[base]
                 else:
                     meanVals[pos] = line[base] 
                 if maxVals.__contains__(pos): 
                     if line[base] > maxVals[pos]: 
                         maxVals[pos] = line[base] 
                 else: 
                     maxVals[pos] = line[base] 
                 if minVals.__contains__(pos): 
                     if line[base] < minVals[pos]: 
                         minVals[pos] = line[base] 
                 else: 
                     minVals[pos] = line[base] 
                 #for now, store the counts of each quality score encountered 
                 #for a given position 
                 if modeVals.__contains__(pos): 
                     if modeVals[pos].__contains__(line[base]): 
                         modeVals[pos][line[base]]+=1 
                     else: 
                         modeVals[pos][line[base]]=1 
                 else: 
                     modeVals[pos] = dict() 
                     modeVals[pos][line[base]]=1
    f.close() 
    #print "finished reading qual file \n"
    #print "positions: " + str(positions) +  "\n" 
    #print "meanVals: " + str(meanVals) +"\n" 
    #print "minVals: " + str(minVals) + "\n"
    #print "maxVals: " + str(maxVals) + "\n" 
    #print "modeVals: " + str(modeVals) + "\n" 
    #find mean quality score for each position 
    for position in positions: 
         posCount = positions[position] 
         meanVals[position] = meanVals[position]/float(posCount) 
    #print "meanVals adjusted : " + str(meanVals) +"\n"
    #find the mode of the quality scores for each position
    filteredMode = dict()
    for position in positions: 
         filteredMode[position]= [] 
         posCounts = modeVals[position] 
         maxCount = max(posCounts.values())
         for key in posCounts.keys(): 
             if posCounts[key] == maxCount: 
                 filteredMode[position].append(key) 
    #output the min, max, mean, mode in a .csv file 
    #print "filteredMode: " + str(filteredMode) + "\n" 
    toStrip = sys.argv[1].find(".qual") 
    name = sys.argv[1][0:toStrip]+"qualHist.csv"
    plotname= sys.argv[1][0:toStrip]+"CharacterizationQualHist.png"
    fout = open(name,'w') 
    for position in positions:
         fout.write(str(position)+"," + str(minVals[position])+"," + str(maxVals[position])+"," + str(meanVals[position]))
         for i in filteredMode[position]:
             fout.write("," + str(i)) 
         fout.write("\n") 
    fout.close()
    #plot the quality histogram (if corresponding flag is passed as a command line argument) and 
    #save to a .png file. 

    if plothist: 
        fig=plt.figure() 
        ax=fig.add_subplot(111) 
        ax.bar(meanVals.keys(),meanVals.values()) 
        plt.setp(ax.get_xticklabels(),fontsize=18)
        plt.setp(ax.get_yticklabels(),fontsize=18)        
        ax.set_xlabel('Base Position Along Read',fontsize=20) 
        ax.set_ylabel('Mean Quality Score',fontsize=20) 
        ax.set_title('Mean Quality Score Across Base Position \n Dataset ' + str(sys.argv[1].split('/')[-1]),fontsize=20)
        plt.grid(True) 
        plt.savefig(plotname,bbox_inches=0) 
        #plt.show() 
        

if __name__ == '__main__': 
    main() 
         
     
         
                     
