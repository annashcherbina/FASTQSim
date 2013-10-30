import sys 
import os 
def main(): 
    if len(sys.argv)<2: 
        raise SyntaxError("specify the fasta file for which the read length distribution is to be calculated.")
    readLengths = dict() 
    f = open(sys.argv[1],'r') 
    plothist=False
    if len(sys.argv) > 2: 
        if sys.argv[2] == '-plothistogram': 
            plothist=True 
            import matplotlib.pyplot as plt 
    counter = 0 
    num_reads_total =0 
    
    first = True
    current_read_name = "" 
    current_read = "" 
    

    for line in f: 
        if line.startswith('>'):
            if first == True: 
                first = False 
                current_read_name = line 
            else: 
                if readLengths.__contains__(counter) == False: 
                    readLengths[counter] = 1 
                else: 
                    readLengths[counter]+=1 
                num_reads_total+=1 
                counter=0 
            current_read_name = line
            current_read = "" 
        else: 
            current_read+=line 
            line = line.strip('\n') 
            line = line.strip('\r') 
            counter+= len(line)  
 
    f.close() 
    toStrip = sys.argv[1].find(".fasta") 
    name =  sys.argv[1][0:toStrip]+"readHist.csv" 
    plotname=sys.argv[1][0:toStrip]+"CharacterizationReadHist.png" 
    fout = open(name,'w')
    readKeys = readLengths.keys() 
    readKeys.sort() 
    for key in readKeys: 
        fout.write(str(key) +"," +str(readLengths[key])+"\n") 
    fout.close() 

    #plot the read lengths histogram (if corresponding flag is passed as command line argument) and 
    #save to a .png file. 
    if plothist: 
        fig=plt.figure() 
        ax=fig.add_subplot(111) 
        ax.bar(readLengths.keys(),readLengths.values()) 
        plt.setp(ax.get_xticklabels(),fontsize=18)
        plt.setp(ax.get_yticklabels(),fontsize=18)
        ax.set_xlabel('Read Length (in bases)',fontsize=20) 
        ax.set_ylabel('Number of Reads',fontsize=20) 
        ax.set_title('Read Length Histogram \n Dataset ' + str(sys.argv[1].split('/')[-1]),fontsize=20)
        plt.grid(True)
        plt.savefig(plotname, bbox_inches=0)
        
        
    print "total reads: " + str(num_reads_total) 

 
if __name__ =='__main__':
    main() 
