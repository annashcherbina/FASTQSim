import sys
input = sys.argv[1] 
import math 
outputfasta =  input.replace('fastq','fasta') 
outputfasta = outputfasta.replace('fq','fasta') 
outputqual = input.replace('fastq','qual') 
outputqual = outputqual.replace('fq','qual') 
print outputqual+"\n" 
print outputfasta+"\n" 

fin = open(input,'r') 
foutfasta = open(outputfasta,'w') 
foutqual = open(outputqual,'w') 

seqName = None
qual=False
curseq=""
curqual=""
first=True
for line in fin:
    #print str(line)+'\n'
    if line.startswith('@'): #mke sure it is not a qual field that starts with @
        #print "case 1\n"
        if curseq!="" and curqual!="" and seqName!=None: 
            foutfasta.write(seqName.replace('@','>')+"\n")
            foutfasta.write(curseq+"\n")
            foutqual.write(seqName.replace('@','>')+"\n")
            
            #curqual=[str(10*math.log(1 + 10 ** ((ord(i) - 64.0) / 10.0)) / math.log (10)) for i in curqual]
            curqual=str([ord(i)-33 for i in curqual]).replace(']','').replace('[','').replace(',','') 
            foutqual.write(curqual+'\n')
            
            curseq=""
            curqual=""
            seqName = line.rstrip()
            qual=False 
        elif first==True: 
            seqName = line.rstrip() 
            first=False 
        else: 
            #print "case1, quality score \n"
            curqual=curqual+line.rstrip() 

            
    elif line.startswith('+') and len(line.rstrip())==1:
        #print "case 2\n"
        qual=True 
    elif qual==False:
        #print "case 3\n" 
        curseq=curseq+line.rstrip() 
    else:
        #print "case 4\n" 
        curqual=curqual+line.rstrip() 
foutfasta.write(seqName.replace('@','>')+"\n")
foutfasta.write(curseq+"\n")
foutqual.write(seqName.replace('@','>')+"\n")
#curqual=[str(10*math.log(1 + 10 ** ((ord(i) - 64.0) / 10.0)) / math.log (10)) for i in curqual]
curqual=str([ord(i)-33 for i in curqual]).replace(']','').replace('[','').replace(',','') 
foutqual.write(curqual+'\n')
fin.close()
foutfasta.close()
foutqual.close()

