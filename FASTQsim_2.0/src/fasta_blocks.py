#chops input fasta file into blocks of a specified number of sequences
import sys 
fasta_name = sys.argv[1] 
fastaf = open(sys.argv[1])
blocksize = int(sys.argv[2])
current_block=0
num_added=0
fout=open(fasta_name+"_"+str(current_block),"w") 
line = fastaf.readline() 
while line!="":
    if line.startswith('>')==False: 
        fout.write(line) 
    elif num_added < blocksize: 
        fout.write(line)         
    else: 
        current_block+=1 
        fout=open(fasta_name+"_"+str(current_block),"w") 
        fout.write(line)
        num_added=0 
    line=fastaf.readline() 
    num_added+=1     
print str(current_block)
              
