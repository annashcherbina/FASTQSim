# @author        Anna Shcherbina 
#
#License:          GNU GPL license (http://www.gnu.org/licenses/gpl.html)  
#
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.



For additional information about the algorithms used in this software, see the accompanying paper "Platform-Independent Data Characterizationa and In Silico Read Generation for NGS Datasets"

Dependencies: 
------------
Python version 2.7
       This can be downloaded here http://www.python.org/download/releases/2.7.5/
matplotlib (Python library) version 1.1.0
	This can be downloaded here http://matplotlib.org/downloads.html
BLAST (included) version 2.2.26+
      The included binary is compiled for Feora v. 17, and will not work on MAC or Windows operating systems. 
      The appropriate version of BLAST for your system can be download here ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
scipy (Python library) version 0.12.0
      To download, follow instructions here: http://www.scipy.org/install.html
numpy (Python library) version 1.7.1
      To download, follow instructions here: http://www.scipy.org/install.html
Java Blast Parser (included) version 1.0 
samtools (included) version 0.1.19

The software provides two types of capabilities: 
-----------------------------------------------
1. Next Generation Sequencing (NGS) datasets characterization. Characterization refers to the following: 
   a. Histogram of read length.
   b. Min/Max/Mean quality scores for each base position in a read. 
   c. Insertion characterization. This refers to calculating the probability of an insertion at each position in a sequence, the probability of insertions of a particular size, and the probability of a repeat insertion of a particular size. 
   d. Deletion characterization, same parameters as insertions. 
   e. Single base pair mutations. This refers to calculating the probability of a mutation at each position in a sequence, the probability of mutating from one base to another       (i.e. c--> g versus c--> a) 
   f. Primer histogram. By computing the frequency of each primer, "pileups" of low complexity reads can be detected. 
   g. Read alignment quality histogram. Read alignment refers to the percent of a read from the data that aligns to a matching hit from the BLAST database. This value depends on the data quality and the instrument used: the value is around 98% for Roche 454 data, but around 90% for PacBio data.   
2. "Spiking" of reads from a fasta source into a background fastq dataset. The spiking is performed in a manner such that the characteristics of the spiked reads are identical    to those of the background reads, so a spiked read is indistinguishable from a background read. 

Executing the software: 
-----------------------
There are two main scripts in the top-level "FASTQSim" directory: FASTQcharacterize.sh and FASTQspike.sh These should be run to perform the characterization and spiking tasks described above. Helper python files are located in the "src" directory. 
An example is provided in the "example.sh" script in the top-level directory. 

1. FASTQcharacterize.sh 
The script can be run in three ways: 
    a. Using an externally generated alignment file (BAM or SAM format) to perform characterization 

    sh FASTQcharacterize.sh -input Example5.fq -reference homo_sapien.gdna -sam Example5.sam 

    b. Using an externally generated BLAST alignment file to perform characterization 

    sh FASTQcharacterize.sh -input Example5.fq -reference homo_sapien.gdna -blastresults Example5.blast 

    c. Aligning a subset of the FASTQ file to a reference with the BLAST algorithm and performing characterization on the resulting alignment 

    sh FASTQcharacterize.sh -input Example5.fq -reference homo_sapien.gdna  
    

Inputs: 
	a. The fastq file to be characterized. (-input ) 
	b. A FASTA database to be used for BLAST. (-reference)  The FASTA database should mirror the source of the host/background in the input FASTQ file. In the provided example, "Example5.fq", there is a human host, so the associated database is "homo_sapiens.gdna". By default, the characterization program will select 5000 reads from the input file, perform a BLAST search against the referene database, and parse the output to determine the input file characteristics. Two-threaded BLAST with a maximum of 10 hits per input sequence is used. These parameters are selected for efficiency. For a more thorough characterization of the data the "-full" flag should be used. If this flag is supplied, BLAST will be performed on the entire input FASTQ file, not just 5000 sequences. 
	c. Alternatively, a BLAST results file may be used as an input (-blastresults). If the BLAST file is provided for the input FASTQ, a reference sequence does not need to be supplied. The software will use the supplied BLAST result to characterize the data. 
	d. BAM/SAM alignment file (-bam / -sam) 
	If a BAM or SAM file is specified, this file will be used to calculate characterization statistics for the dataset. The algorithm converts a BAM file to SAM format and calls the "calmd" command from samtools to identify SNPs in the aligned reads. 
	e.If the -plothistogram flag is set, graphs will be generated for the output characterization files. 

Outputs: 
All outputs will be stored in the "FASTQsim" directory. The following files will be generated (using the provided Example5.fq input file as an example): 

a. Example5summary.csv -- This is a list of all output characterization files: 

Example5readHist.csv
Example5qualHist.csv
Example5Usage.csv
Example5primerCheck.csv
Example5posCount.csv
Example5delCount.csv
Example5delsByRead.csv
Example5delSize.csv
Example5insertCount.csv
Example5insertsByRead.csv
Example5insertSize.csv
Example5mutationCount.csv
Example5mutationType.csv


b. Example5readHist.csv
Read length, number of reads with this length in the input fastq file
5,144
6,321
7,530
8,692
9,676
10,697
11,651
12,569
13,611
14,633
15,580
16,601
... (truncated)


c. Example5qualHist.csv 
Read base position, minimum observed quality score, maximum observed quality score, mean observed quality score, mode observed quality score 
0,4,37,29.6080046091,34 (base 0, minimum score 4, maximum score 37, mean score 29.608, mode score 34) 
1,4,37,30.0692716194,34
2,4,37,29.5020540696,34
3,4,37,28.915623454,34
... (truncated) 


d. Example5Usage.csv 
Read Name, Read length, Number of bases aligned with BLAST target, first base aligned to target, last base aligned to target

9HKOK:1139:2115,222,220,1,220
9HKOK:688:1123,161,161,1,161
9HKOK:2070:515,208,207,1,207
9HKOK:659:1833,167,167,1,167
9HKOK:2263:2608,90
9HKOK:1859:1922,215,213,1,213
9HKOK:1999:1449,15 --> No alignment was found, read length was 15 and other fields are empty
9HKOK:2280:1572,212,212,1,212
... (truncated) 

e. Example5primerCheck.csv
The purpose of this file is to detect "pile ups" of low-complexity reads, which are likely to have the same group of starting bases. 
First four bases of read, Number of times this set of bases was observed 
-aaa,1
-aac,1
-aat,1
-aca,2
-acc,2
-act,1
-aga,3
-agc,1
-agg,1
-agt,1

f. Example5posCount.csv  
Base position on a read, Number of reads that include this base position (i.e. if a read is 250 bases long, it will have base positions 1-250). 
1,1873
2,1886
3,1908
4,1917
5,1917
6,1921
... (truncated) 

g. Example5delCount.csv
Base position, Number of deletions observed at base position 
Dividing this number by the position Count (above) gives the probability of deletion at a given base. 
3,1
4,2
5,1
6,2
7,4
8,2
9,2
10,3
11,4
12,4
... (truncated) 


h. Example5delsByRead.csv
Number of deletions in each read, order is the same as order of reads in input fastq file. 
0,0,1,0,0,1,1,4,1,0,1,0,2,1,0,1,0,0,1,0,0,0,1,0,0,...(truncated) 

i. Example5delSize.csv
Deletion size, Number of deletions observed of this size, Number of the observed deletions that were repeat deletions
1,1382,1065
2,45,17
3,7,4
4,5,5
5,2,2
6,7,4
7,2,2
8,1,1
9,1,0
12,1,0
14,1,0
15,1,0
16,1,0

j. Example5insertCount.csv
Base position, Number of insertions observed at base position 
Dividing this number by the position Count (above) gives the probability of insertion at a given base. 
3,6
4,10
5,12
6,9
7,15
8,16
9,8
...(truncated)

k. Example5insertsByRead.csv
Number of insertions in each read, order is the same as order of reads in input fastq file. 
5,3,0,1,1,2,0,4,1,3,2,1,...(truncated) 

l. Example5insertSize.csv
Insertion size, Number of insertions observed of this size, Number of the observed insertions that were repeat insertions
1,4084,3042
2,71,15
3,4,3
6,1,0

m. Example5mutationCount.csv
Base position, Number of single-base mutations observed at base position
Dividing this number by the position count (above) gives the probability of mutation at a given base. 
3,4
4,8
5,8
6,5
7,3
8,11
9,2
10,6
...(truncated)

n. Example5mutationType.csv

Starting base, Final base 1, Number of times observed, Final base 2, Number of times observed, Final base 3, number of times observed, Final base "n", number of times observed

a,c,56,t,88,g,167,n,0
c,a,63,t,109,g,45,n,0
t,a,80,c,183,g,51,n,0
g,a,145,c,69,t,73,n,0
n,a,0,c,0,t,1,g,1

2. FASTQspike.sh
In the provided example, the background file is "Example5.fq", the same that was characterized in the above step. An E. coli genome is spiked in 
at a coverage of 20x. This genome is circular. A plasmid, EU496103_with_xis_AAv.fasta is spiked in at a coverage of 200x. This plasmid is also circular. See the script example.sh for how to perform 
this spiking operation. 

Inputs: 
a. Background FASTQ file to be spiked (Example5.fq) 

It may be the case that only a fraction of the background file should be spiked with the in silico reads, inputs  b and c allow for this capabilitly.  
b. Number of blocks to split the background file into (set to "1" if entire background file should be spiked)
c. The block number in the background file to use for in silico spiking (set to "0" if entire background should be used) 

d. List of background characterization files. This is typically the summary.csv file generated by the characterization script (i.e. Example5summary.csv) 

Multiple read sources can be used for spiking. Provide inputs e,f,g for each source of spiked reads. 
e. coverage for spiked reads (20x in example.sh)  
f. fasta file for spiked reads (ecoli.fasta) 
g. Are the spiked reads from a circular genome (True/False)
h If the optional -plothistogram flag is set, statistical plots will be generated comparing the spiked reads with the background reads.  

Outputs: Stored in the "FASTQsim" directory 

a. FASTQ file containing the spiked reads generated from the sources (example/Example5_Spiked.fq) 
b. FASTQ file containing the full dataset consisting of the background FASTQ file and the spiked reads inserted at random positions (example/Example5_Full.fq) 



