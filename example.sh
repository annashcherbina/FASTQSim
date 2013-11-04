#inputfastq=Example5.fq
#inputref=homo_sapiens.gdna


#characterization

#using -sam input 
#sh FASTQcharacterize.sh -input $inputfastq -reference $inputref -sam Example5.sam 

#using BLAST alignment tool
#sh FASTQcharacterize.sh -input $inputfastq -reference $inputref 


#spiking 

#sh FASTQspike.sh $inputfastq  Example5summary.csv 0.5 ecoli.fasta True 1 EU496103_with_xis_AAV.fasta True -o SpikedResultExample -threads 10

#if you want to plot figures: 

#sh FASTQcharacterize.sh -input $inputfastq -reference $inputref -plothistogram 
#sh FASTQspike.sh $inputfastq  Example5summary.csv 20 ecoli.fasta True 200 EU496103_with_xis_AAV.fasta True -plothistogram -o SpikedResultExample -threads 10

