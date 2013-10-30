inputfastq=Example5.fq
inputref=homo_sapiens.gdna


#using -sam input 
sh FASTQcharacterize.sh -input $inputfastq -reference $inputref -sam Example5.sam 

#sh FASTQcharacterize.sh -input lambda_aligned_filtered_subreads_5k.fastq -reference cchf_r2_falsepos.fasta -sam lambda_bwasw.sam 
#sh FASTQspike.sh lambda_aligned_filtered_subreads_5k.fastq 1 0 lambda_aligned_filtered_subreads_5ksummary.csv 200 lambda.fasta False -o lambda_test

#original approach 

#sh FASTQcharacterize.sh -input $inputfastq -reference $inputref 
#sh FASTQspike.sh $inputfastq 1 0 Example5summary.csv 0.5 ecoli.fasta True 1 EU496103_with_xis_AAV.fasta True -o SpikedResultExample

#if you want to plot figures: 

#sh FASTQcharacterize.sh -input $inputfastq -reference $inputref -plothistogram 
#sh FASTQspike.sh $inputfastq 1 0 Example5summary.csv 20 ecoli.fasta True 200 EU496103_with_xis_AAV.fasta True -plothistogram -o SpikedResultExample
