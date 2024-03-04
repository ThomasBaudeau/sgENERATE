#!/bin/bash

var1=../Tool/sgENERATE/result
var2=../../datasets/
count=0
tt= ls |  wc -l 

conda activate snakemake;
for f in *;
	do cp $f $var1;
	count=$(($count+1))
	mv $var1/$f $var1/final_COV_agregate.fastq ;
	echo $f > $var1/file_name.txt
	cd ../Tool/sgENERATE/;
	sgENERATE --real True --mode extraction ||exit
	cat result/non_COV_cano.fasta >> finalresult/all_ncsgRNA.fasta
	rm -rf result;
	mkdir result;
	rm -rf Periscope;
	rm -rf Periscope_mult;
	cd $var2;
	echo echo $count / $tt
done ||exit 
echo "FIN" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "FIN" ;
echo "NADINE" | mail -r bws@univ-lille.fr antoine.limasset@gmail.com -s "NADINE" ;
