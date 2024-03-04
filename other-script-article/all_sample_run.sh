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
	cd ../Tool/sgENERATE/;
	snakemake -c1 all --use-conda ||(echo "Faut tout relancer" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "Plantage" ; exit 1);
	rm -rf result;
	mkdir result;
	mv benchmarks/PeriscopeMult_COV.txt benchmarks/PeriscopeMult_COV_$f.txt;
	mv benchmarks/Periscope_COV.txt benchmarks/Periscope_COV_$f.txt;
	tail -n+2 Periscope/COV_periscope_counts.csv >> finalresult/peri/COV_periscope_counts.csv;
	tail -n+2 Periscope/COV_periscope_novel_counts.csv >> finalresult/peri/COV_periscope_novel_counts.csv;
	tail -n+2 Periscope_mult/COV_periscope_counts.csv >> finalresult/peri2/COV_periscope_counts.csv;
	tail -n+2 Periscope_mult/COV_periscope_novel_counts.csv >> finalresult/peri2/COV_periscope_novel_counts.csv;
	rm -rf Periscope;
	rm -rf Periscope_mult;
	cd $var2;
	echo echo $count / $tt
done ||exit 1
echo "FIN" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "FIN" ;
echo "NADINE" | mail -r bws@univ-lille.fr antoine.limasset@gmail.com -s "NADINE" ;
