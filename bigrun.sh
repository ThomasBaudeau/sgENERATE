#!/bin/bash

var1=../Tool/sgENERATE/result
var2=../../datasets/


conda activate snakemake;
for f in *;
	do cp $f $var1;
	mv $var1/$f 'final_COV_agregate.fastq.gz'  ;
	cd ../Tool/sgENERATE/;
	snakemake -c16 all --use-conda ||(echo "Faut tout relancer" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "Plantage" ; exit);
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
done
echo "FIN" | mail -r bws@univ-lille.fr thomas.baudeau@univ-lille.fr -s "FIN" ;
echo "NADINE" | mail -r bws@univ-lille.fr antoine.limasset@gmail.com -s "NADINE" ;
