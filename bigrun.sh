#!/bin/bash

var1=../Tool/sgENERATE/result
var2=../../datasets/
for f in *;
	do mv $f 'final_COV_agregate.fastq.gz';
	mv 'final_COV_agregate.fastq.gz' $var1 ;
	conda activate Snakemake;
	cd ../Tool/sgENERATE/;
	snakemake -c8 all --use-conda;
	rm -r result;
	mkdir result;
	mv benchmarks/PeriscopeMult_COV.txt benchmarks/PeriscopeMult_COV_$f.txt;
	mv benchmarks/Periscope_COV.txt benchmarks/Periscope_COV_$f.txt;
	tail -n+2 Periscope/COV_periscope_counts.csv >> finalresult/peri/COV_periscope_counts.csv;
	tail -n+2 Periscope/COV_periscope_novel_counts.csv >> finalresult/peri/COV_periscope_novel_counts.csv;
	tail -n+2 Periscope_mult/COV_periscope_counts.csv >> finalresult/peri2/COV_periscope_counts.csv;
	tail -n+2 Periscope_mult/COV_periscope_novel_counts.csv >> finalresult/peri2/COV_periscope_novel_counts.csv;
	rm -r Periscope;
	rm -r Periscope_mult;
	cd $var2;



done

