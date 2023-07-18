#!/bin/bash

var1=/home/baudeau/Documents/BenchMapL/ProjetSuede/sGenerate/result/final_COV_agregate.fastq


cat result.txt | while read l; 
	do echo $l; 
done
