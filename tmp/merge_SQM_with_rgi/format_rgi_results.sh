#!/bin/bash

sed 's/ # /	/g' ../rgi_results/rgi_ORF_chicken.txt > ../rgi_results/temp.txt
sed 's/ID=.*					//g' ../rgi_results/temp.txt > ../rgi_results/temp2.txt && mv ../rgi_results/temp2.txt ../rgi_results/temp.txt
sed 's/	Contig//g' ../rgi_results/temp.txt > ../rgi_results/rgi_ORF_chicken_formatted.txt
rm ../rgi_results/temp.txt
