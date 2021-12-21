#!/bin/bash
#BSUB -n 40                                 # number of threads
#BSUB -W 180                               # estimated time to run
#BSUB -R "rusage[mem=5000, scratch=2000]"   # memory and disk space needed
#BSUB -e results_human/error.log                          # error file
#BSUB -o results_human/out.log                            # output file
#BSUB -u ssundar@student.ethz.ch                   # specify your email address
#BSUB -B                                    # send email when job starts
#BSUB -N                                    # send email when job ends

#installed package as in https://github.com/arpcard/rgi. If installation does not does not work, try downgrading to python version 3.6. This fixes the issue

#load CARD data locally  with the following commands
#wget https://card.mcmaster.ca/latest/data
#tar -xvf data ./card.json
#rgi load --card_json /path/to/card.json --local

rgi main --input_sequence data/01.SqueezeHuman.fasta  --output_file results_human/rgi_annotation_human --input_type contig --local --alignment_tool DIAMOND --num_threads 40 --clean


#takes around 90 min to run the human contig data on 40 threads and 5 GB of memory per thread. Max memory used in the job though was 12 GB so in reality we don't need this much memory.

#takes around 10 min to run the chicken orf data on 40 threads and 2 GB of memeory. Also was generous here with memory allocation, probably do not need this much memory. 
