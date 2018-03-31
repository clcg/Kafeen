#!/bin/bash

# Submit a qsub job for each chromosome

for chr in {1..22} MT X Y 
do
  qsub -N create_1KG_VCF.chr$chr -v CHR=$chr create_1KG_VCF.qsub.sh
done
