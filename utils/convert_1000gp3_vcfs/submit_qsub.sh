#!/bin/bash
#
#$ -S /bin/bash
#$ -cwd
#$ -j y

# This is a helper script to launch separate qsub jobs for
# separate chromosomes.

# Insert your shell commands below...
ruby generate_1KG_AC-AN-AF_VCFs.rb $CHR
