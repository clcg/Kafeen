#!/bin/bash

##
# Remove Duplicates in a VCF
#
# A duplicate variant is when multiple records have the same
# CHROM, POS, REF, and ALT. This script will pick *one* of the
# duplicate variants and discard the rest. The record that is
# picked is the one that comes first in sorting order.
#
# Input:
#  vcf[.gz] or bcf[.gz]
# Output:
#  vcf
#
# Example usage:
#   ./remove_duplicates.sh dups.[vcf|bcf][.gz] > nodups.vcf
##

myvcf="$1"

# Print header
bcftools view -h -O v "$myvcf"

# Sort and remove duplicates per chromosome
for i in {1..22} MT X Y 
do
  bcftools view -H -r "$i" -O v "$myvcf" | sort -u -k2,2n -k4,4d -k5,5d
done
