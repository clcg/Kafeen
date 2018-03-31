#!/bin/bash

# Pare down the INFO field for ExAC to just ACs, ANs, and AFs for every population.
# The AFs are not in the original VCF, so they are calculated using AC/AN.
#
# Requirements:
#   - bcftools must be in your $PATH
#
# Input:
#   - Original ExAC .vcf/vcf.gz/bcf/bcf.gz file
# Output:
#   - Altered ExAC .vcf.gz file (with only allele counts + freqs)
#   - Tabix index file
#
# Example usage:
# ./get_ExAC_AC-AN-AF_tags.sh exac.sample.vcf exac.sample.NEW.vcf.gz

echo "Creating new ExAC file..."
bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/AC_Adj\t%INFO/AN_Adj\t%INFO/AC_AFR\t%INFO/AN_AFR\t%INFO/AC_AMR\t%INFO/AN_AMR\t%INFO/AC_EAS\t%INFO/AN_EAS\t%INFO/AC_FIN\t%INFO/AN_FIN\t%INFO/AC_NFE\t%INFO/AN_NFE\t%INFO/AC_OTH\t%INFO/AN_OTH\t%INFO/AC_SAS\t%INFO/AN_SAS\n' \
  $1 \
  | awk -f includes/lib.awk -f includes/print_new_vcf.awk \
  | bcftools convert -O z -o $2
echo "Done! Output written to $2"

echo "Creating tabix index file..."
bcftools index -ft $2
echo "Done! Output written to $2.tbi"
