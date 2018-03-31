#!/bin/bash

# Automatically retrieve the latest version of ClinVar from the
# MacArthur Lab git repository (https://github.com/macarthur-lab/clinvar)
#
# Example usage:
#   ./get_newest_ClinVar_from_MacArthur.sh

# Stop on first error
set -e

tsv="hg19_ClinVar_$(date +%Y%m%d).tsv"
vcf="hg19_ClinVar_$(date +%Y%m%d).MORL.LA-norm.vcf"

echo "Retrieving..."
curl -o "$tsv" 'https://raw.githubusercontent.com/macarthur-lab/clinvar/master/output/clinvar.tsv'

echo "Pre-formatting..." # Remove leading/trailing commas and semi-colons
sed -i -e $'s/\t[,;]/\t/g' -e $'s/[,;]\t/\t/g' "$tsv"
echo "Done! Output written to $tsv"

# Check for tab2vcf
if [ -z "$(which tab2vcf 2> /dev/null)" ]; then
  echo "Skipping VCF conversion... tab2vcf is not in your \$PATH"
  echo "Done! Output files:"
  echo "- $tsv"
  exit
fi

echo "Converting to VCF..."
tab2vcf \
  --source 'ClinVar' \
  --reference 'GRCh37.p13' \
  --prefix 'CLINVAR_' \
  --info-tag-map 'info_tag_map.txt' \
  "$tsv" > "$vcf"

if [ -n "$(which bgzip 2> /dev/null)" ] && [ -n "$(which tabix 2> /dev/null)" ]; then
  # Compess/index with bgzip/tabix
  echo "Compressing and indexing $vcf..."
  bgzip -f "$vcf"
  tabix -fp vcf "$vcf.gz"
  echo "Done! Output files:"
  echo "- $tsv"
  echo "- $vcf.gz"
  echo "- $vcf.gz.tbi"
else
  # Skip compression/indexing
  echo "Skipping VCF compression and indexing...  bgzip and/or tabix are not in your \$PATH"
  echo "Done! Output files:"
  echo "- $tsv"
  echo "- $vcf"
fi
