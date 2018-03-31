#!/bin/bash

##
# Concatenate and Normalize 1KG VCFs
#
# This is an example script to show you can concatenate all
# your VCFs into one BCF, and then left-align and normalize
# the BCF.
#
# The GRCh37 FASTA file can be downloaded with the following command:
#   wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.*
#
# Example usage:
#   ./concat_and_normalize_1KG_VCFs.sh
##

set -e

# Compile list of VCFs to concatenate
vcf_list=''
for chr in {1..22} MT X Y
do
  vcf_list="$vcf_list hg19_1KG.chr$chr.vcf.gz"
done

echo "Concatenating chr1-22|MT|X|Y VCFs..."
bcftools concat -Ob -o hg19_1000G_phase3v5a.MORL.bcf.gz $vcf_list
echo "Done! Output written to hg19_1000G_phase3v5a.MORL.bcf.gz"

echo "Indexing hg19_1000G_phase3v5a.MORL.bcf.gz..."
bcftools index -fc hg19_1000G_phase3v5a.MORL.bcf.gz
echo "Done! Index written to hg19_1000G_phase3v5a.MORL.bcf.gz.csi"

echo "Left aligning and normalizing BCF..."
bcftools norm -m- -f human_g1k_v37.fasta -Ob -o hg19_1000G_phase3v5a.MORL.LA-norm.bcf.gz hg19_1000G_phase3v5a.MORL.bcf.gz
echo "Done! Output written to hg19_1000G_phase3v5a.MORL.LA-norm.bcf.gz"

echo "Indexing hg19_1000G_phase3v5a.MORL.LA-norm.bcf.gz..."
bcftools index -fc hg19_1000G_phase3v5a.MORL.LA-norm.bcf.gz
echo "Done! Index written to hg19_1000G_phase3v5a.MORL.LA-norm.bcf.gz.csi"
