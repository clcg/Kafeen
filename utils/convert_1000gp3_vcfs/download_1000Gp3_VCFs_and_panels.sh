#!/bin/bash

# Download panel files, all VCFs, and all indexes from the 1000 Genomes FTP server
#
# Example usage:
#   ./download_1000Gp3_VCFs_and_panels.sh

wget -A "*.genotypes.vcf.gz*,*.panel" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*
