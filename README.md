# Kafeen

## Description
Kafeen is a comprehensive open-source pipeline for the collection, annotation and classification of Genetic Variants. Its output is configurable and is meant to be tailored to your own local database needs

## System Requirements
### Hardware Requirements
Kafeen was developed and tested on the [University of Iowa Argon High Performance Computing Cluster](https://wiki.uiowa.edu/display/hpcdocs/Argon+Cluster), utilizing the following computing resources:
- 4 2.4GHz Intel Broadwell processor cores
- 38 GB memory

### Software Requirements
#### OS Requirements
Kafeen was developed and tested in a CentOS 7.4 environemnt

#### Package Requirements
- Ruby 2.0+
   - If your Ruby version is not up to date, it highly recommended to use [RVM](https://github.com/rvm/rvm) to safely install another version on your machine. No root permissions required.
- `bcftools` version 1.3+ must be in your `$PATH` (download [here](https://github.com/samtools/bcftools/releases))
- `bedtools` version 2.0+ must be in your `$PATH` (download [here](https://github.com/arq5x/bedtools2/releases))
- `tabix` and `bgzip` must be in your `$PATH` (download/compile them as part of the [HTSlib package](https://github.com/samtools/htslib/releases))
- Java 1.7+ is required if running ASAP

## Installation Guide
### Code Installation
The latest release of Kafeen is available for download from this GitHub repository.

### Data Formatting
Configuration of the Kafeen pipeline is open and flexible to the user's needs and decisions.
These annotation sources and formatting methods were used in development of Kafeen and the Deafness Variation Database:
#### 1000 Genomes
- Source: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- Formatting: utils/convert_1000gp3_vcfs
#### ClinVar
- Source: https://github.com/macarthur-lab/clinvar/blob/master/output/b37/single/clinvar_alleles.single.b37.tsv.gz
- Formatting: utils/convert_macarthur_clinvar_tsv_to_vcf
#### EVS ESP6500SI
- Source: http://evs.gs.washington.edu/EVS/#tabs-7
- Formatting: utils/get_esp6500si_ac_an_af_tags
#### ExAC
- Source: ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3/
- Formatting: utils/get_exac_ac_an_af_tags
#### dbNSFP
- Source: https://sites.google.com/site/jpopgen/dbNSFP
- Formatting: utils/convert_dbnsfp_to_bcf
#### dbSNP
- Source: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF
- Formatting: utils/convert_dbsnp
#### HGMD
- Source: Genome Trax license required (https://www.qiagenbioinformatics.com/products/genome-trax/)
- Formatting: utils/hgmd2vcf

### Configuration
The Kafeen pipeline configuration (config/config.yml) allows for easy inclusion and exclusion of annotation sources and fields.
At minimum, the configuration must specify paths to the included annotation sources.

### Estimated Installation Time
From scratch, the full download, formatting, and configuration of Kafeen and associated data inputs can take 1 hour to complete.

## Running Kafeen
To run the Kafeen pipeline, two files are required as input:
- a text file containing list of genes according to HGNC guidelines, one per line (GJB2\nPTPRQ\n, eg.)
- properly-configured YAML configuration file

Running the Kafeen pipeline from the command line:
- `ruby kafeen.rb -c <path_to_config> -o <output_prefix> <gene_list>`


## Authors

- Sean Ephraim
- Brad Crone
