# Convert 1000G phase 3 VCF

This script will pare down the original 1KG VCF files to only contain AC, AN, and AF tags in the INFO column for all populations. AC and AN do not exist for sub-populations in the original VCFs, so they are calculated using BCFtools. The CHROM, POS, ID, REF, ALT, QUAL, and FILTER columns will remain exactly the same. Any columns after the INFO column will be dropped.

## Output

A new VCF file will be created for each specified chromosome. All output files will be written to the current directory unless you edit the `OUTPUT_DIR` constant within the script itself. All output file names will be in the format *hg19\_1KG.chr1.vcf[.gz][.tbi]*, *hg19\_1KG.chr2.vcf[.gz][.tbi]*, *hg19\_1KG.chr3.vcf[.gz][.tbi]*, etc.

The INFO column will only contain the following tags:

- **1KG\_ALL\_AC** - Total number of alternate alleles in called genotypes
- **1KG\_ALL\_AN** - Total number of alleles in called genotypes
- **1KG\_ALL\_AF** - Allele frequency in the all populations calculated from AC and AN, in the range (0,1)
- **1KG\_AFR\_AC** - African/African American alternate allele counts
- **1KG\_AFR\_AN** - African/African American allele counts
- **1KG\_AFR\_AF** - Allele frequency in the African/African American populations calculated from AC and AN, in the range (0,1)
- **1KG\_AMR\_AC** - American alternate allele counts
- **1KG\_AMR\_AN** - American allele counts
- **1KG\_AMR\_AF** - Allele frequency in the American populations calculated from AC and AN, in the range (0,1)
- **1KG\_EAS\_AC** - East Asian alternate allele counts
- **1KG\_EAS\_AN** - East Asian allele counts
- **1KG\_EAS\_AF** - Allele frequency in the East Asian populations calculated from AC and AN, in the range (0,1)
- **1KG\_EUR\_AC** - European alternate allele counts
- **1KG\_EUR\_AN** - European allele counts
- **1KG\_EUR\_AF** - Allele frequency in the European populations calculated from AC and AN, in the range (0,1)
- **1KG\_SAS\_AC** - South Asian alternate allele counts
- **1KG\_SAS\_AN** - South Asian allele counts
- **1KG\_SAS\_AF** - Allele frequency in the South Asian populations calculated from AC and AN, in the range (0,1)

## Usage

### Downloading 1KG

First, you will need to download the original VCF files (25 total), index files (25 total) and panel files (2 total) from the 1KG server [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). **These files must be in the same directory as the scripts.**

The following command will download of necessary files:

    wget -A "*.genotypes.vcf.gz*,*.panel" ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/*

You can also run:

    ./download_1000Gp3_VCFs_and_panels.sh

Then your ready to run the conversion script on the files. See the example below on how to run the script.

This script will work for chr1-22/MT/X/Y VCFs and will create a VCF for each input chromosome (i.e. 25 output VCFs if you run all chromosomes). BCFtools is used to recalculate AC and AN using a population-specific subset of samples.

### Converting 1KG VCFs

#### Example 1: Running the script on all chromosomes [default]

    ./generate_1KG_AC-AN-AF_VCFs.rb
    
**Tip:** If you need to run the script on all chromosomes then it is recommended to submit one `qsub` job per chromosome (25 total jobs) using the `submit_multiple_qsubs.sh` script.

    ./submit_multiple_qsubs.sh
    
This will automatically launch 25 separate `qsub` jobs, each running the conversion script on a single chromosome.
    
#### Example 2: Running the script on all chromosome 4 only

    ./generate_1KG_AC-AN-AF_VCFs.rb 4
    
#### Example 3: Running the script on all chromosomes 9, 20, X, and MT only

    ./generate_1KG_AC-AN-AF_VCFs.rb 9,20,X,MT
    
### Combining all output VCFs

The `concat_and_normalize_1KG_VCFs.sh` script will concatenate all output VCFs (i.e. *hg19\_1KG.chr1.vcf.gz*, *hg19\_1KG.chr2.vcf.gz*, *hg19\_1KG.chr3.vcf.gz*, etc.) into a single BCF file. Then it will attempt to normalize and left-align the BCF file.

**In order for the normalization and left-alignment to work, you must have the *human\_g1k\_v37.fasta* and *human\_g1k\_v37.fasta.fai* files in the same directory as the script itself.** You can download these two files with the following command:

    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.*

Then uncompress the FASTA file:

    gunzip human_g1k_v37.fasta.gz

Now you can run the following script to concatenate, normalize, and left-align the output:

    ./concat_and_normalize_1KG_VCFs.sh

**Note:** If for some reason the normalization / left-alignmentment process fails, the non-normalized / non-left-aligned BCF will still remain.

## Requirements

`bcftools` must be in your `$PATH` (download it [here](https://github.com/samtools/bcftools/releases))

## Author

Sean Ephraim
