# Get ESP6500SI AC, AN, and AF tags

This script will convert the downloadable ESP6500 VCF file into a more readable format that separates out the allele counts and frequencies in a way that is similar to ExAC and 1000 Genomes.

In the original VCF file, the alternate and reference allele counts are both stored into the AC (alternate allele count) tags (i.e. EA\_AC and AA\_AC). No AN (total allele count) tag exists, so you must add up the alternate and reference alleles yourself to get this number. Additionally, no AF (alternate allele frequency) tag exists, so you must calculate this number as well using AC/AN.

The point of this script is to split out these fields in a way that makes sense and is easy to work with for annotation purposes.

## Output

The INFO field will be pared down to only the following tags:

- **EVS\_ALL\_AC** - Total number of alternate alleles in called genotypes
- **EVS\_ALL\_AN** - Total number of alleles in called genotypes
- **EVS\_ALL\_AF** - Estimated allele frequency in all populations in the range (0,1)
- **EVS\_EA\_AC** - Number of alternate alleles in called genotypes in the European American population
- **EVS\_EA\_AN** - Number of alleles in called genotypes in the European American population
- **EVS\_EA\_AF** - Allele frequency in the European American populations calculated from AC and AN, in the range (0,1)
- **EVS\_AA\_AC** - Number of alternate alleles in called genotypes in the African American population
- **EVS\_AA\_AN** - Number of alleles in called genotypes in the African American population
- **EVS\_AA\_AF** - Allele frequency in the African American populations calculated from AC and AN, in the range (0,1)

## Usage

First, you will need to download the original VCF files from the EVS server [here](http://evs.gs.washington.edu/evs_bulk_data/ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.vcf.tar.gz), which are are split up by chromosome. You can run this script on each individual file or you can use `bcftools concat` to first combine the files into a single file (as documented [here](https://samtools.github.io/bcftools/bcftools.html#concat)) before running the script. See the example below on how to run the script.

### Example

This script comes with a sample file that you can use to see what the output will look like. You can run the script like follows:

    ./get_ESP6500SI_AC_AN_AF_tags.sh esp6500SI.sample.vcf esp6500SI.sample.NEW.vcf.gz

#### Input file (argument #1):
- Original ESP6500SI .vcf/vcf.gz/bcf/bcf.gz file

#### Output file (argument #2):
- Altered ESP6500SI .vcf.gz file (with .tbi index file) with easier to read allele counts and freqs

## Requirements

- *bcftools* must be in your $PATH (download it [here](https://github.com/samtools/bcftools/releases))

## Author

Sean Ephraim
