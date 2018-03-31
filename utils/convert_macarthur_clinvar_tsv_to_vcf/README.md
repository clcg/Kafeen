# Convert MacArthur ClinVar TSV to VCF

This script will (1.) download the [latest release of MacArthur Lab's version of ClinVar](https://github.com/macarthur-lab/clinvar/tree/master/output) in TSV format and (2.) convert it to a VCF file.

## Output

The output will be a VCF file (compressed and indexed if possible). The MacArthur Lab splits and left-aligns their TSV file, so in turn, the output VCF will also be split and left-aligned. There is no need for further normalization with BCFtools.

The INFO column will contain the following tags:

- **CLINVAR\_CLNALLELE** - Tells you which allele, REF or ALT, is the one to which the annotations (e.g. pathogenic assertion) refer
- **CLINVAR\_VID** - Variation ID; unique identifier for the set of sequence changes that were interpreted; access online at ncbi.nlm.nih.gov/clinvar/variation/{VID}
- **CLINVAR\_HGNC** - HGNC gene symbol
- **CLINVAR\_CLNSIG** - Clinical significance (e.g. "Pathogenic", "Likely pathogenic", etc.)
- **CLINVAR\_REVSTAT** - Clinical review status
- **CLINVAR\_HGVS\_C** - HGVS cDNA name
- **CLINVAR\_HGVS\_P** - HGVS protein name
- **CLINVAR\_SUBMITTERS** - The names of clinical review submitters
- **CLINVAR\_DISEASE** - Variant disease name(s)
- **CLINVAR\_PMID** - Related PubMed IDs
- **CLINVAR\_PATHOGENIC** - Has this variant been asserted 'Pathogenic' or 'Likely pathogenic' by any submitter for any phenotype? 1 - Yes, 0 - No
- **CLINVAR\_CONFLICTED** - Has this variant ever been asserted 'Pathogenic' or 'Likely pathogenic' by any submitter for any phenotype and also been asserted 'Benign' or 'Likely benign' by any submitter for any phenotype? 1 - Yes, 0 - No; Note that having one assertion of pathogenic and one of uncertain significance does not count as conflicted for this column

## Usage

Simply run:

    ./get_newest_ClinVar_from_MacArthur.sh
   
## What is *info\_tag\_map.txt*?

This is a map file that `get_newest_ClinVar_from_MacArthur.sh` automatically looks for when converting the original TSV file to a VCF file. It contains 3 tab-separated columns:

- **Column 1**: Names of the columns in the original TSV file
- **Column 2**: Descriptions to use for the corresponding INFO tag in the output VCF file
- **Column 3**: Preferred names of the respective INFO tags in the output VCF file

For example:

     symbol      HGNC gene symbol    HGNC
     hgvs_c      HGVS cDNA name      HGVS_C
     hgvs_p      HGVS protein name   HGVS_P
     all_pmids   PubMed IDs          PMID

All INFO tags will automatically be prefixed with the `CLINVAR_` string when the VCF file is written. If a 

## Handling illegal VCF characters

The VCF format restricts all spaces, semi-colons, and equals-signs from being present in an any INFO field. Additionally, commas are reserved for separating allele-specific values. Therefore, these characters are [URL-encoded](http://www.w3schools.com/tags/ref_urlencode.asp) before they are written to the output VCF.

| Illegal VCF character | URL-encoding |
| --------------------- | ------------ |
| *space*               | %20          |
| ,                     | %2C          |
| ;                     | %3B          |
| =                     | %3D          |

## Requirements

- `tab2vcf` must be in your `$PATH` (download it [here](https://github.com/sephraim/bin4matics)). VCF conversion will be skipped if this is not found.
- `bgzip` and `tabix` must both be in your `$PATH` (download with the HTSlib package [here](https://github.com/samtools/htslib/releases/)). If one or both of these are not found, then compression and indexing will be skipped.

## Author

Sean Ephraim
