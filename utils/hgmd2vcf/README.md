# HGMD2VCF

## Description

Extract HGMD data from a Genome Trax MySQL database and put it into an hg19 VCF file.

## Output

### Files

This script will produce a single VCF output file. However, if `bgzip` and `tabix` are in your `$PATH` then this script will automatically compress the VCF file into a *.vcf.gz* file and create a *.vcf.gz.tbi* index file. This script requires one argument, which is the preferred name of the output VCF file. See below for an example.

### VCF columns

The Genome Trax database contains a table called `ngs_feature` which includes a column called `description`. This is where most of the annotation fields are stored. The entire `description` column gets translated into the INFO column of the output VCF. Additionally, some values are also used for other VCF columns. See below for details.

#### Column mapping

| Output VCF column | Genome Trax column / field                 |
| ----------------- | ------------------------------------------ |
| CHROM             | `chromosome` column                        |
| POS               | `feature_start` column                     |
| ID                | `rsid` field from the `description` column |
| REF               | `ref` field from the `description` column  |
| ALT               | `alt` field from the `description` column  |
| QUAL              | No translation; always a `.`               |
| FILTER            | No translation; always a `.`               |
| INFO              | All fields from the `description` column   |

### INFO tags

The fields in the `description` column of Genome Trax are used to create the INFO tags of the output VCF file. Take a look at the manual included with your version of Genome Trax to see what columns exist. To create the INFO tag name, the original field name is (1.) converted to all-caps, and (2.) prefixed with `HGMD_`.

As an example, the `variantType` field in the `description` column will become a VCF INFO tag called `HGMD_VARIANTTYPE`.

While all of the fields in the `description` column get converted to INFO tags in the output VCF, the following is a list of notable INFO tags:

- **HGMD_VARIANTTYPE** - Clinical significance value
- **HGMD_DISEASE** - Disease description
- **HGMD_CONFIDENCE** - Clinical significance confidence value
- **HGMD_PMID** - PubMed IDs for this variant

If you are unsure of how the script converts some `description` fields, you can view all INFO tags in the output VCF by running:

    bcftools view -h hgmd.vcf.gz | grep '^##INFO='

## Example usage

This script only requires one argument, which is the name you would like to call your output VCF file. For example:

    ruby hgmd2vcf.rb hgmd.vcf

In this case, the output will be a VCF file called *hgmd.vcf*. However, if `bgzip` and `tabix` are in your `$PATH` then the VCF will automatically be compressed and indexed, and there will be two output files: *hgmd.vcf.gz* and *hgmd.vcf.gz.tbi*.
 
## Requirements

This script requires the `mysql2` Ruby gem. Install it with:
    
    gem install mysql2
    
Although it is not required, this script will attempt to compress and index the output VCF with `bgzip` and `tabix`, respectively. If these utilities are not in your `$PATH` then the script will skip this step.

## Genome Trax database specs

### Login credentials

The script already contains the database credentials to automatically log into the Genome Trax database on the `crick.healthcare.uiowa.edu` server and dump the latest version of HGMD in hg19 coordinates.

There is a section at the top of the code where you can edit these database credentials if you need to. You can also define which genome assembly (i.e. hg19 or hg38) you would like to get records for.

### Schema

A new manual is provided with every new release of Genome Trax and contains the details for the database schema. This repository contains a previous version of the manual for reference purposes, but it is recommended to use the manual that is provided with the release you are using.

The database contains a table called `ngs_feature`, which includes all of the annotation data for all of the data sources. Most of the annotation data is stored in the `description` column, much like the INFO column of a VCF file. Although it works similarly to the INFO column, pipes (`|`) are used instead of equals-signs (`=`), and the values may contain characters that are considered illegal in the an INFO column. See below for details. Additionally, the `description` column uses `N/A` for empty values, which the script will convert into the VCF-preferred period (`.`) for empty values.

## Handling illegal VCF characters

The VCF format restricts all spaces, semi-colons, and equals-signs from being present in an any INFO field. Additionally, commas are reserved for separating allele-specific values. Therefore, these characters are [URL-encoded](http://www.w3schools.com/tags/ref_urlencode.asp) before they are written to the output VCF.

| Illegal VCF character | URL-encoding |   
| --------------------- | ------------ |
| *space*               | %20          |
| ,                     | %2C          |
| ;                     | %3B          |
| =                     | %3D          |

## Author

Sean Ephraim
