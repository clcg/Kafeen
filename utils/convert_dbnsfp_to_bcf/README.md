# Convert dbNSFP v3.x to hg19 BCF

This script will combine all dbNSFP v3.x chromosome files into a single TSV file and then convert it to a single BCF file.

## Output

### Files

The following output files will be created:

- *dbNSFPv3.X.bcf.gz*
- *dbNSFPv3.X.bcf.gz.tbi*

Additionally, the following VCF files are also created and used as intermediate files to eventually create the final BCF file. These intermediate files are preserved, and it is up to the user whether or not to delete them:

- *dbNSFPv3.X.vcf.gz*
- *dbNSFPv3.X.vcf.gz.tbi*
- *dbNSFPv3.X.unsorted.vcf*

Lastly, a tab-separated file is also preserved, which is a concatenation of all the dbNSFP chromosome files. This is another intermediate file, and it is up to the user whether or not to delete it:

- *dbNSFP3.X.txt*

### INFO tags

The column headers of the dbNSFP chromosome files are used in order to create the INFO tags of the output BCF file. Take a look at the README file included with your dbNSFP download to see what columns exist. To create the INFO tag name, the original column name is (1.) converted to all-caps, (2.) prefixed with `DBNSFP_`, and (3.) sanitized for illegal tag name characters (more on this below).

As an example, the `SIFT_score` and `SIFT_pred` columns will become `DBNSFP_SIFT_SCORE` and `DBNSFP_SIFT_PRED` INFO tags, respectively.

As another example, the `GERP++_RS` and `GERP++_RS_rankscore` columns will become `DBNSFP_GERP_RS` and `DBNSFP_GERP_RS_RANKSCORE` INFO tags, respectively. Notice that the `+` character is removed in the INFO tags because BCFtools cannot properly query tags that contain this character.

If you are unsure of how the script converted some headers, you can view all INFO in the output BCF by running:

    bcftools view -h dbNSFP3.X.bcf.gz | grep '^##INFO='

## Usage

You will first need to download dbNSFP v3.x [here](https://sites.google.com/site/jpopgen/dbNSFP). Then run the conversion script on the entire directory (see below).

### Example

The conversion script requires only one argument, which is the path to the dbNSFP v3.x directory containing all the chromosome annotation files (e.g. *dbNSFP3.2a_variant.chr1*, *dbNSFP3.2a_variant.chr2*, etc.). For example, if you downloaded dbNSFP v3.2a to the directory path `/path/to/dbNSFPv3.2a/`, then you just need to pass this path to the conversion script: 

    ./convert_dbNSFP3_to_hg19_BCF.sh /path/to/dbNSFPv3.2a/
   
All output files will be written to the user's current directory and will be prefixed with `dbNSFPv3.2a`.

## Tips

### Get a job

This script will take a *very* long time to complete, and therefore it is highly recommended to submit it as a `qsub` job.

### Subset before querying

Due to the fact that dbNSFP is so large, simple queries can take a very long time to complete. In order to speed up queries, it is recommended to first subset the dbNSFP BCF to your regions of interest before performing queries. Use the `bcftools view -r ...` option to subset the BCF.

## A warning about intermediate files

This script creates several intermediate files in order to create the final BCF file. The user should be aware of the following in regards to these files:

1. dbNSFP is *very* large, and therefore the intermediate files are *very* large. You will need upwards of ~250 GB of free space (in addition to the size of dbNSFP itself) in order to run the script. The final BCF file will only be around ~20 GB.
1. If the script is interrupted, it will try to pick up where it left off when resubmitted. The script will automatically skip a step if the resulting intermediate file already exists. However, it is **very important** that all intermediate files be complete. If you know that an intermediate file is not complete, then delete it and let the script redo that step. For example, if your logs show that the script failed while trying to create *dbNSFPv3.X.vcf.gz* (perhaps because `vcf-sort` was not in your `$PATH`) then you should delete that file and rerun the script. The script will immediately attempt to create *dbNSFPv3.X.vcf.gz* again.

## Handling illegal VCF characters

The VCF format restricts all spaces, semi-colons, and equals-signs from being present in an any INFO field. Additionally, commas are reserved for separating allele-specific values. Therefore, these characters are [URL-encoded](http://www.w3schools.com/tags/ref_urlencode.asp) before they are written to the output VCF.

| Illegal VCF character | URL-encoding |
| --------------------- | ------------ |
| *space*               | %20          |
| ,                     | %2C          |
| ;                     | %3B          |
| =                     | %3D          |

Likewise, some additional characters are prohibited from being used in INFO tag names. If an INFO tag contains an illegal character, then BCFtools will not be able to query it. In short:

- Any of the following characters are removed from INFO tags: `#+,=;`
- Any of the following characters are replaced with underscores in INFO tags: `()-`

## Requirements

The following utilities must be in your `$PATH` in order for this script to work:

- `tab2vcf` (download it [here](https://github.com/sephraim/bin4matics))
- `bgzip` and `tabix` (download with the HTSlib package [here](https://github.com/samtools/htslib/releases/))
- `bcftools` (download it [here](https://github.com/samtools/bcftools/releases))
- `vcf-sort` (download it with the VCFtools package [here](https://vcftools.github.io/index.html))

## Author

Sean Ephraim
