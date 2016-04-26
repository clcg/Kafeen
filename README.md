# Kafeen

******PLEASE NOTE THAT THIS README IS STILL VERY MUCH A WORK IN PROGRESS******

## VCF output fields
The following INFO fields can be found in the output VCF file.

### ExAC fields:

- **ExAC\_ALL\_AC** - Total number of alternate alleles in called genotypes
- **ExAC\_ALL\_AN** - Total number of alleles in called genotypes
- **ExAC\_ALL\_AF** - Allele frequency in the all populations calculated from AC and AN, in the range (0,1)
- **ExAC\_AFR\_AC** - African/African American alternate allele counts
- **ExAC\_AFR\_AN** - African/African American allele counts
- **ExAC\_AFR\_AF** - Allele frequency in the African/African American populations calculated from AC and AN, in the range (0,1)
- **ExAC\_AMR\_AC** - American alternate allele counts
- **ExAC\_AMR\_AN** - American allele counts
- **ExAC\_AMR\_AF** - Allele frequency in the American populations calculated from AC and AN, in the range (0,1)
- **ExAC\_EAS\_AC** - East Asian alternate allele counts
- **ExAC\_EAS\_AN** - East Asian allele counts
- **ExAC\_EAS\_AF** - Allele frequency in the East Asian populations calculated from AC and AN, in the range (0,1)
- **ExAC\_FIN\_AC** - Finnish alternate allele counts
- **ExAC\_FIN\_AN** - Finnish allele counts
- **ExAC\_FIN\_AF** - Allele frequency in the Finnish populations calculated from AC and AN, in the range (0,1)
- **ExAC\_NFE\_AC** - Non-Finnish European alternate allele counts
- **ExAC\_NFE\_AN** - Non-Finnish European allele counts
- **ExAC\_NFE\_AF** - Allele frequency in the Non-Finnish populations calculated from AC and AN, in the range (0,1)
- **ExAC\_OTH\_AC** - Other alternate allele counts
- **ExAC\_OTH\_AN** - Other allele counts
- **ExAC\_OTH\_AF** - Allele frequency in the Other populations calculated from AC and AN, in the range (0,1)
- **ExAC\_SAS\_AC** - South Asian alternate allele counts
- **ExAC\_SAS\_AN** - South Asian allele counts
- **ExAC\_SAS\_AF** - Allele frequency in the South Asian populations calculated from AC and AN, in the range (0,1)

### EVS fields

- **EVS\_ALL\_AC** - Total number of alternate alleles in called genotypes
- **EVS\_ALL\_AN** - Total number of alleles in called genotypes
- **EVS\_ALL\_AF** - Estimated allele frequency in all populations in the range (0,1)
- **EVS\_EA\_AC** - Number of alternate alleles in called genotypes in the European American population
- **EVS\_EA\_AN** - Number of alleles in called genotypes in the European American population
- **EVS\_EA\_AF** - Allele frequency in the European American populations calculated from AC and AN, in the range (0,1)
- **EVS\_AA\_AC** - Number of alternate alleles in called genotypes in the African American population
- **EVS\_AA\_AN** - Number of alleles in called genotypes in the African American population
- **EVS\_AA\_AF** - Allele frequency in the African American populations calculated from AC and AN, in the range (0,1)

### 1000 Genomes fields

- **1KG\_ALL\_AC** - Total number of alternate alleles in called genotypes
- **1KG\_ALL\_AN** - Total number of alleles in called genotypes
- **1KG\_ALL\_AF** - Estimated allele frequency in ALL populations in the range (0,1)
- **1KG\_AFR\_AC** - Number of alternate alleles in called genotypes in the AFR populations
- **1KG\_AFR\_AN** - Number of alleles in called genotypes in the AFR populations
- **1KG\_AFR\_AF** - Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)
- **1KG\_AMR\_AC** - Number of alternate alleles in called genotypes in the AMR populations
- **1KG\_AMR\_AN** - Number of alleles in called genotypes in the AMR populations
- **1KG\_AMR\_AF** - Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)
- **1KG\_EAS\_AC** - Number of alternate alleles in called genotypes in the EAS populations
- **1KG\_EAS\_AN** - Number of alleles in called genotypes in the EAS populations
- **1KG\_EAS\_AF** - Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)
- **1KG\_EUR\_AC** - Number of alternate alleles in called genotypes in the EUR populations
- **1KG\_EUR\_AN** - Number of alleles in called genotypes in the EUR populations
- **1KG\_EUR\_AF** - Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)
- **1KG\_SAS\_AC** - Number of alternate alleles in called genotypes in the SAS populations
- **1KG\_SAS\_AN** - Number of alleles in called genotypes in the SAS populations
- **1KG\_SAS\_AF** - Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)

## Under the hood

Merging data from multiple databases often requires subjecting the data to normalization and complex decision-making processes. This section is meant to expose the inner workings of the pipeline and reveal some of the complex decisions made in order to allow the data to be merged.

### Normalizing clinical significance nomenclature in HGMD

The clinical significance nomenclature used in ClinVar is very similar to that of this pipeline and is relatively easily to map. However HGMD's nomenclature model is very different. The following table shows how HGMD's nomenclature is mapped to the preferred nomenclature of this pipeline.

| HGMD clinical significance  | Meaning         | Normalized clinical significance | 
| --------------------------- | --------------- | -------------------------------- | 
| DM  | Disease-causing mutation                | Pathogenic                       |
| DM? | Variant is mentioned in the publication | Unknown significance             |
| DP  | Disease-associated polymorphism         | Benign                           |
| DFP | Disease-associated polymorphism with additional supporting functional evidence | Benign |
| FP  | In vitro/laboratory or in vivo functional polymorphism | Benign            |
| FTV | Frameshift / truncating variant         | Benign                           |
| CNV | Copy number variation                   | Benign                           |
| R   | Removed from HGMD                       | Unknown significance             |

### Handling confidence levels in HGMD

HGMD does not have an equivalent classification for ClinVar's "Likely pathogenic" or "Likely benign" categories. However HGMD supplies a confidence level (either "Low" or "High") with their clinical significance value. If variant's clinical significance has a "Low" confidence level, then the clinical significance value will automatically be prefixed with the word "Likely" to denote uncertainty. On the other hand if a variant has a "High" confidence level then nothing will change in the final output.

#### Example
If a variant has a clinical significance of "Benign" in HGMD and also has an associated "Low" confidence level, then it will automatically be converted to "Likely benign" in the final output. Conversely, if a variant has a clinical significance of "Benign" with an associated "High" confidence level then the clinical significance will remain as "Benign" in the final output.

### Handling discrepancies between ClinVar submissions

If multiple ClinVar submissions are available for a single variant, this pipeline will use the *most pathogenic* submission. It will also add a comment to notify the user that not all ClinVar submitters agree.

#### ClinVar clinical significance categories

The following are all possible classifications for clinical significance in ClinVar (order from *most pathogenic* to *least pathogenic*:

- Pathogenic
- Likely pathogenic
- Uncertain significance
- Likely benign
- Benign

#### Example
If a variant has the three ClinVar classifications "Benign", "Likely benign", and "Pathogenic", the most pathogenic classification will be used in the final output. In this case the *most pathogenic* classification is "Pathogenic" and will therefore be used in the final output.

### Handling discrepancies between ClinVar and HGMD

Once the *most pathogenic* value has been selected from ClinVar, the pathogenicities from ClinVar and HGMD are compared. There are several possible outcomes for this comparison, depending on the situation. The possible situations are as follows:

1. **Both sources totally agree**
  - This is the simplest case in which no extra interpretation is necessary. The final output will simply be whatever both sources say.
1. **Both sources totally disagree**
  - This happens when one source says "Pathogenic" (or "Likely pathogenic"), and the other source says "Benign" (or "Likely benign"). If this happens then the variant is automatically listed with "Unknown significance" in the final output. A comment is also made to notify the user of what the original sources say.
1. **Both sources mostly agree**
  - This happens when (a.) one source says "Pathogenic" and the other says "Likely pathogenic", or (b.) one source says "Benign" and the other says "Likely benign". If this happens then the variant will be classified as the lesser of the two (i.e. "Likely pathogenic" or "Likely benign").
1. **One source is certain, the other classifies as VUS**
  - This happens when one source calls a variant "Pathogenic", "Likely pathogenic", "Benign" or "Likely benign", and the other source calls it a VUS. In this case the variant will use the non-VUS call in the final output. This is the same as when a one source contains a variant and the other does not.

### Labeling a variant as "Benign*" (i.e. "benign star")

An asterisk (*) can be added to a benign label when the following two criteria are met:

1. The variant is benign because it has a MAF that is greater than or equal to the maximum allowed MAF.
1. The variant is currently reported as pathogenic in ClinVar, HGMD, or both. In ClinVar, it must be reported as "Pathogenic", not "Likely pathogenic". In HGMD, it must be reported as "DM" (i.e. pathogenic) with "High" confidence (not "Low" confidence).

If the first criteria is met but the second is not, the variant will simply be labeled benign. If the "benign star" option is disabled, then variants that meet these two criteria will simply be labled benign (with no star).
