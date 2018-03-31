# Create dbNSFP3 hg19 BCF file
#
# Convert a directory of dbNSFP3 variant annotation files into a 
# single BCF file (with .csi index file). Output will be written
# in the present working directory.
#
# Example usage:
#   ./create_dbNSFP3_hg19_BCF.sh /path/to/dbNSFP3.x/

DBNSFP_DIR=$(echo "$1" | sed 's/\/*$//') # e.g. /path/to/dbNSFP3.0 (all trailing slashes will be removed)
SOURCE=$(echo "$DBNSFP_DIR" | sed 's/.*\///g') # e.g dbNSFP3.0

if [ ! -f "$SOURCE.txt" ]; then
  echo "- Combining all variant files from $DBNSFP_DIR into $SOURCE.txt..."
  # Remove illegal/undesired characters
  # - Convert hashes and plus-signs to nothing
  # - Convert parentheses and hyphens to underscores
  # - Remove trailing underscores
  cat <(head -1 $DBNSFP_DIR/*_variant.chr1 | sed -e 's/[#+,=;]//g' -e 's/[()-]/_/g' -e 's/_\(\s\|$\)/\1/g') \
      <(tail -qn +2 $DBNSFP_DIR/*_variant.chr1 \
                    $DBNSFP_DIR/*_variant.chr2 \
                    $DBNSFP_DIR/*_variant.chr3 \
                    $DBNSFP_DIR/*_variant.chr4 \
                    $DBNSFP_DIR/*_variant.chr5 \
                    $DBNSFP_DIR/*_variant.chr6 \
                    $DBNSFP_DIR/*_variant.chr7 \
                    $DBNSFP_DIR/*_variant.chr8 \
                    $DBNSFP_DIR/*_variant.chr9 \
                    $DBNSFP_DIR/*_variant.chr10 \
                    $DBNSFP_DIR/*_variant.chr11 \
                    $DBNSFP_DIR/*_variant.chr12 \
                    $DBNSFP_DIR/*_variant.chr13 \
                    $DBNSFP_DIR/*_variant.chr14 \
                    $DBNSFP_DIR/*_variant.chr15 \
                    $DBNSFP_DIR/*_variant.chr16 \
                    $DBNSFP_DIR/*_variant.chr17 \
                    $DBNSFP_DIR/*_variant.chr18 \
                    $DBNSFP_DIR/*_variant.chr19 \
                    $DBNSFP_DIR/*_variant.chr20 \
                    $DBNSFP_DIR/*_variant.chr21 \
                    $DBNSFP_DIR/*_variant.chr22 \
                    $DBNSFP_DIR/*_variant.chrM \
                    $DBNSFP_DIR/*_variant.chrX \
                    $DBNSFP_DIR/*_variant.chrY) \
                    > $SOURCE.txt
  echo "   \--> Done! $SOURCE.txt has been created"
else
  echo "- $SOURCE.txt exists... Skipping combining of all files"
fi

# Create VCF (unsorted)
if [ ! -f "$SOURCE.unsorted.vcf" ]; then
  echo "- Creating VCF file ($SOURCE.unsorted.vcf)..."
  tab2vcf --chr 8 \
          --pos 9 \
          --ref 3 \
          --alt 4 \
          --source $SOURCE \
          --reference hg19 \
          --prefix DBNSFP_ \
          $SOURCE.txt \
          > $SOURCE.unsorted.vcf
  echo "   \--> Done! $SOURCE.unsorted.vcf has been created"
else
  echo "- $SOURCE.unsorted.vcf already exists... Skipping VCF creation"
fi

# Create sorted VCF and index it
if [ ! -f "$SOURCE.vcf.gz" ]; then
  echo "- Creating sorted VCF ($SOURCE.vcf.gz)..."
  vcf-sort -c "$SOURCE.unsorted.vcf" 2> /dev/null | bgzip -c > "$SOURCE.vcf.gz"
  echo "   \--> Done! $SOURCE.vcf.gz has been created"
  echo "- Creating index file for $SOURCE.vcf.gz ($SOURCE.vcf.gz.tbi)..."
  tabix -fp vcf "$SOURCE.vcf.gz"
  echo "   \--> Done! $SOURCE.vcf.gz.tbi has been created"
else
  echo "- $SOURCE.vcf.gz already exists... Skipping sorting"
fi

# Convert to BCF and index it
if [ ! -f "$SOURCE.bcf.gz" ]; then
  echo "- Converting to BCF file ($SOURCE.bcf.gz)..."
  bcftools view -O b -o $SOURCE.bcf.gz $SOURCE.vcf.gz
  echo "   \--> Done! $SOURCE.bcf.gz has been created"
  echo "- Creating index file for $SOURCE.bcf.gz ($SOURCE.bcf.gz.csi)..."
  bcftools index -fc $SOURCE.bcf.gz
  echo "   \--> Done! $SOURCE.bcf.gz.csi has been created"
else
  echo "- $SOURCE.bcf.gz already exists... Skipping BCF creation"
fi

# Done
echo "- Final output written to:"
echo "  * $SOURCE.bcf.gz"
echo "  * $SOURCE.bcf.gz.tbi"
echo "  * $SOURCE.vcf.gz       (delete if unneeded)"
echo "  * $SOURCE.vcf.gz.tbi   (delete if unneeded)"
echo "  * $SOURCE.unsorted.vcf (delete if unneeded)"
echo "  * $SOURCE.txt          (delete if unneeded)"
