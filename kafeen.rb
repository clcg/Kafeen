#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')
require 'logger'
require 'open3'

cmd = Command.new(log_level: 'info')
log = Logger.new(STDOUT)

if TEST_MODE
  # Skip collecting variants if test mode not enabled...
  # Set VCF and BED files that will be used for subsequent steps
  
  # Compress VCF if not compressed already
  if F_IN.match(/.+\.vcf$/i)
    stdout, stderr = Open3.capture3("bcftools view -O z -o '#{F_IN}.gz' '#{F_IN}'")
    unless stderr.empty?
      log.error("bcftools was not able to compress #{F_IN}...") 
      log.error("bcftools error is: #{stderr}") 
      abort
    end
  end

  # Create new index for VCF
  stdout, stderr = Open3.capture3("bcftools index --force --tbi '#{F_IN}.gz'")
  unless stderr.empty?
    log.error("bcftools was not able to index #{F_IN}...") 
    log.error("bcftools error is: #{stderr}") 
    abort
  end

  vcf_file = F_IN
  bed_file = F_BED
  merged_bed_file = F_BED
else
  # Convert gene names to gene regions
  cmd.genes2regions(genes_file: F_IN,
                    ref_file: CONFIG['gene_regions_ref_file'],
                    out_file_prefix: FILE_PREFIX)
  
  # Get variants from gene regions
  cmd.regions2variants(bed_file: cmd.genes2regions_merged_result,
                       vcf_files: CONFIG['annotation_files'],
                       out_file_prefix: FILE_PREFIX)

  # Set VCF and BED files that will be used for subsequent steps
  vcf_file = cmd.regions2variants_result
  bed_file = cmd.genes2regions_result
  merged_bed_file = cmd.genes2regions_merged_result
end

# Add gene symbol to each record in VCF
cmd.add_genes(bed_file: bed_file,
              vcf_file: vcf_file,
              out_file_prefix: FILE_PREFIX)

# Annotate with dbNSFP
cmd.add_predictions(dbnsfp_file: CONFIG['annotation_files']['dbnsfp'],
                    vcf_file: cmd.add_genes_result,
                    bed_file: merged_bed_file,
                    out_file_prefix: FILE_PREFIX,
                    clinical_labels: CONFIG['clinical_labels'])

# Add HGVS notation (using ASAP)
cmd.add_asap(vcf_file: cmd.add_predictions_result,
             out_file_prefix: FILE_PREFIX,
             asap_path: CONFIG['third_party']['asap']['path'],
             ref_flat: CONFIG['third_party']['asap']['ref_flat'],
             ref_seq_ali: CONFIG['third_party']['asap']['ref_seq_ali'],
             fasta: CONFIG['third_party']['asap']['fasta'])

# Add final pathogenicity
cmd.finalize_pathogenicity(vcf_file: cmd.add_predictions_result,
                           out_file_prefix: FILE_PREFIX,
                           clinical_labels: CONFIG['clinical_labels'],
                           enable_benign_star: CONFIG['enable_benign_star'])

# TODO Re-header

# Run tests
if TEST_MODE
  cmd.test(vcf_file: cmd.add_predictions_result,
           assertion_tags: CONFIG['test']['assertion_tags'],
           out_file_prefix: FILE_PREFIX)
end
