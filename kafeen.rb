#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')
require 'logger'
require 'open3'

cmd = Command.new(log_level: 'info')
log = Logger.new(STDOUT)

# Validate config file
annotator = CONFIG['third_party']['annotator'].downcase
valid_annotators = ['asap', 'vep', 'both']
if not valid_annotators.include?(annotator)
  log.error("Invalid annotator: #{annotator}. Valid options: #{valid_annotators.join(', ')}")
  abort
end

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

# Add gene symbol to each record in VCF
cmd.add_genes(bed_file: bed_file,
              vcf_file: vcf_file,
              out_file_prefix: FILE_PREFIX)

# Add HGVS notation (using ASAP and/or VEP, as specified in config)
if ['asap', 'both'].include?(annotator)
  cmd.add_asap(vcf_file: cmd.add_genes_result,
               out_file_prefix: FILE_PREFIX,
               asap_path: CONFIG['third_party']['asap']['path'],
               ref_flat: CONFIG['third_party']['asap']['ref_flat'],
               ref_seq_ali: CONFIG['third_party']['asap']['ref_seq_ali'],
               fasta: CONFIG['third_party']['asap']['fasta'])
end
if ['vep', 'both'].include?(annotator)
  cmd.add_vep(vcf_file: cmd.add_genes_result,
             out_file_prefix: FILE_PREFIX,
             vep_path: CONFIG['third_party']['vep']['path'],
             vep_cache_path: CONFIG['third_party']['vep']['cache_path'],
             vep_config_path: CONFIG['third_party']['vep']['config_path'])
end

include_dbnsfp = CONFIG['annotation_files']['dbnsfp']['include']
if include_dbnsfp == false
  # Exclude dbNSFP
  log.info("dbNSFP inclusion set to false in config... Skipping dbNSFP annotation")

  # Calculate prediction totals without adding dbNSFP
  cmd.add_predictions(annotator: annotator,
                      vcf_file: cmd.add_genes_result,
                      out_file_prefix: FILE_PREFIX,
                      clinical_labels: CONFIG['clinical_labels'])
else
  # Include dbNSFP
  if include_dbnsfp == true
    log.info("dbNSFP inclusion set to true in config... Including dbNSFP annotation")
  elsif include_dbnsfp.nil?
    log.warning("dbNSFP inclusion not specified in config... Including dbNSFP annotation by default")
  else
    log.warning("dbNSFP inclusion expected true / false in config, but got '#{include_dbnsfp}' instead... Including dbNSFP annotation by default")
  end   

  # Annotate with dbNSFP and calculate prediction totals
  cmd.add_predictions(annotator: annotator,
                      dbnsfp_file: CONFIG['annotation_files']['dbnsfp'],
                      bed_file: merged_bed_file,
                      vcf_file: cmd.add_genes_result,
                      out_file_prefix: FILE_PREFIX,
                      clinical_labels: CONFIG['clinical_labels'])
end

# Add final pathogenicity
cmd.finalize_pathogenicity(vcf_file: cmd.add_predictions_result,
                           out_file_prefix: FILE_PREFIX,
                           clinical_labels: CONFIG['clinical_labels'],
                           enable_benign_star: CONFIG['enable_benign_star'])

# Clean up meta-info in final VCF
cmd.cleanup_meta_info(vcf_file: cmd.finalize_pathogenicity_result,
                      out_file_prefix: FILE_PREFIX)
