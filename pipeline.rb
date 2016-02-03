#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')

cmd = Command.new(log_level: 'info')

# Convert gene names to gene regions
cmd.genes2regions(genes_file: F_IN,
                  ref_file: CONFIG['gene_regions_ref_file'],
                  out_file_prefix: FILE_PREFIX)

# Get variants from gene regions
cmd.regions2variants(bed_file: cmd.genes2regions_result,
                     vcf_files: CONFIG['annotation_files'],
                     out_file_prefix: FILE_PREFIX)

# Add gene symbol to each record in VCF
cmd.add_genes(bed_file: cmd.genes2regions_result,
             vcf_file: cmd.regions2variants_result,
             out_file_prefix: FILE_PREFIX)

# Annotate with dbNSFP
cmd.add_predictions(dbnsfp_file: CONFIG['annotation_files']['dbnsfp'],
                   vcf_file: cmd.add_genes_result,
                   bed_file: cmd.genes2regions_result,
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
