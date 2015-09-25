#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')

F_IN = ARGV[0]
FILE_PREFIX = F_IN.sub(/#{File.extname(F_IN)}/, '')

cmd = Command.new

# Convert gene names to gene regions
cmd.genes2regions(genes_file: F_IN,
                  ref_file: CONFIG['gene_regions_ref_file'],
                  out_file_prefix: FILE_PREFIX)

# Get variants from gene regions
cmd.regions2variants(bed_file: cmd.genes2regions_result,
                     vcf_files: CONFIG['annotation_files'],
                     out_file_prefix: FILE_PREFIX)

# Add gene symbol to each record in VCF
cmd.addgenes(bed_file: cmd.genes2regions_result,
             vcf_file: cmd.regions2variants_result,
             out_file_prefix: FILE_PREFIX)

# Annotate with dbNSFP
cmd.addpredictions(dbnsfp_file: CONFIG['annotation_files']['dbnsfp'],
                   vcf_file: cmd.addgenes_result,
                   bed_file: cmd.genes2regions_result,
                   out_file_prefix: FILE_PREFIX,
                   clinical_labels: CONFIG['clinical_labels'])

# TODO Add HGVS notation

# TODO Add final pathogenicity
