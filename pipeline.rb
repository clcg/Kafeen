#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')

F_IN = ARGV[0]
FILE_PREFIX = F_IN.sub(/#{File.extname(F_IN)}/, '')

comm = Command.new

# Convert gene names to gene regions
comm.genes2regions(genes_file: F_IN, ref_file: CONFIG['gene_regions_ref_file'], out_file_prefix: FILE_PREFIX)

# TODO Get variants from gene regions
comm.regions2variants(bed_file: comm.genes2regions_result, vcf_files: CONFIG['annotation_files'], out_file_prefix: FILE_PREFIX)
                   
# TODO Add expert-curated variants that are
#      missing in the original list of variants
