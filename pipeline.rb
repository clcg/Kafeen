#!/usr/bin/env ruby

require_relative File.join('lib', 'boot.rb')

F_IN = ARGV[0]
FILE_PREFIX = F_IN.sub(/#{File.extname(F_IN)}/, '')

comm = Command.new

# Convert gene names to gene regions
comm.genes2regions(genes_file: F_IN, ref_file: CONFIG['data_files']['gene_regions'])
puts comm.genes2regions_result

# TODO Get variants from gene regions
                   
# TODO Add expert-curated variants that are
#      missing in the original list of variants
