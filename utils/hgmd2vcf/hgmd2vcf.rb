##
# HGMD2VCF
#
# Extract HGMD data from a Genome Trax MySQL database and
# put it into a VCF file.
#
# Edit the database login credentials below to determine
# which version of HGMD should be dumped.
# 
# This script requires the mysql2 gem. Install it with:
#   gem install mysql2
# 
# Example usage:
#   ruby hgmd2vcf.rb <OUTPUT_FILE>
#   ruby hgmd2vcf.rb hgmd.vcf
##

############### EDIT #################
#**** DATABASE LOGIN CREDENTIALS *****
  HOST = 'crick.healthcare.uiowa.edu'
  #DATABASE = 'genometrax2014r1'
  #DATABASE = 'genometrax2014r3'
  DATABASE = 'genometrax2015r2'
  USERNAME = 'hgmd'
  PASSWORD = '*gaattc#'
#********** GENOME ASSEMBLY **********
  BUILD = 'hg19'
  #BUILD = 'hg38'
#*********** DBSNP VERSION ***********
  #DBSNP_BUILD_ID = 138 # This is for the VCF meta-info (use for versions prior to 2015r2)
  DBSNP_BUILD_ID = 142  # This is for the VCF meta-info (use for version 2015r2)
######################################

require 'mysql2'
require 'uri'

##
# Get Field
#
# Get specific field from description column.
##
def get_field(field, description)
  result = description.match(/(^|;)#{field}\|([^;]*)/)[2]
  result = '.' if result == 'N/A'
  return result
end

OUT_FILE_NAME = ARGV[0]
F_OUT = File.open(OUT_FILE_NAME, 'w')

NGS_ONTOLOGY_NO = 6 # This is HGMD's ontology # in Genome Trax

puts "Connecting to database..."
CLIENT = Mysql2::Client.new(:host     => HOST,
                            :database => DATABASE,
                            :username => USERNAME,
                            :password => PASSWORD)

puts "Retrieving data..."
results = CLIENT.query("
  SELECT *
  FROM ngs_feature
  WHERE genome = '#{BUILD}'
  AND ngs_ontology_no = #{NGS_ONTOLOGY_NO}
")

results = results.to_a
if results.empty?
  abort("ERROR: No results found")
end

puts "Printing results..."
header = ""
results.each_with_index do |row, row_num|
  # Print header
  if row_num == 0
    # Begin adding meta-info
    TIME = Time.now.strftime("%Y%m%d")
    header = "##fileformat=VCFv4.2\n" +
             "##fileDate=#{TIME}\n" +
             "##source=#{DATABASE}\n" +
             "##dbSNP_BUILD_ID=#{DBSNP_BUILD_ID}\n" +
             "##reference=#{BUILD}\n"

    # Add INFO tags
    fields = row['description'].split(';')
    fields.each do |f|
      f.gsub!(/\|.*$/, '')
      header += "##INFO=<ID=HGMD_#{f.upcase},Number=.,Type=String,Description=\"NA\">\n"
    end

    # Add contigs
    ('1'..'22').to_a.push('X').push('Y').each {|chr| header += "##contig=<ID=#{chr}>\n" }

    # Add column header
    header += ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'].join("\t")
    header += "\n"
  end

  output = {}
  output['chr']    = row['chromosome'].sub(/^chr/, '')
  output['pos']    = row['feature_start']
  output['id']     = get_field('rsid', row['description'])
  output['ref']    = get_field('ref', row['description'])
  output['alt']    = get_field('alt', row['description'])
  output['qual']   = '.'
  output['filter'] = '.'
  # Format INFO field
  # - Make all INFO tags uppercase 
  # - Replace all 'N/A' with '.'
  output['info'] = row['description'].gsub(/([^;]*)\|/) {|s| "HGMD_#{$1.upcase}|"}.gsub(/([,|])N\/A([,;]|$)/, '\1.\2')

  # Remove illegal VCF INFO characters (i.e. spaces, equals-signs, and commas)
  # - Semi-colons can be left alone here, as they already delimit fields
  output['info'] = URI.escape(output['info'], ',= ')

  # Convert "key|value" notation to "key=value" notation
  # - Note that not encoding the string will result in 'invalid byte sequence' error
  # - Source: https://robots.thoughtbot.com/fight-back-utf-8-invalid-byte-sequences
  output['info'] = output['info'].encode('UTF-8', 'binary', invalid: :replace, undef: :replace, replace: '').gsub('|', '=')

  F_OUT.puts output.values.join("\t")
end
F_OUT.close

puts "Sorting chromosomal positions..."
TMP_FILE_NAME = ".#{OUT_FILE_NAME}.tmp"
`sort -k2,2n -k4,4d -k5,5d #{OUT_FILE_NAME} > #{TMP_FILE_NAME}`

puts "Beginning to print output..."
# Write header
File.open(OUT_FILE_NAME, 'w') {|f| f.write(header) }
# Print records
('1'..'22').to_a.push('X').push('Y').each do |chr|
  puts "- Chromosome #{chr}"
  `grep '^#{chr}\t' #{TMP_FILE_NAME} >> #{OUT_FILE_NAME}`
end
`rm -f #{TMP_FILE_NAME}`

begin 
  puts "Compressing output..."
  `bgzip -f #{OUT_FILE_NAME}`
  
  puts "Creating index..."
  `tabix -p vcf #{OUT_FILE_NAME}.gz`
  
  puts "Output file written to #{OUT_FILE_NAME}.gz"
  puts "Index file written to #{OUT_FILE_NAME}.gz.tbi"
rescue Exception => err
  puts "Skipping compression/indexing step - bgzip and/or tabix could not be found in your $PATH"
  puts "Output file written to #{OUT_FILE_NAME}"
end
