# A library of available pipeline commands
#
# @author Sean Ephraim

require 'logger'
require 'open3'

class Command
  # Result file paths
  attr_reader :genes2regions_result
  attr_reader :regions2variants_result
  attr_reader :addgenes_result
  attr_reader :addpredictions_result

  @@log = Logger.new(STDERR)
  @@log.level = Logger::INFO

  ##
  # Genes 2 Regions
  ##
  def genes2regions(genes_file:, ref_file:, out_file_prefix:)
    # Gene region reference file
    # Reference columns:
    #   chr, start, stop, gene_symbol
    f_regions = File.open(ref_file, 'r')

    # Set output file
    @genes2regions_result = "#{out_file_prefix}.gene_regions.bed"
    f_out = File.open(@genes2regions_result, 'w')
    
    File.open(genes_file, 'r').each_line do |gene|
      gene.chomp!
    
      # Get gene region
      result = f_regions.grep(/([^a-zA-Z0-9-]|^)#{gene}([^a-zA-Z0-9-]|$)/)
    
      # Print result
      if !result.empty?
        f_out.puts result
      end
    
      f_regions.rewind # reset file pointer
    end
    f_regions.close
    f_out.close
    @@log.info("Gene regions written to #{@genes2regions_result}")
  end

  ##
  # Get a list of all variants within specified regions
  ##
  def regions2variants(bed_file:, vcf_files:, out_file_prefix:, keep_tmp_files: false)
    tmp_vcfs = {}
    File.open(bed_file).each do |region|
      chr,pos_start,pos_end,gene = region.chomp.split("\t")
      chr.sub!('chr', '')

      # Query all VCF files for variants
      vcf_files.each do |key, vcf|
        next if key == 'dbnsfp' # DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
        tmp_source_vcf = "#{out_file_prefix}.#{vcf['source']}.tmp.vcf.gz"
        if vcf['fields'].nil?
          # Remove all INFO tags
          fields = 'INFO'
        else
          # Keep only the following INFO tags (indicated by ^)
          fields = "^" + vcf['fields'].map { |f| "INFO/#{f}" }.join(',')
        end

        # Query...
        @@log.info("Querying #{vcf['source']}...")
        stdout, stderr = Open3.capture3(
          "bcftools annotate \
             --remove '#{fields}' \
             --regions-file '#{bed_file}' \
             --exclude 'TYPE=\"other\"' \
             --output #{tmp_source_vcf} \
             --output-type z \
             #{vcf['filename']}"
        )

        # Did bcftools return an error?
        # Try again and don't remove any INFO tags this time
        if !stderr.empty?
          @@log.warn("bcftools returned an error for #{vcf['source']}. Trying another query method...")
          stdout, stderr = Open3.capture3(
            "bcftools annotate \
               --regions-file '#{bed_file}' \
               --exclude 'TYPE=\"other\"' \
               --output #{tmp_source_vcf} \
               --output-type z \
               #{vcf['filename']}"
          )
        end

        # Index the results file...
        if !stderr.empty?
          # ERROR
          @@log.error("bcftools was not able to query #{vcf['source']}. Please check that file name and INFO tags are set correctly in your config file.")
        else
          # SUCCESS -- now create index file
          @@log.info("Successfully queried #{vcf['source']}")
          @@log.info("Creating index file for #{vcf['source']}...")
          `bcftools index --force --tbi #{tmp_source_vcf}`

          # Store tmp file name (filename) and the original VCF that the data came from (parent)
          tmp_vcfs[key] = {'filename' => tmp_source_vcf, 'parent' => vcf['filename']}
        end
      end
    end

    # Construct list of VCFs to merge
    files_to_merge = []
    tmp_vcfs.each do |key, tmp_vcf|
      next if key == 'dbnsfp' # DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
      files_to_merge << tmp_vcf['filename']
    end
    files_to_merge << tmp_vcfs['dbsnp']['filename']

    # Merge VCFs...
    @regions2variants_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Merging results...")
    `bcftools merge \
       --merge all \
       --output #{@regions2variants_result} \
       --output-type z \
       #{files_to_merge.join(' ')}`
    @@log.info("Done merging results")

    @@log.info("Creating new index file for #{@regions2variants_result}...")
    `bcftools index --force --tbi #{@regions2variants_result}`
    @@log.info("Done creating index file")

    # Remove tmp files
    if !keep_tmp_files
      @@log.info("Removing temp files...")
      tmp_vcfs.each do |key, tmp_vcf|
        File.unlink(tmp_vcf['filename']) if File.exist?(tmp_vcf['filename'])
        File.unlink("#{tmp_vcf['filename']}.tbi") if File.exist?("#{tmp_vcf['filename']}.tbi")
      end
      @@log.info("Done removing temp files")
    end
  end

  ##
  # Take genes from BED file and add to VCF file
  ##
  def addgenes(bed_file:, vcf_file:, out_file_prefix:)
    # Prepare header file
    header_file = "#{out_file_prefix}.header.tmp.txt"
    header_line = '##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC gene symbol">'
    File.open(header_file, 'w') {|f| f.write(header_line) }

    # Prepare BED file using bgzip and tabix
    `bgzip -c #{bed_file} > #{bed_file}.tmp.gz`
    `tabix -fp bed #{bed_file}.tmp.gz`

    # Add genes to VCF file
    tmp_output = "#{out_file_prefix}.tmp.vcf.gz"
    @@log.info("Adding gene annotations to #{tmp_output}")
    `bcftools annotate \
       --annotations #{bed_file}.tmp.gz \
       --columns CHROM,FROM,TO,GENE \
       --header-lines #{header_file} \
       --output #{tmp_output} \
       --output-type z \
       #{vcf_file}`
    @@log.info("Genes added to #{tmp_output}")

    # Move tmp output to be the new output file
    @addgenes_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Moving output to #{@addgenes_result}...")
    File.rename(tmp_output, @addgenes_result)

    # Create index
    @@log.info("Creating index for #{@addgenes_result}...")
    `bcftools index --force --tbi #{@addgenes_result}`
    @@log.info("Done creating index file")

    # Remove tmp files
    @@log.info("Removing all tmp files...")
    File.unlink(header_file) if File.exist?(header_file)
    File.unlink("#{bed_file}.tmp.gz") if File.exist?("#{bed_file}.tmp.gz")
    File.unlink("#{bed_file}.tmp.gz.tbi") if File.exist?("#{bed_file}.tmp.gz.tbi")
  end

  ##
  # Add predictions from dbNSFP
  ##
  def addpredictions(dbnsfp_file:, vcf_file:, bed_file:, out_file_prefix:)
    # Get only regions of interest from dbNSFP
    @@log.info("Subsetting dbNSFP for faster annotation...")
    dbnsfp_subset_file = "#{out_file_prefix}.dbNSFP_subset.tmp.bcf.gz"
    `bcftools view \
       --regions-file #{bed_file} \
       --output-type b \
       --output-file #{dbnsfp_subset_file} \
       #{dbnsfp_file['filename']}`
    @@log.info("dbNSFP subset written to #{dbnsfp_subset_file}")
    
    # Create index
    @@log.info("Creating index file for #{dbnsfp_subset_file}...")
    `bcftools index --force --csi #{dbnsfp_subset_file}`
    @@log.info("Done creating index file")

    # TODO Add dbNSFP predictions
    tmp_output = "#{out_file_prefix}.tmp.vcf.gz"
    @@log.info("Adding dbNSFP pedictions to #{tmp_output}...")
    `bcftools annotate \
       --annotations #{dbnsfp_subset_file} \
       --columns #{dbnsfp_file['fields'].map { |f| "INFO/#{f}" }.join(',')} \
       --output #{tmp_output} \
       --output-type z \
       #{vcf_file}`
#       .each_line do |vcf_row|
#         if vcf_row.match(/^#/)
#           # TODO Print header
#         else
#           total_num_preds = 0
#           num_path_preds = 0
#           dbnsfp_file['fields'].each do |field|
#             # Skip field if: 
#             # 1. It's not a prediction
#             # 2. The prediction is unknown
#             next if !field.match(/_PRED$/i)
#             next if preds.all? {|pred| pred == '.' || pred == 'U'}
#
#             preds = vcf_row.match(/[\t;]#{field}=([^;\t]*)/)[1].split(';')
#             if field.match(/^SIFT_PRED$/) && preds.include('D')
#               # SIFT prediction
#             elsif field.match(/^POLYPHEN2_PRED$/)
#               # TODO Analyze all other predictions too
#             end
#           end
#         end
#       end

    # TODO Add *nominal* predictions for convervation scores (GERP++ and phyloP)
    # TODO Add final prediction

    #@addpredictions_result = "#{out_file_prefix}.vcf.gz"
    #@@log.info("Predictions added to #{@addpredictions_result}")
  end
end
