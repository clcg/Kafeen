# A library of available pipeline commands
#
# @author Sean Ephraim and Rob Marini

require 'logger'
require 'open3'
require 'uri'
require 'fileutils'
require 'find'
require_relative 'core_extensions'
require_relative 'logic_tracker'

# Monkey-patch the String class
String.include CoreExtensions::String::Vcf
String.include CoreExtensions::String::Colorize

class Command
  # Result file paths
  attr_reader :genes2regions_result
  attr_reader :genes2regions_merged_result
  attr_reader :genes2regions_result_intmed
  attr_reader :genes2regions_merged_result_intmed
  attr_reader :regions2variants_result
  attr_reader :regions2variants_result_intmed
  attr_reader :add_genes_result
  attr_reader :add_genes_result_intmed
  attr_reader :add_predictions_result
  attr_reader :add_predictions_intmed
  attr_reader :add_asap_result
  attr_reader :add_asap_result_intmed
  attr_reader :add_vep_result
  attr_reader :add_vep_result_intmed
  attr_reader :finalize_pathogenicity_result
  attr_reader :finalize_pathogenicity_intmed

  ##
  # Initialize
  ##
  def initialize(out_file_prefix: "NA_PREFIX", log_level: 'info', log_out: STDOUT)
    # Set logger
    @@log = Logger.new(log_out)
    if log_level.upcase == 'UNKNOWN'
      @@log.level = Logger::UNKNOWN
    elsif log_level.upcase == 'FATAL'
      @@log.level = Logger::FATAL
    elsif log_level.upcase == 'ERROR'
      @@log.level = Logger::ERROR
    elsif log_level.upcase == 'WARN'
      @@log.level = Logger::WARN
    elsif log_level.upcase == 'INFO'
      @@log.level = Logger::INFO
    elsif log_level.upcase == 'DEBUG'
      @@log.level = Logger::DEBUG
      @@log.debug("Debugging messages enabled")
    else
      @@log.level = Logger::DEBUG
      @@log.error("#{log_level} is not a valid log level")
      @@log.warn("Debugging messages enabled by default")
    end
    
    @tracker_delim = "|" #reused later when creating clean/new LogicTracker objects
    @out_file_prefix = out_file_prefix
    @classify_trackers_list = @out_file_prefix + ".classify_track_list.tsv" #output file for writing a list of all tracks possible for decoding classify trails during review
    @classify_predictions = LogicTracker.new(@tracker_delim)
    
  end

  ##
  # Genes to Regions
  #
  # Get genomic regions for each HGNC gene symbol
  ##
  def genes2regions(genes_file:, ref_file:, out_file_prefix:, intmed_output: false)
    # Gene region reference file
    # Reference columns:
    #   chr, start, stop, gene_symbol
    f_regions = File.open(ref_file)

    # Set output file
    @genes2regions_result = "#{out_file_prefix}.gene_regions.bed"
    f_out = File.open(@genes2regions_result, 'w')

    File.open(genes_file).each_line do |gene|
      gene.chomp!
      @@log.debug("Retrieving gene region for #{gene}...")
    
      # Get gene region
      result = f_regions.grep(/([^a-zA-Z0-9-]|^)#{Regexp.escape(URI.escape(gene, ';,= '))}([^a-zA-Z0-9-]|$)/)
    
      # Print result
      if !result.empty?
        f_out.puts result
      else
        @@log.error("No matching HGNC symbol for #{gene}")
      end
    
      f_regions.rewind # reset file pointer
    end
    f_regions.close
    f_out.close

    # Create a temporary merged regions file
    @@log.info("Merging regions to use for querying...")
    @genes2regions_merged_result = "#{out_file_prefix}.gene_regions.merged.bed"
    tmp_sorted_bed = "#{out_file_prefix}.gene_regions.sorted.tmp.bed"
    `sort -k1,1 -k2,2n #{@genes2regions_result} > #{tmp_sorted_bed}`
    `bedtools merge -i #{tmp_sorted_bed} > #{@genes2regions_merged_result}`

    # Overwrite original .bed file with sorted .bed file
    
    File.rename(tmp_sorted_bed, @genes2regions_result)

    @@log.info("Gene regions written to #{@genes2regions_result}")
    @@log.info("Merged regions written to #{@genes2regions_merged_result}")
    
    if intmed_output
      @genes2regions_merged_result_intmed = "#{out_file_prefix}.intmed-gene_regions.merged.bed"
      @genes2regions_result_intmed = "#{out_file_prefix}.intmed-gene_regions.result.bed"
      FileUtils.copy_file(@genes2regions_merged_result,@genes2regions_merged_result_intmed)
      @@log.info("An intermediate file of genes2regions written to #{@genes2regions_merged_result_intmed}")
      FileUtils.copy_file(@genes2regions_result,@genes2regions_result_intmed)
      @@log.info("An intermediate file of genes2regions written to #{@genes2regions_result_intmed}")
    end
    
  end

  ##
  # Regions to Variants
  #
  # Get a list of all variants within specified regions
  ##
  def regions2variants(bed_file:, vcf_files:, out_file_prefix:, keep_tmp_files: false, intmed_output: false)
    
    # Note: Custom VCF tags dependent on the vcf sources are added to output in this function
    #  For more information, see the _add_max_mafs function
    
    if intmed_output
      keep_tmp_files = true;
    end
    
    tmp_vcfs = {}
    # Query all VCF files for variants
    vcf_files.each do |key, vcf|
      if vcf['include'] == false # RJM: config.yml vcf source flag check end
          @@log.info("Excluded #{vcf['source']} explicitly in config.yml. Skipping #{vcf['source']} annotation source")
      else # RJM: checks flag in config.yml if the vcf source should be included or not in the kafeen run. true allows the source to be considered as normal, false ignores the source
        if vcf['include'] == true
          @@log.info("Included #{vcf['source']} explicitly by a valid include tag")
        elsif !vcf.has_key?('include')
          @@log.warning("Included #{vcf['source']} implicitly. #{vcf['source']} did not have an include tag in config.yml")
        else
          @@log.warning("Included #{vcf['source']} implicitly... Incompatible include input within config.yml. Expected true or false, config.yml provided: (no value).\n\tPlease review the config.yml and indicated whether or not you would like to include #{vcf['source']} (include: true) or not (include: false)")
        end
        next if key == 'dbnsfp' # DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
        tmp_source_vcf = "#{out_file_prefix}.#{vcf['source']}.tmp.vcf.gz"
        if vcf['fields'].nil?
          # Remove all INFO tags
          fields = 'INFO'
        elsif vcf['fields'] == '*'
          # Keep all INFO tags
          fields = '*'
        else
          # Keep only the following INFO tags (indicated by ^)
          fields = "^" + vcf['fields'].map { |f| "INFO/#{f}" }.join(',')
        end
  
        # Remove ID column as well (except for dbSNP)
        if key != 'dbsnp'
          if fields != '*'
            fields = "ID,#{fields}"
          else
            fields = "ID"
          end
        end
  
        # Query...
        @@log.info("Querying #{vcf['source']}...")
        if fields != '*'
          stdout, stderr = Open3.capture3(
            "bcftools annotate \
               --remove '#{fields}' \
               --regions-file '#{bed_file}' \
               --include 'REF ~ \"^[ACGT]\\+$\" && ALT ~ \"^[ACGT,]\\+$\"' \
               --output-type z \
               --output #{tmp_source_vcf} \
               #{vcf['filename']}"
          )
  
          stderr = stderr.gsub(/^Lines.*$/, '').strip # remove unnecessary "error"
          if !stderr.empty?
            @@log.debug("bcftools returned an error for #{vcf['source']}. Trying another query method...")
            @@log.debug("bcftools error is: #{stderr}")
          end
        end
  
        # Did bcftools return an error?
        # Try again and don't remove any INFO tags this time
        if fields == '*' || !stderr.empty?
          stdout, stderr = Open3.capture3(
            "bcftools view \
               --regions-file '#{bed_file}' \
               --include 'REF ~ \"^[ACGT]\\+$\" && ALT ~ \"^[ACGT,]\\+$\"' \
               --output-type z \
               --output-file #{tmp_source_vcf} \
               #{vcf['filename']}"
          )
          stderr = stderr.gsub(/^Lines.*$/, '').strip # remove unnecessary "error"
        end
        
        # Index the results file...
        if !stderr.empty?
          # ERROR
          @@log.error("bcftools was not able to query #{vcf['source']}. Please check that file name and INFO tags are set correctly in your config file.")
          @@log.debug("bcftools error is: #{stderr}")
        else
          # SUCCESS -- create index file
          @@log.info("Successfully queried #{vcf['source']}")
          @@log.info("Creating index file for #{vcf['source']}...")
          `bcftools index --force --tbi #{tmp_source_vcf}`
          @@log.info("Done creating index file")
        end
        
        # get <source>_MAX_MAF and <source>_MAX_MAF_SOURCE
        if !stderr.empty?
          # ERROR
          @@log.error("bcftools was not able to query #{vcf['source']}. Please check that file name and INFO tags are set correctly in your config file.")
          @@log.error("not able to add MAX_MAF information for this source.")
          @@log.debug("bcftools error is: #{stderr}")
        else
          # SUCCESS
          # Store tmp file name (filename) and the original VCF that the data came from (parent)
          max_mafs_output_vcf = "#{out_file_prefix}.#{vcf['source']}.tmp.max_mafs.vcf.gz"
          tmp_source_vcf_max_mafs = _add_max_mafs(tmp_source_vcf, max_mafs_output_vcf, tag_prefix = vcf['field_prefix'])
          tmp_vcfs[tag_prefix] = {'filename' => tmp_source_vcf_max_mafs, 'parent' => vcf['filename']}
        end
          
      end # RJM: config.yml vcf source flag check end
    end # do end

    # Construct list of VCFs to merge
    files_to_merge = []
    tmp_vcfs.each do |key, tmp_vcf|
      next if key == 'dbnsfp' # DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
      files_to_merge << tmp_vcf['filename']
    end
    
    @regions2variants_result = "#{out_file_prefix}.vcf.gz"
    tmp_merged = "#{out_file_prefix}.merged.tmp.vcf.gz"
    if files_to_merge.length == 1
      # Skip merging if only 1 VCF provided
      File.rename(files_to_merge.join, tmp_merged)
    else
      # Merge VCFs if more than 1 provided
      @@log.info("Merging results...")
      `bcftools merge \
         --merge none \
         --output-type z \
         --output #{tmp_merged} \
         #{files_to_merge.join(' ')}`
      @@log.info("Done merging results")
    end
    
    # Create index
    @@log.info("Creating index file for #{tmp_merged}...")
    `bcftools index --force --tbi #{tmp_merged}`
    @@log.info("Done creating index file")
    
    # get overall MAX_MAF and MAX_MAF_SOURCE
    @regions2variants_result = _add_max_mafs(tmp_merged, @regions2variants_result)

    # Remove tmp files
    if !keep_tmp_files
      @@log.info("Removing temp files...")
      tmp_vcfs.each do |key, tmp_vcf|
        File.unlink(tmp_vcf['filename']) if File.exist?(tmp_vcf['filename'])
        File.unlink("#{tmp_vcf['filename']}.tbi") if File.exist?("#{tmp_vcf['filename']}.tbi")
      end
      
      File.unlink(tmp_merged) if File.exist?(tmp_merged)
      File.unlink("#{tmp_merged}.tbi") if File.exist?("#{tmp_merged}.tbi")
      
      #gather all other *.tmp* files in output directory
      tmp_file_paths = []
      Find.find(File.dirname(out_file_prefix)) do |path|
        tmp_file_paths << path if path.include?".tmp"
      end
      
      # delete found *.tmp* files in output directory
      tmp_file_paths.each do |tmp_file|
        File.unlink(tmp_file) if File.exist?(tmp_file)
      end
      
      @@log.info("Done removing temp files")
    end
    
    if intmed_output
      @regions2variants_result_intmed = "#{out_file_prefix}.intmed-regions2variants.vcf.gz"
      FileUtils.copy_file(@regions2variants_result,@regions2variants_result_intmed)
      @regions2variants_result_intmed_tbi = "#{@regions2variants_result_intmed}.tbi"
      FileUtils.copy_file("#{@regions2variants_result}.tbi",@regions2variants_result_intmed_tbi)
      @@log.info("An intermediate file of regions2variants written to #{@regions2variants_result_intmed} and its .tbi")
    end
    
  end

  ##
  # Add Max MAF information to vcf
  #
  # expects the out_file to have the .vcf.gz tail to it
  ##
  def _add_max_mafs(tmp_source_vcf, out_file, tag_prefix = false)
    
    # Note: Custom VCF tags dependent on the vcf sources are added to output in this function
    #  These tags follow the following convention: <source>_MAX_MAF, <source>_MAX_MAF_SOURCE
    #  Additionally, 2 tags that describe the overall max maf (MAX_MAF) and its source (MAX_MAF_SOURCE) are added to output
    # Set custom MAX MAF VCF tags to be added to output

    @max_maf_tag = "MAX_MAF"
    @max_maf_source_tag = "MAX_MAF_SOURCE"
    
    if(!tag_prefix)
      tag_prefix = "OVERALL"
    end
    
    source_max_maf_tag = "#{tag_prefix}_#{@max_maf_tag}" # INFO tag for this source's max maf
    source_max_maf_source_tag = "#{tag_prefix}_#{@max_maf_source_tag}" # INFO tag for this source's max maf's source (aka. pop)
    
    tmp_source_vcf_max_mafs = out_file.gsub(".gz","")
    File.open(tmp_source_vcf_max_mafs, 'w') do |f|
      @@log.info("Adding #{tag_prefix}'s max maf and source to #{tmp_source_vcf}")
      `bcftools view \
        --output-type v \
        #{tmp_source_vcf}`
      .each_line do |vcf_row|
        vcf_row.chomp!
        if vcf_row.match(/^##/)
          #Print meta-info
          f.puts vcf_row
        elsif vcf_row.match(/^#[^#]/)
          # Add new <source>_MAX_MAF and and <source>_MAX_MAF_SOURCE tags to meta-info
          f.puts "##INFO=<ID=#{source_max_maf_tag},Number=.,Type=Float,Description=\"Max MAF for #{tag_prefix}, in range (0,1)\">"
          f.puts "##INFO=<ID=#{source_max_maf_source_tag},Number=.,Type=String,Description=\"Max MAF population for #{tag_prefix}\">"
          
          # Print header
          f.puts vcf_row
        else
          # data row
          
          vcf_cols = vcf_row.split("\t")
          @@log.debug("Processing: #{vcf_cols[0]}\t#{vcf_cols[1]}\t#{vcf_cols[3]}\t#{vcf_cols[4]}")
          
          # Initialize max maf fields
          max_maf = {}
          max_maf[:maf] = ''
          max_maf[:source] = ''
            
          vcf_infos = vcf_cols[7].split(";")
          vcf_infos.each do |info_tag|
            if (info_tag.include? "_AF=")
              tag_val = info_tag.split("=")
              if(max_maf[:maf].empty?)
                # first AF found in tags and is set to max maf
                max_maf[:source] = tag_val[0]
                max_maf[:maf] = tag_val[1]
              elsif(tag_val[1].to_f > max_maf[:maf].to_f)
                # new max maf, update both :maf and :source
                max_maf[:source] = tag_val[0]
                max_maf[:maf] = tag_val[1]
              elsif (tag_val[1].to_f == max_maf[:maf].to_f)
                # equal max maf, update only the :source
                max_maf[:source] << ",#{tag_val[0]}"
              end # if(max_maf[:maf].empty?)
            end # if (info_tag.include? "AF")
            
          end # vcf_infos.each do |info_tag|
          
          # Update INFO column
          if(max_maf[:maf].empty?)
            max_maf[:maf] = '.'
          end
          
          if(max_maf[:source].empty?)
            max_maf[:source] = '.'
          end
          
          if(vcf_cols[7] != '.') #only adds the MAF_MAF tags if the row had info fields (was an issue with dbSNP
            vcf_cols[7] = [
              vcf_cols[7],
                "#{source_max_maf_tag}=#{max_maf[:maf]}",
                "#{source_max_maf_source_tag}=#{max_maf[:source]}"
                ].map { |e| e.force_encoding('UTF-8') } .join(";")
          end
          
          # Print updated VCF row
          f.puts vcf_cols.join("\t")
          
          @@log.debug("- This variant's #{tag_prefix}'s max maf is #{max_maf[:maf]} from #{max_maf[:source]}")
          @@log.debug("------------------------------------------------------")
          
        end # vcf_row.match(/^##/)
      end #.each_line do |vcf_row|
    end #File.open(vcf['filename'], 'w') do |f|
    
    @@log.info("Compressing #{tmp_source_vcf_max_mafs}")
    `bcftools view \
      --output-type z \
      --output-file #{out_file} \
      #{tmp_source_vcf_max_mafs}
    `
    @@log.info("Compressed output written to #{out_file}")
    
    @@log.info("Creating index file for #{out_file}...")
    `bcftools index --force --tbi #{out_file}`
    @@log.info("Done creating index file")
    
    File.unlink(tmp_source_vcf_max_mafs) if File.exist?(tmp_source_vcf_max_mafs)
    
    @@log.info("Done adding max maf and source to #{tmp_source_vcf}")
    return out_file
    
  end
  
  ##
  # Take genes from BED file and add to VCF file
  #
  #
  #
  ##
  def add_genes(bed_file:, vcf_file:, out_file_prefix:, intmed_output: false)
    tmp_output_file = "#{out_file_prefix}.tmp.vcf.gz"
    @@log.info("Adding gene annotations to #{tmp_output_file}")
    # Prepare header file
    header_file = "#{out_file_prefix}.header.tmp.txt"
    header_line = '##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC gene symbol">'
    File.open(header_file, 'w') {|f| f.write(header_line) }

    # Prepare BED file using bgzip and tabix
    @@log.debug("Creating tmp compressed BED file #{tmp_output_file}...")
    `sort -u -k1,1 -k2,2n -k3,3n #{bed_file} | bgzip -c > #{bed_file}.tmp.gz`
    @@log.debug("Compressed tmp BED file written to #{tmp_output_file}.tmp.gz...")
    @@log.debug("Creating tmp index file for #{tmp_output_file}.tmp.gz...")
    `tabix -fp bed #{bed_file}.tmp.gz`
    @@log.debug("Done creating tmp index file...")

    # Add genes to VCF file
    `bcftools annotate \
       --annotations #{bed_file}.tmp.gz \
       --columns CHROM,FROM,TO,GENE \
       --header-lines #{header_file} \
       --output #{tmp_output_file} \
       --output-type z \
       #{vcf_file}`
    @@log.info("Genes added to #{tmp_output_file}")

    # Move tmp output to be the new output file
    @add_genes_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Moving output to #{@add_genes_result}...")
    File.rename(tmp_output_file, @add_genes_result)

    # Create index
    @@log.info("Creating index for #{@add_genes_result}...")
    `bcftools index --force --tbi #{@add_genes_result}`
    @@log.info("Done creating index file")

    # Remove tmp files
    @@log.info("Removing all tmp files...")
    File.unlink(header_file) if File.exist?(header_file)
    File.unlink("#{bed_file}.tmp.gz") if File.exist?("#{bed_file}.tmp.gz")
    File.unlink("#{bed_file}.tmp.gz.tbi") if File.exist?("#{bed_file}.tmp.gz.tbi")
    
    if intmed_output
      @add_genes_result_intmed = "#{out_file_prefix}.intmed-add_genes.vcf.gz"
      FileUtils.copy_file(@add_genes_result,@add_genes_result_intmed)
      @add_genes_result_intmed_tbi = "#{@add_genes_result_inmed}.tbi"
      FileUtils.copy_file("#{@add_genes_result}.tbi",@add_genes_result_intmed_tbi)
      @@log.info("An intermediate file of add_genes written to #{@add_genes_result_intmed} and its .tbi")
    end
    
  end

  ##
  # Add Predictions
  #
  # Add predictions from dbNSFP, and tally up prediction totals from dbNSFP & VEP.
  # Omitting dbNSFP from input will tally up totals from VEP only.
  ##
  def add_predictions(annotator:, vcf_file:, out_file_prefix:, clinical_labels:, dbnsfp_file: nil, bed_file: nil, intmed_output: false)
    # Set custom VCF tags to be added to output
    @dbnsfp_gerp_pred_tag = "DBNSFP_GERP_PRED"
    @dbnsfp_phylop20way_mammalian_pred_tag = "DBNSFP_PHYLOP20WAY_MAMMALIAN_PRED"
    @num_path_preds_tag = "NUM_PATH_PREDS"
    @total_num_preds_tag = "TOTAL_NUM_PREDS"
    @final_pred_tag = "FINAL_PRED"
    @classify_track_part_tag = "CLASSIFY_TRAIL_PARTIAL"
    
    # Annotate with dbNSFP
    if !dbnsfp_file.nil?
      if bed_file.nil?
        # ERROR: Missing BED file for dbNSFP subsetting
        @@log.error("Failed to add dbNSFP because BED file was not supplied to add_predictions()")
      elsif !File.file?(bed_file)
        # ERROR: BED file supplied for dbNSFP subsetting does not exist
        @@log.error("Failed to add dbNSFP because BED file '#{bed_file}' does not exist")
      else
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

        # Add dbNSFP predictions
        tmp_output_file = "#{out_file_prefix}.tmp.bcf"
        @@log.info("Adding dbNSFP predictions to #{tmp_output_file}...")
        `bcftools annotate \
           --annotations #{dbnsfp_subset_file} \
           --columns #{dbnsfp_file['fields'].map { |f| "INFO/#{f}" }.join(',')} \
           --output-type u \
           --output #{tmp_output_file} \
           #{vcf_file}`

        # Change working VCF file to this file
        File.rename(tmp_output_file, vcf_file)

        @@log.info("Removing tmp dbNSFP files...")
        File.unlink(dbnsfp_subset_file) if File.exist?(dbnsfp_subset_file)
        File.unlink("#{dbnsfp_subset_file}.csi") if File.exist?("#{dbnsfp_subset_file}.csi")
        @@log.info("Done removing tmp dbNSFP files")
      end
    end

    # Select all *_PRED tags (as well as GERP++ and phyloP tags) from VCF header
    pred_fields = `bcftools view -h #{vcf_file}`.scan(/^##INFO=<ID=([^,]*)/).flatten
                   .select { |e| e.match(/(?:_PRED$|^DBNSFP_GERP_RS$|^DBNSFP_PHYLOP20WAY_MAMMALIAN$)/i) }

    # Tally up prediction totals
    tmp_output_file = "#{out_file_prefix}.tmp.vcf"
    File.open(tmp_output_file, 'w') do |f|
      `bcftools view \
         --output-type v \
         #{vcf_file}`
      .each_line do |vcf_row|
           vcf_row.chomp!
           #vcf header
           if vcf_row.match(/^##/)
             # Print meta-info
             f.puts vcf_row
           elsif vcf_row.match(/^#[^#]/)
             # Add custom INFO tags
             f.puts "##INFO=<ID=#{@num_path_preds_tag},Number=.,Type=String,Description=\"Number of pathogenic predictions from dbNSFP\">"
             f.puts "##INFO=<ID=#{@total_num_preds_tag},Number=.,Type=String,Description=\"Total number of prediction scores available from dbNSFP\">"
             f.puts "##INFO=<ID=#{@final_pred_tag},Number=.,Type=String,Description=\"Final prediction consensus based on majority vote of prediction scores\">"
             f.puts "##INFO=<ID=#{@classify_track_part_tag},Number=.,Type=String,Description=\"Classify tracking partial string. To be read in and used in other functions. \">"
             
             # Add tags for any dbNSFP predictions computed on-the-fly
             if !dbnsfp_file.nil?
               # Add GERP++ prediction tag to meta-info
               if dbnsfp_file['fields'].any?{ |e| e == 'DBNSFP_GERP_RS' }
                 f.puts "##INFO=<ID=#{@dbnsfp_gerp_pred_tag},Number=.,Type=String,Description=\"NA\">"
               end
               # Add phyloP20way mammalian prediction tag to meta-info
               if dbnsfp_file['fields'].any?{ |e| e == 'DBNSFP_PHYLOP20WAY_MAMMALIAN' }
                 f.puts "##INFO=<ID=#{@dbnsfp_phylop20way_mammalian_pred_tag},Number=.,Type=String,Description=\"NA\">"
               end
             end

             # Print header (i.e. "CHROM  POS  ID ...")
             f.puts vcf_row
           #vcf data
           else           
             vcf_cols = vcf_row.split("\t")
  
             # Analyze each *_PRED field (as well as GERP++ and phyloP)
             # Tally up pathogenic predictions
             output = {}
             output[:total_num_preds] = 0
             output[:num_path_preds] = 0
             output[:classify_track_part] = "";

             #create classification tracker partial path and create a clean LogicTracker instance
             @classify_predictions.clear_trail() #start with a clean classification trail
               
             pred_fields.each do |field|
               # Get all predictions for this algorithm
               match = vcf_row.get_vcf_field(field)
  
               # No data for this algorithm -- skip it
               next if match.empty?
  
               # Get all predictions for this algorithm
               preds = match.split(/[^a-zA-Z0-9.-]+/)
  
               if field == 'DBNSFP_SIFT_PRED' and ['vep','both'].include?(annotator)
#                 vep_sift_pred = vcf_row.get_vcf_field(@vep_sift_pred_tag) # 06/24/19 - RJM
                 vep_sift_pred = vcf_row.get_vcf_field(@vep_tags['vep_sift_pred']['vcf_tag'])
                 if !vep_sift_pred.empty?
                    preds = vep_sift_pred.split(/[^a-zA-Z0-9.-]+/)
                 end
               elsif field == 'DBNSFP_POLYPHEN2_HDIV_PRED' and ['vep','both'].include?(annotator)
#                 vep_poly_pred = vcf_row.get_vcf_field(@vep_polyphen_pred_tag) # 06/24/19 - RJM
                 vep_poly_pred = vcf_row.get_vcf_field(@vep_tags['vep_polyphen_pred']['vcf_tag'])
                 if !vep_poly_pred.empty?
                    preds = vep_poly_pred.split(/[^a-zA-Z0-9.-]+/)
                 end
               end
  
               # No data for this algorithm -- skip it
               next if preds.all? { |pred| pred == '.' || pred == 'U' }
                 
               if field == 'DBNSFP_SIFT_PRED'
                 # SIFT prediction
                 # See if there is a VEP override
                 output[:num_path_preds] += 1 if preds.include?('D') || preds.include?('D-LC') # <-- "Damaging"
                 output[:total_num_preds] += 1
               elsif field == 'DBNSFP_POLYPHEN2_HDIV_PRED'
                 # Polyphen2 (HDIV) prediction
                 # See if there is a VEP override
                 output[:num_path_preds] += 1 if preds.include?('D') || preds.include?('P') # <-- "Deleterious" or "Possibly damaging"
                 output[:total_num_preds] += 1
               elsif field == 'DBNSFP_LRT_PRED'
                 # LRT prediction
                 output[:num_path_preds] += 1 if preds.include?('D') # <-- "Deleterious"
                 output[:total_num_preds] += 1
               elsif field == 'DBNSFP_MUTATIONTASTER_PRED'
                 # MutationTaster prediction
                 output[:num_path_preds] += 1 if preds.include?('D') || preds.include?('A') # <-- "Disease-causing" or "Disease-causing (automatic)"
                 output[:total_num_preds] += 1
               elsif field == 'DBNSFP_GERP_RS'
                 # GERP++ prediction
                 if preds.any? { |pred| pred.to_f > 0.0 }
                   # Conserved
                   output[:num_path_preds] += 1
                   vcf_cols[7] = [vcf_cols[7], "#{@dbnsfp_gerp_pred_tag}=C"].join(";")
                 else
                   # Non-conserved
                   vcf_cols[7] = [vcf_cols[7], "#{@dbnsfp_gerp_pred_tag}=N"].join(";")
                 end
                 output[:total_num_preds] += 1
               elsif field == 'DBNSFP_PHYLOP20WAY_MAMMALIAN'
                 # phyloP20way mammalian prediction
                 if preds.any? { |pred| pred.to_f >= 0.95 }
                   # Conserved
                   output[:num_path_preds] += 1
                   vcf_cols[7] = [vcf_cols[7], "#{@dbnsfp_phylop20way_mammalian_pred_tag}=C"].join(";")
                 else
                   # Non-conserved
                   vcf_cols[7] = [vcf_cols[7], "#{@dbnsfp_phylop20way_mammalian_pred_tag}=N"].join(";")
                 end
                 output[:total_num_preds] += 1
               end
             end
  
             # Add final prediction
             if output[:total_num_preds] == 0
               # No predictions available
               @classify_predictions.add_track("PredNA","No predictions available.")
               output[:final_pred] = '.'
               output[:num_path_preds] = '.'
               output[:total_num_preds] = '.'
             elsif output[:total_num_preds] >= 5
               @classify_predictions.add_track("PredAvail","Predictions available.")
               @classify_predictions.add_track("5+preds","5 or greater predictions available.")
               path_score = output[:num_path_preds].to_f/output[:total_num_preds].to_f
               if path_score >= 0.6
                 # Predicted pathogenic
                 @classify_predictions.add_track("path_score_patho>=.6","Pathogenicity score >= 0.6 --> user mapping for predicted pathogenic = #{clinical_labels['pred_pathogenic']}") #needs testing
                 output[:final_pred] = clinical_labels['pred_pathogenic']
               elsif path_score <= 0.4
                 # Predicted benign
                 @classify_predictions.add_track("path_score_patho<.6","Pathogenicity score < 0.6.")
                 @classify_predictions.add_track("path_score_benign","Pathogenicity score <= 0.4 --> user mapping for predicted benign = #{clinical_labels['pred_pathogenic']}") #needs testing
                 output[:final_pred] = clinical_labels['pred_benign']
               else
                 # Predicted unknown (benign predictions approx. equal to pathogenic)
                 @classify_predictions.add_track("path_score_patho<.6","Pathogenicity score < 0.6.")
                 @classify_predictions.add_track("path_score_patho>.4","Pathogenicity score > 0.4.")
                 @classify_predictions.add_track("path_score_VUS","number of benign predictions approx. equal to pathogenic --> user mapping for predicted VUS = #{clinical_labels['unknown']}") #needs testing
                 output[:final_pred] = clinical_labels['unknown']
               end
             else
               # Predicted unknown (not enough predictions)
               @classify_predictions.add_track("PredAvail","Predictions available.")
               @classify_predictions.add_track("NotEnough","Not enough predictions to make a decision. --> user mapping for predicted VUS = #{clinical_labels['unknown']}")
               output[:final_pred] = clinical_labels['unknown']
             end
             
             output[:classify_track_part] = @classify_predictions.get_trail() #output now retains the partial compacted trail string delimited by the @tracker_delim
             @classify_predictions.clear_trail() #clear classification trail
               
             # Remove illegal characters and set all empty values to '.'
             output.each do |k, v|
               next if output[k] == '.'
               if v.to_s.strip.empty?
                 output[k] = '.'
               else
                 output[k] = URI.escape(output[k].to_s, ';,= ')
               end
             end
  
             # Update INFO column
             if output[:total_num_preds] != 0 && output[:total_num_preds] != '.'
               vcf_cols[7] = [
                 vcf_cols[7], 
                 "#{@num_path_preds_tag}=#{output[:num_path_preds]}",
                 "#{@total_num_preds_tag}=#{output[:total_num_preds]}", 
                 "#{@final_pred_tag}=#{output[:final_pred]}"
               ].join(";")
             end
             
             #Update INFO column with a tracker
             vcf_cols[7] = [
                 vcf_cols[7], 
                 "#{@classify_track_part_tag}=#{output[:classify_track_part]}"
               ].join(";")
  
             # Print updated VCF row
             f.puts vcf_cols.join("\t")
           end # <-- end else
      end # <-- End .each_line do |vcf_row|
    end # <-- End File.open(tmp_output_file, ...)
    
    @@log.info("Total prediction counts added to #{tmp_output_file}")

    @add_predictions_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Compressing #{tmp_output_file}")
    # Compress the output file
    `bcftools view \
       --output-type z \
       --output-file #{@add_predictions_result} \
       #{tmp_output_file}`
    @@log.info("Compressed output written to #{@add_predictions_result}...")

    # Index output file
    @@log.info("Indexing #{@add_predictions_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@add_predictions_result}`
    @@log.info("Done creating index file")

    @@log.info("Removing tmp files...")
    File.unlink(tmp_output_file) if File.exist?(tmp_output_file)
    @@log.info("Done removing tmp files")
    
    if intmed_output
      @add_predictions_result_intmed = "#{out_file_prefix}.intmed-add_predictions.vcf.gz"
      FileUtils.copy_file(@add_predictions_result,@add_predictions_result_intmed)
      @add_predictions_result_intmed_tbi = "#{@add_predictions_result_inmed}.tbi"
      FileUtils.copy_file("#{@add_predictions_result}.tbi",@add_predictions_result_intmed_tbi)
      @@log.info("An intermediate file of add_predictions written to #{@add_predictions_result_intmed} and its .tbi")
    end
    
  end

  def _get_asap_variant(chr, pos, ref, alt)
    pos = pos.to_i
    if ref == '-' || alt == '-'
      # Already in correct format
      # -- DO NOTHING
    elsif ref.length == 1 && alt.length == 1
      # Substitution (already in correct format)
      # -- DO NOTHING
    elsif ref.length == 1 && alt.length > 1
      # INSERTION
      # Remove redundant nucleotides from beginning of sequences
      ref.length.times do |i|
        if alt.match(/^#{ref[0..-1-i]}/)
          alt = alt.sub(/^#{ref[0..-1-i]}/, '')
          ref = ref.sub(/^#{ref[0..-1-i]}/, '')
          break
        end
      end
      ref = '-' if ref.length == 1
      pos += 1
    elsif ref.length > 1 && alt.length == 1
      # DELETION
      # Remove redundant nucleotides from beginning of sequences
      alt.length.times do |i|
        if ref.match(/^#{alt[0..-1-i]}/)
          ref = ref.sub(/^#{alt[0..-1-i]}/, '')
          alt = alt.sub(/^#{alt[0..-1-i]}/, '')
          break
        end
      end
      alt = '-' if alt.length == 1
      pos += 1
    else
      # DELETION/INSERTION or no match
      # -- DO NOTHING
    end
    ref = '-' if ref.empty?
    alt = '-' if alt.empty?
    return "chr#{chr}:#{pos}:#{ref}>#{alt}"
  end

  ##
  # Add ASAP
  #
  ##
  
  def add_asap(vcf_file:, out_file_prefix:, asap_path:, ref_flat:, ref_seq_ali:, fasta:, intmed_output: false)
    # Set INFO tags
    @asap_variant_tag = 'ASAP_VARIANT'
    @asap_hgvs_c_tag = 'ASAP_HGVS_C'
    @asap_hgvs_p_tag = 'ASAP_HGVS_P'
    @asap_locale_tag = 'ASAP_LOCALE'
    @asap_function_tag = 'ASAP_FUNCTION'

    # Print ASAP input file
    tmp_asap_input  = "#{out_file_prefix}.asap.in.tmp.txt"
    File.open(tmp_asap_input, 'w') do |f|
      @@log.info("Preparing ASAP input file...")
      `bcftools query \
         --format '%CHROM\\t%POS\\t%REF\\t%ALT\\n' \
         #{vcf_file}`
      .each_line do |vcf_row|
        chr,pos,ref,alt = vcf_row.chomp.split("\t")
        f.puts [_get_asap_variant(chr, pos, ref, alt), chr, pos, ref, alt].join("\t")
      end
    end
    # ^ End ASAP conversion

    # Run ASAP
    tmp_asap_output = "#{out_file_prefix}.asap.out.tmp.txt"
    asap_error_log = "#{out_file_prefix}.asap.errors.log"
    validate = true
    @@log.info("Running ASAP...")
    File.open(tmp_asap_output, 'w') do |f|
      `java -Xmx2048m -jar #{asap_path} \
         --refFlat #{ref_flat} \
         --refSeqAli #{ref_seq_ali} \
         --fasta #{fasta} \
         --positions-file #{tmp_asap_input} 2> #{asap_error_log}`
      .each_line do |line|
        line.chomp!

        # Validate first line to see if ASAP is running okay
        if validate
          validate = false
          if line == "Use: manditory*"
            break
          end
        end

        # Parse ASAP output
        fields = line.split("\t", -1)

        # Print: CHROM, POS, REF, ALT, ASAP_variant, HGVS_c, HGVS_p, locale, impact
        if !line.include?("ERROR_")
          f.puts [fields[1..4], fields[0], fields[6..7], fields[10], fields[12]].flatten.join("\t")
        else
          # If error, only print: CHROM, POS, REF, ALT, ASAP_variant
          f.puts [fields[1..4], fields[0], '', '', '', ''].flatten.join("\t")
        end
      end
    end # <-- End writing output file

    # Check output
    if File.zero?(asap_error_log)
      # No errors occurred
      @@log.info("ASAP output written to #{tmp_asap_output}")
      File.unlink(asap_error_log) if File.exist?(asap_error_log)

      # Compress output
      @@log.info("Compressing #{tmp_asap_output}...")
      `bgzip -f #{tmp_asap_output}`
      @@log.info("Compressed output written to #{tmp_asap_output}")

      # Index output
      @@log.info("Indexing #{tmp_asap_output}.gz...")
      `tabix -f -s1 -b2 -e2 #{tmp_asap_output}.gz`
      @@log.info("Done creating index file")

      # Create VCF header file
      @@log.info("Creating VCF header file...")
      tmp_vcf_header = "#{out_file_prefix}.header.tmp.txt"
      header = [
        "##INFO=<ID=#{@asap_variant_tag},Number=1,Type=String,Description=\"ASAP-style variant notation\">",
        "##INFO=<ID=#{@asap_hgvs_c_tag},Number=1,Type=String,Description=\"HGVS cDNA annotation according to ASAP\">",
        "##INFO=<ID=#{@asap_hgvs_p_tag},Number=1,Type=String,Description=\"HGVS protein annotation according to ASAP\">",
        "##INFO=<ID=#{@asap_locale_tag},Number=1,Type=String,Description=\"Variant locale (e.g. exon, intron, etc.) according to ASAP\">",
        "##INFO=<ID=#{@asap_function_tag},Number=1,Type=String,Description=\"Variant function (e.g. missense, frameshift, etc.) according to ASAP\">",
      ].join("\n") + "\n"
      File.open(tmp_vcf_header, 'w') { |f| f.write(header) }
      @@log.info("VCF header file written to #{tmp_vcf_header}")

      # Annotate VCF with ASAP output
      @add_asap_result = "#{out_file_prefix}.vcf.gz"
      tmp_output_vcf = "#{out_file_prefix}.tmp.vcf.gz"
      @@log.info("Annotating VCF with ASAP output...")
      `bcftools annotate \
         --annotations #{tmp_asap_output}.gz \
         --columns CHROM,POS,REF,ALT,#{@asap_variant_tag},#{@asap_hgvs_c_tag},#{@asap_hgvs_p_tag},#{@asap_locale_tag},#{@asap_function_tag} \
         --header-lines #{tmp_vcf_header} \
         --output #{tmp_output_vcf} \
         --output-type z \
         #{vcf_file}`

      # Rename tmp output to final output
      File.rename(tmp_output_vcf, @add_asap_result)
      @@log.info("Output written to #{@add_asap_result}")

      # Index output file
      @@log.info("Indexing #{@add_asap_result}...")
      `bcftools index  \
         --force \
         --tbi \
         #{@add_asap_result}`
      @@log.info("Done creating index file")
    else
      # Errors occurred... skipping annotation
      @@log.error("There was a problem running ASAP -- refer to #{asap_error_log}")
    end

    # Remove tmp files
    @@log.info("Removing tmp files...")
    File.unlink(tmp_asap_input) if File.exist?(tmp_asap_input)
    File.unlink("#{tmp_asap_output}.gz") if File.exist?("#{tmp_asap_output}.gz")
    File.unlink("#{tmp_asap_output}.gz.tbi") if File.exist?("#{tmp_asap_output}.gz.tbi")
    File.unlink(tmp_vcf_header) if File.exist?(tmp_vcf_header)
    File.unlink(tmp_output_vcf) if File.exist?(tmp_output_vcf)
    @@log.info("Done removing tmp files")
    
    if intmed_output
      @add_asap_result_intmed = "#{out_file_prefix}.intmed-add_asap.vcf.gz"
      FileUtils.copy_file(@add_asap_result,@add_asap_result_intmed)
      @add_asap_result_intmed_tbi = "#{@add_asap_result_inmed}.tbi"
      FileUtils.copy_file("#{@add_asap_result}.tbi",@add_asap_result_intmed_tbi)
      @@log.info("An intermediate file of add_asap written to #{@add_asap_result_intmed} and its .tbi")
    end
    
  end

  ##
  # Add VEP
  #
  ##

  def add_vep(vcf_file:, out_file_prefix:, vep_path:, vep_cache_path:, vep_config_path:, vep_tags:, intmed_output: false)

#    # Tag names as they will appear in the end DVD. (moved to yaml config)
    # 06/11/19 - RJM:
    # removing hardcoding of variable, set in config with boot.rb required fields

    @add_vep_result = "#{out_file_prefix}.vcf.gz"
    @vep_tags = vep_tags

    # Run VEP
    # Make tsv containing relevant fields in vep output,
    # to annotate kafeen vcf with using bcftools
    tmp_vep_output = "#{out_file_prefix}.vep.out.tmp.txt"
    vep_stats_log = "#{out_file_prefix}.vep.stats.html"
    vep_error_log = "#{out_file_prefix}.vep.errors.log"
    vep_header = nil
    @@log.info("Running VEP...")

    # Header Indicies for fields used in scoring or additional processing from VEP.  
    # These should be set once from the VEP header line.
    
    # 06/11/19 - RJM:
    # removing hardcoding of variable, set in config with boot.rb required fields
    csq_indices = {}
    vep_tags_ol = []
    vep_tags.keys.each do |tag|
      vep_tags_ol.push(tag)
      unless vep_tags[tag]['csq_head'].nil? #equivalent to if not true
        csq_indices[tag] = -1
      end
    end

    # Intersting thread - it seemed like VEP was waiting to finish before STDOUT processing began
    # The delay is not bad for everyday runtime, but makes debugging a pain.
    # http://stackoverflow.com/questions/1154846/continuously-read-from-stdout-of-external-process-in-ruby
    # Good flowchart: http://stackoverflow.com/questions/6338908/ruby-difference-between-exec-system-and-x-or-backticks
    
    File.open(tmp_vep_output, 'w') do |f|
      IO.popen("#{vep_path} -config #{vep_config_path} -i #{vcf_file} -stats_file #{vep_stats_log} -o STDOUT 2> #{vep_error_log}") do |vep|
        vep.each do |line|
            line.chomp!
            
            begin
            
            # 06/11/19 - RJM:
            # removing hardcoding of variable, set in config with boot.rb required fields
            pick_csq = ""
            csq_picks = nil
              
            pick_score = -1

            if line.start_with?("#")
              if line.include? "<ID=CSQ"
                vep_header = line.split("Format: ")[1].split("\">")[0].split("|")
                
              # 06/11/19 - RJM:
              # removing hardcoding of variable, set in config with boot.rb required fields
                vep_tags.keys.each do |tag|
                  unless vep_tags[tag]['csq_head'].nil?
                    csq_indices[tag] = vep_header.index(vep_tags[tag]['csq_head'])
                  end
                end
              end
              next
            end

            all_features = []

            fields = line.split("\t")
            fields.each_with_index do |field, i|
              if field.include? "CSQ="
                  
                csq = field.split("CSQ=")[1]
                csq.split(",").each do |transcript|
                  vep_vals = {};
                  score = 0
                  # Init variables to avoid nil errors
                  
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  vep_tags.keys.each do |tag|
                    vep_vals[tag] = ""
                  end

                  # Must use -1 limit or trailing empty fields will be excluded
                  vep_fields = transcript.split("|", -1)
                  vep_tags_ol.each do |tag|
                    if ! csq_indices[tag].nil?
                      vep_vals[tag] = vep_fields[csq_indices[tag]]
                    end
                  end
                  # specifically looking at feature vs HGVS nomenclature because
                  # the HGVS nomenclature may not be set for upstream/downstream
                  # variants or other secnearios.
                  
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if vep_vals['vep_feature'].start_with?("X")
                       score = score + 100000
                  elsif vep_vals['vep_feature'].start_with?("N")
                       score = score + 300000
                  else
                       score = score + 200000
                  end
                  
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if vep_vals['vep_impact'] == "HIGH"
                       score = score + 40000
                  elsif vep_vals['vep_impact'] == "MODERATE"
                       score = score + 30000
                  elsif vep_vals['vep_impact'] == "LOW"
                       score = score + 20000
                  elsif vep_vals['vep_impact'] == "MODIFIER"
                       score = score + 10000
                  end
                  vep_vals['vep_hgvs_c'] = vep_fields[csq_indices['vep_hgvs_c']]
                  vep_vals['vep_hgvs_p'] = vep_fields[csq_indices['vep_hgvs_p']]
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if !(vep_vals['vep_hgvs_p'].to_s == '')
                       score = score + 2000
                  elsif !(vep_vals['vep_hgvs_c'].to_s == '')
                       score = score + 1000
                  end
                  # 08/07/17 - BC:
                  # Update SIFT/PolyPhen scoring from raw scores to prediction
                  
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  vep_sift_raw = vep_fields[csq_indices['vep_sift_score']]
                  if !vep_sift_raw.nil? && !vep_sift_raw.empty?
                        score = score + 100
                        vep_vals['vep_sift_score'] = vep_sift_raw.gsub(/[{()}]/, " ").split(" ")[1]
                        vep_sift_pred_parts = vep_sift_raw.split("(")[0].split("_")
                        vep_vals['vep_sift_pred'] = vep_sift_pred_parts[0].upcase[0,1]
                        if vep_sift_pred_parts.length > 1
                            vep_sift_confidence = vep_sift_pred_parts[1].upcase[0,1] + "C"
                            vep_vals['vep_sift_pred'] = vep_vals['vep_sift_pred'] + "-" + vep_sift_confidence
                        end
                        if vep_vals['vep_sift_pred'] == "D" || vep_vals['vep_sift_pred'] == "D-LC" 
                          score = score + 100
                        end
                  end
                  
                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  vep_polyphen_raw = vep_fields[csq_indices['vep_polyphen_score']]
                  if !vep_polyphen_raw.nil? && !vep_polyphen_raw.empty?
                      score = score + 100
                      vep_vals['vep_polyphen_score'] = vep_polyphen_raw.gsub(/[{()}]/, " ").split(" ")[1]
                      vep_polyphen_pred_parts = vep_polyphen_raw.split("(")[0].split("_")
                      # 08/31/17 - BC:
                      # Match VEP PolyPhen prediction on full VEP text
                      # (probably_damaging, possibly_damaging, benign, unknown)
                      if vep_polyphen_pred_parts[0] == "probably"
                          vep_vals['vep_polyphen_pred'] = "D"
                          score = score + 100
                      elsif vep_polyphen_pred_parts[0] == "possibly"
                          vep_vals['vep_polyphen_pred'] = "P"
                      else
                          # Matches on benign and unknown
                          vep_vals['vep_polyphen_pred'] = vep_polyphen_pred_parts[0].upcase[0,1]
                      end
                  end

                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if !vep_vals['vep_canonical'].empty?
                      score = score + 10
                  end

                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if !vep_vals['vep_pick'].empty?
                      score = score+1
                  end

                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  if score > pick_score
                    pick_csq = transcript
                    pick_score = score
                    csq_picks = vep_vals
                  elsif csq_picks.nil?
                    # init incase first entry has a score of 0
                    pick_csq = transcript
                    pick_score = score
                    csq_picks = vep_vals
                  end

                  # 06/11/19 - RJM:
                  # removing hardcoding of variable, set in config with boot.rb required fields
                  feature_info_arr = []
                  vep_tags_ol.each do |vep_val_tag| #vep_tags_ol keeps a specific order, vep_vals.keys does not
                    feature_info_arr.push(vep_vals[vep_val_tag])
                  end
                  feature_info = feature_info_arr.join("|")
                  
                  all_features.push(feature_info)
                end
                # Have final csq to use - log it. 
                
                # 06/11/19 - RJM:
                # removing hardcoding of variable, set in config with boot.rb required fields
                vep_tags_ol.each_with_index do |tag_ol, i| #vep_tags_ol keeps features columns in static order
                  if tag_ol == "vep_all_features"
                    csq_picks[tag_ol] = all_features.join(",")
                  end
                  
                  if i < (vep_tags_ol.count - 1)
                    f.print(csq_picks[tag_ol] + "\t")
                  else
                    f.print(csq_picks[tag_ol])
                  end
                end
            
              # Not a CSQ info field, so spit out chrom, pos, ref, alt if that is the index the logic is at
              elsif [0, 1, 3, 4].include?(i)
                  f.print(field + "\t")
              end
              
            # Finished with line field interations
            end
            f.print("\n")
        rescue => error
          p line
          puts error.backtrace
          raise
          exit
        end
        
        # Finished processing line
        end
      # Finish VEP popen
      end
    # Finished output file open
    end

    #Annotate kafeen vcf with info from vep
    @@log.info("Compressing #{tmp_vep_output}...")
    `bgzip -f #{tmp_vep_output}`
    @@log.info("Compressed VEP output written to #{tmp_vep_output}.gz")

    @@log.info("Indexing #{tmp_vep_output}.gz...")
    `tabix -f -s1 -b2 -e2 #{tmp_vep_output}.gz`
    @@log.info("Done creating VEP index file")

    tmp_vep_output_vcf = "#{out_file_prefix}.vep.vcf.gz"
    tmp_vep_header = "#{out_file_prefix}.vep_header.tmp.txt"

    @@log.info("Creating VEP VCF header file...")
    
    # 06/11/19 - RJM:
    # removing hardcoding of variable, set in config with boot.rb required fields
    header = [] #converted to '\n' delimited string later
    bcftools_annotate_cols = [] #converted to ',' delimited string later
    vep_tags_ol.each do |tag_ol|
      header.push("##INFO=<ID=#{vep_tags[tag_ol]['vcf_tag']},Number=#{vep_tags[tag_ol]['vcf_number']},Type=#{vep_tags[tag_ol]['vcf_type']},Description=\"#{vep_tags[tag_ol]['vcf_description']}\">")
      bcftools_annotate_cols.push(vep_tags[tag_ol]['vcf_tag'])
    end
    header = header.join("\n") + "\n"
    bcftools_annotate_cols = bcftools_annotate_cols.join(',')

    File.open(tmp_vep_header, 'w') { |f| f.write(header) }
    @@log.info("VEP VCF header file written to #{tmp_vep_header}")

    # 06/11/19 - RJM:
    # removing hardcoding of variable, set in config with boot.rb required fields
    `bcftools annotate \
       --annotations #{tmp_vep_output}.gz \
       --columns CHROM,POS,REF,ALT,#{bcftools_annotate_cols} \
       --header-lines #{tmp_vep_header} \
       --output #{tmp_vep_output_vcf} \
       --output-type z \
       #{vcf_file}`
    
    @@log.debug("bcftools annotate command to create temp_vep_output_vcf
      bcftools annotate 
       --annotations #{tmp_vep_output}.gz 
       --columns CHROM,POS,REF,ALT,#{bcftools_annotate_cols} 
       --header-lines #{tmp_vep_header} 
       --output #{tmp_vep_output_vcf} 
       --output-type z 
       #{vcf_file}")

    # Rename tmp output to final output
    File.rename(tmp_vep_output_vcf, @add_vep_result)
    @@log.info("Output written to #{@add_vep_result}")

    # Index output file
    @@log.info("Indexing #{@add_vep_result}...")
    `bcftools index \
       --force \
       --tbi \
       #{@add_vep_result}`
    @@log.info("Done creating index file")

    # Remove tmp files
    if intmed_output
      @add_vep_result_intmed = "#{out_file_prefix}.intmed-add_vep.vcf.gz"
      FileUtils.copy_file(@add_vep_result,@add_vep_result_intmed)
      @add_vep_result_intmed_tbi = "#{@add_vep_result_inmed}.tbi"
      FileUtils.copy_file("#{@add_vep_result}.tbi",@add_vep_result_intmed_tbi)
      @@log.info("An intermediate file of add_vep written to #{@add_vep_result_intmed} and its .tbi")
    else
      @@log.info("Removing VEP tmp files...")
      File.unlink("#{tmp_vep_output}.gz") if File.exist?("#{tmp_vep_output}.gz")
      File.unlink("#{tmp_vep_output}.gz.tbi") if File.exist?("#{tmp_vep_output}.gz.tbi")
      #File.unlink(vep_error_log) if File.exist?(vep_error_log)
      File.unlink(tmp_vep_header) if File.exist?(tmp_vep_header)
      @@log.info("Done removing VEP tmp files")
    end
    
  end


  ##
  # Finalize Pathogenicity
  ##
  def finalize_pathogenicity(vcf_file:, out_file_prefix:, clinical_labels:, enable_benign_star: false, intmed_output: false, maf_threshold: 0.005)
    # Set custom VCF tags to be added to output
    @final_pathogenicity_tag = "FINAL_PATHOGENICITY"
    @final_diseases_tag = "FINAL_DISEASE"
    @final_pmids_tag = "FINAL_PMID"
    @final_pmids_primary_tag = "FINAL_PMID_PRIMARY"
    @final_comments_tag = "FINAL_COMMENTS"
    @final_pathogenicity_source_tag = "FINAL_PATHOGENICITY_SOURCE"
    @final_pathogenicity_reason_tag = "FINAL_PATHOGENICITY_REASON"
    @clinvar_hgmd_conflict_tag = "CLINVAR_HGMD_CONFLICTED"
    @classify_tracking_tag = "CLASSIFY_TRAIL"


    @@log.debug("Finalizing pathogenicity...")
    tmp_output_file = "#{out_file_prefix}.tmp.vcf"

    # Initialize final pathogenicity fields
    final = {}
    final[:pathogenicity] = '.'
    final[:diseases] = '.'
    final[:source] = '.'
    final[:pmids] = '.'
    final[:pmids_primary] = '.'
    final[:reason] = '.'
    final[:comments] = '.'

    # Set HGMD pathogenicity dictionary
    hgmd_pathogenicity_map = {
      'DM'  => clinical_labels['pathogenic'], # Disease mutation
      'DM?' => clinical_labels['unknown'],    # Possible disease mutation
      'DP'  => clinical_labels['benign'],     # Disease-associated polymorphism
      'DFP' => clinical_labels['benign'],     # Disease-associated polymorphism with additional supporting functional evidence
      'FP'  => clinical_labels['benign'],     # In vitro/laboratory or in vivo functional polymorphism
      'FTV' => clinical_labels['benign'],     # Frameshift / truncating variant 
      'CNV' => clinical_labels['benign'],     # Copy number variation
      'R'   => clinical_labels['unknown'],    # Removed from HGMD
      ''    => ''                             # Nothing
    }

    File.open(tmp_output_file, 'w') do |f|
      @@log.info("Adding final pathogenicity to #{tmp_output_file}...")
      `bcftools view \
         --output-type v \
         #{vcf_file}`
      .each_line do |vcf_row|
        vcf_row.chomp!
        if vcf_row.match(/^##/)
          # Print meta-info
          f.puts vcf_row
        elsif vcf_row.match(/^#[^#]/)
          # Add new tags to meta-info
          f.puts "##INFO=<ID=#{@final_pathogenicity_tag},Number=.,Type=String,Description=\"Final curated pathogenicity\">"
          f.puts "##INFO=<ID=#{@final_diseases_tag},Number=.,Type=String,Description=\"Final curated disease\">"
          f.puts "##INFO=<ID=#{@final_pathogenicity_source_tag},Number=.,Type=String,Description=\"Source for final pathogenicity\">"
          f.puts "##INFO=<ID=#{@final_pmids_tag},Number=.,Type=String,Description=\"PubMed IDs\">"
          f.puts "##INFO=<ID=#{@final_pmids_primary_tag},Number=.,Type=String,Description=\"PubMed IDs that factored into the DVD classification\">"
          f.puts "##INFO=<ID=#{@final_pathogenicity_reason_tag},Number=.,Type=String,Description=\"Brief reason for final pathogenicity\">"
          f.puts "##INFO=<ID=#{@final_comments_tag},Number=.,Type=String,Description=\"Additional comments from curator\">"
          f.puts "##INFO=<ID=#{@clinvar_hgmd_conflict_tag},Number=.,Type=String,Description=\"ClinVar and HGMD disagree (0 - No, 1 - Yes)\">"
          f.puts "##INFO=<ID=#{@classify_tracking_tag},Number=.,Type=String,Description=\"A string of numbers delimited by #{@tracker_delim}. Represents the logical path a variant took to its classification beginning with the 1st number and ending with the last. Each number is a unique logical decision. Refer to FILE_PREFIX.classify_track_list.tsv for decoding the string.\">"
  
          # Print header
          f.puts vcf_row
        else
          
          @classify_predictions.clear_trail() #clears any previous trails so that each variant begins with a clean trail
          
          vcf_cols = vcf_row.split("\t")
          @@log.debug("Processing: #{vcf_cols[0]}\t#{vcf_cols[1]}\t#{vcf_cols[3]}\t#{vcf_cols[4]}")
  
          # Initialize final pathogenicity fields
          final = {}
          final[:pathogenicity] = clinical_labels['unknown']
          final[:diseases] = '.'
          final[:source] = "."
          final[:pmids] = '.' # all pmids associated with the variant (concatenation of mancur, clinvar, and hgmd)
          final[:pmids_primary] = '.' # decision driving pmids
          final[:clinvar_hgmd_conflict] = '.'
          final[:reason] = '.'   # <- This field is for internal use only
          final[:comments] = '.' # <- Comments are for public and internal use
          final[:classify_trail] = vcf_row.get_vcf_field(@classify_track_part_tag) # internal use only, begins with the partial trail from add_predictions
  
          # Finalize pathogenicity
          if vcf_cols[7].scan(/(?:^|[\t;])CURATED_PATHOGENICITY=([^;\t]*)/).flatten.any? { |p| p != '.'  } == true
            @@log.debug("- Pathogenicity is based on expert curation")
            @classify_predictions.add_track("mancur","Pathogenicity is based on expert curation. --> pathogenicity = [Manually Curated Pathogenicity].")
            # ^Check for expert-curated pathogenicity
            final[:pathogenicity] = vcf_cols[7].get_vcf_field('CURATED_PATHOGENICITY')
            final[:diseases] = vcf_cols[7].get_vcf_field('CURATED_DISEASE')
              
            final[:pmids_primary] = vcf_cols[7].get_vcf_field('CURATED_PMID')
#              # rob edit - 4/15/2019
#            if final[:pmids] == '.'
#              final[:pmids] = vcf_cols[7].get_vcf_field('CURATED_PMID')
#            else
#              final[:pmids] = (final[:pmids] + ',' + vcf_cols[7].get_vcf_field('CURATED_PMID')).split(/\D+/).uniq.join(',')
#            end
            
            final[:source] = "Expert-curated"
            final[:reason] = "This variant has been expertly curated"
            final[:comments] = vcf_cols[7].get_vcf_field('CURATED_COMMENTS')
          elsif vcf_cols[7].scan(/[^;\t]*_AF=([^;\t]*)/).flatten.any? { |af| af.to_f >= maf_threshold } == true
            @@log.debug("- Pathogenicity is based on MAF (>= #{maf_threshold} in at least one population)")
            @classify_predictions.add_track("not_mancur","Pathogenicity is not based on expert curation.")
            @classify_predictions.add_track("MAF>=#{maf_threshold}","Pathogenicity is based on MAF (>=#{maf_threshold} in at least 1 population).")
            # ^Check if max MAF >= maf_threshold
            final[:pathogenicity] = clinical_labels['benign']
            final[:diseases] = "."
            final[:pmids_primary] = "."
            final[:source] = "MAF"
            final[:reason] = "MAF >= #{maf_threshold}"
            final[:comments] = "This variant contains a MAF in at least one population that meets or exceeds our maximum cutoff of #{maf_threshold}."
            # Convert to "Benign*" if previously reported pathogenic
            if enable_benign_star == true
              # Is it pathogenic in ClinVar and/or HGMD (with high confidence)?
              pathogenic_in_clinvar = !vcf_cols[7].match(/(?:^|[\t;])CLINVAR_CLNSIG=[^;]*(?<![-_a-zA-Z])Pathogenic(?![-_a-zA-Z])/i).nil?
              pathogenic_in_hgmd = (!vcf_cols[7].match(/(?:^|[\t;])HGMD_VARIANTTYPE=DM(?:[;\t]|$)/i).nil? && !vcf_cols[7].match(/(?:^|[\t;])HGMD_CONFIDENCE=High(?:[;\t]|$)/i).nil?)
              
              #clarifying pathogenic hgmd decisions for logic tracking....needs testing
              if(!vcf_cols[7].match(/(?:^|[\t;])HGMD_VARIANTTYPE=DM(?:[;\t]|$)/i).nil?)
                @classify_predictions.add_track("HGMD_patho","Reported pathogenic in HGMD with Low or High confidence.")
              end
              
              if(!vcf_cols[7].match(/(?:^|[\t;])HGMD_CONFIDENCE=High(?:[;\t]|$)/i).nil?)
                @classify_predictions.add_track("HGMD_high","Reported in HGMD with High confidence.")
              else
                @classify_predictions.add_track("HGMD_low","Reported in HGMD with Low confidence.")
              end
              
              if pathogenic_in_clinvar && pathogenic_in_hgmd
                # ^Pathogenic in ClinVar *and* HGMD
                @@log.debug("- Reported pathogenic in ClinVar and HGMD (with high confidence)")
                @classify_predictions.add_track("clinvarHGMD_patho","Reported pathogenic in ClinVar and HGMD (with high confidence) --> pathogenicity = user mapping for Benign + * = #{clinical_labels['benign']}*") #needs to be tested
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in ClinVar/HGMD"
                final[:comments] += " Additionally this variant has been reported pathogenic in both ClinVar and the literature provided in PubMed."
              elsif pathogenic_in_clinvar && !pathogenic_in_hgmd
                # ^Pathogenic in ClinVar (not HGMD)
                @@log.debug("- Reported pathogenic in ClinVar but not HGMD")
                @classify_predictions.add_track("clinvar_notHGMD","Reported pathogenic in ClinVar but not HGMD --> pathogenicity = user mapping for Benign + * = #{clinical_labels['benign']}*") #needs to be tested
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in ClinVar"
                final[:comments] += " Additionally this variant has been reported pathogenic in ClinVar."
              elsif !pathogenic_in_clinvar && pathogenic_in_hgmd
                # ^Pathogenic in HGMD (not ClinVar)
                @@log.debug("- Reported pathogenic in HGMD (with high confidence) but not ClinVar")
                @classify_predictions.add_track("HGMD_notclinvar","Reported pathogenic in HGMD (wit high confidence) but not ClinVar. --> pathogenicity = user mapping for Benign + * = #{clinical_labels['benign']}*") #needs to be tested
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in HGMD"
                final[:comments] += " Additionally this variant has been reported pathogenic in the literature provided in PubMed."
              else
                @classify_predictions.add_track("MAF>=#{maf_threshold},notClinVar/HGMD","MAF>=#{maf_threshold} and not found pathogenic in either HGMD or ClinVar. --> pathogenicity = user mapping for Benign = #{clinical_labels['benign']}") #needs to be tested
              end
            else
              @classify_predictions.add_track("benign*_disabeled","MAF based --> pathogenicity = user mapping for Benign = #{clinical_labels['benign']}") #needs to be tested
            end
          elsif !(vcf_cols[7].match(/(?:^|[\t;])(?:CLINVAR_CLNSIG|HGMD_VARIANTTYPE)=(?:[^;\t]*)/)).nil?
            @classify_predictions.add_track("MAF<#{maf_threshold}","MAF is not >= #{maf_threshold} in any population.")
            @classify_predictions.add_track("not_mancur","Pathogenicity is not based on expert curation.")
            # ^Check for HGMD / ClinVar pathogenicity
  
            # Get ClinVar pathogenicities
            clinvar = {}
            clinvar[:all_pathogenicities] = vcf_cols[7].get_vcf_field('CLINVAR_CLNSIG')
  
            # Get ClinVar diseases, and...
            # 1.) Remove duplicates
            # 2.) Remove meaningless values (e.g. 'not specified', 'not provided', or 'AllHighlyPenetrant')
            clinvar[:diseases] = vcf_cols[7].get_vcf_field('CLINVAR_DISEASE').split(/[|;]/).uniq.delete_if { |e| e.match(/^(?:not specified|not provided|AllHighlyPenetrant)$/i) }.join('; ')
  
            # Get ClinVar PMIDs (remove spaces)
            clinvar[:pmids] = vcf_cols[7].get_vcf_field('CLINVAR_PMID').gsub(' ', '')
              
#            #rob edit - 4/15/2019
#            if final[:pmids] == '.'
#              final[:pmids] = clinvar[:pmids]
#            else
#              final[:pmids] = (final[:pmids] + ',' + clinvar[:pmids]).split(/\D+/).uniq.join(',')
#            end
            
            # Get ClinVar submission conflicts
            clinvar[:conflicted] = vcf_cols[7].get_vcf_field('CLINVAR_CONFLICTED')
  
            # Translate ClinVar pathogenicity
            if !clinvar[:all_pathogenicities].empty?
              if !(clinvar[:all_pathogenicities].match(/(?:^|[,|;])Pathogenic(?:[,|;]|$)/i)).nil?
                # Pathogenic
                clinvar[:worst_pathogenicity] = clinical_labels['pathogenic']
              elsif !(clinvar[:all_pathogenicities].match(/(?:^|[,|;])Likely[-_ ]pathogenic(?:[,|;]|$)/i)).nil?
                # Likely pathogenic
                clinvar[:worst_pathogenicity] = clinical_labels['likely_pathogenic']
              elsif !(clinvar[:all_pathogenicities].match(/(?:^|[,|;])Likely[-_ ](?:benign|non[-_ ]?pathogenic)(?:[,|;]|$)/i)).nil?
                # Likely benign
                clinvar[:worst_pathogenicity] = clinical_labels['likely_benign']
              elsif !(clinvar[:all_pathogenicities].match(/(?:^|[,|;])(?:Benign|Non[-_ ]?pathogenic)(?:[,|;]|$)/i)).nil?
                # Benign
                clinvar[:worst_pathogenicity] = clinical_labels['benign']
              else
                # Unknown significance
                clinvar[:worst_pathogenicity] = clinical_labels['unknown']
              end
            else
              # Unknown significance
              clinvar[:worst_pathogenicity] = ""
            end
  
            # HGMD fields
            hgmd = {}
            hgmd[:pathogenicity] = hgmd_pathogenicity_map[vcf_cols[7].get_vcf_field('HGMD_VARIANTTYPE')]
            hgmd[:diseases] = vcf_cols[7].get_vcf_field('HGMD_DISEASE')
            hgmd[:pmids] = vcf_cols[7].get_vcf_field('HGMD_PMID').gsub(' ', '')
            hgmd[:confidence] = vcf_cols[7].get_vcf_field('HGMD_CONFIDENCE')
            # Low confidence --> Convert pathogenicity to "Likely ..."
            if hgmd[:confidence] == 'Low'
              @classify_predictions.add_track("HGMD_low","Reported in HGMD with Low confidence.") #needs to be tested
              if hgmd[:pathogenicity] == clinical_labels['pathogenic']
                hgmd[:pathogenicity] = clinical_labels['likely_pathogenic']
              elsif hgmd[:pathogenicity] == clinical_labels['benign']
                hgmd[:pathogenicity] = clinical_labels['likely_benign']
              end
            end
            
#            #rob edit - 4/15/2019
#            if final[:pmids] == '.'
#              final[:pmids] = hgmd[:pmids]
#            else
#              final[:pmids] = (final[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
#            end
  
            # Finalize pathogenicity fields...
            if !clinvar[:worst_pathogenicity].empty? && (hgmd[:pathogenicity].empty? || hgmd[:pathogenicity] == clinical_labels['unknown'])
              # ^Only found in ClinVar
              @@log.debug("- Pathogenicity is based on ClinVar submissions only")
              @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
              @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
              @classify_predictions.add_track("onlyClinvar","Pathogenicity based on ClinVar only. --> pathogenicity = [ClinVar]")
              final[:pathogenicity] = clinvar[:worst_pathogenicity]
              final[:diseases] = clinvar[:diseases]
              final[:source] = "ClinVar"
              final[:pmids_primary] = clinvar[:pmids]
              final[:reason] = "Found in ClinVar but not in HGMD"
              final[:comments] = "Pathogenicity is based on ClinVar submissions."
              # Add notes about submission conflicts (if any)
              if clinvar[:conflicted] != '0'
                final[:comments] += " Please be aware that not all submitters agree with this pathogenicity. "
              else
                final[:comments] += " All submitters agree with this pathogenicity."
              end
            elsif !hgmd[:pathogenicity].empty? && (clinvar[:worst_pathogenicity].empty? || clinvar[:worst_pathogenicity] == clinical_labels['unknown'])
              # ^Only found in HGMD
              @@log.debug("- Pathogenicity is based on HGMD only")
              @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
              @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
              @classify_predictions.add_track("onlyHGMD","Pathogenicity based on HGMD only. --> pathogenicity = [HGMD]") #needs to be tested
              final[:pathogenicity] = hgmd[:pathogenicity]
              final[:diseases] = hgmd[:diseases]
              final[:source] = "HGMD"
              final[:pmids_primary] = hgmd[:pmids]
              final[:reason] = "Found in HGMD but not in ClinVar"
              final[:comments] = "Pathogenicity is based on the literature provided in PubMed."
            elsif !clinvar[:worst_pathogenicity].empty? && !hgmd[:pathogenicity].empty?
              # ^Found in ClinVar and HGMD
              @@log.debug("- Pathogenicity is based on ClinVar and HGMD")
              if clinvar[:worst_pathogenicity] == hgmd[:pathogenicity]
                # ClinVar and HGMD totally agree
                @@log.debug("- ClinVar/HGMD totally agree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                @classify_predictions.add_track("clinvarHGMD_agree","Pathogenicity based on agreeing ClinVar and HGMD. --> pathogenicity = [ClinVar/HGMD]")
                final[:pathogenicity] = clinvar[:worst_pathogenicity]
                final[:diseases] = hgmd[:diseases]
                final[:source] = "ClinVar/HGMD"
                final[:pmids_primary] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar/HGMD agree"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed."
                final[:clinvar_hgmd_conflict] = 0
              elsif (clinvar[:worst_pathogenicity] == clinical_labels['pathogenic'] && hgmd[:pathogenicity] == clinical_labels['likely_pathogenic']) || (clinvar[:worst_pathogenicity] == clinical_labels['likely_pathogenic'] && hgmd[:pathogenicity] == clinical_labels['pathogenic'])
                # ClinVar says "Pathogenic", and HGMD says "Likely pathogenic" *OR* vice versa
                @@log.debug("- ClinVar/HGMD mostly agree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                @classify_predictions.add_track("clinvarHGMD_mostlyagree_patho","ClinVar/HGMD mostly agree pathogenic (ie. ClinVar = Pathogenic, HGMD = Likely Pathogenic, or vice versa). --> pathogenicity = user mapping for likely_pathogenic = #{clinical_labels['likely_pathogenic']}") #needs to be tested
                final[:pathogenicity] = clinical_labels['likely_pathogenic']
                final[:diseases] = hgmd[:diseases]
                final[:source] = "ClinVar/HGMD mostly agree"
                final[:pmids_primary] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed. It is important to note that while ClinVar calls this variant '#{clinvar[:worst_pathogenicity]}', the consensus of the literature is that the variant is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 0
              elsif (clinvar[:worst_pathogenicity] == clinical_labels['benign'] && hgmd[:pathogenicity] == clinical_labels['likely_benign']) || (clinvar[:worst_pathogenicity] == clinical_labels['likely_benign'] && hgmd[:pathogenicity] == clinical_labels['benign'])
                # ClinVar says "Benign", and HGMD says "Likely benign" *OR* vice versa
                @@log.debug("- ClinVar/HGMD mostly agree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                @classify_predictions.add_track("clinvarHGMD_mostlyagree_benign","ClinVar/HGMD mostly agree benign (ie. ClinVar = Benign, HGMD = Likely benign or vice versa). --> pathogenicity = user mapping for likely_benign = #{clinical_labels['likely_benign']}") #needs to be tested
                final[:pathogenicity] = clinical_labels['likely_benign']
                final[:diseases] = '.'
                final[:source] = "ClinVar/HGMD mostly agree"
                final[:pmids_primary] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed. It is important to note that while ClinVar calls this variant '#{clinvar[:worst_pathogenicity]}', the consensus of the literature is that the variant is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 0
              else
                # ClinVar and HGMD totally disagree
                @@log.debug("- ClinVar/HGMD totally disagree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                @classify_predictions.add_track("clinvarHGMD_totallyDisagree","ClinVar and HGMD totally disagree. --> pathogenicity = user mapping for VUS = #{clinical_labels['unknown']}") #needs to be tested
                final[:pathogenicity] = clinical_labels['unknown']
                final[:pmids_primary] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:source] = "ClinVar/HGMD"
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity cannot accurately be determined at this time due to conflicts between ClinVar and the consensus of the literature provided in PubMed. ClinVar determines this variant to be '#{clinvar[:worst_pathogenicity]}' while the consensus of the literature seems to be that it is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 1
              end
            else
              # ^Not found in ClinVar or HGMD
              @@log.debug("- SOMETHING WENT WRONG")
              @classify_predictions.add_track("debugError_clinvarHGMD_foundbutnotfound","Not found in ClinVar nor HGMD. SOMETHING WENT WRONG. --> pathogenicity = user mapping for VUS = #{clinical_labels['unknown']}") #needs to be tested
              abort if @@log.level == Logger::DEBUG
              final[:pathogenicity] = clinical_labels['unknown']
              final[:source] = "."
              final[:reason] = "Not enough information"
              final[:comments] = "This variant is a VUS because it does not have enough information."
            end
          elsif !(match = vcf_cols[7].get_vcf_field(@final_pred_tag)).empty? && match != '.'
            @@log.debug("- Pathogenicity is based on predictions from dbNSFP")
            @classify_predictions.add_track("MAF<#{maf_threshold}","MAF is not >=#{maf_threshold} in any population.")
            @classify_predictions.add_track("not_mancur","Pathogenicity is not based on expert curation.")
            @classify_predictions.add_track("notIn_HGMD/ClinVar","Neither ClinVar nor HGMD provided classification for this variant.")
            # ^Use dbNSFP prediction
            if match == clinical_labels['unknown']
              # ^Not enough prediction data
              @classify_predictions.add_track("dbNSFP_VUS","Classified based on dbNSFP, but dbNSFP does not have enough prediction data. --> pathogencity = user mapping for VUS  = #{clinical_labels['unknown']}") #needs to be tested
              final[:pathogenicity] = clinical_labels['unknown']
              final[:diseases] = '.'
              final[:source] = '.'
              final[:pmids_primary] = '.'
              final[:reason] = "Not enough information"
              final[:comments] = "This variant is a VUS because it does not have enough information."
              @@log.debug("- Not enough predictions")
            else
              # ^Set final pathogenicity as predicted pathogenicity
              # these predicted pathogenicity classifications are from add_predictions fx
              @classify_predictions.add_track("dbNSFP_pred","Classified based on dbNSFP. --> pathogenicity = [dbNSFP].")
              final[:pathogenicity] = match
              final[:diseases] = '.'
              final[:source] = "dbNSFP"
              final[:pmids_primary] = '.'
              if !(num_path_preds = vcf_cols[7].get_vcf_field(@num_path_preds_tag)).empty?
                # Get pathogenic prediction fraction
                if num_path_preds != '.'
                  total_num_preds = vcf_cols[7].get_vcf_field(@total_num_preds_tag)
                  final[:reason] = "#{num_path_preds}/#{total_num_preds} pathogenic"
                  final[:comments] = "Pathogenicity is based on prediction data only. #{num_path_preds} out of #{total_num_preds} predictions were pathogenic."
                  @@log.debug("- #{num_path_preds}/#{total_num_preds} pathogenic predictions")
                end
              else
                # Could not find prediction numbers (ideally this should never happen)
                final[:reason] = "Pathogenicity is based on prediction data only."
                final[:comments] = "Pathogenicity is based on prediction data only."
                @@log.error("- PREDICTION COUNTS WERE NOT FOUND")
                @classify_predictions.add_track("NA_pred_counts","Could not find prediction numbers (ideally this should never happen). An error.")
              end
            end
          else
            # Unknown significance
            @@log.debug("- This variant is a VUS because it does not have enough information.")
            @classify_predictions.add_track("not_mancur","Pathogenicity is not based on expert curation.")
            @classify_predictions.add_track("MAF<#{maf_threshold}","MAF is not >=#{maf_threshold} in any population.")
            @classify_predictions.add_track("not_enough_info","This variant does not have enough information. --> pathogenicity = user mapping for VUS = #{clinical_labels['unknown']}") #needs to be tested
            
            final[:pathogenicity] = clinical_labels['unknown']
            final[:reason] = "Not enough information"
            final[:comments] = "This variant is a VUS because it does not have enough information."
          end
          
          final[:classify_trail] += @tracker_delim + @classify_predictions.get_trail()
          @classify_predictions.clear_trail()
          
          # compound all pmid information into final[:pmids] tag, assumption that final[:pmids] is previously instantiated (as it should be)
          
          # manually curated
          if vcf_cols[7].scan(/(?:^|[\t;])CURATED_PATHOGENICITY=([^;\t]*)/).flatten.any? { |p| p != '.'  } == true
            # rob edit - 4/16/2019
            if final[:pmids] == '.'
              final[:pmids] = vcf_cols[7].get_vcf_field('CURATED_PMID')
            else
              final[:pmids] = (final[:pmids] + ',' + vcf_cols[7].get_vcf_field('CURATED_PMID')).split(/\D+/).uniq.join(',')
            end
          end
          
          # clinvar
          if !(vcf_cols[7].match(/(?:^|[\t;])(?:CLINVAR_CLNSIG)=(?:[^;\t]*)/)).nil?
            clinvar = {}
            clinvar[:pmids] = vcf_cols[7].get_vcf_field('CLINVAR_PMID').gsub(' ', '')
            #rob edit - 4/16/2019
            if final[:pmids] == '.'
              final[:pmids] = clinvar[:pmids]
            else
              final[:pmids] = (final[:pmids] + ',' + clinvar[:pmids]).split(/\D+/).uniq.join(',')
            end
          end
          
          # hgmd
          if !(vcf_cols[7].match(/(?:^|[\t;])(?:HGMD_VARIANTTYPE)=(?:[^;\t]*)/)).nil?
            hgmd = {}
            hgmd[:pmids] = vcf_cols[7].get_vcf_field('HGMD_PMID').gsub(' ', '')
            #rob edit - 4/16/2019
            if final[:pmids] == '.'
              final[:pmids] = hgmd[:pmids]
            else
              final[:pmids] = (final[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
            end
          end
#          if final[:pmids] == '.'
#            final[:pmids] = vcf_cols[7].get_vcf_field('CURATED_PMID')
#          else
#            final[:pmids] = (final[:pmids] + ',' + vcf_cols[7].get_vcf_field('CURATED_PMID')).split(/\D+/).uniq.join(',')
#          end
          
#          if !(vcf_cols[7].match(/(?:^|[\t;])(?:CLINVAR_CLNSIG)=(?:[^;\t]*)/)).nil? && !clinvar[:pmids].empty?
#            if final[:pmids] == '.'
#              final[:pmids] = clinvar[:pmids]
#            else
#              final[:pmids] = (final[:pmids] + ',' + clinvar[:pmids]).split(/\D+/).uniq.join(',')
#            end
#          end
          
#          if !(vcf_cols[7].match(/(?:^|[\t;])(?:HGMD_VARIANTTYPE)=(?:[^;\t]*)/)).nil? && !hgmd[pmids].empty?
#            if final[:pmids] == '.'
#              final[:pmids] = hgmd[:pmids]
#            else
#              final[:pmids] = (final[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
#            end
#          end
          
          # Remove illegal characters and set all empty values to '.'
          final.each do |k, v|
            next if final[k] == '.'
            if v.to_s.strip.empty?
              final[k] = '.'
            else
              final[k] = URI.escape(final[k].to_s, ';,= ')
            end
          end
  
          # Update INFO column
          vcf_cols[7] = [
            vcf_cols[7],
            "#{@final_pathogenicity_tag}=#{final[:pathogenicity]}",
            "#{@final_pmids_tag}=#{final[:pmids]}",
            "#{@final_pmids_primary_tag}=#{final[:pmids_primary]}",
            "#{@final_comments_tag}=#{final[:comments]}",
            "#{@final_pathogenicity_source_tag}=#{final[:source]}",
            "#{@final_pathogenicity_reason_tag}=#{final[:reason]}",
            "#{@final_diseases_tag}=#{final[:diseases]}",
            "#{@clinvar_hgmd_conflict_tag}=#{final[:clinvar_hgmd_conflict]}",
            "#{@classify_tracking_tag}=#{final[:classify_trail]}"
          ].map { |e| e.force_encoding('UTF-8') }.join(";")
          
          # Print updated VCF row
          f.puts vcf_cols.join("\t")
  
          @@log.debug("- This variant has been labeled \"#{final[:pathogenicity]}\"")
          @@log.debug("------------------------------------------------------")
        end
      end # <-- End parsing bcftools result
    end # <-- End printing output file

    @@log.info("Final pathogenicity added to #{tmp_output_file}")
    
    @@log.info("Writing classification tracking key to #{@classify_trackers_list}...")
    File.open(@classify_trackers_list, 'w') do |cf|
      cf.puts(@classify_predictions.list_trackers())
    end
    @@log.info("Classification tracking key written to #{@classify_trackers_list}")
    
    @finalize_pathogenicity_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Compressing #{tmp_output_file}...")
    # Compress the output file
    `bcftools view \
       --output-type z \
       --output-file #{@finalize_pathogenicity_result} \
       #{tmp_output_file}`
    @@log.info("Compressed output written to #{@finalize_pathogenicity_result}")

    # Index output file
    @@log.info("Indexing #{@finalize_pathogenicity_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@finalize_pathogenicity_result}`
    @@log.info("Done creating index file")

    # Remove tmp files  
    @@log.info("Removing tmp files...")
    File.unlink(tmp_output_file) if File.exist?(tmp_output_file)
    @@log.info("Done removing tmp files")
    
    if intmed_output
      @finalize_pathogenicity_result_intmed = "#{out_file_prefix}.intmed-finalize_pathogenicity.vcf.gz"
      FileUtils.copy_file(@finalize_pathogenicity_result,@finalize_pathogenicity_result_intmed)
      @finalize_pathogenicity_result_intmed_tbi = "#{@finalize_pathogenicity_result_inmed}.tbi"
      FileUtils.copy_file("#{@finalize_pathogenicity_result}.tbi",@finalize_pathogenicity_result_intmed_tbi)
      @@log.info("An intermediate file of finalize_pathogenicity written to #{@finalize_pathogenicity_result_intmed} and its .tbi")
    end
    
  end

  ##
  # Cleanup Meta-info
  #
  # Spruce up the meta-information (i.e. lines that start with ## at the
  # top of the VCF) in the VCF. This will:
  #   - Update fileformat to current VCF version
  #   - Update the fileDate to the current date
  #   - Keep all ##INFO and ##reference lines
  #   - Filter out non-chromosomal contigs (i.e. not 1-22|X|Y|M/MT)
  #   - Remove all other meta-info not necessary for final output
  ##
  def cleanup_meta_info(vcf_file:, out_file_prefix:, intmed_output: false)
    @@log.info("Cleaning up VCF meta-info...")

    tmp_header_file = "#{out_file_prefix}.header.tmp.txt"
    tmp_output_file = "#{out_file_prefix}.tmp.vcf.gz"
    @cleanup_header_result = "#{out_file_prefix}.vcf.gz"

    # Create new header file
    File.open(tmp_header_file, 'w') do |f|
      # Print VCF version & current date
      f.puts "##fileformat=VCFv4.3"
      f.puts "##fileDate=#{Time.now.strftime("%Y%m%d")}"

      # Remove unnecessary/inaccurate header lines
      `bcftools view \
         --header-only \
         #{vcf_file}`
      .each_line do |line|
        line.chomp!
        if line.match(/^#(?:#reference=|#INFO=|CHROM)/) ||
               line.match(/^##contig=<ID=(?:chr)?(?:[XYM]|MT|[1-9]?|1[0-9]|2[0-2])>$/)
          f.puts line
        end
      end
    end
    
    # Add new header to tmp file
    `bcftools reheader \
       --header #{tmp_header_file} \
       --output #{tmp_output_file} \
       #{vcf_file}`

    # Index output file
    @@log.info("Indexing #{@cleanup_header_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@cleanup_header_result}`
    @@log.info("Done creating index file")

    # Remove tmp files
    @@log.info("Removing tmp files...")
    File.unlink(tmp_header_file) if File.exist?(tmp_header_file)
    @@log.info("Done removing tmp files")

    # Move tmp file to be final output file
    File.rename(tmp_output_file, @cleanup_header_result)
    @@log.info("Output written to #{@cleanup_header_result}")
  end

end
