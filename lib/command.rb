# A library of available pipeline commands
#
# @author Sean Ephraim

require 'logger'
require 'open3'
require 'uri'
require_relative 'core_extensions'

# Monkey-patch the String class
String.include CoreExtensions::String::Vcf
String.include CoreExtensions::String::Colorize

class Command
  # Result file paths
  attr_reader :genes2regions_result
  attr_reader :genes2regions_merged_result
  attr_reader :regions2variants_result
  attr_reader :add_genes_result
  attr_reader :add_predictions_result
  attr_reader :add_asap_result
  attr_reader :finalize_pathogenicity_result

  ##
  # Initialize
  ##
  def initialize(log_level: 'info', log_out: STDOUT)
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
  end

  ##
  # Genes to Regions
  #
  # Get genomic regions for each HGNC gene symbol
  ##
  def genes2regions(genes_file:, ref_file:, out_file_prefix:)
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
  end

  ##
  # Regions to Variants
  #
  # Get a list of all variants within specified regions
  ##
  def regions2variants(bed_file:, vcf_files:, out_file_prefix:, keep_tmp_files: false)
    tmp_vcfs = {}
    # Query all VCF files for variants
    vcf_files.each do |key, vcf|
      if vcf['include'] or vcf['include'].nil? or !vcf.has_key? :include # RJM: checks flag in config.yml if the vcf source should be included or not in the kafeen run. true allows the source to be considered as normal, false ignores the source
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
  
          # Store tmp file name (filename) and the original VCF that the data came from (parent)
          tmp_vcfs[key] = {'filename' => tmp_source_vcf, 'parent' => vcf['filename']}
        end
      else if !vcf['include']
        @@log.info("EXCLUDED #{vcf['source']}...Source flagged for exclusion in config.yaml")
        @@log.info("Skipping annotation source: #{vcf['source']}.")
      else
        @@log.error("EXCLUDED #{vcf['source']}...Incompatible include input within config.yaml.")
        @@log.error("Please review the config.yaml and indicated whether or not you would like to include #{vcf['source']} (include: true) or not (include: false)")
      end # RJM: config.yml vcf source flag check end
    end

    # Construct list of VCFs to merge
    files_to_merge = []
    tmp_vcfs.each do |key, tmp_vcf|
      next if key == 'dbnsfp' # DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
      files_to_merge << tmp_vcf['filename']
    end
#    files_to_merge << tmp_vcfs['dbsnp']['filename']

    # Merge VCFs...
    @regions2variants_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Merging results...")
    `bcftools merge \
       --merge none \
       --output-type z \
       --output #{@regions2variants_result} \
       #{files_to_merge.join(' ')}`
    @@log.info("Done merging results")

    # Create index
    @@log.info("Creating index file for #{@regions2variants_result}...")
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
  def add_genes(bed_file:, vcf_file:, out_file_prefix:)
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
  end

  ##
  # Add predictions from dbNSFP
  ##
  def add_predictions(dbnsfp_file:, vcf_file:, bed_file:, out_file_prefix:, clinical_labels:)
    # Set custom VCF tags to be added to output
    @gerp_pred_tag = "DBNSFP_GERP_PRED"
    @phylop20way_mammalian_pred_tag = "DBNSFP_PHYLOP20WAY_MAMMALIAN_PRED"
    @num_path_preds_tag = "NUM_PATH_PREDS"
    @total_num_preds_tag = "TOTAL_NUM_PREDS"
    @final_pred_tag = "FINAL_PRED"

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
    tmp_output_file = "#{out_file_prefix}.tmp.vcf"
    @@log.info("Adding dbNSFP pedictions to #{tmp_output_file}...")
    File.open(tmp_output_file, 'w') do |f|
    `bcftools annotate \
       --annotations #{dbnsfp_subset_file} \
       --columns #{dbnsfp_file['fields'].map { |f| "INFO/#{f}" }.join(',')} \
       --output-type v \
       #{vcf_file}`
    .each_line do |vcf_row|
         vcf_row.chomp!
         if vcf_row.match(/^##/)
           # Print meta-info
           f.puts vcf_row
         elsif vcf_row.match(/^#[^#]/)
           # Add custom INFO tags
           f.puts "##INFO=<ID=#{@num_path_preds_tag},Number=.,Type=String,Description=\"Number of pathogenic predictions from dbNSFP\">"
           f.puts "##INFO=<ID=#{@total_num_preds_tag},Number=.,Type=String,Description=\"Total number of prediction scores available from dbNSFP\">"
           f.puts "##INFO=<ID=#{@final_pred_tag},Number=.,Type=String,Description=\"Final prediction consensus based on majority vote of prediction scores\">"

           # Add GERP++ prediction tag to meta-info
           if dbnsfp_file['fields'].any?{ |e| e == 'DBNSFP_GERP_RS' }
             f.puts "##INFO=<ID=#{@gerp_pred_tag},Number=.,Type=String,Description=\"NA\">"
           end
           # Add phyloP20way mammalian prediction tag to meta-info
           if dbnsfp_file['fields'].any?{ |e| e == 'DBNSFP_PHYLOP20WAY_MAMMALIAN' }
             f.puts "##INFO=<ID=#{@phylop20way_mammalian_pred_tag},Number=.,Type=String,Description=\"NA\">"
           end
           # Print header (i.e. "CHROM  POS  ID ...")
           f.puts vcf_row
         else
           vcf_cols = vcf_row.split("\t")

           # Analyze each *_PRED field (as well as GERP++ and phyloP)
           # Tally up pathogenic predictions
           output = {}
           output[:total_num_preds] = 0
           output[:num_path_preds] = 0
           dbnsfp_file['fields'].select { |e| e.match(/(?:_PRED$|^DBNSFP_GERP_RS$|^DBNSFP_PHYLOP20WAY_MAMMALIAN$)/i) }.each do |field|
             # Get all predictions for this algorithm
             match = vcf_row.get_vcf_field(field)

             # No data for this algorithm -- skip it
             next if match.empty?

             # Get all predictions for this algorithm
             preds = match.split(/[^a-zA-Z0-9.-]+/)

             # No data for this algorithm -- skip it
             next if preds.all? { |pred| pred == '.' || pred == 'U' }
               
             if field == 'DBNSFP_SIFT_PRED'
               # SIFT prediction
               output[:num_path_preds] += 1 if preds.include?('D') # <-- "Damaging"
               output[:total_num_preds] += 1
             elsif field == 'DBNSFP_POLYPHEN2_HDIV_PRED'
               # Polyphen2 (HDIV) prediction
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
                 vcf_cols[7] = [vcf_cols[7], "#{@gerp_pred_tag}=C"].join(";")
               else
                 # Non-conserved
                 vcf_cols[7] = [vcf_cols[7], "#{@gerp_pred_tag}=N"].join(";")
               end
               output[:total_num_preds] += 1
             elsif field == 'DBNSFP_PHYLOP20WAY_MAMMALIAN'
               # phyloP20way mammalian prediction
               # TODO: The cutoff should be changed from 0.95 to 1.0
               if preds.any? { |pred| pred.to_f >= 0.95 }
                 # Conserved
                 output[:num_path_preds] += 1
                 vcf_cols[7] = [vcf_cols[7], "#{@phylop20way_mammalian_pred_tag}=C"].join(";")
               else
                 # Non-conserved
                 vcf_cols[7] = [vcf_cols[7], "#{@phylop20way_mammalian_pred_tag}=N"].join(";")
               end
               output[:total_num_preds] += 1
             end
           end

           # Add final prediction
           if output[:total_num_preds] == 0
             # No predictions available
             output[:final_pred] = '.'
             output[:num_path_preds] = '.'
             output[:total_num_preds] = '.'
           elsif output[:total_num_preds] >= 5
             path_score = output[:num_path_preds].to_f/output[:total_num_preds].to_f
             if path_score >= 0.6
               # Predicted pathogenic
               output[:final_pred] = clinical_labels['pred_pathogenic']
             elsif path_score <= 0.4
               # Predicted benign
               output[:final_pred] = clinical_labels['pred_benign']
             else
               # Predicted unknown (benign predictions approx. equal to pathogenic)
               output[:final_pred] = clinical_labels['unknown']
             end
           else
             # Predicted unknown (not enough predictions)
             output[:final_pred] = clinical_labels['unknown']
           end

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
               "#{@final_pred_tag}=#{output[:final_pred]}",
             ].join(";")
           end

           # Print updated VCF row
           f.puts vcf_cols.join("\t")
         end
       end
    end # <-- End printing output file
    
    @@log.info("Predictions added to #{tmp_output_file}")

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
    File.unlink(dbnsfp_subset_file) if File.exist?(dbnsfp_subset_file)
    File.unlink("#{dbnsfp_subset_file}.csi") if File.exist?("#{dbnsfp_subset_file}.csi")
    @@log.info("Done removing tmp files")
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

  def add_asap(vcf_file:, out_file_prefix:, asap_path:, ref_flat:, ref_seq_ali:, fasta:)
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
        #if !fields[5].include?("ERROR_")
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
  end


  def add_vep(vcf_file:, out_file_prefix:, vep_path:, vep_cache_path:)
    @vep_consequence_tag = 'VEP_CONSEQUENCE'
    @vep_exon_tag = 'VEP_EXON'
    @vep_intron_tag = 'VEP_INTRON'
    @vep_hgvs_c_tag = 'VEP_HGVS_C'
    @vep_hgvs_p_tag = 'VEP_HGVS_P'

    @add_vep_result = "#{out_file_prefix}.vcf.gz"

    # Run VEP
    # Make tsv containing relevant fields in vep output,
    # to annotate kafeen vcf with using bcftools
    tmp_vep_output = "#{out_file_prefix}.vep.out.tmp.txt"
    vep_error_log = "#{out_file_prefix}.vep.errors.log"
    vep_header = nil
    @@log.info("Running VEP...")
    File.open(tmp_vep_output, 'w') do |f|
      `#{vep_path} \
        -species human -refseq -format vcf -force -fork 4 -port 3337 -cache \
        -dir #{vep_cache_path} \
        -i #{vcf_file} -o STDOUT -vcf -hgvs -numbers -pick 2> #{vep_error_log}`
      .each_line do |line|
        line.chomp!
        #@@log.info(line)
        if line.start_with?("#")
          if line.include? "<ID=CSQ"
            vep_header = line.split("Format: ")[1].split("\">")[0].split("|")
          end
          next
        end
        fields = line.split("\t")
        fields.each_with_index do |field, i|
          if field.include? "CSQ="
            csq = field.split("CSQ=")[1]
            #csq.split(",").each do |transcript|
            #transcript.split("|").each_with_index do |subfield, j|
            csq.split("|").each_with_index do |subfield, j|
              categ = vep_header[j]
              if ["Consequence", "HGVSc", "HGVSp", "EXON", "INTRON"].include?(categ)
                  f.print(subfield + "\t")
              end
            end
          # chrom, pos, ref, alt
          elsif [0, 1, 3, 4].include?(i)
              f.print(field + "\t")
          end
        end
      f.print("\n")
      end
    end

    #Annotate kafeen vcf with info from vep
    @@log.info("Compressing #{tmp_vep_output}...")
    `bgzip -f #{tmp_vep_output}`
    @@log.info("Compressed VEP output written to #{tmp_vep_output}")

    @@log.info("Indexing #{tmp_vep_output}.gz...")
    `tabix -f -s1 -b2 -e2 #{tmp_vep_output}.gz`
    @@log.info("Done creating VEP index file")

    tmp_vep_output_vcf = "#{out_file_prefix}.vep.vcf.gz"
    tmp_vep_header = "#{out_file_prefix}.vep_header.tmp.txt"

    @@log.info("Creating VEP VCF header file...")
    header = [
      "##INFO=<ID=#{@vep_consequence_tag},Number=1,Type=String,Description=\"Predicted variant consequence.\">",
      "##INFO=<ID=#{@vep_exon_tag},Number=1,Type=String,Description=\"Affected exon numbering. Format is Number\/Total.\">",
      "##INFO=<ID=#{@vep_intron_tag},Number=1,Type=String,Description=\"Affected intron numbering. Format is Number\/Total.\">",
      "##INFO=<ID=#{@vep_hgvs_c_tag},Number=1,Type=String,Description=\"HGVS coding sequence names based on Ensembl stable identifiers.\">",
      "##INFO=<ID=#{@vep_hgvs_p_tag},Number=1,Type=String,Description=\"HGVS protein sequence names based on Ensembl stable identifiers.\">",
    ].join("\n") + "\n"

    File.open(tmp_vep_header, 'w') { |f| f.write(header) }
    @@log.info("VEP VCF header file written to #{tmp_vep_header}.")

    `bcftools annotate \
       --annotations #{tmp_vep_output}.gz \
       --columns CHROM,POS,REF,ALT,#{@vep_consequence_tag},#{@vep_exon_tag},#{@vep_intron_tag},#{@vep_hgvs_c_tag},#{@vep_hgvs_p_tag} \
       --header-lines #{tmp_vep_header} \
       --output #{tmp_vep_output_vcf} \
       --output-type z \
       #{vcf_file}`

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
    @@log.info("Removing VEP tmp files...")
    File.unlink("#{tmp_vep_output}.gz") if File.exist?("#{tmp_vep_output}.gz")
    File.unlink("#{tmp_vep_output}.gz.tbi") if File.exist?("#{tmp_vep_output}.gz.tbi")
    File.unlink(vep_error_log) if File.exist?(vep_error_log)
    File.unlink(tmp_vep_header) if File.exist?(tmp_vep_header)
    @@log.info("Done removing VEP tmp files")
  end


  ##
  # Finalize Pathogenicity
  ##
  def finalize_pathogenicity(vcf_file:, out_file_prefix:, clinical_labels:, enable_benign_star: false)
    # Set custom VCF tags to be added to output
    @final_pathogenicity_tag = "FINAL_PATHOGENICITY"
    @final_diseases_tag = "FINAL_DISEASE"
    @final_pmids_tag = "FINAL_PMID"
    @final_comments_tag = "FINAL_COMMENTS"
    @final_pathogenicity_source_tag = "FINAL_PATHOGENICITY_SOURCE"
    @final_pathogenicity_reason_tag = "FINAL_PATHOGENICITY_REASON"
    @clinvar_hgmd_conflict_tag = "CLINVAR_HGMD_CONFLICTED"


    @@log.debug("Finalizing pathogenicity...")
    tmp_output_file = "#{out_file_prefix}.tmp.vcf"

    # Initialize final pathogenicity fields
    final = {}
    final[:pathogenicity] = '.'
    final[:diseases] = '.'
    final[:source] = '.'
    final[:pmids] = '.'
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
          f.puts "##INFO=<ID=#{@final_pathogenicity_reason_tag},Number=.,Type=String,Description=\"Brief reason for final pathogenicity\">"
          f.puts "##INFO=<ID=#{@final_comments_tag},Number=.,Type=String,Description=\"Additional comments from curator\">"
          f.puts "##INFO=<ID=#{@clinvar_hgmd_conflict_tag},Number=.,Type=String,Description=\"ClinVar and HGMD disagree (0 - No, 1 - Yes)\">"
  
          # Print header
          f.puts vcf_row
        else
          vcf_cols = vcf_row.split("\t")
          @@log.debug("Processing: #{vcf_cols[0]}\t#{vcf_cols[1]}\t#{vcf_cols[3]}\t#{vcf_cols[4]}")
  
          # Initialize final pathogenicity fields
          final = {}
          final[:pathogenicity] = clinical_labels['unknown']
          final[:diseases] = '.'
          final[:source] = "."
          final[:pmids] = '.'
          final[:clinvar_hgmd_conflict] = '.'
          final[:reason] = '.'   # <- This field is for internal use only
          final[:comments] = '.' # <- Comments are for public and internal use
  
          # Finalize pathogenicity
          if vcf_cols[7].scan(/(?:^|[\t;])CURATED_PATHOGENICITY=([^;\t]*)/).flatten.any? { |p| p != '.'  } == true
            @@log.debug("- Pathogenicity is based on expert curation")
            # ^Check for expert-curated pathogenicity
            final[:pathogenicity] = vcf_cols[7].get_vcf_field('CURATED_PATHOGENICITY')
            final[:diseases] = vcf_cols[7].get_vcf_field('CURATED_DISEASE')
            final[:pmids] = vcf_cols[7].get_vcf_field('CURATED_PMID')
            final[:source] = "Expert-curated"
            final[:reason] = "This variant has been expertly curated"
            final[:comments] = vcf_cols[7].get_vcf_field('CURATED_COMMENTS')
          elsif vcf_cols[7].scan(/[^;\t]*_AF=([^;\t]*)/).flatten.any? { |af| af.to_f >= 0.005 } == true
            @@log.debug("- Pathogenicity is based on MAF (>=0.005 in at least one population)")
            # ^Check if max MAF >= 0.005
            final[:pathogenicity] = clinical_labels['benign']
            final[:diseases] = "."
            final[:pmids] = "."
            final[:source] = "MAF"
            final[:reason] = "MAF >= 0.005"
            final[:comments] = "This variant contains a MAF in at least one population that meets or exceeds our maximum cutoff of 0.005."
            # Convert to "Benign*" if previously reported pathogenic
            if enable_benign_star == true
              # Is it pathogenic in ClinVar and/or HGMD (with high confidence)?
              pathogenic_in_clinvar = !vcf_cols[7].match(/(?:^|[\t;])CLINVAR_CLNSIG=[^;]*(?<![-_a-zA-Z])Pathogenic(?![-_a-zA-Z])/i).nil?
              pathogenic_in_hgmd = (!vcf_cols[7].match(/(?:^|[\t;])HGMD_VARIANTTYPE=DM(?:[;\t]|$)/i).nil? && !vcf_cols[7].match(/(?:^|[\t;])HGMD_CONFIDENCE=High(?:[;\t]|$)/i).nil?)
              if pathogenic_in_clinvar && pathogenic_in_hgmd
                # ^Pathogenic in ClinVar *and* HGMD
                @@log.debug("- Reported pathogenic in ClinVar and HGMD (with high confidence)")
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in ClinVar/HGMD"
                final[:comments] += " Additionally this variant has been reported pathogenic in both ClinVar and the literature provided in PubMed."
              elsif pathogenic_in_clinvar && !pathogenic_in_hgmd
                # ^Pathogenic in ClinVar (not HGMD)
                @@log.debug("- Reported pathogenic in ClinVar but not HGMD")
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in ClinVar"
                final[:comments] += " Additionally this variant has been reported pathogenic in ClinVar."
              elsif !pathogenic_in_clinvar && pathogenic_in_hgmd
                # ^Pathogenic in HGMD (not ClinVar)
                @@log.debug("- Reported pathogenic in HGMD (with high confidence) but not ClinVar")
                final[:pathogenicity] += '*'
                final[:reason] += " & pathogenic in HGMD"
                final[:comments] += " Additionally this variant has been reported pathogenic in the literature provided in PubMed."
              end
            end
          elsif !(vcf_cols[7].match(/(?:^|[\t;])(?:CLINVAR_CLNSIG|HGMD_VARIANTTYPE)=(?:[^;\t]*)/)).nil?
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
              if hgmd[:pathogenicity] == clinical_labels['pathogenic']
                hgmd[:pathogenicity] = clinical_labels['likely_pathogenic']
              elsif hgmd[:pathogenicity] == clinical_labels['benign']
                hgmd[:pathogenicity] = clinical_labels['likely_benign']
              end
            end
  
            # Finalize pathogenicity fields...
            if !clinvar[:worst_pathogenicity].empty? && (hgmd[:pathogenicity].empty? || hgmd[:pathogenicity] == clinical_labels['unknown'])
              # ^Only found in ClinVar
              @@log.debug("- Pathogenicity is based on ClinVar submissions only")
              @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
              @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
              final[:pathogenicity] = clinvar[:worst_pathogenicity]
              final[:diseases] = clinvar[:diseases]
              final[:source] = "ClinVar"
              final[:pmids] = clinvar[:pmids]
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
              final[:pathogenicity] = hgmd[:pathogenicity]
              final[:diseases] = hgmd[:diseases]
              final[:source] = "HGMD"
              final[:pmids] = hgmd[:pmids]
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
                final[:pathogenicity] = clinvar[:worst_pathogenicity]
                final[:diseases] = hgmd[:diseases]
                final[:source] = "ClinVar/HGMD"
                final[:pmids] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar/HGMD agree"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed."
                final[:clinvar_hgmd_conflict] = 0
              elsif (clinvar[:worst_pathogenicity] == clinical_labels['pathogenic'] && hgmd[:pathogenicity] == clinical_labels['likely_pathogenic']) || (clinvar[:worst_pathogenicity] == clinical_labels['likely_pathogenic'] && hgmd[:pathogenicity] == clinical_labels['pathogenic'])
                # ClinVar says "Pathogenic", and HGMD says "Likely pathogenic" *OR* vice versa
                @@log.debug("- ClinVar/HGMD mostly agree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                final[:pathogenicity] = clinical_labels['likely_pathogenic']
                final[:diseases] = hgmd[:diseases]
                final[:source] = "ClinVar/HGMD mostly agree"
                final[:pmids] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed. It is important to note that while ClinVar calls this variant '#{clinvar[:worst_pathogenicity]}', the consensus of the literature is that the variant is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 0
              elsif (clinvar[:worst_pathogenicity] == clinical_labels['benign'] && hgmd[:pathogenicity] == clinical_labels['likely_benign']) || (clinvar[:worst_pathogenicity] == clinical_labels['likely_benign'] && hgmd[:pathogenicity] == clinical_labels['benign'])
                # ClinVar says "Benign", and HGMD says "Likely benign" *OR* vice versa
                @@log.debug("- ClinVar/HGMD mostly agree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                final[:pathogenicity] = clinical_labels['likely_benign']
                final[:diseases] = '.'
                final[:source] = "ClinVar/HGMD mostly agree"
                final[:pmids] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity is based on ClinVar submissions and the literature provided in PubMed. It is important to note that while ClinVar calls this variant '#{clinvar[:worst_pathogenicity]}', the consensus of the literature is that the variant is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 0
              else
                # ClinVar and HGMD totally disagree
                @@log.debug("- ClinVar/HGMD totally disagree")
                @@log.debug("  * ClinVar says: #{clinvar[:worst_pathogenicity]}")
                @@log.debug("  * HGMD says: #{hgmd[:pathogenicity]}")
                final[:pathogenicity] = clinical_labels['unknown']
                final[:pmids] = (clinvar[:pmids] + ',' + hgmd[:pmids]).split(/\D+/).uniq.join(',')
                final[:source] = "ClinVar/HGMD"
                final[:reason] = "ClinVar says '#{clinvar[:worst_pathogenicity]}'/HGMD says '#{hgmd[:pathogenicity]}'"
                final[:comments] = "Pathogenicity cannot accurately be determined at this time due to conflicts between ClinVar and the consensus of the literature provided in PubMed. ClinVar determines this variant to be '#{clinvar[:worst_pathogenicity]}' while the consensus of the literature seems to be that it is '#{hgmd[:pathogenicity]}'"
                final[:clinvar_hgmd_conflict] = 1
              end
            else
              # ^Not found in ClinVar or HGMD
              @@log.debug("- SOMETHING WENT WRONG")
              abort if @@log.level == Logger::DEBUG
              final[:pathogenicity] = clinical_labels['unknown']
              final[:source] = "."
              final[:reason] = "Not enough information"
              final[:comments] = "This variant is a VUS because it does not have enough information."
            end
          elsif !(match = vcf_cols[7].get_vcf_field(@final_pred_tag)).empty? && match != '.'
            @@log.debug("- Pathogenicity is based on predictions from dbNSFP")
            # ^Use dbNSFP prediction
            if match == clinical_labels['unknown']
              # ^Not enough prediction data
              final[:pathogenicity] = clinical_labels['unknown']
              final[:diseases] = '.'
              final[:source] = '.'
              final[:pmids] = '.'
              final[:reason] = "Not enough information"
              final[:comments] = "This variant is a VUS because it does not have enough information."
              @@log.debug("- Not enough predictions")
            else
              # ^Set final pathogenicity as predicted pathogenicity
              final[:pathogenicity] = match
              final[:diseases] = '.'
              final[:source] = "dbNSFP"
              final[:pmids] = '.'
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
              end
            end
          else
            # Unknown significance
            @@log.debug("- This variant is a VUS because it does not have enough information.")
            final[:pathogenicity] = clinical_labels['unknown']
            final[:reason] = "Not enough information"
            final[:comments] = "This variant is a VUS because it does not have enough information."
          end
  
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
            "#{@final_comments_tag}=#{final[:comments]}",
            "#{@final_pathogenicity_source_tag}=#{final[:source]}",
            "#{@final_pathogenicity_reason_tag}=#{final[:reason]}",
            "#{@final_diseases_tag}=#{final[:diseases]}",
            "#{@clinvar_hgmd_conflict_tag}=#{final[:clinvar_hgmd_conflict]}",
          ].join(";")
          
          # Print updated VCF row
          f.puts vcf_cols.join("\t")
  
          @@log.debug("- This variant has been labeled \"#{final[:pathogenicity]}\"")
          @@log.debug("------------------------------------------------------")
        end
      end # <-- End parsing bcftools result
    end # <-- End printing output file

    @@log.info("Final pathogenicity added to #{tmp_output_file}")

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
  end

  ##
  # Cleanup header
  #
  # NOTE: Do not use yet
  ##
  def cleanup_header(vcf_file:, out_file_prefix:)
    @@log.info("Cleaning up VCF header...")

    tmp_header_file = "#{out_file_prefix}.header.tmp.txt"
    tmp_output_file = "#{out_file_prefix}.tmp.vcf.gz"

    # Create header file
    File.open(tmp_header_file, 'w') do |f|
      `bcftools view \
         --header-only \
         #{vcf_file}`
      .each_line do |line|
        line.chomp!
        if line.match(/^#(?:#fileformat|#INFO|CHROM)=/)
          f.puts line
        end
      end
      f.puts ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'].join("\t")
    end
    
    # Re-header
#    `bcftools reheader \
#       --header #{tmp_header_file} \
#       --output-type u \
#       #{vcf_file} \
#     | bcftools view \
#         --output-type z \
#         --output-file #{tmp_output_file}`

puts "Beginning 'bcftools reheader' on #{vcf_file}..."
    `bcftools reheader \
       --header #{tmp_header_file} \
       --output #{tmp_output_file} \
       #{vcf_file}`

    # NOTE: This is causing issues right now
puts "Beginning 'bcftools view' conversion on #{tmp_output_file}..."
    @cleanup_header_result = "#{out_file_prefix}.vcf.gz"
    `bcftools view \
       --output-type z \
       --output-file #{@cleanup_header_result} \
       #{tmp_output_file}`

    # Move tmp file to output file
    #File.rename(tmp_output_file, @cleanup_header_result)
    @@log.info("Output written to #{@cleanup_header_result}")

    # Index output file
    @@log.info("Indexing #{@cleanup_header_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@cleanup_header_result}`
    @@log.info("Done creating index file")

#    # Remove tmp files
#    @@log.info("Removing tmp files...")
#    File.unlink(tmp_header_file) if File.exist?(tmp_header_file)
#    File.unlink(tmp_output_file) if File.exist?(tmp_output_file)
#    @@log.info("Done removing tmp files")
  end

  ##
  # Test
  ##
  def test(vcf_file:, out_file_prefix:, assertion_tags:)
    # Create bcftools query string
    fields = []
    assertion_tags.each do |a|
      fields << "%INFO/#{a}"
      fields << "%INFO/ASSERT_#{a}"
    end
    
    # Produce output that will be tested
    stdout, stderr = Open3.capture3(
      "bcftools query \
         --format '%CHROM\\t%POS\\t%REF\\t%ALT\t#{fields.join('\t')}\\n' \
         #{vcf_file}"
    )

    # Did bcftools produce an error?
    unless stderr.empty?
      @@log.error("bcftools was not able to compress #{F_IN}...") 
      @@log.error("bcftools error is: #{stderr}") 
      abort
    end

    total_num_tests = 0  # Total number of tests
    num_tests_passed = 0 # Number of passed tests

    # Go forth... do testing
    puts "BEGIN ASSERTION TESTING".blue
    stdout.each_line do |line|
      fields = line.chomp.split("\t")
      puts "- " + "TESTING VARIANT: ".blue + fields[0..3].join("\t")

      j = 0 # Set assertion_tags iterator

      # Compare field pairs
      # ...Start from fields[4] in order to skip CHROM, POS, REF, ALT
      (4..(fields.length-1)).step(2) do |i|
        total_num_tests += 1
        if fields[i] == fields[i+1]
          # PASS - fields are the same
          puts "  - " + "PASS: ".green + "[#{assertion_tags[j]}] #{URI.unescape(fields[i])} == #{URI.unescape(fields[i+1])} [ASSERT_#{assertion_tags[j]}]"
          num_tests_passed += 1
        else
          # FAIL - fields are not the same
          puts "  - " + "FAIL: ".red + "[#{assertion_tags[j]}] #{URI.unescape(fields[i])} != #{URI.unescape(fields[i+1])} [ASSERT_#{assertion_tags[j]}]"
        end
        j += 1
      end
    end

    # Print verdict (including num. tests passed)
    puts "END ASSERTION TESTING".blue
    if num_tests_passed == total_num_tests
      # PASS
      puts "- " + "VERDICT: ".blue + "PASS :)".green
      puts "  - " + "#{num_tests_passed}/#{total_num_tests} tests passed".green
    else
      # FAIL
      puts "- " + "VERDICT: ".blue + "FAIL :(".red
      puts "  - " + "#{num_tests_passed}/#{total_num_tests} tests passed".red
    end
  end
end
