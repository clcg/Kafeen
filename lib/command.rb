# A library of available pipeline commands
#
# @author Sean Ephraim

require 'logger'
require 'open3'
require 'uri'

class Command
  # Result file paths
  attr_reader :genes2regions_result
  attr_reader :regions2variants_result
  attr_reader :addgenes_result
  attr_reader :addpredictions_result
  attr_reader :finalize_pathogenicity_result

  ##
  # Initialize
  ##
  def initialize(log_level: 'INFO', log_out: STDOUT)
    # Set custom VCF tags that will be added
    @gerp_pred_tag = "GERP++_PRED"
    @phylop20way_mammalian_pred_tag = "PHYLOP20WAY_MAMMALIAN_PRED"
    @num_path_preds_tag = "NUM_PATH_PREDS"
    @total_num_preds_tag = "TOTAL_NUM_PREDS"
    @final_pred_tag = "FINAL_PRED"
    @final_pathogenicity_tag = "FINAL_PATHOGENICITY"
    @final_pathogenicity_comments_tag = "FINAL_NOTES"
    @final_pathogenicity_source_tag = "FINAL_PATHOGENICITY_SOURCE"
    @final_pathogenicity_reason_tag = "FINAL_PATHOGENICITY_REASON"
    @final_disease_tag = "FINAL_DISEASE"

    # Set logger
    @@log = Logger.new(log_out)
    if log_level == 'UNKNOWN'
      @@log.level = Logger::UNKNOWN
    elsif log_level == 'FATAL'
      @@log.level = Logger::FATAL
    elsif log_level == 'ERROR'
      @@log.level = Logger::ERROR
    elsif log_level == 'WARN'
      @@log.level = Logger::WARN
    elsif log_level == 'DEBUG'
      @@log.level = Logger::DEBUG
    else
      @@log.level = Logger::INFO
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
  # Regions to Variants
  #
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
        elsif vcf['fields'] == '*'
          # Keep all INFO tags
          fields = '*'
        else
          # Keep only the following INFO tags (indicated by ^)
          fields = "^" + vcf['fields'].map { |f| "INFO/#{f}" }.join(',')
        end

        # Query...
        @@log.info("Querying #{vcf['source']}...")
        if fields != '*'
          stdout, stderr = Open3.capture3(
            "bcftools annotate \
               --remove '#{fields}' \
               --regions-file '#{bed_file}' \
               --exclude 'TYPE=\"other\"' \
               #{vcf['filename']} \
             | bcftools norm \
                 --multiallelics '-' \
                 --output-type z \
                 --output #{tmp_source_vcf}"
          )
          stderr = stderr.sub(/^Lines total\/modified\/skipped.*/, '').strip # remove unnecessary "error"
          if !stderr.empty?
            @@log.warn("bcftools returned an error for #{vcf['source']}. Trying another query method...")
          end
        end

        # Did bcftools return an error?
        # Try again and don't remove any INFO tags this time
        if fields == '*' || !stderr.empty?
          stdout, stderr = Open3.capture3(
            "bcftools view \
               --regions-file '#{bed_file}' \
               --exclude 'TYPE=\"other\"' \
               #{vcf['filename']} \
             | bcftools norm \
                 --multiallelics '-' \
                 --output-type z \
                 --output #{tmp_source_vcf}"
          )
          stderr = stderr.sub(/^Lines total\/modified\/skipped.*/, '').strip # remove unnecessary "error"
        end

        # Index the results file...
        if !stderr.empty?
          # ERROR
          @@log.error("bcftools was not able to query #{vcf['source']}. Please check that file name and INFO tags are set correctly in your config file.")
        else
          # SUCCESS -- create index file
          @@log.info("Successfully queried #{vcf['source']}")
          @@log.info("Creating index file for #{vcf['source']}...")
          `bcftools index --force --tbi #{tmp_source_vcf}`
          @@log.info("Done creating index file")

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
  def addgenes(bed_file:, vcf_file:, out_file_prefix:)
    # Prepare header file
    header_file = "#{out_file_prefix}.header.tmp.txt"
    header_line = '##INFO=<ID=GENE,Number=1,Type=String,Description="HGNC gene symbol">'
    File.open(header_file, 'w') {|f| f.write(header_line) }

    # Prepare BED file using bgzip and tabix
    `bgzip -c #{bed_file} > #{bed_file}.tmp.gz`
    `tabix -fp bed #{bed_file}.tmp.gz`

    # Add genes to VCF file
    tmp_output_file = "#{out_file_prefix}.tmp.vcf.gz"
    @@log.info("Adding gene annotations to #{tmp_output_file}")
    `bcftools annotate \
       --annotations #{bed_file}.tmp.gz \
       --columns CHROM,FROM,TO,GENE \
       --header-lines #{header_file} \
       --output #{tmp_output_file} \
       --output-type z \
       #{vcf_file}`
    @@log.info("Genes added to #{tmp_output_file}")

    # Move tmp output to be the new output file
    @addgenes_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Moving output to #{@addgenes_result}...")
    File.rename(tmp_output_file, @addgenes_result)

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
  def addpredictions(dbnsfp_file:, vcf_file:, bed_file:, out_file_prefix:, clinical_labels:)
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
    f_tmp_output_file = File.open(tmp_output_file, 'w')
    @@log.info("Adding dbNSFP pedictions to #{tmp_output_file}...")
    `bcftools annotate \
       --annotations #{dbnsfp_subset_file} \
       --columns #{dbnsfp_file['fields'].map { |f| "INFO/#{f}" }.join(',')} \
       --output-type v \
       #{vcf_file}`.each_line do |vcf_row|
         vcf_row.chomp!
         if vcf_row.match(/^##/)
           # Print meta-info
           f_tmp_output_file.puts vcf_row
         elsif vcf_row.match(/^#[^#]/)
           f_tmp_output_file.puts "##INFO=<ID=#{@num_path_preds_tag},Number=.,Type=String,Description=\"Number of pathogenic predictions from dbNSFP\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@total_num_preds_tag},Number=.,Type=String,Description=\"Total number of prediction scores available from dbNSFP\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@final_pred_tag},Number=.,Type=String,Description=\"Final prediction consensus based on majority vote of prediction scores\">"

           # Add GERP++ prediction tag to meta-info
           if dbnsfp_file['fields'].any?{ |e| e == 'GERP++_RS' }
             f_tmp_output_file.puts "##INFO=<ID=#{@gerp_pred_tag},Number=.,Type=String,Description=\"NA\">"
           end
           # Add phyloP20way mammalian prediction tag to meta-info
           if dbnsfp_file['fields'].any?{ |e| e == 'PHYLOP20WAY_MAMMALIAN' }
             f_tmp_output_file.puts "##INFO=<ID=#{@phylop20way_mammalian_pred_tag},Number=.,Type=String,Description=\"NA\">"
           end
           # Print header
           f_tmp_output_file.puts vcf_row
         else
           vcf_cols = vcf_row.split("\t")

           # Analyze each *_PRED field (as well as GERP++ and phyloP)
           # Tally up pathogenic predictions
           total_num_preds = 0
           num_path_preds = 0
           dbnsfp_file['fields'].select { |e| e.match(/(?:_PRED$|^GERP\+\+_RS$|^PHYLOP20WAY_MAMMALIAN$)/i) }.each do |field|
             # Get all predictions for this algorithm
             match = vcf_row.match(/(?:^|[\t;])#{Regexp.escape(field)}=([^;\t]*)/)

             # No data for this algorithm -- skip it
             next if match.nil?

             # Get all predictions for this algorithm
             preds = match[1].split(/[^a-zA-Z0-9.-]+/)

             # No data for this algorithm -- skip it
             next if preds.all? { |pred| pred == '.' || pred == 'U' }
               
             if field == 'SIFT_PRED'
               # SIFT prediction
               num_path_preds += 1 if preds.include?('D') # <-- Damaging
               total_num_preds += 1
             elsif field == 'POLYPHEN2_HDIV_PRED'
               # Polyphen2 (HDIV) prediction
               num_path_preds += 1 if preds.include?('D') || preds.include?('P') # <-- "Deleterious" or "Possibly damaging"
               total_num_preds += 1
             elsif field == 'LRT_PRED'
               # LRT prediction
               num_path_preds += 1 if preds.include?('D') # <-- "Deleterious"
               total_num_preds += 1
             elsif field == 'MUTATIONTASTER_PRED'
               # MutationTaster prediction
               num_path_preds += 1 if preds.include?('D') || preds.include?('A') # <-- "Disease-causing" or "Disease-causing (automatic)"
               total_num_preds += 1
             elsif field == 'GERP++_RS'
               # GERP++ prediction
               if preds.any? { |pred| pred.to_f > 0.0 }
                 # Conserved
                 num_path_preds += 1
                 vcf_cols[7] = [vcf_cols[7], "#{@gerp_pred_tag}=C"].join(";")
               else
                 # Non-conserved
                 vcf_cols[7] = [vcf_cols[7], "#{@gerp_pred_tag}=N"].join(";")
               end
               total_num_preds += 1
             elsif field == 'PHYLOP20WAY_MAMMALIAN'
               # phyloP20way mammalian prediction
               if preds.any? { |pred| pred.to_f >= 0.95 }
                 # Conserved
                 num_path_preds += 1
                 vcf_cols[7] = [vcf_cols[7], "#{@phylop20way_mammalian_pred_tag}=C"].join(";")
               else
                 # Non-conserved
                 vcf_cols[7] = [vcf_cols[7], "#{@phylop20way_mammalian_pred_tag}=N"].join(";")
               end
               total_num_preds += 1
             end
           end

           # Add final prediction
           if total_num_preds > 5
             path_score = num_path_preds.to_f/total_num_preds.to_f
             if path_score >= 0.6
               # Predicted pathogenic
               final_pred = URI.escape(clinical_labels['pred_pathogenic'])
             elsif path_score <= 0.4
               # Predicted benign
               final_pred = URI.escape(clinical_labels['pred_benign'])
             else
               # Predicted unknown (benign predictions approx. equal to pathogenic)
               final_pred = URI.escape(clinical_labels['unknown'])
             end
           else
             # Predicted unknown (not enough predictions)
             final_pred = URI.escape(clinical_labels['unknown'])
           end

           # Update INFO column
           vcf_cols[7] = [vcf_cols[7], "#{@num_path_preds_tag}=#{num_path_preds}", "#{@total_num_preds_tag}=#{total_num_preds}", "#{@final_pred_tag}=#{final_pred}"].join(";")

           # Print updated VCF row
           f_tmp_output_file.puts vcf_cols.join("\t")
         end
       end

    f_tmp_output_file.close
    @@log.info("Predictions added to #{tmp_output_file}")

    @addpredictions_result = "#{out_file_prefix}.vcf.gz"
    @@log.info("Compressing #{tmp_output_file}")
    # Compress the output file
    `bcftools view \
       --output-type z \
       --output-file #{@addpredictions_result} \
       #{tmp_output_file}`
    @@log.info("Compressed output written to #{@addpredictions_result}...")

    # Index output file
    @@log.info("Indexing #{@addpredictions_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@addpredictions_result}`
    @@log.info("Done creating index file")

    @@log.info("Removing tmp files...")
    File.unlink(tmp_output_file) if File.exist?(tmp_output_file)
    File.unlink(dbnsfp_subset_file) if File.exist?(dbnsfp_subset_file)
    @@log.info("Done removing tmp files")
  end

  ##
  # Finalize Pathogenicity
  ##
  def finalize_pathogenicity(vcf_file:, out_file_prefix:, clinical_labels:)
    tmp_output_file = "#{out_file_prefix}.tmp.vcf"
    f_tmp_output_file = File.open(tmp_output_file, 'w')

    # Set ClinVar pathogenicity dictionary
    clinvar_pathogenicity_map = {
      '2' => clinical_labels['benign'],
      '3' => clinical_labels['likely_benign'],
      '4' => clinical_labels['likely_pathogenic'],
      '5' => clinical_labels['pathogenic'],
      'Benign' => clinical_labels['benign'],
      'Likely benign' => clinical_labels['likely_benign'],
      'Likely pathogenic' => clinical_labels['likely_pathogenic'],
      'Pathogenic' => clinical_labels['pathogenic'],
    }

    # Clinical confidence/review weight
    clinvar_confidence_weight = {
      'not'    => 0,
      'single' => 1,
      'mult'   => 2,
      'exp'    => 3,
      'prof'   => 4,
      'classified by multiple submitters' => 2,
    }

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
    }

    # HGMD confidence weight
    hgmd_confidence_weight = {
      'Low'  => 1,
      'High' => 2,
    }

    @@log.info("Finalizing pathogenicity... #{tmp_output_file}...")
    `bcftools view \
       --output-type v \
       #{vcf_file}`.each_line do |vcf_row|
         vcf_row.chomp!
         if vcf_row.match(/^##/)
           # Print meta-info
           f_tmp_output_file.puts vcf_row
         elsif vcf_row.match(/^#[^#]/)
           # Add final pathogneicity tag to meta-info
           f_tmp_output_file.puts "##INFO=<ID=#{@final_pathogenicity_tag},Number=.,Type=String,Description=\"Final curated pathogenicity\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@final_pathogenicity_comments_tag},Number=.,Type=String,Description=\"Additional comments from curator\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@final_pathogenicity_source_tag},Number=.,Type=String,Description=\"Source for final pathogenicity\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@final_pathogenicity_reason_tag},Number=.,Type=String,Description=\"Reason for final pathogenicity\">"
           f_tmp_output_file.puts "##INFO=<ID=#{@final_disease_tag},Number=.,Type=String,Description=\"Final curated disease\">"

           # Print header
           f_tmp_output_file.puts vcf_row
         else
           vcf_cols = vcf_row.split("\t")
           # TODO Finalize here
           if (vcf_cols[7].match(/(?:^|[\t;])(?:MORL|CURATED)_PATHOGENICITY=([^;\t]*)/)).nil?
             # Check for curated pathogenicity
             final[:pathogenicity] = vcf_cols[7].scan(/(?:^|[\t;])(?:MORL|CURATED)_PATHOGENICITY=([^;\t]*)/)
             final[:disease] = vcf_cols[7].scan(/(?:^|[\t;])(?:MORL|CURATED)_PATHOGENICITY=([^;\t]*)/)
             final[:pmid] = vcf_cols[7].scan(/(?:^|[\t;])(?:MORL|CURATED)_PMID=([^;\t]*)/)
             final[:source] = "Expert-curated"
             final[:reason] = URI.escape("This variant has been expertly curated.")
             final[:comments] = vcf_cols[7].scan(/(?:^|[\t;])(?:MORL|CURATED)_COMMENTS=([^;\t]*)/)
           elsif vcf_cols[7].scan(/[^;\t]*_AF=([^;\t]*)/).flatten.any? { |af| af.to_f >= 0.005 } == true
             # Check if max MAF >= 0.005
             final[:pathogenicity] = URI.escape(clinical_labels['benign'])
             final[:disease] = "."
             final[:pmid] = "."
             final[:source] = "MAF"
             final[:reason] = "MAF>=0.005"
             final[:comments] = URI.escape("This variant contains a MAF greater than 0.005 in at least one population.")
           elsif (vcf_cols[7].match(/(?:^|[\t;])(?:CLNSIG|(?:HGMD_)?VARIANTTYPE)=([^;\t]*)/)).nil?
             # TODO Check for HGMD / ClinVar
             
             # ClinVar fields
             clinvar[:pathogenicities] = vcf_cols[7].scan(/(?:^|[\t;])CLNSIG=([^;\t]*)/).flatten[0].to_s.split(/[,|]/)
             clinvar[:diseases] = vcf_cols[7].scan(/(?:^|[\t;])CLNDBN=([^;\t]*)/).flatten[0].to_s.split(/[,|]/)
             clinvar[:pmid] = "."
             clinvar[:confidence] = vcf_cols[7].scan(/(?:^|[\t;])CLNREVSTAT=([^;\t]*)/).flatten[0].to_s.split(/[,|]/)

             ## TODO Parse ClinVar result
             #if clinvar[:pathogenicities].empty?
             #  
             #end

             # HGMD fields
             hgmd[:pathogenicity] = vcf_cols[7].scan(/(?:^|[\t;])(?:HGMD_)?VARIANTTYPE=([^;\t]*)/).flatten[0]
             hgmd[:disease] = vcf_cols[7].scan(/(?:^|[\t;])(?:HGMD_)?DISEASE=([^;\t]*)/).flatten[0]
             hgmd[:pmid] = vcf_cols[7].scan(/(?:^|[\t;])(?:HGMD_)?PMID=([^;\t]*)/).flatten[0]
             hgmd[:confidence] = vcf_cols[7].scan(/(?:^|[\t;])(?:HGMD_)?CONFIDENCE=([^;\t]*)/).flatten[0]

             if !clinvar[:pathogenicity].nil? && hgmd[:pathogenicity].nil?
               # Only found in ClinVar
               final[:pathogenicity] = URI.escape(clinical_labels['unknown']) # TODO <-- change this
               final[:source] = "ClinVar"
               final[:pmid] = "."
               final[:reason] = URI.escape("Found in ClinVar, not in HGMD")
               final[:comments] = "."
             elsif !hgmd[:pathogenicity].nil? && clinvar[:pathogenicity].nil?
               # Only found in HGMD
               final[:pathogenicity] = URI.escape(clinical_labels['unknown']) # TODO <-- change this
               final[:source] = "HGMD"
               final[:pmid] = vcf_cols[7].scan(/(?:^|[\t;])PMID=([^;\t]*)/).flatten[0]
               final[:reason] = URI.escape("Found in HGMD, not in ClinVar")
               final[:comments] = "."
             elsif !clinvar[:pathogenicity].nil? && !hgmd[:pathogenicity].nil?
               # TODO Found in ClinVar and HGMD
             else
               # Not found in ClinVar or HGMD
               final[:pathogenicity] = URI.escape(clinical_labels['unknown'])
               final[:source] = "."
               final[:reason] = URI.escape("Not enough information")
             end
           elsif (match = vcf_cols[7].match(/(?:^|[\t;])#{@final_pred_tag}=([^;\t]*)/)).nil?
             # TODO Use dbNSFP prediction
           else
             # Unknown significance
             final[:pathogenicity] = URI.escape(clinical_labels['unknown'])
             final[:source] = "."
             final[:reason] = URI.escape("Not enough information")
           end
         end
       end

    f_tmp_output_file.close

    @@log.info("Predictions added to #{tmp_output_file}")

    @finalize_pathogenicity_result = "#{out_file_prefix}.vcf.gz}"
    @@log.info("Compressing #{tmp_output_file}")
    # Compress the output file
    `bcftools view \
       --output-type z \
       --output-file #{@finalize_pathogenicity_result} \
       #{tmp_output_file}`
    @@log.info("Compressed output written to #{@finalize_pathogenicity_result}...")

    # Index output file
    @@log.info("Indexing #{@finalize_pathogenicity_result}...")
    `bcftools index  \
       --force \
       --tbi \
       #{@finalize_pathogenicity_result}`
    @@log.info("Done creating index file")

    @@log.info("Removing tmp files...")
    File.unlink(tmp_output_file) if File.exist?(tmp_output_file)
    @@log.info("Done removing tmp files")
  end
end
