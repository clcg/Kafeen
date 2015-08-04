# A library of available pipeline commands
#
# @author Sean Ephraim
class Command
  attr_reader :genes2regions_result
  attr_reader :regions2variants_result
  attr_reader :addgenes_result

  ##
  # Input columns:
  #   gene_symbol
  #
  # Example usage:
  #   ruby get_gene_regions.rb mygenes.txt > myregions.txt
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
  end

  ##
  # Get a list of all variants within specified regions
  #
  # Input file format (tab-separated):
  #   chr  pos_start  pos_end  gene
  ##
  def regions2variants(bed_file:, vcf_files:, out_file_prefix:, keep_tmp_files: false)
    tmp_vcfs = {}
    File.open(bed_file).each do |region|
      chr,pos_start,pos_end,gene = region.chomp.split("\t")
      chr.sub!('chr', '')

      # Query all VCF files for variants
      vcf_files.each do |key, vcf|
        tmp_source_vcf = "#{out_file_prefix}.#{vcf['source']}.tmp.vcf.gz"
        if vcf['fields'].nil?
          # Remove all INFO tags
          fields = 'INFO'
        else
          # Keep only the following INFO tags (indicated by ^)
          fields = "^" + vcf['fields'].map { |f| "INFO/#{f}" }.join(',')
        end

        # Query...
        `bcftools annotate \
           --remove '#{fields}' \
           --regions-file '#{bed_file}' \
           #{vcf['filename']} | \
         bcftools filter \
           --exclude 'ALT !~ "^[actgACGT,]*$"' \
           --output #{tmp_source_vcf} \
           --output-type z`
        
        # Index the results file...
        `bcftools index --force --tbi #{tmp_source_vcf}`

        # Store tmp file name (filename) and the original VCF that the data came from (parent)
        tmp_vcfs[key] = {'filename' => tmp_source_vcf, 'parent' => vcf['filename']}
      end
    end

    # Construct list of VCFs to merge
    # NOTE: DO NOT MERGE dbNSFP - ONLY ANNOTATE WITH IT
    files_to_merge = []
    tmp_vcfs.each do |key, tmp_vcf|
      if key != 'dbnsfp'
        files_to_merge << tmp_vcf['filename']
      end
    end
    files_to_merge << tmp_vcfs['dbsnp']['filename']

    # Merge VCFs...
    @regions2variants_result = "#{out_file_prefix}.vcf.gz"
    `bcftools merge \
       --merge all \
       --output #{@regions2variants_result} \
       --output-type z \
       #{files_to_merge.join(' ')}`

    # Remove tmp files
    if !keep_tmp_files
      tmp_vcfs.each do |key, tmp_vcf|
        File.unlink(tmp_vcf['filename']) if File.exist?(tmp_vcf['filename'])
        File.unlink("#{tmp_vcf['filename']}.tbi") if File.exist?("#{tmp_vcf['filename']}.tbi")
      end
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
    `tabix -p bed #{bed_file}.tmp.gz`

    # Add genes to VCF file
    @addgenes_result = "#{out_file_prefix}.vcf.gz"
    `bcftools annotate \
       --annotations #{bed_file}.tmp.gz \
       --columns CHROM,FROM,TO,GENE \
       --header-lines #{header_file} \
       --output #{@addgenes_result} \
       --output-type z \
       #{vcf_file}`

    # Remove tmp files
    File.unlink(header_file) if File.exist?(header_file)
    File.unlink("#{bed_file}.tmp.gz") if File.exist?("#{bed_file}.tmp.gz")
    File.unlink("#{bed_file}.tmp.gz.tbi") if File.exist?("#{bed_file}.tmp.gz.tbi")
  end
end
