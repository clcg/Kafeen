# A library of available pipeline commands
#
# @author Sean Ephraim
class Command
  attr_reader :genes2regions_result

  # Input columns:
  #   gene_symbol
  #
  # Example usage:
  #   ruby get_gene_regions.rb mygenes.txt > myregions.txt
  def genes2regions(genes_file:, ref_file:, out_file: nil)
    # Gene region reference file
    # Reference columns:
    #   chr, start, stop, gene_symbol
    f_regions = File.open(ref_file, 'r')

    # Set output file
    if out_file.nil?
      out_file = genes_file.sub(/#{File.extname(genes_file)}/, '') + '.gene_regions'
      f_out = File.open(out_file, 'w')
    end
    
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
    @genes2regions_result = out_file
  end
end
