module CoreExtensions
  module String
    module Vcf
      ##
      # Get VCF Field
      #
      # Example:
      #   vcf_row.get_vcf_field('CLINVAR_PATHOGENICITY')
      #   vcf_row.get_vcf_field(/CURATED_DISEASE[Ss]/i)
      #
      # @param [String/Regexp] tag Name of the INFO tag
      ##
      def get_vcf_field(tag)
        if !tag.is_a? Regexp
          tag = Regexp.escape(tag)
        end 
        return URI.unescape(self.scan(/(?:^|[\t;])#{tag}=([^;\t]*)/).flatten[0].to_s)
      end 
    end
  end
end
