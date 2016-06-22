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
    module Colorize
      ##
      # Colorize
      #
      # Print colored text to stdout.
      #
      # @param [Integer] color_code Numeric color code
      ##
      def colorize(color_code)
        "\e[#{color_code}m#{self}\e[0m"
      end
    
      ##
      # Red
      #
      # Print red text to stdout.
      ##
      def red
        colorize(31)
      end

      ##
      # Green
      #
      # Print green text to stdout.
      ##
      def green
        colorize(32)
      end

      ##
      # Yellow
      #
      # Print yellow text to stdout.
      ##
      def yellow
        colorize(33)
      end
    end
  end
end

module CoreExtensions
  module String
    module Colorize
      # colorization
      def colorize(color_code)
        "\e[#{color_code}m#{self}\e[0m"
      end
    
      def red
        colorize(31)
      end
    
      def green
        colorize(32)
      end
    
      def yellow
        colorize(33)
      end
    
      def blue
        colorize(34)
      end
    
      def pink
        colorize(35)
      end
    
      def light_blue
        colorize(36)
      end
    end
  end
end
