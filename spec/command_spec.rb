require 'command'
require 'pp'

# Prefix for all output files
PREFIX = "/tmp/_kafeen.spec"

# Data for each MAF population
MAF_EXAMPLES = [
  {:key => 'g1000', :tag => '1KG', :example => 'adds MAFs from 1000 Genomes'},
]

describe Command do
  let(:genes_file) {'spec/fixtures/DVD_v8_genes.2016-09-23.txt'}
  let(:ref_file) {'spec/fixtures/hg19_gene_and_mitochondrial_regions.v5.bed'}
  let(:gjb2_genes_file) {'spec/fixtures/DVD_v8_genes.2016-09-23.txt'}
  let(:gjb2_bed_file) {'spec/fixtures/GJB2.gene_regions.merged.bed'}
  let(:gjb2_vcf) {'spec/fixtures/GJB2.vcf.gz'}
  let(:gjb2_blank_vcf) {'spec/fixtures/GJB2.no_INFO.vcf.gz'}

  out_file_prefix = "#{PREFIX}.GJB2"
  f_test_set = "#{out_file_prefix}.test_set.tsv"
  f_kafeen_output = "#{out_file_prefix}.kafeen_output.tsv"

  describe "#genes2regions" do
    context "with options" do
      it "sets two output BED file names" do
        out_file_prefix = "#{PREFIX}.GJB2"
        # Run
        no_output do
          subject.genes2regions(genes_file: gjb2_genes_file, ref_file: gjb2_bed_file, out_file_prefix: out_file_prefix)
        end
        # Test
        expect(subject.genes2regions_result).to eq("#{out_file_prefix}.gene_regions.bed")
        expect(subject.genes2regions_merged_result).to eq("#{out_file_prefix}.gene_regions.merged.bed")
        # Clean up
        remove_created_files(subject.genes2regions_result, subject.genes2regions_merged_result)
      end

      # TODO
      it "creates an output VCF file"

    end
  end

  describe "#regions2variants" do
    # TODO
    it "creates an output VCF file"

    it "adds MAFs from 1000 Genomes" do
      vcf_files = VCF_FILES.select_keys('g1000')

      no_output do
        subject.regions2variants(bed_file: gjb2_bed_file, vcf_files: vcf_files, out_file_prefix: out_file_prefix)
      end

      # Subset Kafeen output for testing
      subset_variants_from_vcf(f_in: subject.regions2variants_result,
                               f_out: f_kafeen_output, tags: '1KG_ALL_AF')

      # Subset GJB2 test set
      subset_variants_from_vcf(f_in: gjb2_vcf, f_out: f_test_set,
                               tags: '1KG_ALL_AF', exclude: '1KG_ALL_AF == "."')

      expect(FileUtils.compare_file(f_test_set, f_kafeen_output)).to be(true)
      remove_created_files(f_test_set, f_kafeen_output)
    end
  end

  describe "#add_genes" do
    it "creates an output VCF file"
  end

  describe "#add_predictions" do
    it "creates an output VCF file"
  end
end
