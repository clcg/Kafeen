require 'command'

describe Command do
  let(:genes_file) {'spec/fixtures/DVD_v8_genes.2016-09-23.txt'}
  let(:ref_file) {'spec/fixtures/hg19_gene_and_mitochondrial_regions.v5.bed'}

  describe "#genes2regions" do
    context "with options" do
      it "sets an output VCF file name" do
        no_output do
          subject.genes2regions(genes_file: genes_file, ref_file: ref_file, out_file_prefix: 'my_genes')
        end
        expect(subject.genes2regions_result).to eq('my_genes.gene_regions.bed')
      end
      it "creates an output VCF file"
    end
  end

  describe "#regions2variants" do
    it "creates an output VCF file"
  end

  describe "#add_genes" do
    it "creates an output VCF file"
  end

  describe "#add_predictions" do
    it "creates an output VCF file"
  end
end
