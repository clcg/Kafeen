boot_fail = false

# Require local libraries
require_relative 'command'
require_relative 'trollop'

# Require built-in gems
require 'yaml'

# Parse options
opts = Trollop::options do
  opt :config, "Configuration override file", :type => :string, :default => nil
  opt :output_prefix, "Prefix to use for output files", :type => :string, :default => nil
#  opt :debug, "Print debugging statements", :default => false
end

# Validate options
Trollop::die :config, "must exist" unless File.exist?(opts[:config]) if opts[:config]

# Set input file
F_IN = ARGV[0]

# Load configuration file
default_config = YAML.load_file(File.join('config', 'config.yml'))
if opts[:config].nil?
  # Use default
  CONFIG = default_config
else
  # Use config override
  CONFIG = default_config.merge(YAML.load_file(opts[:config]))
end

# Set output file prefix
if opts[:output_prefix].nil?
  # Use default
  FILE_PREFIX = F_IN.sub(/#{File.extname(F_IN)}$/, '')
else
  # Use user-specified prefix
  FILE_PREFIX = opts[:output_prefix]
end

# Verify required utilities are in $PATH
utils = ['bcftools', 'tabix', 'bgzip', 'sort', 'java']
utils.each do |util|
  if `which #{util} 2> /dev/null`.empty?
    # Print to stderr
    STDERR.puts "ERROR: #{util} is not in your $PATH, which is required to run this tool"
    boot_fail = true
  end
end
abort if boot_fail
