boot_fail = false

## Configurations
#require_relative File.join('..', 'config', 'config')
#require_relative File.join('..', 'config', 'constants')
## Libraries
require_relative 'command'
#require_relative 'app'
#require_relative 'asap'
#require_relative 'tabix'
#require_relative 'tarball'
#require_relative 'vcfrow'
#require_relative 'trollop'
# Built-in gems
require 'yaml'

# Load configurations
CONFIG = YAML.load_file(File.join('config', 'config.yml'))

# Verify required utilities are in $PATH
utils = ['bcftools', 'tabix', 'bgzip']
utils.each do |util|
  if `which #{util} 2> /dev/null`.empty?
    puts "ERROR: #{util} is not in your $PATH, which is required to run this tool"
    boot_fail = true
  end
end
exit if boot_fail
