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
CONFIG = YAML.load_file(File.join('config', 'config.yml'))
