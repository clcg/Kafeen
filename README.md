# Kafeen

## Description

Kafeen is a data collection pipeline used to aggregate and interpret variants and annotations from various resources.

Please see the [MORL Wiki](https://wiki.uiowa.edu/pages/viewpage.action?pageId=160538797) for full documentation on Kafeen.

## Requirements

- Ruby 2.0+
   - If your Ruby version is not up to date, it highly recommended to use [RVM](https://github.com/rvm/rvm) to safely install another version on your machine. No root permissions required.
- `bcftools` version 1.3+ must be in your `$PATH` (download [here](https://github.com/samtools/bcftools/releases))
- `bedtools` version 2.0+ must be in your `$PATH` (download [here](https://github.com/arq5x/bedtools2/releases))
- `tabix` and `bgzip` must be in your `$PATH` (download/compile them as part of the [HTSlib package](https://github.com/samtools/htslib/releases))
- Java 1.7+ is required if running ASAP

## Author

Sean Ephraim
