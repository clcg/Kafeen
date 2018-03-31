# Print Header
#
# Import a header from a file and print it.
#
# @author Sean Ephraim
# @param file Name of header file
function print_header(file) {
  # Print header
  while (( getline < file) > 0 ) {
     print
  }
  close(file)
}

# Get AC AN AF
#
# Compute the AF for an ExAC variant based on AC/AN.
#
# Example variant:
# 22 24199764 . AC A,ACCC,ACC,CC,ACCCC,ACCCCC ... AC_AFR=97,282,829,1,16,0;...AN_AFR=2466;...
# Example input (e.g. AC_AFR, AN_AFR, output array name):
#   get_ac_an_af("97,282,829,1,16,0", "2466", "afr_ac_an_af")
# Example output (i.e. array[AC, AN, AF]):
#   afr_ac_an_af["97,282,829,1,16,0", "2466", "0.039334955,0.114355231,0.336171938,0.000405515,0.00648824,0"]
#
# @author Sean Ephraim
# @param ac  Original string containing alternate allele counts (comma-separated)
# @param an  Original string containing total allele counts
# @param ac_an_af  Array to store the final AC, AN, and AF values
function get_ac_an_af(ac, an, ac_an_af) {
  len_ac_array = split(ac, ac_array, ",")
  # Calculate AF
  for (i = 1; i <= len_ac_array; i++) {
    if (i == 1) {
      # Add first AF
      if (an == 0) {
        af = 0 # Can't divide by 0
      }
      else {
        af = ac_array[i]/an
      }
    }
    else {
      # Append any subsequent AFs (comma-separated)
      if (an == 0) {
        af = af "," 0 # Can't divide by 0
      }
      else {
        af = af "," ac_array[i]/an
      }
    }
  }

  # Set AC, AN, AF
  ac_an_af[1] = ac
  ac_an_af[2] = an
  ac_an_af[3] = af
}
