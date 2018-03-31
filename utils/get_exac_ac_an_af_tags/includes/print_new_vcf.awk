#!/bin/awk -f

# NOTE: This script requires lib.awk

##
# Input fields:
# $1:  CHROM
# $2:  POS
# $3:  ID
# $4:  REF
# $5:  ALT
# $6:  QUAL
# $7:  FILTER
# $8:  AC_Adj
# $9:  AN_Adj
# $10: AC_AFR
# $11: AN_AFR
# $12: AC_AMR
# $13: AN_AMR
# $14: AC_EAS
# $15: AN_EAS
# $16: AC_FIN
# $17: AN_FIN
# $18: AC_NFE
# $19: AN_NFE
# $20: AC_OTH
# $21: AN_OTH
# $22: AC_SAS
# $23: AN_SAS
##

BEGIN {
  FS="\t"
  OFS="\t"
  # Print the VCF header
  print_header("includes/new_header.txt")
}

# Print variant records
/^[^#]/ {
  # Get AC, AN, and AF for each population
  get_ac_an_af($8, $9, all_ac_an_af)
  get_ac_an_af($10, $11, afr_ac_an_af)
  get_ac_an_af($12, $13, amr_ac_an_af)
  get_ac_an_af($14, $15, eas_ac_an_af)
  get_ac_an_af($16, $17, fin_ac_an_af)
  get_ac_an_af($18, $19, nfe_ac_an_af)
  get_ac_an_af($20, $21, oth_ac_an_af)
  get_ac_an_af($22, $23, sas_ac_an_af)

  # Print record
  print $1, $2, $3, $4, $5, $6, $7, "EXAC_ALL_AC=" all_ac_an_af[1] ";EXAC_ALL_AN=" all_ac_an_af[2] ";EXAC_ALL_AF=" all_ac_an_af[3] \
                                    ";EXAC_AFR_AC=" afr_ac_an_af[1] ";EXAC_AFR_AN=" afr_ac_an_af[2] ";EXAC_AFR_AF=" afr_ac_an_af[3] \
                                    ";EXAC_AMR_AC=" amr_ac_an_af[1] ";EXAC_AMR_AN=" amr_ac_an_af[2] ";EXAC_AMR_AF=" amr_ac_an_af[3] \
                                    ";EXAC_EAS_AC=" eas_ac_an_af[1] ";EXAC_EAS_AN=" eas_ac_an_af[2] ";EXAC_EAS_AF=" eas_ac_an_af[3] \
                                    ";EXAC_FIN_AC=" fin_ac_an_af[1] ";EXAC_FIN_AN=" fin_ac_an_af[2] ";EXAC_FIN_AF=" fin_ac_an_af[3] \
                                    ";EXAC_NFE_AC=" nfe_ac_an_af[1] ";EXAC_NFE_AN=" nfe_ac_an_af[2] ";EXAC_NFE_AF=" nfe_ac_an_af[3] \
                                    ";EXAC_OTH_AC=" oth_ac_an_af[1] ";EXAC_OTH_AN=" oth_ac_an_af[2] ";EXAC_OTH_AF=" oth_ac_an_af[3] \
                                    ";EXAC_SAS_AC=" sas_ac_an_af[1] ";EXAC_SAS_AN=" sas_ac_an_af[2] ";EXAC_SAS_AF=" sas_ac_an_af[3]
}
