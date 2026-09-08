# Shared clumped-isotope temperature calibrations.

# Anderson et al. (2021) unified I-CDES calibration, reported at 90 C:
# D47 = 0.0391 * 10^6 / T(K)^2 + 0.154
anderson_2021_T47_C <- function(D47_iCDES_90) {
  result <- rep(NA_real_, length(D47_iCDES_90))
  valid <- is.finite(D47_iCDES_90) & D47_iCDES_90 > 0.154
  result[valid] <- sqrt(39100 / (D47_iCDES_90[valid] - 0.154)) - 273.15
  result
}
