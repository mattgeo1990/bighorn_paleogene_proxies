# BHB_d18O_alteration_probability.R
# Shared one-dimensional alteration-probability mapping for BHB carbonates.
#
# Interpretation:
#   - the pooled mean d18Ocarb of BHB pedogenic-micrite T47 observations is
#     assigned P(altered) = 0.05;
#   - 20 per mil VSMOW is assigned P(altered) = 0.95;
#   - probabilities are interpolated linearly between those anchors;
#   - values below 20 per mil extrapolate above 0.95 and are capped at 1;
#   - values above the pooled mean are conservatively capped at 0.05.
#
# This is a transparent screening index along the low-d18O direction expected
# for progressive water-rock exchange. It is not a calibrated posterior
# probability and does not include qualitative petrographic priors, D47-D48
# residuals, or temperature plausibility.

calc_d18O_alteration_probability <- function(
    d18Ocarb_vsmow,
    reference_mean_vsmow,
    altered_anchor_vsmow = 20,
    probability_at_mean = 0.05,
    probability_at_altered_anchor = 0.95
) {
  if (
    !is.finite(reference_mean_vsmow) ||
      reference_mean_vsmow <= altered_anchor_vsmow
  ) {
    stop("reference_mean_vsmow must exceed altered_anchor_vsmow.")
  }

  raw_probability <- probability_at_mean +
    (probability_at_altered_anchor - probability_at_mean) *
      (reference_mean_vsmow - d18Ocarb_vsmow) /
      (reference_mean_vsmow - altered_anchor_vsmow)

  dplyr::case_when(
    !is.finite(d18Ocarb_vsmow) ~ NA_real_,
    d18Ocarb_vsmow >= reference_mean_vsmow ~ probability_at_mean,
    TRUE ~ pmin(pmax(raw_probability, probability_at_mean), 1)
  )
}
