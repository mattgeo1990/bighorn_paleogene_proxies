# Bighorn Basin seasonal-temperature evidence database

`BHB_seasonal_temperature_constraints.csv` is the auditable input to
`scripts/08b_synthesize_BHB_seasonal_climate.R`.

The database intentionally distinguishes direct Bighorn Basin proxy evidence,
regional analogs, and climate-model experiments. Values are not treated as
exchangeable measurements. `primary_weight` is a transparent relevance weight,
not a likelihood or an objective measure of study quality.

The resulting distributions are literature-informed Monte Carlo evidence
ensembles, not formal Bayesian posteriors. Where a paper reports only a point
estimate, `sd_c` is an explicit structural-error sensitivity assumption.
Where a paper reports “±” without an uncertainty definition available in the
summary source, the script provisionally treats that value as one standard
deviation and preserves that decision in `uncertainty_basis`.

Metrics:

- `MAAT`: mean annual air temperature (equivalent to MAT here)
- `CMMT`: cold-month mean temperature
- `WMMT`: warm-month mean temperature
- `MART`: mean annual range of temperature, normally WMMT minus CMMT

Key sources:

- Kiehl et al. (2018), *Philosophical Transactions A*,
  <https://doi.org/10.1098/rsta.2017.0085>. Table 3 supplies Bighorn Basin
  observational and PETM model values. The orbital-maximum run is retained as
  an intentionally extreme sensitivity case and excluded from the central
  synthesis.
- Hyland et al. (2018), *Climate of the Past*,
  <https://doi.org/10.5194/cp-14-1391-2018>. Supplies Green River Basin
  carbonate and downscaled GCM context for a dry continental-interior analog.
- Burgener et al. (2026), *Geological Society of America Bulletin*,
  <https://doi.org/10.1130/B38825.1>. Supports keeping dry-soil carbonate and
  wetter leaf assemblage MART estimates separate because they sample different
  environments.
- Zeebe et al. (2017), *Paleoceanography*,
  <https://doi.org/10.1002/2016PA003054>. Supports the orbital interpretation:
  eccentricity maxima amplify precession-scale insolation seasonality, whereas
  exact orbital phase becomes less secure in the older part of the record.
- Sewall and Sloan (2006) is retained as a historical Bighorn Basin GCM
  benchmark and excluded from the primary synthesis because the simulation is
  substantially cooler than the proxy and later-model evidence.

Downstream users should ordinarily plot
`BHB_seasonal_temperature_integrated_draws.csv.gz` or its compact summary, but
should show the source-specific distributions in sensitivity figures.
