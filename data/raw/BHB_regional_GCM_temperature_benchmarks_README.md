# Regional Eocene GCM temperature benchmarks

`BHB_regional_GCM_temperature_benchmarks.csv` contains climate-model values
used only as comparisons with the Bighorn Basin proxy models. These are
equilibrium Eocene climate experiments, not dated observations, so the R
pipeline does not interpolate them through geologic time.

## Numerical sources

- Sewall and Sloan (2006), Table 2: Bighorn Basin regional-model ranges for
  MAT, cold mean month (standardized here as CMMT), and mean annual temperature
  range (standardized as MART). WMMT is an explicitly flagged approximation,
  calculated as `MAT + MART / 2`, solely so the same four metrics can be shown.
- Hyland et al. (2018), Supplement Table S7,
  <https://doi.org/10.5194/cp-14-1391-2018>: RegCM3 monthly temperatures for
  the nearby Green River Basin under LoCO (560 ppm) and HiCO (2240 ppm)
  configurations. MAT is the mean of the 12 monthly values; CMMT and WMMT are
  their minimum and maximum; MART is `WMMT - CMMT`.

## Interpretive sources

Huber and Caballero (2011), <https://doi.org/10.5194/cp-7-603-2011>, showed
that sufficiently strong greenhouse forcing in CCSM3 can broadly reproduce
early Eocene terrestrial MAT and CMMT. Their experiments support comparing the
proxy reconstructions with a range of greenhouse states, but they are not
treated as a dated Bighorn time series.

Burgener et al. (2026), <https://doi.org/10.1130/B38825.1>, argue that
land-cover and proxy-formation environments can generate real local variation
in MART. Accordingly, disagreement among LMA MAT, carbonate formation
temperature, and model seasonality is not interpreted automatically as model
failure; the pipeline preserves the climate quantities and proxy targets as
separate products.

DeepMIP Phase 1 is another appropriate ensemble comparison, but its EECO
simulations are likewise CO2-indexed climate states rather than a GTS2020-dated
series. If a Bighorn paleolocation extraction is added later, it should enter
this table as additional benchmark rows, not as control points in either proxy
temperature model.
