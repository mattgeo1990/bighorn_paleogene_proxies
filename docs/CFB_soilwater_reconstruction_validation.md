# CFB soil-water reconstruction validation

The production reconstruction is implemented in [`scripts/05_reconstruct_CFB_soilwater.R`](../scripts/05_reconstruct_CFB_soilwater.R). The independent validation is implemented in [`scripts/05_validate_CFB_soilwater_reconstruction.R`](../scripts/05_validate_CFB_soilwater_reconstruction.R).

## What is being validated

The validation checks four things:

1. the carbonate–water fractionation equation;
2. the isotope reference-scale and logarithmic notation conventions;
3. the temperature values and temperature uncertainties supplied to the Monte Carlo calculation; and
4. agreement between independent deterministic calculations and the existing Monte Carlo medians.

This is a numerical/procedural validation. It does not determine whether the carbonate was precipitated at equilibrium or whether it is diagenetically altered.

## δ18O soil-water equation

The production script defines:

```text
α18(T) = exp[(18.03 × (1000/TK) − 32.42)/1000]
TK = T°C + 273.15
```

and reconstructs water as:

```text
δ18Owater = [(1000 + δ18Ocarbonate) / α18(T)] − 1000
```

Carbonate δ18O is in VSMOW-compatible conventional delta notation before this calculation. The independent validator reproduces this equation separately and checks the deterministic value against the Monte Carlo median.

## Δ′17O soil-water equation

The production script uses:

```text
Δ′17Owater = Δ′17Ocarbonate
             − 1,000,000 × (θcarbonate − λreference) × ln(α18)
```

with:

```text
θcarbonate = 0.5250
λreference = 0.528
```

Because Δ′17O is correlated with δ′18O, the production Monte Carlo samples carbonate Δ′17O directly rather than independently resampling δ′17O and δ′18O. This avoids generating artificial uncertainty from destroying the covariance embedded in the calculated Δ′17O value.

## Reference scales and units

The validation confirms that:

- carbonate δ18O is treated as a conventional delta value for the fractionation equation;
- Δ′17O is treated as per meg;
- temperature is converted from °C to K only inside the fractionation equation;
- λ and θ are dimensionless slopes;
- `ln(α18)` is dimensionless; and
- the factor 1,000,000 converts the logarithmic fractionation term to per meg.

## Temperature inputs

The production script uses `T_recon_C` and `T_recon_se_C`, which are selected upstream from the project temperature model and measured T47 data according to the existing temperature-input logic. The validator checks that every tested temperature is above absolute zero and that its uncertainty is non-negative.

Temperature uncertainty is sampled as a normal distribution in the Monte Carlo calculation. The validation does not claim that this captures age-model uncertainty, model structural uncertainty, or preservation uncertainty; those remain separate limitations.

## Spot checks

The validator independently recalculates five representative horizons:

- PK95-SC-246
- PK95-SC-279
- PK95-SC-187
- PK95-SC-165
- PK95-SC-185

It writes the detailed results to:

[`data/processed/CFB_soilwater_reconstruction_spot_checks.csv`](../data/processed/CFB_soilwater_reconstruction_spot_checks.csv)

The overall validation results are written to:

[`data/processed/CFB_soilwater_reconstruction_validation_summary.csv`](../data/processed/CFB_soilwater_reconstruction_validation_summary.csv)

The script requires the independent deterministic and Monte Carlo medians to agree within 0.05‰ for δ18Owater and 1 per meg for Δ′17Owater. These are numerical implementation checks, not analytical uncertainty limits.

## Interpretation

Passing this validation means that the stored outputs are numerically consistent with the equations implemented in the production script and that the temperature inputs are physically valid. It does not mean that:

- the carbonate–water system was at equilibrium;
- the reconstructed water was unmodified meteoric water;
- the temperature model is correct;
- all uncertainty sources are represented; or
- diagenetic screening is complete.

Those questions require the separate temperature, preservation, and diagenetic workflows.
