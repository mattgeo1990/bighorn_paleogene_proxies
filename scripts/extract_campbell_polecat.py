"""Extract Campbell et al. (2024) Polecat Bench soil-water isotope profiles.

The legacy project CSV is best reproduced by the equal-weight mean of the four
model cells at 48.316–50.211 N, 272.5–275.0 E and CLM soil layer 5
(center depth 0.366 m). Delta values are calculated from isotope-resolved
volumetric soil water: 1000 * (H2OSOI_H218O / H2OSOI - 1).
"""
from pathlib import Path
import argparse
import csv
import numpy as np
from netCDF4 import Dataset

parser = argparse.ArgumentParser()
parser.add_argument(
    "--nc-dir", type=Path, required=True,
    help="Directory containing the eight Campbell et al. CLM NetCDF files."
)
parser.add_argument(
    "--out-dir", type=Path, default=Path("data/raw"),
    help="Output directory (default: data/raw)."
)
args = parser.parse_args()
NC_DIR = args.nc_dir.expanduser().resolve()
OUT_DIR = args.out_dir.expanduser().resolve()
OUT_DIR.mkdir(exist_ok=True)

FILES = {
    "3x_orbmod": "3x_OrbMod.nc",
    "3x_orbmin": "3x_OrbMin.nc",
    "3x_orbmaxn": "3x_OrbMaxN.nc",
    "3x_orbmaxs": "3x_OrbMaxS.nc",
    "6x_orbmod": "OrbMod.nc",
    "6x_orbmin": "OrbMin.nc",
    "6x_orbmaxn": "OrbMaxN.nc",
    "6x_orbmaxs": "OrbMaxS.nc",
}

monthly_rows = []
depth_rows = []

for experiment, filename in FILES.items():
    with Dataset(NC_DIR / filename) as ds:
        lat = np.asarray(ds["lat"][:])
        lon = np.asarray(ds["lon"][:])
        depths = np.asarray(ds["levgrnd"][:])
        lat_idx = np.array([np.abs(lat - x).argmin() for x in (48.316, 50.211)])
        lon_idx = np.array([np.abs(lon - x).argmin() for x in (272.5, 275.0)])

        h2o = np.asarray(ds["H2OSOI"][:, :, lat_idx, :][:, :, :, lon_idx])
        h218o = np.asarray(ds["H2OSOI_H218O"][:, :, lat_idx, :][:, :, :, lon_idx])
        delta = (h218o / h2o - 1.0) * 1000.0
        box_mean = np.nanmean(delta, axis=(2, 3))

        # Primary extraction: layer center at 0.366 m, identified by matching
        # the legacy project summary across all eight experiments.
        z = int(np.abs(depths - 0.3661).argmin())
        for month, value in enumerate(box_mean[:, z], start=1):
            monthly_rows.append({
                "experiment": experiment,
                "month": month,
                "soil_layer_center_m": float(depths[z]),
                "d18O_soil_water_permil_vsmow": float(value),
            })

        # Full monthly depth profile, retained so other rooting/depth
        # integrations can be evaluated without downloading the NetCDF again.
        for month in range(12):
            for zi, depth in enumerate(depths):
                depth_rows.append({
                    "experiment": experiment,
                    "month": month + 1,
                    "soil_layer": zi + 1,
                    "soil_layer_center_m": float(depth),
                    "d18O_soil_water_permil_vsmow": float(box_mean[month, zi]),
                })

def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

write_csv(OUT_DIR / "campbell2024_polecat_monthly_soilwater.csv", monthly_rows)
write_csv(OUT_DIR / "campbell2024_polecat_monthly_soilwater_all_depths.csv", depth_rows)

print(f"Wrote {len(monthly_rows)} primary monthly values")
print(f"Wrote {len(depth_rows)} monthly depth-profile values")
