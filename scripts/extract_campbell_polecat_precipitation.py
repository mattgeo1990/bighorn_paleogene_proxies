"""Extract Campbell et al. (2024) precipitation isotopes at Polecat Bench.

Reads the official CAM climatology files archived at Zenodo record 7971738.
The extraction uses the same four model cells as the existing soil-water
extraction: 48.316--50.211 N and 272.5--275.0 E.

Monthly regional isotope values are precipitation-flux-weighted across the
four cells. Seasonal and annual isotope values are then precipitation-weighted
across their constituent months. This preserves isotope mass balance; a plain
arithmetic mean of monthly delta values does not.
"""
from pathlib import Path
import argparse
import csv
import calendar

import fsspec
import h5py
import numpy as np


RECORD = "https://zenodo.org/api/records/7971738/files"
FILES = {
    "3x_orbmod": "b.e12.B1850C5CN.f19_g16.iPETM03x.03.cam.h0.2101-2200.climo.nc",
    "3x_orbmin": "b.e12.B1850C5CN.f19_g16.iPETM03x.OrbMin.01.cam.h0.2601-2700.climo.nc",
    "3x_orbmaxn": "b.e12.B1850C5CN.f19_g16.iPETM03x.OrbMaxN.01.cam.h0.2601-2700.climo.nc",
    "3x_orbmaxs": "b.e12.B1850C5CN.f19_g16.iPETM03x.OrbMaxS.01.cam.h0.2601-2700.climo.nc",
    "6x_orbmod": "b.e12.B1850C5CN.f19_g16.iPETM06x.09.cam.h0.2101-2200.climo.nc",
    "6x_orbmin": "b.e12.B1850C5CN.f19_g16.iPETM06x.OrbMin.01.cam.h0.2601-2700.climo.nc",
    "6x_orbmaxn": "b.e12.B1850C5CN.f19_g16.iPETM06x.OrbMaxN.01.cam.h0.2601-2700.climo.nc",
    "6x_orbmaxs": "b.e12.B1850C5CN.f19_g16.iPETM06x.OrbMaxS.01.cam.h0.2601-2700.climo.nc",
}

H2O_COMPONENTS = (
    "PRECRC_H2Or", "PRECRL_H2OR", "PRECSC_H2Os", "PRECSL_H2OS",
)
H218O_COMPONENTS = (
    "PRECRC_H218Or", "PRECRL_H218OR", "PRECSC_H218Os", "PRECSL_H218OS",
)
TARGET_LATS = (48.316, 50.211)
TARGET_LONS = (272.5, 275.0)
SECONDS_PER_DAY = 86400.0


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def weighted_delta(h218o, h2o):
    return 1000.0 * (np.sum(h218o) / np.sum(h2o) - 1.0)


def open_remote(filename):
    url = f"{RECORD}/{filename}/content"
    remote = fsspec.open(url, "rb", block_size=2**20).open()
    return remote, h5py.File(remote, "r")


parser = argparse.ArgumentParser()
parser.add_argument("--out-dir", type=Path, default=Path("data/raw"))
args = parser.parse_args()
args.out_dir.mkdir(parents=True, exist_ok=True)

monthly_rows = []
summary_rows = []

for experiment, filename in FILES.items():
    print(f"Reading {experiment}", flush=True)
    remote, ds = open_remote(filename)
    try:
        lat = np.asarray(ds["lat"][:])
        lon = np.asarray(ds["lon"][:])
        lat_idx = np.array([np.abs(lat - x).argmin() for x in TARGET_LATS])
        lon_idx = np.array([np.abs(lon - x).argmin() for x in TARGET_LONS])

        def box(variable):
            return np.asarray(ds[variable][:, lat_idx, :][:, :, lon_idx], dtype=float)

        h2o = sum(box(name) for name in H2O_COMPONENTS)
        h218o = sum(box(name) for name in H218O_COMPONENTS)
        precip = box("PRECC") + box("PRECL")

        # The isotope-specific H2O sum and standard precipitation field should
        # describe the same total flux. Fail visibly if archive conventions change.
        if not np.allclose(h2o, precip, rtol=2e-4, atol=1e-15, equal_nan=True):
            raise ValueError(f"H2O components do not reproduce PRECC + PRECL: {experiment}")

        for month in range(12):
            days = calendar.monthrange(2001, month + 1)[1]  # non-leap climatology
            monthly_rows.append({
                "experiment": experiment,
                "co2_multiple": int(experiment[0]),
                "month": month + 1,
                "month_name": calendar.month_abbr[month + 1],
                "season": "DJF" if month in (11, 0, 1) else
                          "MAM" if month in (2, 3, 4) else
                          "JJA" if month in (5, 6, 7) else "SON",
                "d18O_precip_permil_vsmow": weighted_delta(h218o[month], h2o[month]),
                "precip_mm_month": float(np.mean(precip[month]) * SECONDS_PER_DAY * days * 1000.0),
                "n_grid_cells": int(h2o[month].size),
            })

        periods = {
            "annual": np.arange(12),
            "DJF": np.array([11, 0, 1]),
            "MAM": np.array([2, 3, 4]),
            "JJA": np.array([5, 6, 7]),
            "SON": np.array([8, 9, 10]),
        }
        for period, months in periods.items():
            month_days = np.array([
                calendar.monthrange(2001, int(month) + 1)[1] for month in months
            ])
            day_weights = month_days[:, None, None]
            summary_rows.append({
                "experiment": experiment,
                "co2_multiple": int(experiment[0]),
                "period": period,
                "d18O_precip_permil_vsmow": weighted_delta(
                    h218o[months] * day_weights, h2o[months] * day_weights
                ),
                "precip_mm_period": float(np.sum(
                    np.mean(precip[months], axis=(1, 2))
                    * SECONDS_PER_DAY * month_days * 1000.0
                )),
                "months": ",".join(str(int(month) + 1) for month in months),
            })
    finally:
        ds.close()
        remote.close()

write_csv(
    args.out_dir / "campbell2024_polecat_monthly_precipitation.csv",
    monthly_rows,
)
write_csv(
    args.out_dir / "campbell2024_polecat_precipitation_summary.csv",
    summary_rows,
)
print(f"Wrote {len(monthly_rows)} monthly rows and {len(summary_rows)} summaries")
