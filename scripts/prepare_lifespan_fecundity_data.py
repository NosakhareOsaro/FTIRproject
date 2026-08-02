"""
prepare_lifespan_fecundity_data.py

Reformat six new phenotype source files into DGRPool-format TSVs
(columns DGRP, sex, value) for the lifespan/fecundity extension:

  Durham et al. 2014 (mated female lifespan and age-specific fecundity):
    S12_lsm_Lifespan_original.tsv        -> Durham_Lifespan_mated.tsv
    S12_lsm_Week1_Fecundity_original.tsv -> Durham_Week1_Fecundity.tsv
    S12_lsm_Week7_Fecundity_original.tsv -> Durham_Week7_Fecundity.tsv

  These three are already in clean DGRP/sex/value format (one row per
  line) and are copied through unchanged. Note: DGRP_042 in the Week7
  file has a negative value (-0.165). This is expected and correct:
  these are least-squares means corrected for body size and block,
  not raw egg counts, so small negative LS means are possible.

  Huang et al. 2020 (mated female lifespan at three temperatures):
    S40_Lifespan_18C_original.tsv -> Huang_Lifespan_18C_female.tsv
    S40_Lifespan_25C_original.tsv -> Huang_Lifespan_25C_female.tsv
    S40_Lifespan_28C_original.tsv -> Huang_Lifespan_28C_female.tsv

  These three are individual fly-level raw data (one row per fly,
  both sexes present, ~72 flies per sex/line/temperature). Filtered
  to sex == 'F', then averaged per DGRP line to produce one female
  line-mean lifespan value per line, matching the DGRPool format
  used elsewhere in this project.
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).parent.parent
RAW_DIR = REPO / "phenotype-data" / "raw"
OUT_DIR = REPO / "phenotype-data"

DURHAM_FILES = {
    "S12_lsm_Lifespan_original.tsv": "Durham_Lifespan_mated.tsv",
    "S12_lsm_Week1_Fecundity_original.tsv": "Durham_Week1_Fecundity.tsv",
    "S12_lsm_Week7_Fecundity_original.tsv": "Durham_Week7_Fecundity.tsv",
}

HUANG_FILES = {
    "S40_Lifespan_18C_original.tsv": "Huang_Lifespan_18C_female.tsv",
    "S40_Lifespan_25C_original.tsv": "Huang_Lifespan_25C_female.tsv",
    "S40_Lifespan_28C_original.tsv": "Huang_Lifespan_28C_female.tsv",
}


def main():
    summary = []

    for src_name, out_name in DURHAM_FILES.items():
        df = pd.read_csv(RAW_DIR / src_name, sep="\t")
        out_path = OUT_DIR / out_name
        df.to_csv(out_path, sep="\t", index=False)
        summary.append({
            "file": out_name,
            "n": len(df),
            "min": df["value"].min(),
            "max": df["value"].max(),
            "mean": df["value"].mean(),
        })

    for src_name, out_name in HUANG_FILES.items():
        df = pd.read_csv(RAW_DIR / src_name, sep="\t")
        female = df[df["sex"] == "F"]
        means = female.groupby("DGRP", as_index=False)["value"].mean()
        means["sex"] = "F"
        means = means[["DGRP", "sex", "value"]]

        out_path = OUT_DIR / out_name
        means.to_csv(out_path, sep="\t", index=False)
        summary.append({
            "file": out_name,
            "n": len(means),
            "min": means["value"].min(),
            "max": means["value"].max(),
            "mean": means["value"].mean(),
        })

    print(f"{'File':<34} {'n':>4} {'min':>12} {'max':>12} {'mean':>12}")
    print("-" * 80)
    for row in summary:
        print(f"{row['file']:<34} {row['n']:>4} {row['min']:>12.4f} {row['max']:>12.4f} {row['mean']:>12.4f}")


if __name__ == "__main__":
    main()
