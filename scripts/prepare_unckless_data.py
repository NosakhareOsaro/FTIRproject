"""
prepare_unckless_data.py

Convert the Unckless et al. 2015 (G3, DOI 10.1534/g3.114.016477) Table S2
nutritional-index supplementary file into per-measure, per-diet DGRPool-
format TSVs: one pooled-diet and two single-diet (high/low glucose) files
per metabolic measure.

Source: phenotype-data/raw/016477_tables2.xlsx, sheet "TableS2.csv".
Male-only data: assayed in pools of 10 adult males per DGRP line
(confirmed from the paper's Materials and Methods).

For each of six measures (glucose, glycerol, glycogen, triglyceride,
protein, meanweight(ug)), three source columns are extracted: "_pooled",
"_high_glucose", and "_low_glucose". Missing values (marked with a
literal "." string in the source file) are converted to NaN and dropped
independently for each column; the set of missing lines differs by diet
condition (not by measure) but is dropped identically across all six
measures within a given diet. DGRP Line integers are reformatted to
DGRP<n> to match our spectral line ID convention (e.g. DGRP21, not the
zero-padded DGRP_021 used by some other DGRPool downloads).

Writes to phenotype-data/, for each of the six measures:
  Unckless_<Measure>_pooled.tsv
  Unckless_<Measure>_highglucose.tsv
  Unckless_<Measure>_lowglucose.tsv

Each TSV has columns: DGRP, sex (hardcoded to M), value.
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).parent.parent
RAW_PATH = REPO / "phenotype-data" / "raw" / "016477_tables2.xlsx"
SHEET_NAME = "TableS2.csv"
OUT_DIR = REPO / "phenotype-data"

# Source column prefix -> output file label
MEASURES = {
    "glucose": "Glucose",
    "glycerol": "Glycerol",
    "glycogen": "Glycogen",
    "triglyceride": "Triglyceride",
    "protein": "Protein",
    "meanweight(ug)": "MeanWeight",
}

# Source column suffix -> output file suffix
DIETS = {
    "pooled": "pooled",
    "high_glucose": "highglucose",
    "low_glucose": "lowglucose",
}


def main():
    df = pd.read_excel(RAW_PATH, sheet_name=SHEET_NAME)

    n_dup = df["DGRP Line"].duplicated().sum()
    if n_dup:
        dup_lines = sorted(int(x) for x in df.loc[df["DGRP Line"].duplicated(keep=False), "DGRP Line"].unique())
        print(f"Note: {n_dup} duplicate DGRP Line row(s) found ({dup_lines}); keeping first occurrence.")
        df = df.drop_duplicates(subset="DGRP Line", keep="first").reset_index(drop=True)

    df["DGRP"] = "DGRP" + df["DGRP Line"].astype(int).astype(str)

    summary = []
    for col_prefix, label in MEASURES.items():
        for diet_suffix, out_suffix in DIETS.items():
            col = f"{col_prefix}_{diet_suffix}"
            values = pd.to_numeric(df[col], errors="coerce")

            out = pd.DataFrame({"DGRP": df["DGRP"], "sex": "M", "value": values})
            out = out.dropna(subset=["value"]).reset_index(drop=True)

            out_path = OUT_DIR / f"Unckless_{label}_{out_suffix}.tsv"
            out.to_csv(out_path, sep="\t", index=False)

            summary.append({
                "file": out_path.name,
                "n": len(out),
                "min": out["value"].min(),
                "max": out["value"].max(),
                "mean": out["value"].mean(),
            })

    print()
    print(f"{'File':<38} {'n':>4} {'min':>12} {'max':>12} {'mean':>12}")
    print("-" * 84)
    for row in summary:
        print(f"{row['file']:<38} {row['n']:>4} {row['min']:>12.4f} {row['max']:>12.4f} {row['mean']:>12.4f}")


if __name__ == "__main__":
    main()
