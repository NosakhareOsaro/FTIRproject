"""
run_dgrpool_phenotype.py

General-purpose runner: takes any DGRPool-format phenotype TSV
(columns DGRP, sex, value) and runs the elastic net LOO-CV pipeline
against DGRP line-mean FTIR spectra. Same model and hyperparameters
as run_regularised_regression.py / run_fecundity_enet.py.

Usage:
  .venv/bin/python scripts/run_dgrpool_phenotype.py <phenotype.tsv> \\
      [--sex F] --study "Morgante 2015" --phenotype "Starvation resistance"

Validation hierarchy (see phenotype-data/README.md):
  S00_EMMeans_starvation.tsv  : our own EMMeans reformatted as a mock
                                 DGRPool TSV. This is the smoke test:
                                 it should reproduce CV R² ≈ 0.673,
                                 confirming the script is correct.
  S24_StarvationRes_summary_mean.tsv (Morgante 2015) and all other
  DGRPool files: genuine external phenotypes. Whatever R² comes out
  is a real result, not expected to match 0.673.

Each run appends one row to results/DGRP/dgrpool_phenotype_summary.csv
(created with a header if absent) so cross-phenotype results
accumulate automatically.
"""

import argparse
import sys
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

# Coordinate descent frequently fails the strict convergence check at p=1,723
# features; this is expected (see run_regularised_regression.py notes) and
# does not affect the selected alpha or final coefficients. Silenced to keep
# LOO-CV output readable across ~100 folds x 30 alphas x 5 l1_ratios.
warnings.filterwarnings("ignore", category=ConvergenceWarning)

REPO = Path(__file__).parent.parent
sys.path.insert(0, str(Path(__file__).parent))
from ftir_loader import load_ftir

L1_RATIOS = [0.5, 0.7, 0.9, 0.95, 1.0]
N_ALPHAS = 30

SUMMARY_PATH = REPO / "results" / "DGRP" / "dgrpool_phenotype_summary.csv"
SUMMARY_COLS = ["phenotype", "study", "n_lines", "cv_r2", "rmse", "spearman_rho", "date_run"]


def normalise_line_id(dgrp_id):
    """DGRP_021 -> DGRP21 ; DGRP229 -> DGRP229 (strip underscore, drop leading zeros)."""
    return "DGRP" + str(int(dgrp_id.replace("DGRP", "").replace("_", "")))


def main():
    parser = argparse.ArgumentParser(
        description="Run elastic net LOO-CV of DGRP line-mean FTIR spectra "
        "against any DGRPool-format phenotype TSV (columns: DGRP, sex, value)."
    )
    parser.add_argument("tsv_path", type=Path, help="Path to the DGRPool phenotype TSV")
    parser.add_argument("--sex", default="F", help="Sex to filter to (default: F)")
    parser.add_argument("--study", required=True, help="Study label, e.g. 'Morgante 2015'")
    parser.add_argument("--phenotype", required=True, help="Phenotype label, e.g. 'Starvation resistance'")
    args = parser.parse_args()

    # ── Load and normalise phenotype ──────────────────────────────────────────
    pheno = pd.read_csv(args.tsv_path, sep="\t")
    pheno = pheno[pheno["sex"] == args.sex].copy()
    pheno["line"] = pheno["DGRP"].apply(normalise_line_id)
    pheno_map = pheno.set_index("line")["value"].to_dict()

    # ── Load spectral data, filter to sex, average per line ──────────────────
    meta, spectra = load_ftir(REPO / "FTIR-data" / "DGRPFTIR.dat")
    sex_mask = meta["Sex"] == args.sex
    meta = meta[sex_mask].reset_index(drop=True)
    spectra = spectra[sex_mask].reset_index(drop=True)

    wavenumbers = np.array(spectra.columns.tolist())

    spec_df = spectra.copy()
    spec_df["Genot."] = meta["Genot."].values
    X_line_df = spec_df.groupby("Genot.").mean()

    spectral_lines = set(X_line_df.index)
    pheno_lines = set(pheno_map)
    overlap = sorted(spectral_lines & pheno_lines)
    missing_from_pheno = sorted(spectral_lines - pheno_lines)
    n_lines = len(overlap)

    # ── Header block ───────────────────────────────────────────────────────────
    print("=" * 68)
    print(f"Phenotype : {args.phenotype}")
    print(f"Study     : {args.study}")
    print(f"File      : {args.tsv_path}")
    print(f"Sex       : {args.sex}")
    print("-" * 68)
    print(f"Spectral lines ({args.sex})         : {len(spectral_lines)}")
    print(f"Phenotype lines ({args.sex})        : {len(pheno_lines)}")
    print(f"Overlap                       : {n_lines}")
    if missing_from_pheno:
        print(f"Spectral lines with no phenotype value ({len(missing_from_pheno)}): {missing_from_pheno}")
    print("=" * 68)
    print()

    if n_lines < 10:
        print(f"ERROR: only {n_lines} overlapping lines — too few for LOO-CV. Aborting.")
        sys.exit(1)

    X_line = X_line_df.loc[overlap, wavenumbers.tolist()].values
    y_line = np.array([pheno_map[ln] for ln in overlap])

    print(f"Analysis matrix : {n_lines} lines x {X_line.shape[1]} wavenumbers")
    print(f"Phenotype range : {y_line.min():.3f} - {y_line.max():.3f}  "
          f"(mean {y_line.mean():.3f}, SD {y_line.std():.3f})")
    print()

    # ── LOO-CV elastic net ────────────────────────────────────────────────────
    print(f"Running elastic net LOO-CV (l1_ratio={L1_RATIOS}, alphas={N_ALPHAS}, cv=3) ...")

    loo = LeaveOneOut()
    y_pred = np.zeros(n_lines)

    for train_idx, test_idx in loo.split(X_line):
        sc = StandardScaler()
        X_tr = sc.fit_transform(X_line[train_idx])
        X_te = sc.transform(X_line[test_idx])
        mdl = ElasticNetCV(
            cv=3, l1_ratio=L1_RATIOS,
            alphas=N_ALPHAS, max_iter=5000, tol=0.01,
        )
        mdl.fit(X_tr, y_line[train_idx])
        y_pred[test_idx] = mdl.predict(X_te)

    r2 = r2_score(y_line, y_pred)
    rmse = np.sqrt(mean_squared_error(y_line, y_pred))
    rho = spearmanr(y_line, y_pred).statistic
    pred_sd = y_pred.std()
    true_sd = y_line.std()

    print()
    print("-" * 68)
    print("RESULTS")
    print("-" * 68)
    print(f"  n lines     : {n_lines}")
    print(f"  CV R²       : {r2:+.3f}")
    print(f"  RMSE        : {rmse:.4f}")
    print(f"  Spearman ρ  : {rho:+.3f}")

    collapsed = pred_sd < 0.2 * true_sd
    if collapsed:
        print()
        print("  NOTE: predictions collapse toward the training mean "
              f"(pred SD={pred_sd:.3f} vs true SD={true_sd:.3f}).")
        print("        No reliable spectral signal detected for this phenotype.")
        print("        Spearman ρ may be a LOO mean-shift artefact rather than a")
        print("        genuine rank correlation — treat it with caution.")
    print("-" * 68)
    print()

    # ── Append to summary CSV ─────────────────────────────────────────────────
    row = pd.DataFrame([{
        "phenotype": args.phenotype,
        "study": args.study,
        "n_lines": n_lines,
        "cv_r2": r2,
        "rmse": rmse,
        "spearman_rho": rho,
        "date_run": date.today().isoformat(),
    }])

    SUMMARY_PATH.parent.mkdir(parents=True, exist_ok=True)
    write_header = not SUMMARY_PATH.exists()
    row.to_csv(SUMMARY_PATH, mode="a", header=write_header, index=False, columns=SUMMARY_COLS)
    print(f"Appended result to {SUMMARY_PATH.relative_to(REPO)}")


if __name__ == "__main__":
    main()
