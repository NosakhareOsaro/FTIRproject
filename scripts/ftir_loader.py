"""
ftir_loader.py

Load any of the six FTIR .dat files. Returns a (meta, spectra) tuple:
  meta    — normalised metadata, columns vary by file
  spectra — absorbance matrix, columns are wavenumber integers (3900 → 456 cm⁻¹)

Both DataFrames share the same 0-based integer index (row = one fly spectrum).
"""

from pathlib import Path
import warnings
import pandas as pd

# ── Known metadata column names ───────────────────────────────────────────────
_META_COLS = {"Genot.", "Sex", "Age", "Diet", "StoTime"}

# ── Normalisation maps ────────────────────────────────────────────────────────

# Diet: two coding schemes across experiments, mapped to descriptive strings
_DIET_MAP = {
    "SY":  "standard_yeast",  # standard yeast diet (DGRP, SexGeno, Age, Diet1 control)
    "LI":  "high_fat",        # lipid-supplemented (Diet1 treatment)
    "D05": "yeast_5pct",      # 5 % yeast dietary restriction (Diet2, Diet-Age)
    "D10": "yeast_10pct",     # 10 % yeast dietary restriction
    "D20": "yeast_20pct",     # 20 % yeast (ad-libitum equivalent)
}

# Age: leading-zero numeric codes → compact descriptive strings
_AGE_MAP = {
    "09D": "9d",       # 9-day-old
    "04W": "4w",       # 4-week-old
    "06W": "6w",       # 6-week-old
    "03W": "3w",       # 3-week-old (Diet-Age only)
    "UNK": "unknown",  # age not recorded (DGRP, SexGeno, Diet1)
}

# Sex: F/M are unambiguous.
# 'S' appears only in Diet1FTIR.dat.  The pre-print states diet experiments used
# female flies, which makes 'S' unlikely to mean "unknown sex" — but the exact
# meaning has not been confirmed.
# TODO: confirm what 'S' represents with Dr Rita Ibrahim before using Diet1 sex labels.
_SEX_MAP = {
    "F": "F",
    "M": "M",
    "S": "S",  # unresolved — left as-is pending confirmation with data owner
}

# Genot. is left unnormalised because values are experiment-specific:
#   DGRP lines  : DGRP100, DGRP405, …  (108 inbred lines)
#   Mito-nuclear: AA1/AA2/AA3, AB1/…, BA1/…, BB1/…  (SexGeno experiment)
#   Lab stocks  : WDH (likely wDah/Dahomey w1118 background — confirm with data owner)
#                 DPM (identity unknown — confirm with data owner)
# TODO: confirm WDH and DPM stock identities with Dr Rita Ibrahim.

# Verified from all six files: 3900, 3898, …, 456 step -2
_EXPECTED_WNS = tuple(range(3900, 454, -2))  # 1723 wavenumbers


def load_ftir(path):
    """
    Load one FTIR .dat file.

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    meta : pd.DataFrame
        Normalised metadata. Columns: subset of [Genot., Sex, Age, Diet, StoTime].
        Sex, Age, Diet are mapped through the module-level normalisation maps.
        Genot. and StoTime are returned as-is.
    spectra : pd.DataFrame
        Absorbance matrix, shape (n_spectra, 1723).
        Columns are wavenumber integers (cm⁻¹), descending 3900 → 456.

    Raises
    ------
    ValueError
        If the wavenumber axis does not match the expected 3900 → 456 range.
    """
    path = Path(path)
    df = pd.read_csv(path, sep="\t")

    meta_cols = [c for c in df.columns if c in _META_COLS]
    wn_cols   = [c for c in df.columns if c not in _META_COLS]

    try:
        wns_int = [int(c) for c in wn_cols]
    except ValueError as exc:
        raise ValueError(
            f"{path.name}: unexpected non-numeric column in spectrum region: {exc}"
        ) from exc

    if tuple(wns_int) != _EXPECTED_WNS:
        raise ValueError(
            f"{path.name}: wavenumber axis differs from expected.\n"
            f"  Got   : {wns_int[0]} → {wns_int[-1]}  ({len(wns_int)} columns)\n"
            f"  Expect: {_EXPECTED_WNS[0]} → {_EXPECTED_WNS[-1]}  ({len(_EXPECTED_WNS)} columns)"
        )

    meta    = df[meta_cols].copy().reset_index(drop=True)
    spectra = df[wn_cols].astype(float).copy().reset_index(drop=True)
    spectra.columns = wns_int

    # ── Normalise metadata ────────────────────────────────────────────────────
    for col, mapping in [("Diet", _DIET_MAP), ("Age", _AGE_MAP), ("Sex", _SEX_MAP)]:
        if col not in meta.columns:
            continue
        unmapped = set(meta[col].dropna()) - set(mapping)
        if unmapped:
            warnings.warn(
                f"{path.name}: unknown {col} values (kept as-is): {sorted(unmapped)}"
            )
        meta[col] = meta[col].map(mapping).fillna(meta[col])

    return meta, spectra


def check_wavenumber_consistency(paths):
    """
    Load only the header of each file and verify all share the same WN axis.

    Parameters
    ----------
    paths : dict[str, str | Path]

    Returns
    -------
    bool : True if all files match.
    """
    axes = {}
    for label, path in paths.items():
        header = pd.read_csv(Path(path), sep="\t", nrows=0)
        wn_cols = [c for c in header.columns if c not in _META_COLS]
        axes[label] = tuple(int(c) for c in wn_cols)

    ref_label, ref_axis = next(iter(axes.items()))
    all_match = True
    for label, axis in axes.items():
        if axis != ref_axis:
            print(f"  MISMATCH: {label} vs {ref_label}")
            print(f"    {label:<10}: {axis[0]} → {axis[-1]}  ({len(axis)} cols)")
            print(f"    {ref_label:<10}: {ref_axis[0]} → {ref_axis[-1]}  ({len(ref_axis)} cols)")
            all_match = False

    if all_match:
        print(
            f"  All {len(axes)} files share the same axis: "
            f"{ref_axis[0]} → {ref_axis[-1]}  ({len(ref_axis)} columns)"
        )
    return all_match


# ── Self-check ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    data_dir = Path(__file__).parent.parent / "FTIR-data"

    FILES = {
        "DGRP":     data_dir / "DGRPFTIR.dat",
        "SexGeno":  data_dir / "SexGenoFTIR.dat",
        "Diet2":    data_dir / "Diet2FTIR.dat",
        "Age":      data_dir / "AgeFTIR.dat",
        "Diet-Age": data_dir / "Diet-AgeFTIR.dat",
        "Diet1":    data_dir / "Diet1FTIR.dat",
    }

    print("=" * 60)
    print("Wavenumber axis consistency")
    print("=" * 60)
    check_wavenumber_consistency(FILES)

    for name, path in FILES.items():
        meta, spectra = load_ftir(path)
        print()
        print("=" * 60)
        print(f"FILE : {name}  ({path.name})")
        print(f"Shape: {spectra.shape[0]} spectra × {spectra.shape[1]} wavenumbers")
        print(f"Meta : {list(meta.columns)}")
        for col in meta.columns:
            uvals = sorted(meta[col].dropna().unique().tolist(), key=str)
            print(f"  {col:10s}: {uvals}")
        print(
            f"Absorbance range: "
            f"min={spectra.values.min():.4f}  max={spectra.values.max():.4f}"
        )
