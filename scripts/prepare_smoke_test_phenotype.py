"""
prepare_smoke_test_phenotype.py

Reformats Emmeans.csv (our own starvation-resistance EMMeans, produced by
scripts/run_survival_analysis.R) into a mock DGRPool-format phenotype TSV,
so it can be run through scripts/run_dgrpool_phenotype.py as a smoke test:
since the model is being asked to predict the exact same values it was
already shown (via a different script) to predict with R²=0.673, the
general-purpose runner should reproduce that number almost exactly. If it
doesn't, the generic script has a bug, not the phenotype.

This is the one committed phenotype-data/ file that isn't an external
download or a reformat of one - it's derived entirely from this project's
own output. It went undocumented for a while (built by hand once, early
on, before the general-purpose runner even had other phenotypes to be
tested against) - this script exists so it's reproducible from Emmeans.csv
like everything else.

Output: phenotype-data/S00_EMMeans_starvation.tsv
  Columns: DGRP, sex (hardcoded "F" - EMMeans.csv is female-only), value
  108 rows, matching Emmeans.csv's 108 DGRP lines exactly.
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).parent.parent
EMMEANS_PATH = REPO / "Emmeans.csv"
OUT_PATH = REPO / "phenotype-data" / "S00_EMMeans_starvation.tsv"


def main():
    emmeans = pd.read_csv(EMMEANS_PATH)
    out = pd.DataFrame({
        "DGRP": emmeans["DGRP"],
        "sex": "F",
        "value": emmeans["emmean"],
    })
    out.to_csv(OUT_PATH, sep="\t", index=False)
    print(f"Wrote {len(out)} rows to {OUT_PATH.relative_to(REPO)}")


if __name__ == "__main__":
    main()
