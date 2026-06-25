# phenotype-data

External phenotype data downloaded from DGRPool (dgrpool.epfl.ch) for use as regression targets alongside the FTIR spectra. Add a row here each time a new file is downloaded.

| File | Phenotype | Study | Reference | URL | Downloaded | Notes |
|------|-----------|-------|-----------|-----|------------|-------|
| `S24_StarvationRes_summary_mean.tsv` | Starvation resistance (hours, mean per line) | Study 24 | Morgante et al. 2015 | dgrpool.epfl.ch/phenotypes/2798 | 2026-06-16 | 197 DGRP lines, both sexes (F and M rows); 93 of our 108 FTIR lines present; all 93 have valid female means |
| `S18_LifeFecundity_mean.tsv` | Lifetime fecundity (eggs/female, mean per line) | Study 18 | DGRPool | dgrpool.epfl.ch/studies/18 | 2026-06-25 | 189 DGRP lines, all female; 96 of our 108 FTIR lines present; no NaN values. Elastic net LOO-CV: R²=−0.109 — no spectral signal detected |
