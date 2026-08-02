# FTIR project

MSc dissertation. Can FTIR spectra of *Drosophila melanogaster* (the DGRP
panel) predict continuous phenotype, not just classify it (classification
is what the lab pre-print this builds on did: Ibrahim et al., bioRxiv 2026,
doi:10.64898/2026.03.22.713522).

Short version of where this stands: FTIR predicts starvation resistance
well in our own data (elastic net, LOO-CV R²=0.673). It does not predict
much else. Same phenotype measured by a different lab, no. Lifespan,
fecundity, chill coma recovery, cuticle hydrocarbons, six Unckless
metabolic measures across three diet conditions, no. Roughly thirty
phenotype/diet/temperature tests run in total: two show real signal (both
starvation resistance), two are weak uncorrected candidates, the rest are
nulls. As a final, unrelated check, lipidomics does predict BMI in the
EATRIS-Plus human cohort (R²=0.28). Modest for that field, but real, and
the first positive result in the whole project that isn't starvation
resistance.

Where to look for more:

- `PROJECT_NOTES.md`. The working log. Design decisions and every result,
  in the order things actually happened.
- `REPRODUCE.md`. Every command to reproduce this from a fresh clone,
  raw-data sources included.
- `phenotype-data/README.md`. Every phenotype file used as a target, where
  it came from, what it showed.
- `notebooks/`. Write-ups of the bigger side quests (Unckless, the
  Durham/Huang lifespan and fecundity work, EATRIS-Plus lipidomics).
