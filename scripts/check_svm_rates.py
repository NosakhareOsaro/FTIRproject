"""
Read the per-fold CV values saved by run_dgrp_baseline.py and report
the SVM per-class reclassification rates (TPR = Resistant, TNR = Sensitive).
Target from paper: ~88% resistant, ~84% sensitive.
"""

import pandas as pd
import numpy as np
from pathlib import Path

results_dir = Path(__file__).parent.parent / "results" / "DGRP"
df = pd.read_csv(results_dir / "DGRP_XGBoost_CV_values.csv")

# The Classifier column holds tuple reprs because the dict key is a tuple (name, estimator)
# Find the SVM rows by checking for the string 'Support Vector'
svm_mask = df["Classifier"].astype(str).str.contains("Support Vector")
svm = df[svm_mask]

print(f"SVM folds found: {len(svm)}")
print()

mean_tpr = svm["TPR"].mean()
std_tpr  = svm["TPR"].std()
mean_tnr = svm["TNR"].mean()
std_tnr  = svm["TNR"].std()

print(f"SVM Resistant reclassification (TPR): {mean_tpr:.1%} ± {std_tpr:.1%}")
print(f"SVM Sensitive  reclassification (TNR): {mean_tnr:.1%} ± {std_tnr:.1%}")
print()
print("Paper targets: ~88% resistant, ~84% sensitive")

print()
print("--- Per-fold detail ---")
print(svm[["Fold", "TPR", "TNR", "FPR", "FNR"]].to_string(index=False))
