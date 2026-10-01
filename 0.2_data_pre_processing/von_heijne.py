"""
von Heijne method for signal peptide detection.

- cleavage-site context: [-13, +2]  (15 positions)
- pseudocounts: 1
- background: SwissProt amino-acid composition
- threshold: cross-validated on a validation fold via the PR curve
- 5-fold cross-validation, metrics averaged over the folds
"""

import re
import math
import numpy as np
import pandas as pd
from itertools import product

# =====================================================
# 0. SETTINGS
# =====================================================

CV_FILE = "training_cv.tsv"                 # accession, label, fold
UNIPROT_FILES = ["uniprot_results_pos.tsv", "uniprot_results_neg.tsv"]

UPSTREAM = 13        # positions before the cleavage site
DOWNSTREAM = 2       # positions after
WINDOW = UPSTREAM + DOWNSTREAM
PSEUDO = 1
SCAN = 90            # N-terminal region scanned for the best window

AA = "ACDEFGHIKLMNPQRSTVWY"

# UniProtKB/Swiss-Prot amino-acid composition (%)
SWISSPROT = {
    "A": 8.25, "R": 5.53, "N": 4.06, "D": 5.45, "C": 1.37, "Q": 3.93, "E": 6.75,
    "G": 7.07, "H": 2.27, "I": 5.96, "L": 9.66, "K": 5.84, "M": 2.42, "F": 3.86,
    "P": 4.70, "S": 6.56, "T": 5.34, "W": 1.08, "Y": 2.92, "V": 6.87,
}
BG = {a: SWISSPROT[a] / 100 for a in AA}

# =====================================================
# 1. DATA
# =====================================================

def load_data():
    df = pd.read_csv(CV_FILE, sep="\t", header=None,
                     names=["accession", "label", "fold"], dtype=str)
    df["accession"] = df["accession"].str.strip()
    df["label"] = df["label"].astype(int)
    df["fold"] = df["fold"].astype(int)

    uni = pd.concat([pd.read_csv(f, sep="\t") for f in UNIPROT_FILES])
    uni = uni.drop_duplicates("Entry").rename(columns={"Entry": "accession"})
    df = df.merge(uni[["accession", "Signal peptide", "Sequence"]], on="accession", how="left")

    missing = df["Sequence"].isna().sum()
    if missing:
        raise SystemExit(f"{missing} proteins without sequence: check the UniProt files")

    pat = re.compile(r'^SIGNAL\s+(\d+)\.\.(\d+)')
    df["cleavage"] = df["Signal peptide"].apply(
        lambda s: int(pat.match(str(s).strip()).group(2)) if pat.match(str(s).strip()) else None
    )
    return df

# =====================================================
# 2. PSWM (log-odds), built on positives only
# =====================================================

def build_pswm(seqs, sites):
    counts = np.full((WINDOW, 20), float(PSEUDO))
    idx = {a: i for i, a in enumerate(AA)}

    used = 0
    for seq, cs in zip(seqs, sites):
        cs = int(cs)
        if cs < UPSTREAM or len(seq) < cs + DOWNSTREAM:
            continue
        window = seq[cs - UPSTREAM: cs + DOWNSTREAM]
        if len(window) != WINDOW:
            continue
        for p, a in enumerate(window):
            if a in idx:
                counts[p, idx[a]] += 1
        used += 1

    freq = counts / counts.sum(axis=1, keepdims=True)
    pswm = np.array([[math.log(freq[p, idx[a]] / BG[a]) for a in AA] for p in range(WINDOW)])
    return pswm, used

def score_sequence(seq, pswm):
    idx = {a: i for i, a in enumerate(AA)}
    best = -np.inf
    last = min(len(seq) - WINDOW, SCAN)
    for start in range(0, max(last, 0) + 1):
        window = seq[start: start + WINDOW]
        s = sum(pswm[p, idx[a]] for p, a in enumerate(window) if a in idx)
        best = max(best, s)
    return best if best > -np.inf else 0.0

# =====================================================
# 3. METRICS
# =====================================================

def metrics(y_true, y_pred):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    tp = int(((y_true == 1) & (y_pred == 1)).sum())
    tn = int(((y_true == 0) & (y_pred == 0)).sum())
    fp = int(((y_true == 0) & (y_pred == 1)).sum())
    fn = int(((y_true == 1) & (y_pred == 0)).sum())

    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    acc = (tp + tn) / len(y_true)
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / den) if den else 0.0
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "Precision": prec, "Recall": rec, "Accuracy": acc, "F1": f1, "MCC": mcc}

def best_threshold(scores, labels):
    """Threshold maximising F1 along the PR curve."""
    best_t, best_f1 = None, -1
    for t in np.unique(scores):
        f1 = metrics(labels, (scores >= t).astype(int))["F1"]
        if f1 > best_f1:
            best_f1, best_t = f1, t
    return best_t, best_f1

# =====================================================
# 4. 5-FOLD CROSS-VALIDATION
# =====================================================

def main():
    df = load_data()
    folds = sorted(df["fold"].unique())
    rows = []

    for test in folds:
        val = folds[(folds.index(test) + 1) % len(folds)]
        train = [f for f in folds if f not in (test, val)]

        tr = df[df["fold"].isin(train) & (df["label"] == 1) & df["cleavage"].notna()]
        pswm, used = build_pswm(tr["Sequence"], tr["cleavage"])

        dv = df[df["fold"] == val]
        vs = np.array([score_sequence(s, pswm) for s in dv["Sequence"]])
        thr, vf1 = best_threshold(vs, dv["label"].values)

        dt = df[df["fold"] == test]
        ts = np.array([score_sequence(s, pswm) for s in dt["Sequence"]])
        m = metrics(dt["label"].values, (ts >= thr).astype(int))
        m.update({"fold": test, "threshold": thr, "n_train_SP": used, "val_F1": vf1})
        rows.append(m)
        print(f"test fold {test} | train {train} | val {val} | "
              f"thr={thr:.3f} | MCC={m['MCC']:.3f} F1={m['F1']:.3f}")

    res = pd.DataFrame(rows).set_index("fold")
    cols = ["Precision", "Recall", "Accuracy", "F1", "MCC"]

    summary = pd.DataFrame({
        "mean": res[cols].mean(),
        "stderr": res[cols].std(ddof=1) / math.sqrt(len(res)),
    })

    print("\nPER-FOLD\n", res[["threshold"] + cols].round(3).to_string())
    print("\nAVERAGE OVER THE 5 FOLDS\n", summary.round(3).to_string())

    res.to_csv("cv_per_fold.tsv", sep="\t")
    summary.to_csv("cv_summary.tsv", sep="\t")

    # final model: PSWM on all the training set, threshold = mean of the folds
    allpos = df[(df["label"] == 1) & df["cleavage"].notna()]
    pswm, used = build_pswm(allpos["Sequence"], allpos["cleavage"])
    np.savetxt("pswm_final.tsv", pswm, delimiter="\t", header="\t".join(AA), comments="")
    print(f"\nFinal PSWM on {used} signal peptides -> pswm_final.tsv")
    print("Threshold to use on the benchmark:", round(res['threshold'].mean(), 3))

if __name__ == "__main__":
    main()
