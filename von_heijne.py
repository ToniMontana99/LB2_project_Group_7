import numpy as np
import pandas as pd
from sklearn.metrics import (
    precision_recall_curve,
    confusion_matrix,
    precision_score,
    recall_score,
    accuracy_score,
    matthews_corrcoef,
    f1_score,
)

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")
AA_TO_COL = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
POSITIONS = [-13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2]
WINDOW = len(POSITIONS)
SCAN = 90          # only the N-terminal region is scanned

# UniProtKB/Swiss-Prot amino-acid composition (%)
BACKGROUND_PERCENT = {
    "L": 9.65, "A": 8.26, "G": 7.07, "V": 6.85, "E": 6.71,
    "S": 6.67, "I": 5.90, "K": 5.79, "R": 5.53, "D": 5.46,
    "T": 5.37, "P": 4.75, "N": 4.06, "Q": 3.93, "F": 3.87,
    "Y": 2.92, "M": 2.41, "H": 2.28, "C": 1.39, "W": 1.11,
}
BACKGROUND = {aa: value / 100 for aa, value in BACKGROUND_PERCENT.items()}

# 1 - POSITIVE DATASETS AND CLEAVAGE-SITE CONTEXTS

def get_context(row):
    """Return the 15 residues around the cleavage site, or None if out of range."""
    sequence = str(row["sequence"]).strip().upper()
    cleavage = int(float(row["cleavage_site"]))

    start = cleavage - 13
    end = cleavage + 2

    if start < 0 or end > len(sequence):
        return None

    context = sequence[start:end]
    return context if len(context) == WINDOW else None

training = pd.read_csv("training_annotated.tsv", sep="\t").query("label == 1").copy()
benchmarking = pd.read_csv("benchmarking_annotated.tsv", sep="\t").query("label == 1").copy()

# axis=1 because the function needs both the sequence and the cleavage site of each row
training["cleavage_context"] = training.apply(get_context, axis=1)
benchmarking["cleavage_context"] = benchmarking.apply(get_context, axis=1)

training = training.dropna(subset=["cleavage_context"])
benchmarking = benchmarking.dropna(subset=["cleavage_context"])

training.to_csv("training_positive_contexts.tsv", sep="\t", index=False)
benchmarking.to_csv("benchmarking_positive_contexts.tsv", sep="\t", index=False)

print(f"Positive contexts: {len(training)} (training), {len(benchmarking)} (benchmarking)")

# 2 - FOLD ASSIGNMENT FROM THE 5 SUBSETS

folds = []
for i in range(1, 6):
    for name in (f"positives_subset_{i}.tsv", f"negatives_subset_{i}.tsv"):
        d = pd.read_csv(name, sep="\t", header=None, names=["accession", "label"])
        d["fold"] = i - 1
        folds.append(d) 

folds = pd.concat(folds)
folds["accession"] = folds["accession"].str.strip()

df = pd.read_csv("training_annotated.tsv", sep="\t")
df = df.dropna(subset=["label", "sequence"])
df["accession"] = df["accession"].str.strip()
df = df.merge(folds[["accession", "fold"]], on="accession", how="left")

missing = df["fold"].isna().sum()
if missing:
    raise SystemExit(f"{missing} proteins without a fold: check the subset files")

df["fold"] = df["fold"].astype(int)
df.to_csv("training_5folds.tsv", sep="\t", index=False)

# 3 - PSWM AND SCORING FUNCTIONS

def build_pswm(contexts):
    """Count matrix with pseudocount 1 -> PSPM -> log-odds PSWM."""
    count_matrix = pd.DataFrame(1, index=POSITIONS, columns=AMINO_ACIDS)

    # zip iterates over positions and context residues at the same time
    for context in contexts:
        for position, aa in zip(POSITIONS, context):
            if aa in AA_TO_COL:
                count_matrix.at[position, aa] += 1

    pspm = count_matrix.div(count_matrix.sum(axis=1), axis=0)
    pswm = pspm.copy()
    for aa in pswm.columns:
        pswm[aa] = np.log(pspm[aa] / BACKGROUND[aa])

    return pswm

def von_heijne_score(sequence, weights):
    """Highest window score over the N-terminal region; NaN if no valid window."""
    sequence = str(sequence).strip().upper()[:SCAN]

    if len(sequence) < WINDOW:
        return np.nan

    best_score = -np.inf

    for start in range(len(sequence) - WINDOW + 1):
        window = sequence[start:start + WINDOW]

        score = 0.0
        valid = True
        for i, aa in enumerate(window):
            # non-standard residues invalidate the window
            if aa not in AA_TO_COL:
                valid = False
                break
            score += weights[i, AA_TO_COL[aa]]

        if valid and score > best_score:
            best_score = score

    return best_score if best_score > -np.inf else np.nan

# 4 - 5-FOLD CROSS-VALIDATION

contexts = training.merge(folds[["accession", "fold"]], on="accession", how="left")
results = []

for j in range(5):
    test_fold = j
    validation_fold = (j + 1) % 5
    training_folds = [i for i in range(5) if i not in (test_fold, validation_fold)]

    print(f"Run {j + 1}/5 - train {training_folds}, validation {validation_fold}, "
          f"test {test_fold}", flush=True)

    # PSWM on the positives of the three training folds
    train_contexts = contexts[contexts["fold"].isin(training_folds)]["cleavage_context"]
    pswm = build_pswm(train_contexts)
    pswm.to_csv(f"pswm_run_{j + 1}.tsv", sep="\t")
    weights = pswm.loc[POSITIONS, AMINO_ACIDS].to_numpy()

    # threshold on the validation fold, from the PR curve
    validation = df[df["fold"] == validation_fold].copy()
    validation["score"] = validation["sequence"].apply(
        lambda sequence: von_heijne_score(sequence, weights)
    )
    validation = validation.dropna(subset=["score"])

    precision, recall, thresholds = precision_recall_curve(
        validation["label"], validation["score"]
    )
    f1_curve = (2 * precision[:-1] * recall[:-1]
                / (precision[:-1] + recall[:-1] + 1e-12))
    threshold = thresholds[np.argmax(f1_curve)]

    # evaluation on the test fold
    test = df[df["fold"] == test_fold].copy()
    test["score"] = test["sequence"].apply(
        lambda sequence: von_heijne_score(sequence, weights)
    )
    if test["score"].isna().any():
        raise ValueError(f"Run {j + 1}: some sequences could not be scored")

    test["prediction"] = (test["score"] >= threshold).astype(int)
    test.to_csv(f"test_predictions_run_{j + 1}.tsv", sep="\t", index=False)

    y_true = test["label"]
    y_pred = test["prediction"]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()

    results.append({
        "run": j + 1,
        "threshold": threshold,
        "TN": tn, "FP": fp, "FN": fn, "TP": tp,
        "precision": precision_score(y_true, y_pred, zero_division=0),
        "recall": recall_score(y_true, y_pred, zero_division=0),
        "accuracy": accuracy_score(y_true, y_pred),
        "MCC": matthews_corrcoef(y_true, y_pred),
        "F1": f1_score(y_true, y_pred, zero_division=0),
    })

results = pd.DataFrame(results)
results.to_csv("test_metrics_per_run.tsv", sep="\t", index=False)

# 5 - AVERAGE AND STANDARD ERROR OVER THE 5 RUNS

metrics = ["precision", "recall", "accuracy", "MCC", "F1"]

summary = pd.DataFrame({
    "metric": metrics,
    "mean": [results[m].mean() for m in metrics],
    "standard_error": [results[m].std(ddof=1) / np.sqrt(5) for m in metrics],
})
summary.to_csv("cross_validation_summary.tsv", sep="\t", index=False)

summary["result"] = (summary["mean"].round(4).astype(str)
                     + " ± "
                     + summary["standard_error"].round(4).astype(str))
summary[["metric", "result"]].to_csv(
    "final_cross_validation_results.tsv", sep="\t", index=False
)

# 6 - FINAL MODEL FOR THE BENCHMARK

final_pswm = build_pswm(training["cleavage_context"])
final_pswm.to_csv("pswm_final.tsv", sep="\t")
final_threshold = results["threshold"].mean()

print("\nPER-RUN RESULTS")
print(results.round(4).to_string(index=False))
print("\nCROSS-VALIDATION SUMMARY")
print(summary[["metric", "result"]].to_string(index=False))
print(f"\nFinal PSWM on {len(training)} signal peptides -> pswm_final.tsv")
print(f"Threshold for the benchmark: {final_threshold:.4f}")