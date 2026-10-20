import gzip
import pickle
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import (accuracy_score, f1_score, matthews_corrcoef, precision_score, recall_score)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.svm import SVC


# 0 - SETTINGS AND DATA CHECK 

TRAINING_FILE = "training_5folds.tsv"
AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")

# Kyte-Doolittle hydropathy scale (same values as Bio.SeqUtils.ProtParamData.kd)
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
HYDRO_LENGTH = 40      # N-terminal portion analysed, long enough to contain the whole SP
HYDRO_WINDOW = 5       # sliding window of the profile
NREGION_LENGTH = 5     # first residues used for the n-region charge
CHARGE = {"K": 1, "R": 1, "D": -1, "E": -1}    # residue charge at physiological pH

# Random Forest: fixed parameters (not tuned)
RF_PARAMS = dict(n_estimators = 500, random_state = 0, n_jobs = -1)

# minimal grid (slide "A minimal grid for the inner CV")
C_GRID = [0.25, 0.5, 1, 4, 16, 64]
GAMMA_FACTORS = [0.0625, 0.125, 0.25, 0.5, 1, 2]        # gamma = factor / k  (1/k = gamma='scale' on k features)

def check_data():
    """Quick look at the raw file before anything else."""
    raw = pd.read_csv(TRAINING_FILE, sep = "\t")
    print("Total rows:", len(raw))
    print("Columns:", raw.columns.tolist())
    print("Class distribution (NaN = malformed rows):")
    print(raw["label"].value_counts(dropna = False).to_string(), "\n")


# 1 - THE DATA AND THE 5 PREDEFINED SUBSETS

def load_data():
    data = pd.read_csv(TRAINING_FILE, sep = "\t")
    data = data.dropna(subset = ["label", "sequence", "fold"]).copy()
    data["label"] = data["label"].astype(int)
    data["fold"] = data["fold"].astype(int)
    return data.reset_index(drop = True) # reset the index after removing rows with NaN values


# 2 - MEDIAN SIGNAL-PEPTIDE LENGTH (n)

def median_sp_length(cleavage_sites, y):
    """n = median SP length of the positives in the training set passed in (rounded)."""
    sp_lengths = pd.Series(np.asarray(cleavage_sites)[np.asarray(y) == 1]).dropna() # for every protein create pd series, where y is the label 0, 1 and cleavage site , if positive, report also its lenght 
    if len(sp_lengths) == 0:
        raise ValueError("no positives in this training set: n cannot be computed")
    return int(round(sp_lengths.median()))           # round, do not truncate


# 3 - AMINO-ACID COMPOSITION OF THE FIRST n RESIDUES

#Initialize a character-level vectorizer using a predefined amino acid vocabulary, preserving case sensitivity for protein sequence analysis.
COMPOSITION_VECTORIZER = CountVectorizer(analyzer = "char", lowercase = False,
                                         vocabulary = AMINO_ACIDS)

def aa_composition(sequences, n):
    """20-D composition of the first n residues; non-standard residues ignored; rows sum to 1."""
    counts = COMPOSITION_VECTORIZER.transform(sequences.str[:n]).toarray().astype(float)
    counts[counts.sum(axis=1) == 0] = 1.0             # no standard residue: uniform vector
    return normalize(counts, norm = "l1") # normalization dividing each aa count for the length of the sequence, and this for each protein 


# 4 - FEATURE EXTRACTION (composition + hydrophobicity, D = 22)

def hydrophobicity_profile(sequence, scale=KD, w=HYDRO_WINDOW):
    """
    One value per residue: mean of the scale over the window centred on it.
    Sequence ends and unknown residues are handled with the shrinking window
    (average only over the known residues that exist), as in the course notebook.
    """
    values = np.array([scale.get(aa, np.nan) for aa in sequence], dtype=float)
    known = ~np.isnan(values)
    if not known.any():
        return np.array([np.nan])

    kernel = np.ones(w)
    if len(values) < w:                          # shorter than one window
        return np.array([values[known].mean()])

    total = np.convolve(np.where(known, values, 0.0), kernel, mode = "same")
    count = np.convolve(known.astype(float), kernel, mode = "same")
    with np.errstate(invalid = "ignore", divide = "ignore"):
        return np.where(count > 0, total / count, np.nan)

class SPFeatures(BaseEstimator, TransformerMixin):
    """
    X: DataFrame with the columns 'sequence' and 'cleavage_site'; y: labels (1 = SP).

    fit()       learns n = median SP length of the positives in the training set
    transform() 20 composition features of the first n residues (sum to 1)
                + max, mean and position of the max hydrophobicity (first HYDRO_LENGTH residues)
                + net charge of the n-region (first NREGION_LENGTH residues)
    n is learned like the mean/std of a StandardScaler: on training data only,
    then applied unchanged to validation and test data.
    """

    def __init__(self, hydro_length=HYDRO_LENGTH, hydro_window=HYDRO_WINDOW):
        self.hydro_length = hydro_length
        self.hydro_window = hydro_window

    def fit(self, X, y):
        self.n_ = median_sp_length(X["cleavage_site"], y)      # step 2, training only
        return self

    def transform(self, X):
        sequences = X["sequence"].astype(str).str.strip().str.upper()

        composition = aa_composition(sequences, self.n_)        # step 3

        # 3 hydrophobicity features: max, mean and POSITION of the maximum (1-based)
        # SP h-region peaks around residues 10-15, a TM helix further away (slide 33)
        hydro = []
        for s in sequences:
            p = hydrophobicity_profile(s[:self.hydro_length], KD, self.hydro_window)
            hydro.append([np.nanmax(p), np.nanmean(p), np.nanargmax(p) + 1])
        hydro = np.array(hydro)

        # 1 charge feature: net charge of the n-region (K, R positive; D, E negative)
        charge = np.array([[sum(CHARGE.get(aa, 0) for aa in s[:NREGION_LENGTH])]
                           for s in sequences])

        return np.hstack([composition, hydro, charge])

    def get_feature_names_out(self, input_features=None):
        return np.array([f"AAC_{aa}" for aa in AMINO_ACIDS]
                        + ["KD_max", "KD_mean", "KD_max_position", "nregion_charge"])


# 5 - FEATURE SELECTION: RANDOM-FOREST RANKING

def rf_ranking(X_train, y_train):
    """Feature indices sorted by RF importance, most important first (one RF per training set)."""
    rf = RandomForestClassifier(**RF_PARAMS).fit(X_train, y_train)
    return np.argsort(rf.feature_importances_)[::-1]


# 6 - NESTED CROSS-VALIDATION (k, C, gamma)

def k_grid(D):
    """k in {D/4, D/2, D}, rounded up: D = 22 -> [6, 11, 22]."""
    return sorted({int(np.ceil(D / 4)), int(np.ceil(D / 2)), D})

def metrics(y_true, y_pred):
    return {
        "MCC": matthews_corrcoef(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division = 0),
        "recall": recall_score(y_true, y_pred, zero_division = 0),
        "F1": f1_score(y_true, y_pred, zero_division = 0),
        "accuracy": accuracy_score(y_true, y_pred),
    }

def fit_and_predict(X_df_train, y_train, X_df_test, k, C, gamma):
    """n, RF ranking, scaler and SVM all fitted on the training part; one prediction on the test part."""
    features = SPFeatures().fit(X_df_train, y_train)
    X_train, X_test = features.transform(X_df_train), features.transform(X_df_test)

    top = rf_ranking(X_train, y_train)[:k]
    scaler = StandardScaler().fit(X_train[:, top])
    svm = SVC(kernel = "rbf", C = C, gamma = gamma).fit(scaler.transform(X_train[:, top]), y_train)
    return svm.predict(scaler.transform(X_test[:, top])), features.n_

def choose_hyperparameters(X_df, y, fold):
    """
    Inner CV on the subsets present in `fold`: each one is the validation set once.
    In every inner split, n and the RF ranking are computed once and reused for all
    the grid combinations. Returns the best (k, C, gamma) by mean validation MCC.
    """
    scores = {}
    for val_fold in np.unique(fold):
        fit_mask, val_mask = fold != val_fold, fold == val_fold

        features = SPFeatures().fit(X_df[fit_mask], y[fit_mask])
        X_fit = features.transform(X_df[fit_mask])
        X_val = features.transform(X_df[val_mask])
        y_fit, y_val = y[fit_mask], y[val_mask]

        ranking = rf_ranking(X_fit, y_fit)
        for k in k_grid(X_fit.shape[1]):
            top = ranking[:k]
            scaler = StandardScaler().fit(X_fit[:, top])        # one scaler per k
            Xf, Xv = scaler.transform(X_fit[:, top]), scaler.transform(X_val[:, top])
            for C in C_GRID:
                for factor in GAMMA_FACTORS:
                    svm = SVC(kernel="rbf", C = C, gamma = factor / k).fit(Xf, y_fit)
                    mcc = matthews_corrcoef(y_val, svm.predict(Xv))
                    scores.setdefault((k, C, factor), []).append(mcc)

    mean_scores = {combo: np.mean(v) for combo, v in scores.items()}
    best = max(mean_scores, key=mean_scores.get)
    return best, mean_scores

def warn_if_on_edge(k, C, factor):
    """The slides ask to extend the grid if the best value lies on its edge."""
    edges = []
    if C in (min(C_GRID), max(C_GRID)):
        edges.append(f"C = {C}")
    if factor in (min(GAMMA_FACTORS), max(GAMMA_FACTORS)):
        edges.append(f"gamma = {factor}/k")
    if edges:
        print(f"note: best {', '.join(edges)} on the edge of the grid -> consider extending it")

def nested_cv(data):
    X_df = data[["sequence", "cleavage_site"]]
    y = data["label"].to_numpy()
    fold = data["fold"].to_numpy()

    rows, inner_rows = [], []
    for test_fold in range(5):
        train_mask, test_mask = fold != test_fold, fold == test_fold

        # INNER loop: model selection on the 4 training subsets only
        (k, C, factor), mean_scores = choose_hyperparameters(
            X_df[train_mask], y[train_mask], fold[train_mask])
        for (kk, cc, ff), s in mean_scores.items():
            inner_rows.append({"run": test_fold + 1, "k": kk, "C": cc,
                               "gamma": f"{ff}/k", "mean_inner_MCC": s})

        # refit on all 4 training subsets (n and RF ranking recomputed), test ONCE
        y_pred, n = fit_and_predict(X_df[train_mask], y[train_mask],
                                    X_df[test_mask], k, C, factor / k)
        m = metrics(y[test_mask], y_pred)
        rows.append({"run": test_fold + 1, "test subset": test_fold, "n": n, "k": k,
                     "C": C, "gamma": f"{factor}/k = {factor / k:.4f}", **m})

        print(f"run {test_fold + 1} (test subset {test_fold}): n = {n}  k = {k}  C = {C}  "
              f"gamma = {factor}/k  ->  MCC {m['MCC']:.3f}")
        warn_if_on_edge(k, C, factor)

    return pd.DataFrame(rows), pd.DataFrame(inner_rows)

def summarize(per_fold):
    names = ["MCC", "precision", "recall", "F1", "accuracy"]
    return pd.DataFrame({
        "metric": names,
        "mean": [per_fold[m].mean() for m in names],
        "standard_error": [per_fold[m].std(ddof=1) / np.sqrt(len(per_fold)) for m in names],
    })


# 7 - FINAL MODEL (whole training set, for the benchmark)

def final_model(data):
    """
    Hyper-parameters chosen by CV on all 5 subsets, then one pipeline fitted on the
    whole training set. SelectFromModel with threshold=-inf keeps exactly the top-k
    features by RF importance, i.e. the same selection used in the CV.
    """
    X_df = data[["sequence", "cleavage_site"]]
    y = data["label"].to_numpy()
    (k, C, factor), _ = choose_hyperparameters(X_df, y, data["fold"].to_numpy())

    model = Pipeline([
        ("features", SPFeatures()),
        ("select", SelectFromModel(RandomForestClassifier(**RF_PARAMS),
                                   threshold = -np.inf, max_features = k)),
        ("scale", StandardScaler()),
        ("svm", SVC(kernel = "rbf", C = C, gamma = factor / k)),
    ])
    model.fit(X_df, y)

    with gzip.open("svm_final_model.pkl.gz", "wb") as fh:
        pickle.dump(model, fh)

    selected = model.named_steps["features"].get_feature_names_out()[
        model.named_steps["select"].get_support()]
    print(f"Final model: n = {model.named_steps['features'].n_}  k = {k}  C = {C}  "
          f"gamma = {factor}/k")
    print("Selected features:", ", ".join(selected))
    return model


# 8 - RUN ALL THE STEPS

if __name__ == "__main__":
    check_data()                   # step 0
    data = load_data()             # step 1
    print(f"Training set: {len(data)} proteins  "
          f"({data['label'].sum()} positives, {(data['label'] == 0).sum()} negatives)")
    print(data.groupby("fold")["label"].agg(["count", "sum"])
          .rename(columns={"count": "proteins", "sum": "positives"}).to_string(), "\n")

    per_fold, inner = nested_cv(data)
    summary = summarize(per_fold)

    per_fold.to_csv("svm_cv_per_fold.tsv", sep = "\t", index = False)
    inner.to_csv("svm_inner_scores.tsv", sep = "\t", index = False)
    summary.to_csv("svm_cv_summary.tsv", sep = "\t", index = False)

    print("\nPER-FOLD RESULTS")
    print(per_fold.round(4).to_string(index = False))
    print("\nNESTED CV SUMMARY (mean +- standard error over the 5 outer test folds)")
    for _, r in summary.iterrows():
        print(f"  {r['metric']:10s} {r['mean']:.4f} +- {r['standard_error']:.4f}")

    print()
    final_model(data)
