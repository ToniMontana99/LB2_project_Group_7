import gzip
import pickle
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    precision_score,
    recall_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, normalize
from sklearn.svm import SVC

TRAINING_FILE = "../0.4_Von_Heijne/training_5folds.tsv"
BENCHMARK_FILE = "../0.2_data_pre_processing/0.2.2_train/benchmarking_annotated.tsv"

AMINO_ACIDS = list("ARNDCQEGHILKMFPSTWYV")
AA_TO_COL = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

POSITIONS = [-13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2]
VH_WINDOW = len(POSITIONS)
VH_SCAN = 90

BACKGROUND_PERCENT = {
    "L": 9.65, "A": 8.26, "G": 7.07, "V": 6.85, "E": 6.71,
    "S": 6.67, "I": 5.90, "K": 5.79, "R": 5.53, "D": 5.46,
    "T": 5.37, "P": 4.75, "N": 4.06, "Q": 3.93, "F": 3.87,
    "Y": 2.92, "M": 2.41, "H": 2.28, "C": 1.39, "W": 1.11,
}
BACKGROUND = {aa: v / 100 for aa, v in BACKGROUND_PERCENT.items()}

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2
}
HYDRO_LENGTH = 40
HYDRO_WINDOW = 5
NREGION_LENGTH = 5
CHARGE = {"K": 1, "R": 1, "D": -1, "E": -1}

RF_PARAMS = dict(n_estimators = 500, random_state = 0, n_jobs = -1)
C_GRID = [0.25, 0.5, 1, 4, 16, 64]
GAMMA_FACTORS = [0.0625, 0.125, 0.25, 0.5, 1, 2]


def compute_metrics(y_true, y_pred): 
    return {
        "MCC": matthews_corrcoef(y_true, y_pred),
        "precision": precision_score(y_true, y_pred, zero_division = 0),
        "recall": recall_score(y_true, y_pred, zero_division = 0),
        "F1": f1_score(y_true, y_pred, zero_division = 0),
        "accuracy": accuracy_score(y_true, y_pred),
    }

#VON HEJNE METHOD'S FUNCTIONS

def get_context(row):
    sequence = str(row["sequence"]).strip().upper()
    try:                                            # check for string type before validating numeric format
        cleavage = int(float(row["cleavage_site"]))
    except (ValueError, TypeError):
        return None

    start = cleavage - 13
    end = cleavage + 2
    if start < 0 or end > len(sequence):
        return None
    context = sequence[start:end]
    return context if len(context) == VH_WINDOW else None


def build_pswm(contexts):
    count_matrix = pd.DataFrame(1, index = POSITIONS, columns = AMINO_ACIDS)
    for ctx in contexts:
        for pos, aa in zip(POSITIONS, ctx):
            if aa in AA_TO_COL:
                count_matrix.at[pos, aa] += 1
    pspm = count_matrix.div(count_matrix.sum(axis = 1), axis = 0)
    pswm = pspm.copy()
    for aa in pswm.columns:
        pswm[aa] = np.log(pspm[aa] / BACKGROUND[aa])
    return pswm


def von_heijne_score(sequence, weights):
    sequence = str(sequence).strip().upper()[:VH_SCAN]
    if len(sequence) < VH_WINDOW:
        return np.nan
    best_score = -np.inf
    for start in range(len(sequence) - VH_WINDOW + 1):
        window = sequence[start:start + VH_WINDOW]
        score = 0.0
        valid = True
        for i, aa in enumerate(window):
            if aa not in AA_TO_COL:
                valid = False
                break
            score += weights[i, AA_TO_COL[aa]]
        if valid and score > best_score:
            best_score = score
    return best_score if best_score > -np.inf else np.nan


#SVM METHOD'S FUNCTIONS

COMPOSITION_VECTORIZER = CountVectorizer(
    analyzer = "char", lowercase = False, vocabulary = AMINO_ACIDS
)

def median_sp_length(cleavage_sites, y):
    sp_lengths = pd.Series(np.asarray(cleavage_sites)[np.asarray(y) == 1]).dropna()
    if len(sp_lengths) == 0:
        raise ValueError("Nessun positivo nel training set per calcolare n.")
    return int(round(sp_lengths.median()))


def aa_composition(sequences, n):
    counts = COMPOSITION_VECTORIZER.transform(sequences.str[:n]).toarray().astype(float)
    counts[counts.sum(axis = 1) == 0] = 1.0
    return normalize(counts, norm="l1")


def hydrophobicity_profile(sequence, scale = KD, w = HYDRO_WINDOW):
    values = np.array([scale.get(aa, np.nan) for aa in sequence], dtype = float)
    known = ~np.isnan(values)
    if not known.any():
        return np.array([np.nan])
    kernel = np.ones(w)
    if len(values) < w:
        return np.array([values[known].mean()])
    total = np.convolve(np.where(known, values, 0.0), kernel, mode = "same")
    count = np.convolve(known.astype(float), kernel, mode = "same")
    with np.errstate(invalid = "ignore", divide = "ignore"):
        return np.where(count > 0, total / count, np.nan)


class SPFeatures(BaseEstimator, TransformerMixin):
    def __init__(self, hydro_length = HYDRO_LENGTH, hydro_window = HYDRO_WINDOW):
        self.hydro_length = hydro_length
        self.hydro_window = hydro_window

    def fit(self, X, y):
        self.n_ = median_sp_length(X["cleavage_site"], y)
        return self

    def transform(self, X):
        sequences = X["sequence"].astype(str).str.strip().str.upper()
        composition = aa_composition(sequences, self.n_)
        hydro = []
        for s in sequences:
            p = hydrophobicity_profile(s[:self.hydro_length], KD, self.hydro_window)
            hydro.append([np.nanmax(p), np.nanmean(p), np.nanargmax(p) + 1])
        charge = np.array([[sum(CHARGE.get(aa, 0) for aa in s[:NREGION_LENGTH])] for s in sequences])
        return np.hstack([composition, np.array(hydro), charge])

    def get_feature_names_out(self, input_features = None):
        return np.array(
            [f"AAC_{aa}" for aa in AMINO_ACIDS]
            + ["KD_max", "KD_mean", "KD_max_position", "nregion_charge"]
        )


def rf_ranking(X_train, y_train):
    rf = RandomForestClassifier(**RF_PARAMS).fit(X_train, y_train)
    return np.argsort(rf.feature_importances_)[::-1]


def k_grid(D):
    return sorted({int(np.ceil(D / 4)), int(np.ceil(D / 2)), D})


def choose_hyperparameters(X_df, y, fold):
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
            scaler = StandardScaler().fit(X_fit[:, top])
            Xf, Xv = scaler.transform(X_fit[:, top]), scaler.transform(X_val[:, top])
            for C in C_GRID:
                for factor in GAMMA_FACTORS:
                    svm = SVC(kernel = "rbf", C = C, gamma = factor / k).fit(Xf, y_fit)
                    mcc = matthews_corrcoef(y_val, svm.predict(Xv))
                    scores.setdefault((k, C, factor), []).append(mcc)

    mean_scores = {combo: np.mean(v) for combo, v in scores.items()}
    best = max(mean_scores, key = mean_scores.get)
    return best, mean_scores


#TRAINING AND BENCHMARK TESTING ON BOTH MODELS

if __name__ == "__main__":
    #training data loading
    data = pd.read_csv(TRAINING_FILE, sep = "\t").dropna(subset = ["label", "sequence", "fold"]).copy()
    data["label"] = data["label"].astype(int)
    data["fold"] = data["fold"].astype(int)

    #benchmark data loading
    bench = pd.read_csv(BENCHMARK_FILE, sep="\t").dropna(subset = ["sequence", "label"]).copy()
    bench["label"] = bench["label"].astype(int)
    if "cleavage_site" not in bench.columns:
        bench["cleavage_site"] = np.nan

    print("=" * 65)
    print("1. TRAINING of final VON HEJNE model")
    print("=" * 65)

    #Filter positive entries, extract valid cleavage site windows, and drop invalid rows
    pos_data = data.query("label == 1").copy()
    pos_data["cleavage_context"] = pos_data.apply(get_context, axis = 1)
    pos_contexts = pos_data.dropna(subset=["cleavage_context"])

    #STEP1: training the pswm (without fold 0, used as validation fold) and threshold estimation on validation fold
    val_fold_vh = 4
    train_ctx = pos_contexts[pos_contexts["fold"] != val_fold_vh]["cleavage_context"]
    pswm_calib = build_pswm(train_ctx)
    weights_calib = pswm_calib.loc[POSITIONS, AMINO_ACIDS].to_numpy()

    val_data = data[data["fold"] == val_fold_vh].copy()
    val_data["score"] = val_data["sequence"].apply(lambda s: von_heijne_score(s, weights_calib))
    val_data = val_data.dropna(subset = ["score"])

    precision, recall, thresholds = precision_recall_curve(val_data["label"], val_data["score"])
    f1_curve = (2 * precision[:-1] * recall[:-1]) / (precision[:-1] + recall[:-1] + 1e-12)
    best_vh_threshold = thresholds[np.argmax(f1_curve)]

    #STEP2: Final retraining of the pswm on the positives of the whole training dataset
    final_pswm = build_pswm(pos_contexts["cleavage_context"])
    final_pswm.to_csv("pswm_final.tsv", sep = "\t")
    final_vh_weights = final_pswm.loc[POSITIONS, AMINO_ACIDS].to_numpy()

    print(f"-> Final PSWM built on {len(pos_contexts)} positive contexts.")
    print(f"-> Decision threshold set on fold {val_fold_vh}: {best_vh_threshold:.4f}\n")

    print("=" * 65)
    print("2. TTraining of final SVM model")
    print("=" * 65)

    #STEP1: 5-fold cross-validation to select hyperparameters
    X_df = data[["sequence", "cleavage_site"]]
    y = data["label"].to_numpy()
    (best_k, best_C, best_factor), _ = choose_hyperparameters(X_df, y, data["fold"].to_numpy())

    #STEP2: Final fit of the entire pipeline on the full dataset
    final_svm_pipeline = Pipeline([
        ("features", SPFeatures()),
        ("select", SelectFromModel(RandomForestClassifier(**RF_PARAMS), threshold = -np.inf, max_features = best_k)),
        ("scale", StandardScaler()),
        ("svm", SVC(kernel = "rbf", C = best_C, gamma = best_factor / best_k)),
    ])
    final_svm_pipeline.fit(X_df, y)

    with gzip.open("svm_final_model.pkl.gz", "wb") as fh:
        pickle.dump(final_svm_pipeline, fh)

    feature_extractor = final_svm_pipeline.named_steps["features"]
    feature_selector  = final_svm_pipeline.named_steps["select"]

    all_feature_names = feature_extractor.get_feature_names_out() #retrieve all the 24 features

    # 3. Boolean mask: True for the features kept by the selection, False for the discarded ones
    selected_mask = feature_selector.get_support() #boolean mask(True for important features)

    selected_feat = all_feature_names[selected_mask] #retain True-masked features
    
    print(f"-> Hyperparameters: k = {best_k}, C = {best_C}, gamma = {best_factor}/{best_k}")
    print(f"-> Selected features ({len(selected_feat)}): {', '.join(selected_feat)}\n")

#MODELS TESTING ON BENCHMARK SET
    print("=" * 65)
    print("3. EVALUATION ON THE BENCHMARKING SET")
    print("=" * 65)

    y_bench = bench["label"].to_numpy() 

    # VON HEJNE inference
    bench["vh_score"] = bench["sequence"].apply(lambda s: von_heijne_score(s, final_vh_weights))
    vh_pred = (bench["vh_score"] >= best_vh_threshold).astype(int)
    vh_metrics = compute_metrics(y_bench, vh_pred)

    # SVM inference
    svm_pred = final_svm_pipeline.predict(bench[["sequence", "cleavage_site"]])
    svm_metrics = compute_metrics(y_bench, svm_pred)

    # # Comparative summary table
    comparison = pd.DataFrame({
        "Metric": list(vh_metrics.keys()),
        "von Heijne": list(vh_metrics.values()),
        "SVM": list(svm_metrics.values()),
    })
    print(comparison.round(4).to_string(index = False))

    # Benchmark results 
    bench["pred_von_heijne"] = vh_pred
    bench["pred_svm"] = svm_pred
    bench.to_csv("benchmarking_predictions.tsv", sep = "\t", index = False)
    comparison.to_csv("benchmark_comparison_metrics.tsv", sep = "\t", index = False)
    print(type(final_svm_pipeline))