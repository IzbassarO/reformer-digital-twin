"""ML surrogates of the digital twin with conformal prediction intervals.

Training data: open-loop LHS campaign ``data/lhs_runs/lhs_plantA_v1.csv.gz`` (2000 points, seed 0); test data:
the Saltelli sample ``sobol_saltelli_v1.csv.gz`` (4096 points, different sampling scheme, honest generalisation).
Seven operating inputs (see :mod:`rdt.operate`) plus two optional physics-informed features:

* ``x_eq_CH4``: equilibrium CH4 conversion of the feed at the base outlet temperature (1105.55 K) and the inlet
  pressure (Cantera, gri30), a thermodynamic ceiling that depends on S/C and P;
* ``specific_duty``: combustion heat per unit total feed, ``Q_comb / F_feed`` [MJ/kmol], from firing, load and S/C.

Models per target: HistGradientBoosting, Gaussian process (Matern 5/2 + white noise, on a subset), MLP.
Uncertainty: MAPIE ``CrossConformalRegressor`` (CV+, 5 folds) on the best model per target; GP native intervals
for comparison. Life targets rest on the master curve selected by the run configuration.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.compose import TransformedTargetRegressor

from rdt import operate as op

ROOT = Path(__file__).resolve().parents[2]
MODELS_DIR = ROOT / "models"
METRICS_JSON = ROOT / "data" / "design" / "surrogate_metrics_v1.json"
TRAIN_FILE = op.LHS_DIR / "lhs_plantA_v1.csv.gz"
TEST_FILE = op.LHS_DIR / "sobol_saltelli_v1.csv.gz"

INPUTS = list(op.INPUT_NAMES)
PHYSICS_FEATURES = ["x_eq_CH4", "specific_duty_MJ_per_kmol"]
TARGETS = {"T_wo_max_K": "K", "z_frac_T_wo_max": "-", "T_out_K": "K", "P_out_bar": "bar", "CH4_slip_dry_pct": "%",
           "H2_net_kmol_h": "kmol/h", "duty_W": "W", "log10_t_r_hot": "-"}
MODEL_NAMES = ("hgb", "gp", "mlp")
GP_SUBSET = 1000
SEED = 0


# ---------------------------------------------------------------------------
# Data and features
# ---------------------------------------------------------------------------
def load_data(train_file: Optional[Path] = None, test_file: Optional[Path] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    tr = pd.read_csv(train_file or TRAIN_FILE); te = pd.read_csv(test_file or TEST_FILE)
    tr = tr[tr.converged == True].reset_index(drop=True); te = te[te.converged == True].reset_index(drop=True)  # noqa: E712
    return tr, te


def _feed_composition(base: op.BaseCase, sc: float, load: float) -> Dict[str, float]:
    F = {s: v * load for s, v in base.dry_feed_kmol_h.items()}
    F["H2O"] = sc * base.carbon_kmol_h * load
    return F


def equilibrium_ch4_conversion(sc: np.ndarray, P_bar: np.ndarray, T_K: float = op.TARGET_T_OUT_K,
                               base: Optional[op.BaseCase] = None) -> np.ndarray:
    """Equilibrium CH4 conversion (carbon-equivalent basis, C2+ counted as CH4) at ``T_K`` and ``P_bar``."""
    import cantera as ct
    from rdt import reactor1d as r1

    base = base or op.base_case()
    gas = ct.Solution("gri30.yaml")
    out = np.empty(len(sc))
    for i, (s, P) in enumerate(zip(sc, P_bar)):
        F = r1.convert_higher_alkanes(_feed_composition(base, float(s), 1.0), rule="latham")
        gas.TPX = T_K, float(P) * 1e5, {k: v for k, v in F.items() if v > 0}
        n_c_in = sum(F[k] for k in ("CH4", "CO", "CO2"))
        gas.equilibrate("TP")
        X = gas.mole_fraction_dict()
        # moles conserved per atom: CH4 fraction of carbon at equilibrium
        c_tot = X.get("CH4", 0) + X.get("CO", 0) + X.get("CO2", 0)
        ch4_eq_per_c = X.get("CH4", 0) / c_tot
        ch4_in_per_c = F["CH4"] / n_c_in
        out[i] = 1.0 - ch4_eq_per_c / ch4_in_per_c
    return out


def add_physics_features(df: pd.DataFrame, base: Optional[op.BaseCase] = None) -> pd.DataFrame:
    base = base or op.base_case()
    df = df.copy()
    df["x_eq_CH4"] = equilibrium_ch4_conversion(df.steam_to_carbon.to_numpy(), df.inlet_P.to_numpy(), base=base)
    dry = sum(base.dry_feed_kmol_h.values())
    feed_total = df.feed_per_tube_fraction * (dry + df.steam_to_carbon * base.carbon_kmol_h)
    Q = df.specific_firing_factor * df.feed_per_tube_fraction * base.Q_comb_W / op.base_case().row.n_tubes  # W per tube
    df["specific_duty_MJ_per_kmol"] = Q * 3600 / 1e6 / feed_total
    return df


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
def make_model(name: str, n_features: int, seed: int = SEED) -> Pipeline:
    if name == "hgb":
        est = HistGradientBoostingRegressor(max_iter=600, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=10,
                                            l2_regularization=1e-3, random_state=seed)
    elif name == "gp":
        kernel = ConstantKernel(1.0, (1e-2, 1e2)) * Matern(length_scale=np.ones(n_features), length_scale_bounds=(1e-2, 1e2), nu=2.5) \
            + WhiteKernel(1e-4, (1e-8, 1e-1))
        est = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=1, random_state=seed)
    elif name == "mlp":
        est = TransformedTargetRegressor(
            regressor=MLPRegressor(hidden_layer_sizes=(64, 64), activation="tanh", alpha=1e-4, learning_rate_init=2e-3,
                                   max_iter=3000, early_stopping=True, n_iter_no_change=50, random_state=seed),
            transformer=StandardScaler())
    else:
        raise ValueError(name)
    return Pipeline([("scale", StandardScaler()), ("model", est)])


def _fit(name: str, X: np.ndarray, y: np.ndarray, seed: int = SEED) -> Pipeline:
    m = make_model(name, X.shape[1], seed)
    if name == "gp" and len(X) > GP_SUBSET:
        idx = np.random.default_rng(seed).choice(len(X), GP_SUBSET, replace=False)
        X, y = X[idx], y[idx]
    return m.fit(X, y)


def metrics(y, yhat) -> Dict[str, float]:
    return {"r2": float(r2_score(y, yhat)), "rmse": float(np.sqrt(mean_squared_error(y, yhat))),
            "max_abs_error": float(np.max(np.abs(y - yhat)))}


def predict_time_per_1e5(model, X_like: np.ndarray, n: int = 100_000) -> float:
    rng = np.random.default_rng(1)
    Xr = X_like[rng.integers(0, len(X_like), n)]
    t = time.perf_counter(); model.predict(Xr); return time.perf_counter() - t


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------
def train_all(features: Sequence[str] = INPUTS, with_physics: bool = True, save: bool = True,
              cv_folds: int = 5, seed: int = SEED, verbose: bool = True) -> Dict[str, object]:
    tr, te = load_data()
    base = op.base_case()
    if with_physics:
        tr = add_physics_features(tr, base); te = add_physics_features(te, base)
    feats_base = list(INPUTS); feats_phys = list(INPUTS) + PHYSICS_FEATURES
    result = {"train_file": str(Path(TRAIN_FILE).relative_to(ROOT)), "test_file": str(Path(TEST_FILE).relative_to(ROOT)), "models_dir": str(Path(MODELS_DIR).relative_to(ROOT)),
              "n_train": int(len(tr)), "n_test": int(len(te)), "features_base": feats_base, "features_physics": feats_phys,
              "gp_subset": GP_SUBSET, "seed": seed, "targets": {}}
    Path(MODELS_DIR).mkdir(parents=True, exist_ok=True)
    kf = KFold(cv_folds, shuffle=True, random_state=seed)
    for tgt in TARGETS:
        y_tr, y_te = tr[tgt].to_numpy(float), te[tgt].to_numpy(float)
        rt = {"unit": TARGETS[tgt], "models": {}}
        for fs_name, feats in (("base", feats_base), ("physics", feats_phys)):
            if fs_name == "physics" and not with_physics:
                continue
            X_tr, X_te = tr[feats].to_numpy(float), te[feats].to_numpy(float)
            for name in MODEL_NAMES:
                t0 = time.perf_counter(); m = _fit(name, X_tr, y_tr, seed); fit_t = time.perf_counter() - t0
                te_m = metrics(y_te, m.predict(X_te))
                # 5-fold CV inside the training set (GP on the same subset rule)
                if name == "gp":
                    idx = np.random.default_rng(seed).choice(len(X_tr), min(GP_SUBSET, len(X_tr)), replace=False)
                    cvp = cross_val_predict(make_model(name, X_tr.shape[1], seed), X_tr[idx], y_tr[idx], cv=kf, n_jobs=5)
                    cv_m = metrics(y_tr[idx], cvp)
                else:
                    cvp = cross_val_predict(make_model(name, X_tr.shape[1], seed), X_tr, y_tr, cv=kf, n_jobs=5)
                    cv_m = metrics(y_tr, cvp)
                rt["models"][f"{name}__{fs_name}"] = {"model": name, "features": fs_name, "test": te_m, "cv5": cv_m,
                                                       "fit_time_s": fit_t, "predict_time_per_1e5_s": predict_time_per_1e5(m, X_te)}
                if save:
                    joblib.dump({"pipeline": m, "features": feats, "target": tgt}, MODELS_DIR / f"{tgt}__{name}__{fs_name}.joblib")
                if verbose:
                    print(f"{tgt:18s} {name:4s} {fs_name:8s} test R2 {te_m['r2']:.4f} RMSE {te_m['rmse']:.4g} max {te_m['max_abs_error']:.4g} | cv RMSE {cv_m['rmse']:.4g} | fit {fit_t:.1f}s pred/1e5 {rt['models'][f'{name}__{fs_name}']['predict_time_per_1e5_s']:.3f}s")
        best = min(rt["models"], key=lambda k: rt["models"][k]["test"]["rmse"])
        rt["best"] = best
        # physics-feature effect: best base vs best physics per model
        rt["physics_feature_rmse_change_pct"] = {name: 100 * (rt["models"][f"{name}__physics"]["test"]["rmse"] / rt["models"][f"{name}__base"]["test"]["rmse"] - 1)
                                                for name in MODEL_NAMES if f"{name}__physics" in rt["models"]}
        result["targets"][tgt] = rt
    return {"metrics": result, "train": tr, "test": te}


def conformalize(tgt: str, best_key: str, tr: pd.DataFrame, te: pd.DataFrame, confidence=(0.9, 0.95),
                 seed: int = SEED, save: bool = True) -> Dict[str, object]:
    """MAPIE CV+ conformal intervals for the best model of one target; coverage and width on the test set."""
    from mapie.regression import CrossConformalRegressor

    name, fs_name = best_key.split("__")
    feats = list(INPUTS) + (PHYSICS_FEATURES if fs_name == "physics" else [])
    X_tr, y_tr = tr[feats].to_numpy(float), tr[tgt].to_numpy(float)
    X_te, y_te = te[feats].to_numpy(float), te[tgt].to_numpy(float)
    if name == "gp":
        idx = np.random.default_rng(seed).choice(len(X_tr), min(GP_SUBSET, len(X_tr)), replace=False)
        X_tr, y_tr = X_tr[idx], y_tr[idx]
    t0 = time.perf_counter()
    cc = CrossConformalRegressor(estimator=make_model(name, X_tr.shape[1], seed), confidence_level=list(confidence),
                                 method="plus", cv=5, random_state=seed)
    cc.fit_conformalize(X_tr, y_tr)
    y_pred, intervals = cc.predict_interval(X_te)
    out = {"model": best_key, "method": "MAPIE CrossConformalRegressor CV+ (5 folds)", "fit_time_s": time.perf_counter() - t0, "levels": {}}
    for j, cl in enumerate(confidence):
        lo, hi = intervals[:, 0, j], intervals[:, 1, j]
        out["levels"][str(cl)] = {"picp": float(np.mean((y_te >= lo) & (y_te <= hi))), "mpiw": float(np.mean(hi - lo))}
    if save:
        joblib.dump({"conformal": cc, "features": feats, "target": tgt, "model": best_key}, MODELS_DIR / f"{tgt}__conformal.joblib")
    out["_pred"] = y_pred; out["_intervals"] = intervals; out["_y"] = y_te
    return out


def gp_native_intervals(tgt: str, fs_name: str, tr: pd.DataFrame, te: pd.DataFrame, confidence=(0.9, 0.95),
                        seed: int = SEED) -> Dict[str, object]:
    from scipy.stats import norm
    feats = list(INPUTS) + (PHYSICS_FEATURES if fs_name == "physics" else [])
    m = _fit("gp", tr[feats].to_numpy(float), tr[tgt].to_numpy(float), seed)
    Xs = m.named_steps["scale"].transform(te[feats].to_numpy(float))
    mu, sd = m.named_steps["model"].predict(Xs, return_std=True)
    y = te[tgt].to_numpy(float)
    out = {"levels": {}}
    for cl in confidence:
        z = norm.ppf(0.5 + cl / 2); lo, hi = mu - z * sd, mu + z * sd
        out["levels"][str(cl)] = {"picp": float(np.mean((y >= lo) & (y <= hi))), "mpiw": float(np.mean(hi - lo))}
    return out


def predict_base_case(tgt: str = "T_wo_max_K", model_key: Optional[str] = None) -> float:
    """Surrogate prediction at the Plant A base inputs (uses the best saved model unless given)."""
    if model_key is None:
        model_key = json.loads(METRICS_JSON.read_text())["targets"][tgt]["best"]
    bundle = joblib.load(MODELS_DIR / f"{tgt}__{model_key}.joblib")
    base = op.base_case()
    x = pd.DataFrame([base.inputs()])
    if any(f in PHYSICS_FEATURES for f in bundle["features"]):
        x = add_physics_features(x, base)
    return float(bundle["pipeline"].predict(x[bundle["features"]].to_numpy(float))[0])


def use_config(cfg) -> None:
    """Point the module at a rdt.config.RunConfig (train/test files, metrics JSON, models directory)."""
    global TRAIN_FILE, TEST_FILE, METRICS_JSON, MODELS_DIR
    TRAIN_FILE = Path(cfg.lhs_file); TEST_FILE = Path(cfg.saltelli_file); METRICS_JSON = Path(cfg.surrogate_metrics); MODELS_DIR = Path(cfg.models_dir)
