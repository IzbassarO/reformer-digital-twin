"""RQ2: tube-life cost of flexible and non-steady operation as hourly quasi-steady annual histories.

An hourly history is a set of arrays (load fraction, steam-to-carbon, inlet T, inlet P, excess air, catalyst
activity and either firing or a control mode). Each hour is evaluated as a steady state of the calibrated twin
through the GP surrogates (:mod:`rdt.surrogate`); hours whose inputs fall outside the surrogate training bounds
are evaluated with the physics model (:mod:`rdt.operate`). Control modes invert the outlet-temperature or the
CH4-slip surrogate for the firing factor by vectorised bisection. Creep damage per hour is
``dD = 1 h / t_r(T_wo,max, sigma_hoop)`` with the Yeh (2021) Manaurite XM placeholder Larson-Miller curve and is
accumulated linearly (Robinson). Every 24th hour is re-evaluated with the physics model as a check.

Scope limitation: quasi-steady creep only. Start-up/shutdown thermal fatigue, creep-fatigue interaction,
carburisation and oxidation are outside this module.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import joblib
import numpy as np
import pandas as pd

from rdt import creep
from rdt import operate as op
from rdt import surrogate as sg

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "scenarios"
HOURS = 8760
T_OUT_TARGET = op.TARGET_T_OUT_K          # 1105.55 K
EXCESS_AIR_SCEN = 20.0                    # upper bound of the sampled space (Plant A derived value 21.6 % is outside it)
SURROGATE_TARGETS = ("T_wo_max_K", "T_out_K", "CH4_slip_dry_pct", "H2_net_kmol_h", "sigma_hot_MPa")


# ---------------------------------------------------------------------------
# Surrogate bundle
# ---------------------------------------------------------------------------
class TwinSurrogate:
    """GP surrogates for the quantities needed by the scenarios, with physics fallback."""

    def __init__(self):
        self.metrics = json.loads(sg.METRICS_JSON.read_text())
        self.base = op.base_case()
        self.curve = creep.LarsonMillerCurve.from_yaml()
        self.lo, self.up = op.bounds_arrays()
        self.models: Dict[str, dict] = {}
        for tgt in SURROGATE_TARGETS:
            if tgt == "sigma_hot_MPa":
                self.models[tgt] = self._sigma_model()
            else:
                self.models[tgt] = joblib.load(sg.MODELS_DIR / f"{tgt}__{self.metrics['targets'][tgt]['best']}.joblib")
        self._xeq_cache: Dict[Tuple[float, float], float] = {}
        self.n_physics_calls = 0

    def _sigma_model(self) -> dict:
        path = sg.MODELS_DIR / "sigma_hot_MPa__gp__base.joblib"
        if path.exists():
            return joblib.load(path)
        tr, _ = sg.load_data()
        m = sg._fit("gp", tr[sg.INPUTS].to_numpy(float), tr["sigma_hot_MPa"].to_numpy(float))
        b = {"pipeline": m, "features": list(sg.INPUTS), "target": "sigma_hot_MPa"}
        joblib.dump(b, path)
        return b

    # -- features -------------------------------------------------------------
    def _features(self, X: pd.DataFrame, feats: Sequence[str]) -> np.ndarray:
        X = X.copy()
        if "x_eq_CH4" in feats:
            keys = list(zip(np.round(X.steam_to_carbon.to_numpy(), 6), np.round(X.inlet_P.to_numpy(), 6)))
            new = sorted(set(k for k in keys if k not in self._xeq_cache))
            if new:
                vals = sg.equilibrium_ch4_conversion(np.array([k[0] for k in new]), np.array([k[1] for k in new]), base=self.base)
                self._xeq_cache.update(dict(zip(new, vals)))
            X["x_eq_CH4"] = [self._xeq_cache[k] for k in keys]
        if "specific_duty_MJ_per_kmol" in feats:
            dry = sum(self.base.dry_feed_kmol_h.values())
            feed_total = X.feed_per_tube_fraction * (dry + X.steam_to_carbon * self.base.carbon_kmol_h)
            Q = X.specific_firing_factor * X.feed_per_tube_fraction * self.base.Q_comb_W / self.base.row.n_tubes
            X["specific_duty_MJ_per_kmol"] = Q * 3600 / 1e6 / feed_total
        return X[list(feats)].to_numpy(float)

    def in_bounds(self, X: pd.DataFrame) -> np.ndarray:
        A = X[list(op.INPUT_NAMES)].to_numpy(float)
        return np.all((A >= self.lo - 1e-9) & (A <= self.up + 1e-9), axis=1)

    def predict(self, X: pd.DataFrame, targets: Sequence[str] = SURROGATE_TARGETS, physics_fallback: bool = True,
                force_physics: Optional[np.ndarray] = None) -> pd.DataFrame:
        """Predict the targets for a table of hourly inputs; physics model where out of bounds (or forced)."""
        out = pd.DataFrame(index=X.index)
        for tgt in targets:
            b = self.models[tgt]
            out[tgt] = b["pipeline"].predict(self._features(X, b["features"]))
        use_phys = ~self.in_bounds(X) if physics_fallback else np.zeros(len(X), bool)
        if force_physics is not None:
            use_phys |= np.asarray(force_physics, bool)
        out["evaluated_by"] = np.where(use_phys, "physics", "surrogate")
        if use_phys.any():
            phys = self.physics(X[use_phys])
            for tgt in targets:
                out.loc[use_phys, tgt] = phys[tgt].to_numpy()
        return out

    def physics(self, X: pd.DataFrame) -> pd.DataFrame:
        """Physics model on distinct input rows (cached by rounded inputs), parallel."""
        from joblib import Parallel, delayed
        A = X[list(op.INPUT_NAMES)].round(8)
        keys = [tuple(r) for r in A.to_numpy()]
        distinct = sorted(set(keys))
        self.n_physics_calls += len(distinct)
        res = Parallel(n_jobs=-1)(delayed(_physics_row)(dict(zip(op.INPUT_NAMES, k))) for k in distinct)
        lut = dict(zip(distinct, res))
        rows = [lut[k] for k in keys]
        df = pd.DataFrame(rows, index=X.index)
        return df

    # -- control ------------------------------------------------------------------
    def solve_firing(self, X: pd.DataFrame, mode: str, target: float, lo: float = 0.85, hi: float = 1.15,
                     n_iter: int = 40) -> np.ndarray:
        """Vectorised bisection on the firing factor so that T_out (mode 'hold_T_out') or the dry CH4 slip
        (mode 'hold_CH4_slip') equals ``target`` in every hour. Both are monotone in firing."""
        tgt = {"hold_T_out": "T_out_K", "hold_CH4_slip": "CH4_slip_dry_pct"}[mode]
        sign = 1.0 if mode == "hold_T_out" else -1.0        # T_out increases, slip decreases with firing
        a = np.full(len(X), lo); b = np.full(len(X), hi)
        Xw = X.copy()
        for _ in range(n_iter):
            m = 0.5 * (a + b); Xw["specific_firing_factor"] = m
            y = self.predict(Xw, (tgt,), physics_fallback=False)[tgt].to_numpy()
            too_low = sign * (y - target) < 0
            a = np.where(too_low, m, a); b = np.where(too_low, b, m)
        return 0.5 * (a + b)


def _physics_row(x: Dict[str, float]) -> Dict[str, float]:
    r = op.run_case(x, op.base_case(), creep.LarsonMillerCurve.from_yaml())
    if not r["converged"]:
        return {k: np.nan for k in SURROGATE_TARGETS}
    return {k: r[k] for k in SURROGATE_TARGETS}


# ---------------------------------------------------------------------------
# Histories
# ---------------------------------------------------------------------------
@dataclass
class History:
    name: str
    load: np.ndarray
    steam_to_carbon: np.ndarray
    inlet_T: np.ndarray
    inlet_P: np.ndarray
    excess_air: np.ndarray
    activity: np.ndarray
    control: str = "hold_T_out"                 # or "hold_CH4_slip" or "fixed_firing"
    firing: Optional[np.ndarray] = None         # used when control == "fixed_firing" or as override
    firing_override: Optional[np.ndarray] = None   # multiplicative override on the controlled firing (nan = none)
    force_physics: Optional[np.ndarray] = None
    T_wo_add_K: Optional[np.ndarray] = None     # additive perturbation of T_wo,max (hot band proxy)
    notes: str = ""

    @property
    def n(self) -> int:
        return len(self.load)

    def table(self) -> pd.DataFrame:
        return pd.DataFrame({"feed_per_tube_fraction": self.load, "steam_to_carbon": self.steam_to_carbon,
                             "inlet_T": self.inlet_T, "inlet_P": self.inlet_P, "excess_air": self.excess_air,
                             "catalyst_activity": self.activity,
                             "specific_firing_factor": self.firing if self.firing is not None else np.ones(self.n)})


def base_inputs(n: int = HOURS, excess_air: Optional[float] = None) -> Dict[str, np.ndarray]:
    excess_air = EXCESS_AIR_SCEN if excess_air is None else excess_air
    b = op.base_case().inputs()
    return {"load": np.ones(n), "steam_to_carbon": np.full(n, b["steam_to_carbon"]), "inlet_T": np.full(n, b["inlet_T"]),
            "inlet_P": np.full(n, b["inlet_P"]), "excess_air": np.full(n, excess_air), "activity": np.full(n, b["catalyst_activity"])}


def ramp_limit(x: np.ndarray, max_step: float, x0: Optional[float] = None) -> np.ndarray:
    y = np.empty_like(x); prev = x[0] if x0 is None else x0
    for i, v in enumerate(x):
        prev = prev + np.clip(v - prev, -max_step, max_step); y[i] = prev
    return y


def daily_load_profile(n: int = HOURS, day: float = 1.0, night: float = 0.70, on: int = 6, off: int = 22, ramp: float = 0.10) -> np.ndarray:
    h = np.arange(n) % 24
    target = np.where((h >= on) & (h < off), day, night)
    return ramp_limit(target, ramp, x0=night)


def ou_load_profile(n: int = HOURS, mean: float = 0.85, lo: float = 0.60, hi: float = 1.10, theta: float = 0.05,
                    sigma: float = 0.04, ramp: float = 0.10, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    x = np.empty(n); x[0] = mean
    for i in range(1, n):
        step = theta * (mean - x[i - 1]) + sigma * rng.standard_normal()
        x[i] = np.clip(x[i - 1] + np.clip(step, -ramp, ramp), lo, hi)
    return x


def make_scenarios(seed: int = 0) -> Dict[str, History]:
    S = {}
    def hist(name, control="hold_T_out", **kw):
        d = base_inputs(); d.update({k: v for k, v in kw.items() if k in d})
        h = History(name=name, load=d["load"], steam_to_carbon=d["steam_to_carbon"], inlet_T=d["inlet_T"], inlet_P=d["inlet_P"],
                    excess_air=d["excess_air"], activity=d["activity"], control=control)
        for k, v in kw.items():
            if k not in d: setattr(h, k, v)
        return h
    S["S1_steady"] = hist("S1_steady", notes="load 1.0, S/C 2.89, hold_T_out")
    S["S2_daily"] = hist("S2_daily", load=daily_load_profile(), notes="1.0 06-22 h, 0.70 at night, ramp 0.10/h")
    S["S3_renewable"] = hist("S3_renewable", load=ou_load_profile(seed=seed), notes="bounded OU load [0.60, 1.10], mean 0.85, ramp 0.10/h")
    rng = np.random.default_rng(seed)
    ov = np.full(HOURS, np.nan); starts = rng.choice(HOURS - 2, 24, replace=False)
    for s0 in starts: ov[s0:s0 + 2] = 1.10
    S["S4a_mild_overfire"] = hist("S4a_mild_overfire", firing_override=ov, notes="S1 + 24 events/yr of firing +10 % for 2 h")
    ov = np.full(HOURS, np.nan); starts = rng.choice(HOURS - 1, 6, replace=False); ov[starts] = 1.25
    S["S4b_severe_overfire"] = hist("S4b_severe_overfire", firing_override=ov, force_physics=~np.isnan(ov), notes="S1 + 6 events/yr of firing +25 % for 1 h (physics model)")
    add = np.zeros(HOURS); add[rng.choice(HOURS, 24, replace=False)] = 40.0
    S["S4c_hot_band"] = hist("S4c_hot_band", T_wo_add_K=add, notes="S1 + 40 K added to T_wo,max for 24 h/yr (flame-impingement proxy)")
    S["S6_SC2.5"] = hist("S6_SC2.5", steam_to_carbon=np.full(HOURS, 2.5), notes="S1 at S/C 2.5")
    S["S6_SC3.5"] = hist("S6_SC3.5", steam_to_carbon=np.full(HOURS, 3.5), notes="S1 at S/C 3.5")
    return S


def ageing_campaign(years: int = 4, a0: float = 0.30, a1: float = 0.10, control: str = "hold_CH4_slip",
                    constant_activity: Optional[float] = None) -> List[History]:
    out = []
    for y in range(years):
        d = base_inputs()
        t = (np.arange(HOURS) + y * HOURS) / (years * HOURS)
        act = np.full(HOURS, constant_activity) if constant_activity is not None else a0 + (a1 - a0) * t
        out.append(History(name=f"S5_{control}_{'const' if constant_activity else 'ageing'}_y{y+1}", load=d["load"], steam_to_carbon=d["steam_to_carbon"],
                           inlet_T=d["inlet_T"], inlet_P=d["inlet_P"], excess_air=d["excess_air"], activity=act, control=control,
                           notes=f"year {y+1} of {years}; activity {act[0]:.3f} -> {act[-1]:.3f}; {control}"))
    return out


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------
def evaluate(h: History, tw: TwinSurrogate, slip_target: Optional[float] = None, verify_every: int = 24) -> Dict[str, object]:
    t0 = time.perf_counter()
    X = h.table()
    if h.control == "hold_T_out":
        f = tw.solve_firing(X, "hold_T_out", T_OUT_TARGET)
    elif h.control == "hold_CH4_slip":
        assert slip_target is not None
        f = tw.solve_firing(X, "hold_CH4_slip", slip_target)
    else:
        f = X["specific_firing_factor"].to_numpy(float)
    if h.firing_override is not None:
        m = ~np.isnan(h.firing_override); f = np.where(m, f * h.firing_override, f)
    X["specific_firing_factor"] = f
    pred = tw.predict(X, force_physics=h.force_physics)
    T_wo = pred["T_wo_max_K"].to_numpy(float) + (h.T_wo_add_K if h.T_wo_add_K is not None else 0.0)
    sigma = pred["sigma_hot_MPa"].to_numpy(float)
    t_r = tw.curve.time_to_rupture(T_wo, sigma)
    dD = 1.0 / t_r
    hourly = X.copy()
    hourly["T_out_K"] = pred["T_out_K"]; hourly["CH4_slip_dry_pct"] = pred["CH4_slip_dry_pct"]; hourly["H2_net_kmol_h"] = pred["H2_net_kmol_h"]
    hourly["T_wo_max_K"] = T_wo; hourly["sigma_hot_MPa"] = sigma; hourly["t_r_h"] = t_r; hourly["dD"] = dD; hourly["D_cum"] = np.cumsum(dD)
    hourly["evaluated_by"] = pred["evaluated_by"]
    # verification against physics every verify_every hours
    idx = np.arange(0, h.n, verify_every)
    phys = tw.physics(X.iloc[idx])
    dT = phys["T_wo_max_K"].to_numpy() - pred["T_wo_max_K"].to_numpy()[idx]
    dTo = phys["T_out_K"].to_numpy() - pred["T_out_K"].to_numpy()[idx]
    # average-condition (Jensen) estimate: damage at annual-mean inputs with the same control
    Xm = pd.DataFrame([X[list(op.INPUT_NAMES)].mean()])
    if h.control in ("hold_T_out", "hold_CH4_slip"):
        Xm["specific_firing_factor"] = tw.solve_firing(Xm, h.control, T_OUT_TARGET if h.control == "hold_T_out" else slip_target)
    pm = tw.predict(Xm, physics_fallback=False)
    D_avg = h.n / float(tw.curve.time_to_rupture(float(pm["T_wo_max_K"].iloc[0]) + (float(np.mean(h.T_wo_add_K)) if h.T_wo_add_K is not None else 0.0), float(pm["sigma_hot_MPa"].iloc[0])))
    D = float(dD.sum()); H2 = float(pred["H2_net_kmol_h"].sum())
    return {"name": h.name, "hourly": hourly, "annual_H2_kmol": H2, "annual_damage": D, "years_to_D1": 1.0 / D,
            "damage_per_kmol_H2": D / H2, "D_avg_condition": D_avg, "avg_condition_error": D_avg / D,
            "T_wo_max_mean_K": float(T_wo.mean()), "T_wo_max_max_K": float(T_wo.max()), "firing_mean": float(f.mean()),
            "T_out_mean_K": float(pred["T_out_K"].mean()), "CH4_slip_mean_pct": float(pred["CH4_slip_dry_pct"].mean()),
            "n_physics_hours": int((pred["evaluated_by"] == "physics").sum()),
            "verify_max_abs_dT_wo_K": float(np.nanmax(np.abs(dT))), "verify_max_abs_dT_out_K": float(np.nanmax(np.abs(dTo))),
            "verify_n": int(len(idx)), "wall_time_s": time.perf_counter() - t0, "control": h.control, "notes": h.notes}


def steady_state_closed_form(tw: TwinSurrogate, control: str = "hold_T_out", n_hours: int = HOURS) -> Dict[str, float]:
    """Damage of a steady year from a single evaluation: D = n_hours / t_r."""
    X = pd.DataFrame([{**{k: v[0] for k, v in base_inputs(1).items() if k != "load"}, "feed_per_tube_fraction": 1.0, "specific_firing_factor": 1.0}])
    X = X.rename(columns={"activity": "catalyst_activity"})
    X["specific_firing_factor"] = tw.solve_firing(X, control, T_OUT_TARGET)
    p = tw.predict(X, physics_fallback=False)
    t_r = float(tw.curve.time_to_rupture(float(p["T_wo_max_K"].iloc[0]), float(p["sigma_hot_MPa"].iloc[0])))
    return {"t_r_h": t_r, "D_year": n_hours / t_r, "T_wo_max_K": float(p["T_wo_max_K"].iloc[0]), "firing": float(X["specific_firing_factor"].iloc[0]),
            "CH4_slip_dry_pct": float(p["CH4_slip_dry_pct"].iloc[0]), "H2_year_kmol": n_hours * float(p["H2_net_kmol_h"].iloc[0])}


def use_config(cfg) -> None:
    """Point the module at a rdt.config.RunConfig (base excess air, output names)."""
    global EXCESS_AIR_SCEN, SUMMARY_CSV, SUMMARY_META, TAG
    EXCESS_AIR_SCEN = float(cfg.excess_air_base)
    SUMMARY_CSV = Path(cfg.scenarios_summary); SUMMARY_META = Path(cfg.scenarios_meta); TAG = cfg.tag


SUMMARY_CSV = OUT_DIR / "summary_v1.csv"
SUMMARY_META = OUT_DIR / "summary_v1_meta.json"
TAG = "v1"


def hourly_path(name: str) -> Path:
    return OUT_DIR / (f"hourly_{name}.csv.gz" if TAG == "v1" else f"hourly_{TAG}_{name}.csv.gz")


def run_all_scenarios(tw: Optional["TwinSurrogate"] = None) -> Dict[str, object]:
    """Run S1-S6 (+ S5 campaigns), write the summary CSV / meta JSON and hourly files for the active config."""
    t0 = time.perf_counter()
    tw = tw or TwinSurrogate()
    ss = steady_state_closed_form(tw); slip0 = ss["CH4_slip_dry_pct"]
    results = {}; S = make_scenarios()
    for name, h in S.items():
        results[name] = evaluate(h, tw)
    camp = {}
    for label, kw in (("hold_CH4_slip_ageing", dict(control="hold_CH4_slip")), ("hold_T_out_ageing", dict(control="hold_T_out")),
                      ("hold_T_out_const0.20", dict(control="hold_T_out", constant_activity=0.20))):
        yrs = []
        for h in ageing_campaign(**kw):
            r = evaluate(h, tw, slip_target=slip0); results[h.name] = r; hh = r["hourly"]
            yrs.append({"year": len(yrs) + 1, "activity_start": float(hh.catalyst_activity.iloc[0]), "activity_end": float(hh.catalyst_activity.iloc[-1]),
                        "firing_start": float(hh.specific_firing_factor.iloc[0]), "firing_end": float(hh.specific_firing_factor.iloc[-1]),
                        "T_out_start": float(hh.T_out_K.iloc[0]), "T_out_end": float(hh.T_out_K.iloc[-1]), "T_wo_start": float(hh.T_wo_max_K.iloc[0]), "T_wo_end": float(hh.T_wo_max_K.iloc[-1]),
                        "slip_start": float(hh.CH4_slip_dry_pct.iloc[0]), "slip_end": float(hh.CH4_slip_dry_pct.iloc[-1]), "H2_kmol": r["annual_H2_kmol"], "D": r["annual_damage"]})
        camp[label] = yrs
    D1, H1 = results["S1_steady"]["annual_damage"], results["S1_steady"]["annual_H2_kmol"]
    rows = []
    for name, r in results.items():
        rows.append({"scenario": name, "control": r["control"], "annual_H2_kmol_per_tube": r["annual_H2_kmol"], "annual_damage_D": r["annual_damage"],
                     "years_to_D1_yeh_placeholder": r["years_to_D1"], "life_consumption_per_kmol_H2_rel_S1": (r["annual_damage"] / r["annual_H2_kmol"]) / (D1 / H1),
                     "D_avg_condition": r["D_avg_condition"], "avg_condition_error": r["avg_condition_error"], "T_wo_max_mean_K": r["T_wo_max_mean_K"], "T_wo_max_max_K": r["T_wo_max_max_K"],
                     "firing_mean": r["firing_mean"], "T_out_mean_K": r["T_out_mean_K"], "CH4_slip_mean_pct": r["CH4_slip_mean_pct"], "n_physics_hours": r["n_physics_hours"],
                     "verify_max_abs_dT_wo_K": r["verify_max_abs_dT_wo_K"], "verify_max_abs_dT_out_K": r["verify_max_abs_dT_out_K"], "notes": r["notes"]})
    summ = pd.DataFrame(rows).set_index("scenario")
    dD_base = results["S1_steady"]["hourly"].dD.iloc[0]; ev = {}
    for name in ("S4a_mild_overfire", "S4b_severe_overfire", "S4c_hot_band"):
        hh = results[name]["hourly"]; m = hh.dD > dD_base * 1.0001
        ev[name] = {"event_hours": int(m.sum()), "dD_event_mean": float(hh.dD[m].mean()), "dD_base": float(dD_base), "ratio_event_to_base": float(hh.dD[m].mean() / dD_base),
                    "share_of_annual_damage_pct": float(100 * (hh.dD[m].sum() - m.sum() * dD_base) / hh.dD.sum()), "T_wo_event_mean_K": float(hh.T_wo_max_K[m].mean())}
        summ.loc[name, "event_hours"] = ev[name]["event_hours"]; summ.loc[name, "dD_event_to_base_ratio"] = ev[name]["ratio_event_to_base"]
    OUT_DIR.mkdir(exist_ok=True)
    summ.to_csv(SUMMARY_CSV)
    for name, r in results.items():
        r["hourly"].to_csv(hourly_path(name), index=False, compression="gzip")
    meta = {"created": time.strftime("%Y-%m-%d %H:%M:%S"), "wall_time_s": time.perf_counter() - t0, "excess_air_pct": EXCESS_AIR_SCEN, "T_out_target_K": T_OUT_TARGET,
            "slip_target_pct": slip0, "steady_closed_form": ss, "S4_events": ev, "S5_campaigns": camp, "n_physics_calls": tw.n_physics_calls, "tag": TAG,
            "life_curve": "Yeh 2021 Manaurite XM minimum curve, PLACEHOLDER", "scope": "quasi-steady creep only; start-up/shutdown thermal fatigue excluded"}
    json.dump(meta, open(SUMMARY_META, "w"), indent=2, default=float)
    return {"summary": summ, "meta": meta, "results": results}
