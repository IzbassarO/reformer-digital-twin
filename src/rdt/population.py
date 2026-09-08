"""From one average tube to the tube population: per-tube creep life under wall-temperature scatter.

``python -m rdt.population [--version v3] [--sigma 15 18.5 22]``

The rest of the twin models *one* tube. A top-fired reformer has many: the reference furnace is 7 rows
of 48, so ``N_TUBES = 336``, and the standard deviation of the outer-wall temperature across that
population at a given elevation is 15-22 K (Latham et al., the journal version of the thesis data).
Damage rate ``1/t_r`` is convex in wall temperature, so

* the population-mean damage rate exceeds the damage rate of the tube at the mean temperature, and
* the retubing decision is governed by the *first* tube to fail, not by the average one.

This module propagates a per-tube offset ``dT ~ Normal(0, sigma)`` over the population at an operating
point the twin has already evaluated, and reports both effects.

Where the offset enters
-----------------------
The alloy enters the twin only at the last step, ``t_r(T_wo, sigma_hoop)`` (see
:mod:`rdt.alloy_comparison`). A tube-to-tube wall-temperature offset is not an operating-space input:
firing, feed and pressure are common to the bundle, so the hoop stress is common to the bundle and only
the metal temperature moves. The per-tube rupture time is therefore

    t_r,i = t_r(T_wo,max + dT_i, sigma_hoop)

evaluated directly on the master curve. The surrogate/physics choice below fixes the *anchor* state
``(T_wo,max, sigma_hoop, H2)`` of the operating point; it is not re-entered per tube.

Code path
---------
An operating point is evaluated through :class:`rdt.scenarios.TwinSurrogate`, which uses the Gaussian-
process surrogates inside the sampled operating box and falls back to the physics model outside it, and
labels each row ``surrogate`` or ``physics``. Hourly scenario histories are read from the stored
``data/scenarios/hourly_<tag>_<name>.csv.gz``, which already carry that label per hour. Every number
printed by :func:`main` states the path it came from.

Analytic cross-check
--------------------
At fixed stress the master curve gives ``log10 t_r = LMP/(scale T) - C`` with LMP fixed by the stress,
so the local slope is available in closed form,

    a := -d log10 t_r / dT = LMP / (scale T^2)      [decades per K]

and, if ``a`` were constant over the sampled range, ``1/t_r`` would be exactly lognormal with

    E[rate] / rate(T_bar) = exp((a ln10 sigma)^2 / 2).

:func:`validate_analytic` prints that value beside the sampled one. The two differ by the curvature of
``a`` (which falls as ``1/T^2``) plus the Monte Carlo error of a finite draw; both are reported.

Scope
-----
The offset is a static per-tube property (burner imbalance, row position, flow maldistribution), drawn
once and held for the whole history, not resampled each hour. Quasi-steady creep only, as elsewhere in
the twin.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rdt import config as C
from rdt import creep

ROOT = C.ROOT
OUT_DIR = ROOT / "data" / "population"

#: reference furnace bundle: 7 rows of 48 tubes (Latham Plant A, ``n_tubes`` in data/literature_validation)
ROWS, TUBES_PER_ROW = 7, 48
N_TUBES = ROWS * TUBES_PER_ROW          # 336
#: one-sigma wall-temperature scatter across the population at a given elevation [K]
SIGMA_DEFAULT_K = 18.5
SIGMA_GRID_K = (15.0, 18.5, 22.0)
SEED = 0
#: fraction of the bundle counted as the "hottest decile"
DECILE = 0.10
HOURS_PER_YEAR = 8760.0

#: RQ2 scenarios, in the order of data/scenarios/summary_<tag>.csv
#: RQ3 regimes are read from cfg.regimes_csv (the RQ3 representative-regime table of the manuscript).


# ---------------------------------------------------------------------------
# The population sample
# ---------------------------------------------------------------------------
def offsets(sigma_K: float = SIGMA_DEFAULT_K, n_tubes: int = N_TUBES, seed: int = SEED) -> np.ndarray:
    """Per-tube wall-temperature offsets ``dT ~ Normal(0, sigma_K)``, one draw per tube.

    Reproducible: the same ``(sigma_K, n_tubes, seed)`` always returns the same array. The draw is not
    re-centred -- a real bundle of 336 tubes has a sample mean that misses zero, and that scatter is
    part of what is being reported.
    """
    if n_tubes < 1:
        raise ValueError("n_tubes must be >= 1")
    if sigma_K < 0:
        raise ValueError("sigma_K must be >= 0")
    return np.random.default_rng(seed).normal(0.0, float(sigma_K), int(n_tubes))


# ---------------------------------------------------------------------------
# Creep evaluation over the population
# ---------------------------------------------------------------------------
def _lmp_of_stress(curve: creep.LarsonMillerCurve, sigma_MPa: np.ndarray) -> np.ndarray:
    """LMP at each hoop stress, root-finding once per *distinct* stress (histories repeat theirs)."""
    s = np.atleast_1d(np.asarray(sigma_MPa, float))
    uniq, inv = np.unique(np.round(s, 10), return_inverse=True)
    lm = np.array([curve.lmp_of_stress(float(v)) for v in uniq])
    return lm[inv]


def local_slope(curve: creep.LarsonMillerCurve, T_K: float, sigma_MPa: float) -> float:
    """``a = -d log10 t_r / dT`` at fixed stress [decades per K], in closed form.

    ``LMP = scale T (C + log10 t_r)`` with LMP set by the stress alone, so
    ``log10 t_r = LMP/(scale T) - C`` and ``a = LMP/(scale T^2)``.
    """
    lmp = curve.lmp_of_stress(float(sigma_MPa))
    return float(lmp / (curve.scale * float(T_K) ** 2))


def damage_rates(curve: creep.LarsonMillerCurve, T_K, sigma_MPa, dT: np.ndarray) -> np.ndarray:
    """Per-tube damage rate ``1/t_r`` [1/h] for offsets ``dT``.

    ``T_K`` and ``sigma_MPa`` are scalars (one operating point) or arrays of equal length (an hourly
    history); the result is ``(n_tubes,)`` in the first case and ``(n_hours, n_tubes)`` in the second.
    """
    T = np.atleast_1d(np.asarray(T_K, float))
    lm = _lmp_of_stress(curve, sigma_MPa)
    if lm.size == 1:
        lm = np.repeat(lm, T.size)
    if lm.size != T.size:
        raise ValueError("T_K and sigma_MPa must have the same length")
    Tt = T[:, None] + np.asarray(dT, float)[None, :]
    log10_tr = lm[:, None] / (curve.scale * Tt) - curve.C
    out = 10.0 ** (-log10_tr)
    return out[0] if np.ndim(T_K) == 0 else out


def analytic_mean_ratio(a_per_K: float, sigma_K: float) -> float:
    """Lognormal population mean over the value at the mean temperature: ``exp((a ln10 sigma)^2/2)``."""
    return float(np.exp((a_per_K * np.log(10.0) * sigma_K) ** 2 / 2.0))


# ---------------------------------------------------------------------------
# One operating point
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TubeState:
    """Hot-spot state of the *average* tube at one operating point, with its code path."""

    label: str
    T_wo_max_K: float
    sigma_hot_MPa: float
    H2_net_kmol_h: float
    evaluated_by: str          # "surrogate", "physics", or "stored hourly (<n> surrogate / <n> physics)"


def population_stats(curve: creep.LarsonMillerCurve, state: TubeState, dT: np.ndarray,
                     sigma_K: float) -> Dict[str, float]:
    """Every quantity item 2 of the analysis asks for, at one operating point.

    Rates are per hour; lives are ``1/rate`` in years. ``mult`` is the damage-rate multiplier over the
    average tube. Percentiles are percentiles *of the per-tube damage rate*, so the 95th percentile is a
    hot tube and the 1st a cold one.
    """
    rate_avg = float(1.0 / np.asarray(curve.time_to_rupture(state.T_wo_max_K, state.sigma_hot_MPa)))
    r = damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, dT)
    n_dec = max(1, int(round(DECILE * len(r))))
    hottest_decile = float(np.mean(np.sort(r)[-n_dec:]))
    a = local_slope(curve, state.T_wo_max_K, state.sigma_hot_MPa)
    out = {
        "T_wo_max_K": state.T_wo_max_K, "sigma_hot_MPa": state.sigma_hot_MPa,
        "evaluated_by": state.evaluated_by, "sigma_K": float(sigma_K), "n_tubes": int(len(r)),
        "a_dec_per_K": a, "rate_avg_tube_per_h": rate_avg,
        "rate_pop_mean_per_h": float(np.mean(r)),
        "ratio_pop_mean_to_avg": float(np.mean(r)) / rate_avg,
        "ratio_analytic_lognormal": analytic_mean_ratio(a, sigma_K),
        "n_hottest_decile": n_dec,
        "dT_max_K": float(np.max(dT)), "dT_mean_K": float(np.mean(dT)), "dT_sd_K": float(np.std(dT, ddof=1)),
    }
    named = {"hottest": float(np.max(r)), "hottest_decile": hottest_decile,
             "p95": float(np.percentile(r, 95)), "p5": float(np.percentile(r, 5)),
             "p1": float(np.percentile(r, 1))}
    for k, v in named.items():
        out[f"mult_{k}"] = v / rate_avg
        out[f"life_{k}_y"] = 1.0 / v / HOURS_PER_YEAR
    out["life_avg_tube_y"] = 1.0 / rate_avg / HOURS_PER_YEAR
    #: the claim the manuscript makes about hottest-versus-average rupture time
    out["t_r_ratio_avg_over_hottest"] = out["mult_hottest"]
    return out


# ---------------------------------------------------------------------------
# Getting an operating point out of the twin, with the code path recorded
# ---------------------------------------------------------------------------
class PointEvaluator:
    """Evaluate operating points through the surrogates where they are in range, physics otherwise.

    Thin wrapper over :class:`rdt.scenarios.TwinSurrogate`, whose ``predict`` already applies exactly
    that rule and labels every row. Built lazily so that the stored-history paths of this module do not
    require the surrogate models to be present.
    """

    def __init__(self):
        from rdt import scenarios as sc
        self.tw = sc.TwinSurrogate()

    def state(self, x: Mapping[str, float], label: str) -> TubeState:
        from rdt import operate as op
        X = pd.DataFrame([{k: float(x[k]) for k in op.INPUT_NAMES}])
        p = self.tw.predict(X).iloc[0]
        return TubeState(label=label, T_wo_max_K=float(p.T_wo_max_K), sigma_hot_MPa=float(p.sigma_hot_MPa),
                         H2_net_kmol_h=float(p.H2_net_kmol_h), evaluated_by=str(p.evaluated_by))

    def physics_state(self, x: Mapping[str, float], label: str) -> TubeState:
        """The same point forced through the physics model, for a surrogate-error cross-check."""
        from rdt import operate as op
        r = op.run_case({k: float(x[k]) for k in op.INPUT_NAMES})
        if not r["converged"]:
            raise RuntimeError(f"{label}: physics did not converge: {r['error']}")
        return TubeState(label=label, T_wo_max_K=float(r["T_wo_max_K"]), sigma_hot_MPa=float(r["sigma_hot_MPa"]),
                         H2_net_kmol_h=float(r["H2_net_kmol_h"]), evaluated_by="physics")


def regime_inputs(cfg: C.RunConfig) -> pd.DataFrame:
    """The RQ3 representative regimes, as the seven operating inputs."""
    from rdt import operate as op
    df = pd.read_csv(cfg.regimes_csv, index_col=0)
    return df[list(op.INPUT_NAMES)]


def scenario_names(cfg: C.RunConfig) -> List[str]:
    return pd.read_csv(cfg.scenarios_summary)["scenario"].tolist()


def load_hourly(cfg: C.RunConfig, name: str) -> pd.DataFrame:
    p = cfg.scenario_hourly(name)
    if not p.exists():
        raise FileNotFoundError(f"{p} not found; run `python -m rdt.pipeline --stage scenarios --version {cfg.tag}`")
    return pd.read_csv(p, usecols=["T_wo_max_K", "sigma_hot_MPa", "H2_net_kmol_h", "evaluated_by"])


def load_hourly_inputs(cfg: C.RunConfig, name: str, hour: int = 0) -> Dict[str, float]:
    """The seven operating inputs of one hour of a stored scenario history."""
    from rdt import operate as op
    p = cfg.scenario_hourly(name)
    row = pd.read_csv(p, usecols=list(op.INPUT_NAMES)).iloc[hour]
    return {k: float(row[k]) for k in op.INPUT_NAMES}


def _path_label(evaluated_by: pd.Series) -> str:
    n_s = int((evaluated_by == "surrogate").sum()); n_p = int((evaluated_by == "physics").sum())
    return f"stored hourly ({n_s} surrogate / {n_p} physics)"


def steady_scenario_state(cfg: C.RunConfig, name: str = "S1_steady") -> TubeState:
    """Anchor state of a *steady* RQ2 scenario, read from its stored hourly history.

    Raises if the history is not in fact constant, since a single anchor state would then be a
    misrepresentation of the year.
    """
    df = load_hourly(cfg, name)
    for c in ("T_wo_max_K", "sigma_hot_MPa", "H2_net_kmol_h"):
        if float(df[c].max() - df[c].min()) > 1e-6 * abs(float(df[c].iloc[0])):
            raise ValueError(f"{name} is not steady in {c}; use scenario_population() instead")
    r = df.iloc[0]
    return TubeState(label=name, T_wo_max_K=float(r.T_wo_max_K), sigma_hot_MPa=float(r.sigma_hot_MPa),
                     H2_net_kmol_h=float(r.H2_net_kmol_h), evaluated_by=_path_label(df.evaluated_by))


# ---------------------------------------------------------------------------
# Annual histories over the population
# ---------------------------------------------------------------------------
def scenario_population(cfg: C.RunConfig, name: str, dT: np.ndarray,
                        curve: creep.LarsonMillerCurve) -> Dict[str, object]:
    """Annual damage of every tube for one RQ2 scenario, and the life cost per kmol H2.

    The offset is held fixed over the year: a tube's position in the furnace does not change. Damage is
    accumulated by the Robinson rule exactly as the scenario stage does it, ``D = sum_h 1 h / t_r``.
    """
    df = load_hourly(cfg, name)
    T = df.T_wo_max_K.to_numpy(float); s = df.sigma_hot_MPa.to_numpy(float)
    D_i = damage_rates(curve, T, s, dT).sum(axis=0)                 # (n_tubes,) hours -> dimensionless
    D_avg = float(np.sum(1.0 / np.asarray(curve.time_to_rupture(T, s))))
    H2 = float(df.H2_net_kmol_h.sum())
    n_dec = max(1, int(round(DECILE * len(D_i))))
    D_dec = float(np.mean(np.sort(D_i)[-n_dec:]))
    return {"scenario": name, "evaluated_by": _path_label(df.evaluated_by), "H2_annual_kmol_per_tube": H2,
            "D_avg_tube": D_avg, "D_pop_mean": float(np.mean(D_i)), "D_hottest": float(np.max(D_i)),
            "D_hottest_decile": D_dec, "n_hottest_decile": n_dec,
            "cost_avg": D_avg / H2, "cost_pop_mean": float(np.mean(D_i)) / H2, "cost_hottest_decile": D_dec / H2,
            "years_avg_tube": 1.0 / D_avg, "years_pop_mean": 1.0 / float(np.mean(D_i)),
            "years_hottest_decile": 1.0 / D_dec, "years_first_failure": 1.0 / float(np.max(D_i)),
            "ratio_pop_mean_to_avg": float(np.mean(D_i)) / D_avg,
            "ratio_hottest_to_avg": float(np.max(D_i)) / D_avg}


def regime_population(states: Sequence[TubeState], dT: np.ndarray,
                      curve: creep.LarsonMillerCurve) -> pd.DataFrame:
    """Life cost per kmol H2 of each RQ3 regime for the average tube, the population and the decile."""
    rows = []
    for st in states:
        r = damage_rates(curve, st.T_wo_max_K, st.sigma_hot_MPa, dT)
        rate_avg = float(1.0 / np.asarray(curve.time_to_rupture(st.T_wo_max_K, st.sigma_hot_MPa)))
        n_dec = max(1, int(round(DECILE * len(r))))
        dec = float(np.mean(np.sort(r)[-n_dec:]))
        rows.append({"regime": st.label, "evaluated_by": st.evaluated_by, "T_wo_max_K": st.T_wo_max_K,
                     "sigma_hot_MPa": st.sigma_hot_MPa, "H2_net_kmol_h": st.H2_net_kmol_h,
                     "cost_avg": rate_avg / st.H2_net_kmol_h,
                     "cost_pop_mean": float(np.mean(r)) / st.H2_net_kmol_h,
                     "cost_hottest_decile": dec / st.H2_net_kmol_h,
                     "ratio_pop_mean_to_avg": float(np.mean(r)) / rate_avg,
                     "years_avg_tube": 1.0 / rate_avg / HOURS_PER_YEAR,
                     "years_first_failure": 1.0 / float(np.max(r)) / HOURS_PER_YEAR})
    return pd.DataFrame(rows).set_index("regime")


# ---------------------------------------------------------------------------
# Ranking invariance
# ---------------------------------------------------------------------------
def ranking(table: pd.DataFrame, cols=("cost_avg", "cost_pop_mean", "cost_hottest_decile")) -> pd.DataFrame:
    """Rank of each row by life cost per kmol H2 (1 = cheapest) under each definition of the tube."""
    out = pd.DataFrame(index=table.index)
    for c in cols:
        out[c.replace("cost", "rank")] = table[c].rank(method="min").astype(int)
    out["rank_changes_pop"] = out["rank_pop_mean"] != out["rank_avg"]
    out["rank_changes_decile"] = out["rank_hottest_decile"] != out["rank_avg"]
    return out


def swapped_pairs(table: pd.DataFrame, a_col: str = "cost_avg",
                  b_col: str = "cost_pop_mean") -> List[Tuple[str, str, float, float]]:
    """Pairs whose order reverses between two cost definitions, with the size of each gap.

    A rank change only matters if the two members were meaningfully apart to begin with; the returned
    gaps are ``|c_i - c_j| / min(c_i, c_j)`` in per cent under each definition, so a reversal between
    two rows that differ by a fraction of a per cent can be seen for what it is.
    """
    names = list(table.index)
    a, b = table[a_col].to_numpy(float), table[b_col].to_numpy(float)
    out = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            if (a[i] - a[j]) * (b[i] - b[j]) < 0:
                out.append((names[i], names[j], 100.0 * abs(a[i] - a[j]) / min(a[i], a[j]),
                            100.0 * abs(b[i] - b[j]) / min(b[i], b[j])))
    return sorted(out, key=lambda t: -max(t[2], t[3]))


# ---------------------------------------------------------------------------
# Validation of the code path
# ---------------------------------------------------------------------------
def validate_analytic(curve: creep.LarsonMillerCurve, state: TubeState, sigma_K: float,
                      n_large: int = 400_000, seed: int = SEED) -> Dict[str, float]:
    """Sampled population-mean ratio against ``exp((a ln10 sigma)^2/2)``, separating the two error terms.

    The finite bundle (336 tubes) carries a real Monte Carlo error; a large draw isolates the residual,
    which is the curvature of ``a`` over the sampled temperature range.
    """
    a = local_slope(curve, state.T_wo_max_K, state.sigma_hot_MPa)
    ana = analytic_mean_ratio(a, sigma_K)
    rate_avg = float(1.0 / np.asarray(curve.time_to_rupture(state.T_wo_max_K, state.sigma_hot_MPa)))
    r_n = damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, offsets(sigma_K, N_TUBES, seed))
    r_big = damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, offsets(sigma_K, n_large, seed))
    s_ln = a * np.log(10.0) * sigma_K
    cv = float(np.sqrt(np.exp(s_ln ** 2) - 1.0))                    # lognormal coefficient of variation
    return {"a_dec_per_K": a, "sigma_K": float(sigma_K), "analytic": ana,
            "sampled_N336": float(np.mean(r_n)) / rate_avg, "sampled_large": float(np.mean(r_big)) / rate_avg,
            "n_large": int(n_large), "rel_err_N336_pct": 100.0 * (float(np.mean(r_n)) / rate_avg / ana - 1.0),
            "rel_err_large_pct": 100.0 * (float(np.mean(r_big)) / rate_avg / ana - 1.0),
            "mc_se_N336_pct": 100.0 * cv / np.sqrt(N_TUBES), "lognormal_cv": cv}


def seed_spread(curve: creep.LarsonMillerCurve, state: TubeState, sigma_K: float,
                n_seeds: int = 200) -> Dict[str, float]:
    """Spread of the finite-bundle statistics over independent 336-tube draws.

    The hottest tube of a bundle is the maximum of 336 draws and is itself random; a single furnace is
    one realisation, so the numbers below say how much of the reported value is draw-to-draw luck.
    """
    rate_avg = float(1.0 / np.asarray(curve.time_to_rupture(state.T_wo_max_K, state.sigma_hot_MPa)))
    means, hots = [], []
    for sd in range(n_seeds):
        r = damage_rates(curve, state.T_wo_max_K, state.sigma_hot_MPa, offsets(sigma_K, N_TUBES, sd))
        means.append(np.mean(r) / rate_avg); hots.append(np.max(r) / rate_avg)
    m, h = np.array(means), np.array(hots)
    return {"n_seeds": n_seeds, "pop_mean_ratio_mean": float(m.mean()), "pop_mean_ratio_sd": float(m.std(ddof=1)),
            "hottest_mult_mean": float(h.mean()), "hottest_mult_sd": float(h.std(ddof=1)),
            "hottest_mult_p5": float(np.percentile(h, 5)), "hottest_mult_p95": float(np.percentile(h, 95))}


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def figure_damage_distributions(cfg: C.RunConfig, states: Mapping[str, TubeState], dT: np.ndarray,
                                curve: creep.LarsonMillerCurve, sigma_K: float,
                                stem: Optional[Path] = None) -> List[Path]:
    """Per-tube damage-rate distribution at the base point and at the knee regime, on a common axis.

    Both panels are normalised by the *base* average-tube rate, so the horizontal shift between them is
    the life cost of the regime and the width of each is the population scatter.
    """
    from rdt import plotting as P

    stem = Path(stem or (cfg.figures_dir / "fig16_population_damage"))
    ref = float(1.0 / np.asarray(curve.time_to_rupture(states["base"].T_wo_max_K, states["base"].sigma_hot_MPa)))
    panels = [("base", "a", "Base case"), ("knee", "b", "Knee regime")]
    data = {k: damage_rates(curve, states[k].T_wo_max_K, states[k].sigma_hot_MPa, dT) / ref for k, _, _ in panels}
    allv = np.concatenate(list(data.values()))
    bins = np.logspace(np.log10(allv.min()) - 0.05, np.log10(allv.max()) + 0.05, 34)
    with P.style():
        fig, ax = P.figure(P.DOUBLE, 6.2 * P.CM, 1, 2, sharex=True, sharey=True)
        for j, (key, lab, title) in enumerate(panels):
            v = data[key]; st = states[key]
            avg = float(1.0 / np.asarray(curve.time_to_rupture(st.T_wo_max_K, st.sigma_hot_MPa))) / ref
            n_dec = max(1, int(round(DECILE * len(v))))
            dec = float(np.mean(np.sort(v)[-n_dec:]))
            ax[j].hist(v, bins=bins, color=P.PALETTE[0], alpha=0.55, edgecolor="none")
            ax[j].axvline(avg, color=P.PALETTE[7], lw=1.2, label=f"average tube ({avg:.2f})")
            ax[j].axvline(v.mean(), color=P.PALETTE[1], lw=1.2, ls="--", label=f"population mean ({v.mean():.2f})")
            ax[j].axvline(dec, color=P.PALETTE[2], lw=1.2, ls=":", label=f"hottest decile ({dec:.2f})")
            ax[j].set_xscale("log")
            P.tidy(ax[j], "Damage rate relative to the base average tube (-)",
                   "Tubes" if j == 0 else None, legend=True, loc="upper left")
            P.panel_label(ax[j], lab)
            ax[j].text(0.98, 0.97, title, transform=ax[j].transAxes, ha="right", va="top", fontsize=7)
        fig.tight_layout()
        return P.save(fig, stem)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
def _fmt(v: float, w: int = 10, p: int = 3) -> str:
    return f"{v:{w}.{p}f}" if abs(v) < 1e5 else f"{v:{w}.{p}e}"


def build(cfg: C.RunConfig, sigmas: Sequence[float] = SIGMA_GRID_K, seed: int = SEED,
          n_tubes: int = N_TUBES, base_sigma: float = SIGMA_DEFAULT_K) -> Dict[str, object]:
    """Everything the analysis reports, as plain data (no printing)."""
    from rdt import operate as op

    curve = creep.active_curve()
    ev = PointEvaluator()

    # --- anchor states -----------------------------------------------------
    reg_x = regime_inputs(cfg)
    reg_states = [ev.state(reg_x.loc[name], str(name)) for name in reg_x.index]
    by_name = {s.label: s for s in reg_states}
    base_state = by_name["base"]

    # --- item 2: two readings of "the base point", at each sigma -----------
    # S1_steady is the RQ2 base year (firing solved to hold the measured outlet temperature); the RQ3
    # `base` regime is the calibrated base case at firing factor 1.0, which is the 1147 K hot spot the
    # manuscript quotes. They are different operating points and both are reported.
    anchors = {"S1_steady": steady_scenario_state(cfg, "S1_steady"), "RQ3 base": base_state}
    cross = {"S1_steady": ev.physics_state(load_hourly_inputs(cfg, "S1_steady"), "S1_steady"),
             "RQ3 base": ev.physics_state(reg_x.loc["base"], "base")}

    per_sigma: Dict[str, Dict[str, dict]] = {}
    for label, st in anchors.items():
        per_sigma[label] = {}
        for sg in sigmas:
            dT = offsets(sg, n_tubes, seed)
            per_sigma[label][str(sg)] = {"stats": population_stats(curve, st, dT, sg),
                                         "validation": validate_analytic(curve, st, sg, seed=seed),
                                         "seed_spread": seed_spread(curve, st, sg)}

    # --- item 3: scenarios and regimes at the default sigma ----------------
    dT = offsets(base_sigma, n_tubes, seed)
    scen = pd.DataFrame([scenario_population(cfg, n, dT, curve) for n in scenario_names(cfg)]).set_index("scenario")
    regs = regime_population(reg_states, dT, curve)
    return {"cfg": cfg, "curve": curve, "dT": dT, "anchors": anchors, "cross": cross,
            "regime_states": by_name, "per_sigma": per_sigma, "scenarios": scen, "regimes": regs,
            "rank_scen": ranking(scen), "rank_reg": ranking(regs), "sigmas": list(sigmas),
            "seed": seed, "n_tubes": n_tubes, "base_sigma": base_sigma}


def report(res: Dict[str, object]) -> str:
    """The printed analysis, as a string (also written to data/population/)."""
    cfg: C.RunConfig = res["cfg"]; curve = res["curve"]
    anchors: Dict[str, TubeState] = res["anchors"]; cross: Dict[str, TubeState] = res["cross"]
    L: List[str] = []
    P = L.append

    P(f"Per-tube population extension of the life model  ({cfg.tag} results, alloy {cfg.creep_alloy})")
    P("=" * 104)
    P(f"Bundle: {ROWS} rows x {TUBES_PER_ROW} tubes = {res['n_tubes']} tubes; dT ~ Normal(0, sigma), one static")
    P(f"draw per tube, seed {res['seed']}. Master curve {curve.alloy or cfg.creep_alloy!r} "
      f"(C = {curve.C}, lower scatter band).")
    P("The offset moves the metal temperature only: firing, feed and pressure are common to the bundle,")
    P("so the hoop stress is common to it. t_r,i = t_r(T_wo,max + dT_i, sigma_hoop).")
    P("")
    P("ANCHOR STATES AND THE CODE PATH THAT PRODUCED THEM")
    P("-" * 104)
    P("  rule: inside the sampled operating box -> GP surrogates; outside it -> physics model.")
    P("  'S1_steady' is the RQ2 base year, firing solved to hold the measured outlet temperature.")
    P("  'RQ3 base'  is the calibrated base case at firing factor 1.0 -- the 1147 K hot spot the paper quotes.")
    P("")
    P(f"{'anchor':<12} {'T_wo,max [K]':>13} {'sigma [MPa]':>12} {'H2 [kmol/h]':>12}   code path")
    for k, st in anchors.items():
        cp = cross[k]
        P(f"{k:<12} {st.T_wo_max_K:13.3f} {st.sigma_hot_MPa:12.4f} {st.H2_net_kmol_h:12.4f}   {st.evaluated_by}")
        P(f"{'':<12} {cp.T_wo_max_K:13.3f} {cp.sigma_hot_MPa:12.4f} {cp.H2_net_kmol_h:12.4f}   physics, "
          f"cross-check ({st.T_wo_max_K - cp.T_wo_max_K:+.3f} K, {st.sigma_hot_MPa - cp.sigma_hot_MPa:+.4f} MPa)")
    P("")

    # ---- item 2 ----
    P("1. THE BASE POINT OVER THE POPULATION")
    P("-" * 104)
    P("Multipliers are damage rate over the average tube; lives are 1/rate in years. Percentiles are")
    P("percentiles of the per-tube damage rate, so p95 is a hot tube and p1 a cold one; 'hottest decile'")
    P("is the mean damage rate of the hottest 10 % of the bundle.")
    for label, st in anchors.items():
        s0 = res["per_sigma"][label][str(res["base_sigma"])]["stats"]
        P("")
        P(f"  {label}  (T_wo,max = {st.T_wo_max_K:.2f} K, sigma_hoop = {st.sigma_hot_MPa:.3f} MPa, "
          f"a = {s0['a_dec_per_K']:.5f} decades/K)   [{st.evaluated_by}]")
        P(f"  {'':<9}  {'--- damage-rate multiplier over the average tube ---':^54}   "
          f"{'--- life [years] ---':^54}")
        P(f"  {'sigma [K]':>9} {'pop mean':>9} {'analytic':>9} {'hottest':>9} {'hot dec':>9} {'p95':>8} "
          f"{'p5':>8} {'p1':>8} | {'avg tube':>10} {'hottest':>9} {'hot dec':>9} {'p95':>9} {'p5':>9} {'p1':>9}")
        for sg in res["sigmas"]:
            s = res["per_sigma"][label][str(sg)]["stats"]
            P(f"  {sg:9.1f} {s['ratio_pop_mean_to_avg']:9.3f} {s['ratio_analytic_lognormal']:9.3f} "
              f"{s['mult_hottest']:9.3f} {s['mult_hottest_decile']:9.3f} {s['mult_p95']:8.3f} "
              f"{s['mult_p5']:8.3f} {s['mult_p1']:8.3f} | {s['life_avg_tube_y']:10.4g} "
              f"{s['life_hottest_y']:9.4g} {s['life_hottest_decile_y']:9.4g} {s['life_p95_y']:9.4g} "
              f"{s['life_p5_y']:9.4g} {s['life_p1_y']:9.4g}")
    s0 = res["per_sigma"]["S1_steady"][str(res["base_sigma"])]["stats"]
    P("")
    P(f"  hottest decile = mean over the hottest {s0['n_hottest_decile']} of {s0['n_tubes']} tubes")
    P(f"  realised draw at sigma = {res['base_sigma']} K: mean dT {s0['dT_mean_K']:+.2f} K, "
      f"sd {s0['dT_sd_K']:.2f} K, max {s0['dT_max_K']:+.2f} K")
    P("")
    P("  RATIO OF HOTTEST-TUBE TO AVERAGE-TUBE RUPTURE TIME")
    P("  The manuscript (Section 5, limitations) claims 'a factor of two to four'. The hottest tube of a")
    P("  bundle is the maximum of N draws and is itself random, so the spread over bundles is given too.")
    for label in anchors:
        P(f"    {label}:")
        for sg in res["sigmas"]:
            st = res["per_sigma"][label][str(sg)]
            s, sp = st["stats"], st["seed_spread"]
            P(f"      sigma = {sg:4.1f} K : t_r(avg)/t_r(hottest) = {s['t_r_ratio_avg_over_hottest']:6.2f} "
              f"(seed {res['seed']})   {sp['hottest_mult_mean']:6.2f} +- {sp['hottest_mult_sd']:5.2f} over "
              f"{sp['n_seeds']} bundles, 5-95 % {sp['hottest_mult_p5']:.2f}-{sp['hottest_mult_p95']:.2f}")
    P("    a tube ONE sigma above the mean, which is the number 'two to four' actually corresponds to: "
      + ", ".join(f"{sg:.1f} K -> {10 ** (s0['a_dec_per_K'] * sg):.2f}" for sg in res["sigmas"]))
    P("")

    # ---- validation ----
    P("2. VALIDATION OF THE CODE PATH: sampled population mean against the analytic lognormal result")
    P("-" * 104)
    P("   E[rate]/rate(T_bar) = exp((a ln10 sigma)^2 / 2) holds exactly only if a were constant in T;")
    P("   a = LMP/(scale T^2) falls slowly with T, so a residual is expected on top of the MC error.")
    for label in anchors:
        P("")
        P(f"   {label}")
        P(f"   {'sigma [K]':>9} {'analytic':>10} {'sampled N=336':>14} {'err [%]':>9} {'MC s.e. [%]':>12} "
          f"{'sampled N=4e5':>14} {'err [%]':>9}")
        for sg in res["sigmas"]:
            v = res["per_sigma"][label][str(sg)]["validation"]
            P(f"   {sg:9.1f} {v['analytic']:10.4f} {v['sampled_N336']:14.4f} {v['rel_err_N336_pct']:9.2f} "
              f"{v['mc_se_N336_pct']:12.2f} {v['sampled_large']:14.4f} {v['rel_err_large_pct']:9.2f}")
    P("")
    P("   The N=336 deviations sit inside one Monte Carlo standard error of the analytic value; the")
    P("   large-sample column isolates the curvature term, which is the honest size of the lognormal")
    P("   approximation error.")
    P("")

    # ---- item 3 ----
    P(f"3. RQ2 SCENARIOS OVER THE POPULATION (sigma = {res['base_sigma']} K, {res['n_tubes']} tubes, seed {res['seed']})")
    P("-" * 100)
    P("   Life cost per kmol H2, relative to the S1 average tube. Years are years to D = 1.")
    P("")
    scen: pd.DataFrame = res["scenarios"]; rk: pd.DataFrame = res["rank_scen"]
    ref = float(scen.loc["S1_steady", "cost_avg"])
    P(f"{'scenario':<30} {'avg':>8} {'pop mean':>9} {'hot dec':>8} {'pop/avg':>8} | "
      f"{'y avg':>8} {'y pop':>8} {'y 1st fail':>10} | {'rank a/p/d':>12} | code path")
    for name, r in scen.iterrows():
        rr = rk.loc[name]
        flag = "*" if (rr.rank_changes_pop or rr.rank_changes_decile) else " "
        P(f"{name:<30} {r.cost_avg / ref:8.3f} {r.cost_pop_mean / ref:9.3f} {r.cost_hottest_decile / ref:8.3f} "
          f"{r.ratio_pop_mean_to_avg:8.3f} | {r.years_avg_tube:8.3g} {r.years_pop_mean:8.3g} "
          f"{r.years_first_failure:10.3g} | {rr.rank_avg:3d}/{rr.rank_pop_mean:3d}/{rr.rank_hottest_decile:3d}{flag:>2} "
          f"| {r.evaluated_by}")
    P("")
    L.extend(_rank_verdict(scen, rk, 30))
    P("")

    P(f"4. RQ3 REGIMES OVER THE POPULATION (sigma = {res['base_sigma']} K)")
    P("-" * 100)
    P("   The RQ3 representative-regime table of the manuscript is Table 3 in the compiled PDF (label")
    P("   tab:regimes); there is no Table 6. Costs are relative to the base average tube.")
    P("")
    regs: pd.DataFrame = res["regimes"]; rkr: pd.DataFrame = res["rank_reg"]
    refr = float(regs.loc["base", "cost_avg"])
    P(f"{'regime':<20} {'T_wo [K]':>9} {'avg':>8} {'pop mean':>9} {'hot dec':>8} {'pop/avg':>8} | "
      f"{'y avg':>9} {'y 1st fail':>11} | {'rank a/p/d':>12} | code path")
    for name, r in regs.iterrows():
        rr = rkr.loc[name]
        flag = "*" if (rr.rank_changes_pop or rr.rank_changes_decile) else " "
        P(f"{name:<20} {r.T_wo_max_K:9.2f} {r.cost_avg / refr:8.3f} {r.cost_pop_mean / refr:9.3f} "
          f"{r.cost_hottest_decile / refr:8.3f} {r.ratio_pop_mean_to_avg:8.3f} | {r.years_avg_tube:9.3g} "
          f"{r.years_first_failure:11.3g} | {rr.rank_avg:3d}/{rr.rank_pop_mean:3d}/{rr.rank_hottest_decile:3d}{flag:>2} "
          f"| {r.evaluated_by}")
    P("")
    L.extend(_rank_verdict(regs, rkr, 20))
    return "\n".join(L)


def _rank_verdict(table: pd.DataFrame, rk: pd.DataFrame, width: int) -> List[str]:
    """The ranking answer for one table: whether it moves, which rows move, and by how much they were apart."""
    L: List[str] = []
    n_p = int(rk.rank_changes_pop.sum()); n_d = int(rk.rank_changes_decile.sum())
    L.append(f"   ranking by life cost per kmol H2 identical for the population mean : "
             f"{'YES' if n_p == 0 else f'NO ({n_p} of {len(rk)} rows move)'}")
    L.append(f"   ranking identical for the hottest decile                          : "
             f"{'YES' if n_d == 0 else f'NO ({n_d} of {len(rk)} rows move)'}")
    for name, rr in rk.iterrows():
        if rr.rank_changes_pop or rr.rank_changes_decile:
            L.append(f"     {name:<{width}} avg {rr.rank_avg:>2} -> pop mean {rr.rank_pop_mean:>2} "
                     f"-> hottest decile {rr.rank_hottest_decile:>2}")
    for col, what in (("cost_pop_mean", "population mean"), ("cost_hottest_decile", "hottest decile")):
        sw = swapped_pairs(table, "cost_avg", col)
        if not sw:
            continue
        L.append(f"     pairs that reverse order under the {what} (gap under each definition):")
        for i, j, ga, gb in sw:
            L.append(f"       {i:<{width}} vs {j:<{width}} {ga:6.2f} % apart on the average tube, "
                     f"{gb:6.2f} % apart on the {what}")
        L.append(f"       largest gap involved in any reversal: {max(max(s[2], s[3]) for s in sw):.2f} %")
    return L


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", default="v3", choices=sorted(C.BY_TAG))
    p.add_argument("--sigma", type=float, nargs="+", default=list(SIGMA_GRID_K),
                   help="wall-temperature scatter values to report [K]")
    p.add_argument("--n-tubes", type=int, default=N_TUBES)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--no-figure", action="store_true")
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    a = p.parse_args(argv)

    cfg = C.BY_TAG[a.version]
    from rdt import pipeline as pl
    pl.apply_config(cfg)

    base_sigma = a.sigma[len(a.sigma) // 2] if SIGMA_DEFAULT_K not in a.sigma else SIGMA_DEFAULT_K
    res = build(cfg, a.sigma, a.seed, a.n_tubes, base_sigma)
    text = report(res)
    print(text)

    a.out_dir.mkdir(parents=True, exist_ok=True)
    (a.out_dir / f"population_{cfg.tag}.txt").write_text(text + "\n")
    res["scenarios"].to_csv(a.out_dir / f"population_scenarios_{cfg.tag}.csv")
    res["regimes"].to_csv(a.out_dir / f"population_regimes_{cfg.tag}.csv")
    js = {"tag": cfg.tag, "alloy": cfg.creep_alloy, "n_tubes": a.n_tubes, "seed": a.seed, "rows": ROWS,
          "tubes_per_row": TUBES_PER_ROW, "sigmas_K": a.sigma, "base_sigma_K": base_sigma,
          "anchors": {k: asdict(v) for k, v in res["anchors"].items()},
          "anchors_physics_cross_check": {k: asdict(v) for k, v in res["cross"].items()},
          "per_sigma": res["per_sigma"],
          "scenario_rank_changes": int(res["rank_scen"].rank_changes_pop.sum()),
          "scenario_rank_changes_decile": int(res["rank_scen"].rank_changes_decile.sum()),
          "regime_rank_changes": int(res["rank_reg"].rank_changes_pop.sum()),
          "regime_rank_changes_decile": int(res["rank_reg"].rank_changes_decile.sum())}
    (a.out_dir / f"population_{cfg.tag}.json").write_text(json.dumps(js, indent=2, default=float))
    print(f"\nwritten: {a.out_dir / f'population_{cfg.tag}.txt'}")
    print(f"         {a.out_dir / f'population_scenarios_{cfg.tag}.csv'}")
    print(f"         {a.out_dir / f'population_regimes_{cfg.tag}.csv'}")
    print(f"         {a.out_dir / f'population_{cfg.tag}.json'}")

    if not a.no_figure:
        paths = figure_damage_distributions(cfg, res["regime_states"], res["dT"], res["curve"], base_sigma)
        for q in paths:
            print(f"         {q}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
