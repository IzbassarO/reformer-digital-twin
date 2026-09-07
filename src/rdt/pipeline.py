"""Single entry point for the reproducible pipeline: ``python -m rdt.pipeline --stage all [--version v2]``.

Stages (each writes ``*_<tag>`` files next to the v1 files and never overwrites v1):
``design`` (LHS campaign, closed-loop S/C x load grid, Saltelli/Sobol screening), ``surrogate`` (GP/HGB/MLP
surrogates with CV+ conformal intervals), ``scenarios`` (RQ2 hourly histories), ``pareto`` (RQ3 NSGA-II and the
steam-credit sweep), ``uq`` (RQ4 Monte Carlo and analysis), ``figures`` (journal-style figures from the results).
Runtimes per stage are recorded in ``data/pipeline_<tag>_runtimes.json``. The run configuration also selects the
creep master curve, so ``--version v3`` runs on the Centralloy G 4852 data-sheet curves rather than the placeholder.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

from rdt import config as C

STAGES = ("design", "surrogate", "scenarios", "pareto", "uq", "figures")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def apply_config(cfg: C.RunConfig) -> None:
    from rdt import operate, optimize, scenarios, steam_credit, surrogate, uq
    from rdt import creep
    for m in (creep, operate, surrogate, scenarios, optimize, steam_credit, uq):
        m.use_config(cfg)


def stage_design(cfg: C.RunConfig) -> dict:
    import numpy as np, pandas as pd
    from SALib.analyze import sobol
    from SALib.sample import saltelli
    from rdt import operate as op
    out = op.run_lhs_campaign(n=cfg.n_lhs, seed=cfg.lhs_seed, name=cfg.lhs_name, ranges_csv=cfg.ranges_csv)
    _log(f"LHS {out['meta']['n_converged']}/{out['meta']['n_samples']} in {out['meta']['wall_time_s']:.0f} s")
    grid = op.run_closed_loop_grid(n_sc=cfg.grid_n, n_load=cfg.grid_n, name=cfg.grid_name, ranges_csv=cfg.ranges_csv)
    _log(f"grid {int((grid.converged == True).sum())}/{len(grid)}")
    ranges = op.load_ranges(cfg.ranges_csv); base = op.base_case(); b = op.base_run(base)
    problem = {"num_vars": len(op.INPUT_NAMES), "names": list(op.INPUT_NAMES), "bounds": [[float(ranges.loc[k, "min"]), float(ranges.loc[k, "max"])] for k in op.INPUT_NAMES]}
    X = saltelli.sample(problem, cfg.sobol_N, calc_second_order=True)
    t = time.perf_counter(); df = op.run_batch(pd.DataFrame(X, columns=op.INPUT_NAMES), b["life_consumption_rate_per_h"], -1); dt = time.perf_counter() - t
    df.to_csv(cfg.saltelli_file, index=False, compression="gzip")
    res = {"problem": problem, "N": cfg.sobol_N, "n_runs": int(len(df)), "n_converged": int((df.converged == True).sum()), "wall_time_s": dt, "git_commit": op.git_hash(), "tag": cfg.tag, "outputs": {}}
    for y in ("T_wo_max_K", "log10_t_r_hot", "CH4_slip_dry_pct"):
        Y = pd.to_numeric(df[y], errors="coerce").to_numpy(float); Y = np.where(np.isfinite(Y), Y, np.nanmedian(Y))
        Si = sobol.analyze(problem, Y, calc_second_order=True, print_to_console=False)
        res["outputs"][y] = {k: dict(zip(op.INPUT_NAMES, map(float, Si[k]))) for k in ("S1", "S1_conf", "ST", "ST_conf")}
    json.dump(res, open(cfg.sobol_json, "w"), indent=2)
    _log(f"Sobol {res['n_converged']}/{res['n_runs']} in {dt:.0f} s")
    return {"lhs_converged": out["meta"]["convergence_rate"], "grid_converged": float((grid.converged == True).mean()), "sobol_converged": res["n_converged"] / res["n_runs"]}


def stage_surrogate(cfg: C.RunConfig) -> dict:
    import numpy as np
    from rdt import operate as op
    from rdt import surrogate as s
    out = s.train_all(verbose=False); M, tr, te = out["metrics"], out["train"], out["test"]
    M["conformal"] = {}; M["gp_native_intervals"] = {}
    for tgt in s.TARGETS:
        c = s.conformalize(tgt, M["targets"][tgt]["best"], tr, te); M["conformal"][tgt] = {k: v for k, v in c.items() if not k.startswith("_")}
    for tgt in ("T_wo_max_K", "CH4_slip_dry_pct", "log10_t_r_hot"):
        fs = "physics" if M["targets"][tgt]["models"]["gp__physics"]["test"]["rmse"] < M["targets"][tgt]["models"]["gp__base"]["test"]["rmse"] else "base"
        M["gp_native_intervals"][tgt] = {"features": fs, **s.gp_native_intervals(tgt, fs, tr, te)}
    b = op.base_case(); t = time.perf_counter(); r = op.run_case(b.inputs(), b); phys_t = time.perf_counter() - t
    best = M["targets"]["T_wo_max_K"]["best"]; sur_t = M["targets"]["T_wo_max_K"]["models"][best]["predict_time_per_1e5_s"] / 1e5
    M["speed"] = {"physics_model_s_per_run": phys_t, "best_T_wo_surrogate_s_per_run": sur_t, "speedup_factor": phys_t / sur_t, "physics_base_case_T_wo_max_K": r["T_wo_max_K"]}
    M["tag"] = cfg.tag
    json.dump(M, open(s.METRICS_JSON, "w"), indent=2)
    M["speed"]["surrogate_base_case_T_wo_max_K"] = s.predict_base_case("T_wo_max_K"); json.dump(M, open(s.METRICS_JSON, "w"), indent=2)
    _log(f"surrogates: best T_wo RMSE {M['targets']['T_wo_max_K']['models'][best]['test']['rmse']:.3f} K")
    return {"best_T_wo_rmse_K": M["targets"]["T_wo_max_K"]["models"][best]["test"]["rmse"]}


def stage_scenarios(cfg: C.RunConfig) -> dict:
    from rdt import scenarios as sc
    out = sc.run_all_scenarios()
    _log(f"scenarios: {len(out['results'])} scenario-years, physics calls {out['meta']['n_physics_calls']}")
    return {"n_scenarios": len(out["results"])}


def stage_pareto(cfg: C.RunConfig) -> dict:
    from rdt import optimize as ox
    from rdt import steam_credit as scr
    out = ox.run_all(pop=cfg.nsga_pop, gens=cfg.nsga_gens)
    _log(f"pareto: {out['n_feasible_phys']}/{out['n_pareto']} feasible under physics")
    ox.run_variant_total_heat(pop=cfg.nsga_pop, gens=cfg.nsga_gens)
    r = scr.run_steam_credit()
    _log("steam credit done")
    return {"n_feasible": out["n_feasible_phys"], "steam_credit_alphas": r["alphas"]}


def stage_uq(cfg: C.RunConfig) -> dict:
    from rdt import uq
    uq.run_campaign(n=cfg.n_mc, n_group=cfg.n_mc_group, log=_log)
    S = uq.analyse(log=_log)
    return {"N": S["N"], "convergence": S["convergence"]}


def stage_figures(cfg: C.RunConfig) -> dict:
    from rdt import figures
    made = figures.make_all(cfg)
    _log(f"figures: {len(made)} written")
    return {"n_figures": len(made)}


def run(stages, cfg: C.RunConfig) -> dict:
    apply_config(cfg)
    runtimes = json.load(open(cfg.runtimes_json)) if cfg.runtimes_json.exists() else {}
    for st in stages:
        _log(f"=== stage {st} ({cfg.tag}) ===")
        t0 = time.perf_counter()
        info = globals()[f"stage_{st}"](cfg)
        runtimes[st] = {"wall_time_s": time.perf_counter() - t0, "finished": time.strftime("%Y-%m-%d %H:%M:%S"), **info}
        json.dump(runtimes, open(cfg.runtimes_json, "w"), indent=2, default=float)
        _log(f"=== stage {st} done in {runtimes[st]['wall_time_s']:.0f} s ===")
    return runtimes


def main(argv=None):
    ap = argparse.ArgumentParser(description="Reformer digital twin pipeline")
    ap.add_argument("--stage", default="all", help="one of " + ", ".join(STAGES) + " or 'all' or a comma list")
    ap.add_argument("--version", default="v2", choices=sorted(C.BY_TAG))
    a = ap.parse_args(argv)
    cfg = C.BY_TAG[a.version]
    stages = list(STAGES) if a.stage == "all" else [s.strip() for s in a.stage.split(",")]
    for s in stages:
        if s not in STAGES:
            sys.exit(f"unknown stage {s}")
    run(stages, cfg)


if __name__ == "__main__":
    main()
