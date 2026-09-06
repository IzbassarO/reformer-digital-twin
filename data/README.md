# Data inventory

Generated 2026-09-06 17:03 by `rdt.report.write_data_readme` at commit 3275df4c. Sizes in MB; commit = last commit touching the file (`uncommitted` = not yet in git).

| file | size (MB) | commit | produced by |
|---|---|---|---|
| data/LICENSE.md | 0.000 | d984aac | — |
| data/creep_derived/creep_config.yaml | 0.000 | uncommitted | rdt.creep (Yeh digitisation, NIMS ingestion template/config) |
| data/creep_derived/nims_creep_data_schema.yaml | 0.001 | 311241c | rdt.creep (Yeh digitisation, NIMS ingestion template/config) |
| data/creep_derived/nims_transcription_template.csv | 0.000 | uncommitted | rdt.creep (Yeh digitisation, NIMS ingestion template/config) |
| data/creep_derived/yeh2021_manaurite_xm.yaml | 0.006 | 311241c | rdt.creep (Yeh digitisation, NIMS ingestion template/config) |
| data/design/operating_space.yaml | 0.003 | bd43a52 | hand-written |
| data/design/parameter_ranges.csv | 0.001 | bd43a52 | hand-written operating-space definition |
| data/design/parameter_ranges_v2.csv | 0.001 | uncommitted | hand-written operating-space definition |
| data/design/sobol_v1.json | 0.005 | bd43a52 | rdt.pipeline stage design (rdt.operate) |
| data/design/sobol_v2.json | 0.005 | uncommitted | rdt.pipeline stage design (rdt.operate) |
| data/design/surrogate_metrics_v1.json | 0.033 | fd9d85d | rdt.pipeline stage surrogate (rdt.surrogate) |
| data/design/surrogate_metrics_v2.json | 0.031 | uncommitted | rdt.pipeline stage surrogate (rdt.surrogate) |
| data/design/v1_vs_v2_diff.md | 0.013 | uncommitted | rdt.report.write_diff_report |
| data/kinetics/xu_froment_1989.yaml | 0.010 | 9d782fa | hand transcription of Xu & Froment 1989 Table 5-7 |
| data/lhs_runs/grid_sc_load_v1.csv.gz | 0.052 | bd43a52 | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/grid_sc_load_v2.csv.gz | 0.052 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v1.csv.gz | 0.572 | bd43a52 | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v1_failed.csv | 0.000 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v1_meta.json | 0.002 | bd43a52 | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v2.csv.gz | 0.572 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v2_failed.csv | 0.000 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/lhs_plantA_v2_meta.json | 0.002 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/sobol_saltelli_v1.csv.gz | 0.837 | bd43a52 | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/lhs_runs/sobol_saltelli_v2.csv.gz | 0.837 | uncommitted | rdt.pipeline stage design (rdt.operate: LHS, closed-loop grid, Saltelli) |
| data/literature_validation/digitized/xf1989_fig3_p_t.csv | 0.018 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_tidy.csv | 0.224 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_x_CH4.csv | 0.018 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_x_CO2.csv | 0.016 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_x_CO2_T_gas.csv | 0.014 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_x_CO2_T_wall_inner.csv | 0.010 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/digitized/xf1989_fig3_x_CO2_T_wall_outer.csv | 0.005 | ca57cf9 | manual WebPlotDigitizer export + validate_xf1989 tidy step |
| data/literature_validation/latham2008_plant_cases.csv | 0.005 | 2662e30 | scratch build script from thesis Appendix H (see README section) |
| data/literature_validation/latham2011_cases.yaml | 0.017 | 9d782fa | — |
| data/literature_validation/latham_fit.yaml | 0.111 | 311241c | notebooks/04-06 (rdt.latham_cases, rdt.latham_calibration, rdt.creep) |
| data/literature_validation/xf1989_fit.yaml | 0.034 | 234b5a3 | notebooks/03_validation_xf1989.ipynb (rdt.validate_xf1989) |
| data/literature_validation/xu_froment_1989_industrial.yaml | 0.006 | 9d782fa | — |
| data/optimization/pareto_v1.csv | 0.094 | 1846737 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v1_alt_totalheat.csv | 0.096 | 1846737 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v1_alt_totalheat_meta.json | 0.001 | 1846737 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v1_meta.json | 0.003 | 1846737 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v2.csv | 0.093 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v2_alt_totalheat.csv | 0.096 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v2_alt_totalheat_meta.json | 0.001 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/pareto_v2_meta.json | 0.004 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/regimes_v1.csv | 0.002 | 1846737 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/regimes_v2.csv | 0.002 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/steam_credit_fronts_v1.csv | 0.315 | 3275df4 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/steam_credit_fronts_v2.csv | 0.314 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/steam_credit_v1.json | 0.007 | 3275df4 | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/optimization/steam_credit_v2.json | 0.007 | uncommitted | rdt.pipeline stage pareto (rdt.optimize, rdt.steam_credit) |
| data/pipeline_v2_runtimes.json | 0.001 | uncommitted | rdt.pipeline |
| data/scenarios/hourly_S1_steady.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S2_daily.csv.gz | 0.098 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S3_renewable.csv.gz | 0.731 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S4a_mild_overfire.csv.gz | 0.095 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S4b_severe_overfire.csv.gz | 0.093 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S4c_hot_band.csv.gz | 0.093 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_CH4_slip_ageing_y1.csv.gz | 0.617 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_CH4_slip_ageing_y2.csv.gz | 0.622 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_CH4_slip_ageing_y3.csv.gz | 0.628 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_CH4_slip_ageing_y4.csv.gz | 0.636 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_ageing_y1.csv.gz | 0.628 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_ageing_y2.csv.gz | 0.631 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_ageing_y3.csv.gz | 0.637 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_ageing_y4.csv.gz | 0.648 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_const_y1.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_const_y2.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_const_y3.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S5_hold_T_out_const_y4.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S6_SC2.5.csv.gz | 0.092 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_S6_SC3.5.csv.gz | 0.091 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S1_steady.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S2_daily.csv.gz | 0.098 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S3_renewable.csv.gz | 0.735 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S4a_mild_overfire.csv.gz | 0.095 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S4b_severe_overfire.csv.gz | 0.093 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S4c_hot_band.csv.gz | 0.093 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y1.csv.gz | 0.627 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y2.csv.gz | 0.631 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y3.csv.gz | 0.637 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y4.csv.gz | 0.647 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_ageing_y1.csv.gz | 0.638 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_ageing_y2.csv.gz | 0.641 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_ageing_y3.csv.gz | 0.647 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_ageing_y4.csv.gz | 0.657 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_const_y1.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_const_y2.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_const_y3.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S5_hold_T_out_const_y4.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S6_SC2.5.csv.gz | 0.093 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/hourly_v2_S6_SC3.5.csv.gz | 0.092 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/summary_v1.csv | 0.007 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/summary_v1_meta.json | 0.008 | 53f1718 | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/summary_v2.csv | 0.007 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/scenarios/summary_v2_meta.json | 0.008 | uncommitted | rdt.pipeline stage scenarios (rdt.scenarios) |
| data/uq/mc_v1.csv.gz | 13.041 **(>5 MB)** | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_RQ3_knee.csv.gz | 1.106 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_RQ3_knee_lqvar.csv.gz | 1.106 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_RQ3_min_life_iso_H2.csv.gz | 1.108 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_RQ3_min_life_iso_H2_lqvar.csv.gz | 1.107 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S1_base.csv.gz | 1.105 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S1_base_group_creep.csv.gz | 0.126 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S1_base_group_heat_transfer.csv.gz | 0.114 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S1_base_group_kinetics.csv.gz | 0.197 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S1_base_group_measurement.csv.gz | 0.059 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S5_y4_holdslip_end.csv.gz | 1.109 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S6_SC2.5.csv.gz | 1.105 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_S6_SC3.5.csv.gz | 1.106 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_design.csv.gz | 0.743 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_load_0.70.csv.gz | 1.103 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_load_0.70_lqvar.csv.gz | 1.103 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_load_0.85.csv.gz | 1.104 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_load_0.85_lqvar.csv.gz | 1.104 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v1_points.json | 0.002 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2.csv.gz | 6.521 **(>5 MB)** | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_RQ3_knee.csv.gz | 0.554 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_RQ3_knee_lqvar.csv.gz | 0.554 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_RQ3_min_life_iso_H2.csv.gz | 0.555 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_RQ3_min_life_iso_H2_lqvar.csv.gz | 0.554 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S1_base.csv.gz | 0.553 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S1_base_group_creep.csv.gz | 0.127 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S1_base_group_heat_transfer.csv.gz | 0.113 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S1_base_group_kinetics.csv.gz | 0.197 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S1_base_group_measurement.csv.gz | 0.059 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S5_y4_holdslip_end.csv.gz | 0.555 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S6_SC2.5.csv.gz | 0.553 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_S6_SC3.5.csv.gz | 0.554 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_design.csv.gz | 0.372 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_load_0.70.csv.gz | 0.552 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_load_0.70_lqvar.csv.gz | 0.552 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_load_0.85.csv.gz | 0.553 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_load_0.85_lqvar.csv.gz | 0.553 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/mc_v2_points.json | 0.002 | uncommitted | rdt.pipeline stage uq (rdt.uq) |
| data/uq/uncertainty_spec.yaml | 0.004 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/uq_summary_v1.json | 0.027 | 3275df4 | rdt.pipeline stage uq (rdt.uq) |
| data/uq/uq_summary_v2.json | 0.027 | uncommitted | rdt.pipeline stage uq (rdt.uq) |

## Files larger than 5 MB

- data/uq/mc_v1.csv.gz (13.0 MB)
- data/uq/mc_v2.csv.gz (6.5 MB)

## Proposed Zenodo deposit at publication (list only; nothing deleted)

Campaign outputs that are reproducible from `python -m rdt.pipeline --stage all` and are too large or too numerous for the repository:

- data/lhs_runs/lhs_plantA_v1.csv.gz
- data/lhs_runs/lhs_plantA_v1_failed.csv
- data/lhs_runs/lhs_plantA_v1_meta.json
- data/lhs_runs/lhs_plantA_v2.csv.gz
- data/lhs_runs/lhs_plantA_v2_failed.csv
- data/lhs_runs/lhs_plantA_v2_meta.json
- data/lhs_runs/sobol_saltelli_v1.csv.gz
- data/lhs_runs/sobol_saltelli_v2.csv.gz
- data/scenarios/hourly_S1_steady.csv.gz
- data/scenarios/hourly_S2_daily.csv.gz
- data/scenarios/hourly_S3_renewable.csv.gz
- data/scenarios/hourly_S4a_mild_overfire.csv.gz
- data/scenarios/hourly_S4b_severe_overfire.csv.gz
- data/scenarios/hourly_S4c_hot_band.csv.gz
- data/scenarios/hourly_S5_hold_CH4_slip_ageing_y1.csv.gz
- data/scenarios/hourly_S5_hold_CH4_slip_ageing_y2.csv.gz
- data/scenarios/hourly_S5_hold_CH4_slip_ageing_y3.csv.gz
- data/scenarios/hourly_S5_hold_CH4_slip_ageing_y4.csv.gz
- data/scenarios/hourly_S5_hold_T_out_ageing_y1.csv.gz
- data/scenarios/hourly_S5_hold_T_out_ageing_y2.csv.gz
- data/scenarios/hourly_S5_hold_T_out_ageing_y3.csv.gz
- data/scenarios/hourly_S5_hold_T_out_ageing_y4.csv.gz
- data/scenarios/hourly_S5_hold_T_out_const_y1.csv.gz
- data/scenarios/hourly_S5_hold_T_out_const_y2.csv.gz
- data/scenarios/hourly_S5_hold_T_out_const_y3.csv.gz
- data/scenarios/hourly_S5_hold_T_out_const_y4.csv.gz
- data/scenarios/hourly_S6_SC2.5.csv.gz
- data/scenarios/hourly_S6_SC3.5.csv.gz
- data/scenarios/hourly_v2_S1_steady.csv.gz
- data/scenarios/hourly_v2_S2_daily.csv.gz
- data/scenarios/hourly_v2_S3_renewable.csv.gz
- data/scenarios/hourly_v2_S4a_mild_overfire.csv.gz
- data/scenarios/hourly_v2_S4b_severe_overfire.csv.gz
- data/scenarios/hourly_v2_S4c_hot_band.csv.gz
- data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y1.csv.gz
- data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y2.csv.gz
- data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y3.csv.gz
- data/scenarios/hourly_v2_S5_hold_CH4_slip_ageing_y4.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_ageing_y1.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_ageing_y2.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_ageing_y3.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_ageing_y4.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_const_y1.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_const_y2.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_const_y3.csv.gz
- data/scenarios/hourly_v2_S5_hold_T_out_const_y4.csv.gz
- data/scenarios/hourly_v2_S6_SC2.5.csv.gz
- data/scenarios/hourly_v2_S6_SC3.5.csv.gz
- data/uq/mc_v1.csv.gz
- data/uq/mc_v1_RQ3_knee.csv.gz
- data/uq/mc_v1_RQ3_knee_lqvar.csv.gz
- data/uq/mc_v1_RQ3_min_life_iso_H2.csv.gz
- data/uq/mc_v1_RQ3_min_life_iso_H2_lqvar.csv.gz
- data/uq/mc_v1_S1_base.csv.gz
- data/uq/mc_v1_S1_base_group_creep.csv.gz
- data/uq/mc_v1_S1_base_group_heat_transfer.csv.gz
- data/uq/mc_v1_S1_base_group_kinetics.csv.gz
- data/uq/mc_v1_S1_base_group_measurement.csv.gz
- data/uq/mc_v1_S5_y4_holdslip_end.csv.gz
- data/uq/mc_v1_S6_SC2.5.csv.gz
- data/uq/mc_v1_S6_SC3.5.csv.gz
- data/uq/mc_v1_design.csv.gz
- data/uq/mc_v1_load_0.70.csv.gz
- data/uq/mc_v1_load_0.70_lqvar.csv.gz
- data/uq/mc_v1_load_0.85.csv.gz
- data/uq/mc_v1_load_0.85_lqvar.csv.gz
- data/uq/mc_v1_points.json
- data/uq/mc_v2.csv.gz
- data/uq/mc_v2_RQ3_knee.csv.gz
- data/uq/mc_v2_RQ3_knee_lqvar.csv.gz
- data/uq/mc_v2_RQ3_min_life_iso_H2.csv.gz
- data/uq/mc_v2_RQ3_min_life_iso_H2_lqvar.csv.gz
- data/uq/mc_v2_S1_base.csv.gz
- data/uq/mc_v2_S1_base_group_creep.csv.gz
- data/uq/mc_v2_S1_base_group_heat_transfer.csv.gz
- data/uq/mc_v2_S1_base_group_kinetics.csv.gz
- data/uq/mc_v2_S1_base_group_measurement.csv.gz
- data/uq/mc_v2_S5_y4_holdslip_end.csv.gz
- data/uq/mc_v2_S6_SC2.5.csv.gz
- data/uq/mc_v2_S6_SC3.5.csv.gz
- data/uq/mc_v2_design.csv.gz
- data/uq/mc_v2_load_0.70.csv.gz
- data/uq/mc_v2_load_0.70_lqvar.csv.gz
- data/uq/mc_v2_load_0.85.csv.gz
- data/uq/mc_v2_load_0.85_lqvar.csv.gz
- data/uq/mc_v2_points.json
- models/ (trained surrogates, ~400 MB, git-ignored; regenerate with the surrogate stage)
- literature/ (copyrighted PDFs: never deposited)
