# Change report: v1 -> v2 headline numbers

v2 = consolidated release: excess-air range 5-25 % (v1: 5-20 %), base case = calibrated Plant A state (21.6 % excess air; v1 scenarios/optimisation used 20 %), Monte Carlo N = 2048 (v1: 4096). Everything else unchanged. Flag: relative change > 5 % (temperatures: > 5 K; indices, shares and fractions: > 0.05 absolute).

Generated 2026-09-06 17:03 at commit 3275df4c. Runtimes per v2 stage: design 214 s, surrogate 549 s, scenarios 191 s, pareto 172 s, uq 932 s, figures 10 s

| family | item | v1 | v2 | change | flag | explanation |
|---|---|---|---|---|---|---|
| sobol | T_wo_max_K ST steam_to_carbon | 0.308 | 0.284 | -0.02  |  |  |
| sobol | T_wo_max_K ST feed_per_tube_fraction | 0.0246 | 0.0276 | +0.00  |  |  |
| sobol | T_wo_max_K ST specific_firing_factor | 0.395 | 0.353 | -0.04  |  |  |
| sobol | T_wo_max_K ST inlet_T | 0.0606 | 0.0573 | -0.00  |  |  |
| sobol | T_wo_max_K ST inlet_P | 0.0195 | 0.018 | -0.00  |  |  |
| sobol | T_wo_max_K ST catalyst_activity | 0.0833 | 0.0758 | -0.01  |  |  |
| sobol | T_wo_max_K ST excess_air | 0.12 | 0.194 | +0.07  | **changed** | excess air now spans 5-25 % instead of 5-20 %: a wider input range carries more variance, so its total-effect index rises while all others fall slightly |
| sobol | log10_t_r_hot ST steam_to_carbon | 0.268 | 0.251 | -0.02  |  |  |
| sobol | log10_t_r_hot ST feed_per_tube_fraction | 0.0262 | 0.029 | +0.00  |  |  |
| sobol | log10_t_r_hot ST specific_firing_factor | 0.324 | 0.292 | -0.03  |  |  |
| sobol | log10_t_r_hot ST inlet_T | 0.0498 | 0.0476 | -0.00  |  |  |
| sobol | log10_t_r_hot ST inlet_P | 0.175 | 0.164 | -0.01  |  |  |
| sobol | log10_t_r_hot ST catalyst_activity | 0.0701 | 0.0646 | -0.01  |  |  |
| sobol | log10_t_r_hot ST excess_air | 0.0994 | 0.165 | +0.07  | **changed** | excess air now spans 5-25 % instead of 5-20 %: a wider input range carries more variance, so its total-effect index rises while all others fall slightly |
| scenarios | S2_daily life per kmol H2 rel. S1 | 0.985 | 0.991 | +0.6 % |  |  |
| scenarios | S3_renewable life per kmol H2 rel. S1 | 0.974 | 0.983 | +1.0 % |  |  |
| scenarios | S4a_mild_overfire life per kmol H2 rel. S1 | 1.01 | 1.01 | -0.0 % |  |  |
| scenarios | S4b_severe_overfire life per kmol H2 rel. S1 | 1 | 1 | +0.0 % |  |  |
| scenarios | S4c_hot_band life per kmol H2 rel. S1 | 1.03 | 1.03 | +0.0 % |  |  |
| scenarios | S6_SC2.5 life per kmol H2 rel. S1 | 1.15 | 1.16 | +1.1 % |  |  |
| scenarios | S6_SC3.5 life per kmol H2 rel. S1 | 0.822 | 0.834 | +1.5 % |  |  |
| scenarios | S5_hold_CH4_slip_ageing_y4 life per kmol H2 rel. S1 | 2.25 | 2.25 | +0.1 % |  |  |
| scenarios | S5_hold_T_out_ageing_y4 life per kmol H2 rel. S1 | 1.88 | 1.88 | -0.1 % |  |  |
| scenarios | S1 mean T_wo,max | 1.14e+03 | 1.14e+03 | -1.37 K |  |  |
| scenarios | S3 avg-condition error | 0.979 | 0.982 | +0.4 % |  |  |
| pareto | base H2 (kmol/h) | 14.9 | 14.8 | -0.6 % |  |  |
| pareto | base fuel (MJ/kmol) | 157 | 158 | +0.6 % |  |  |
| pareto | base life rel. base | 1 | 1 | +0.0 % |  |  |
| pareto | base T_wo,max | 1.15e+03 | 1.15e+03 | -3.27 K |  |  |
| pareto | base S/C | 2.89 | 2.89 | +0.0 % |  |  |
| pareto | min_life_iso_H2 H2 (kmol/h) | 14.8 | 14.7 | -0.8 % |  |  |
| pareto | min_life_iso_H2 fuel (MJ/kmol) | 128 | 134 | +4.8 % |  |  |
| pareto | min_life_iso_H2 life rel. base | 0.154 | 0.0825 | -46.5 % | **changed** | the minimum-life regime now uses the excess-air headroom above 20 % (leaner, cooler flame), lowering its hot spot by 13 K; the base denominator also lengthened by ~22 % |
| pareto | min_life_iso_H2 T_wo,max | 1.12e+03 | 1.11e+03 | -12.56 K | **changed** | the minimum-life regime now uses excess air above 20 %, which dilutes the flame and lowers the hot-spot temperature |
| pareto | min_life_iso_H2 S/C | 4 | 4 | +0.0 % |  |  |
| pareto | knee H2 (kmol/h) | 16.9 | 17.2 | +2.2 % |  |  |
| pareto | knee fuel (MJ/kmol) | 129 | 131 | +1.3 % |  |  |
| pareto | knee life rel. base | 0.173 | 0.28 | +62.3 % | **changed** | the knee point is selected by distance to the utopia point and moved along the front; its own hot spot is 5 K hotter and the base denominator lengthened by ~22 % |
| pareto | knee T_wo,max | 1.13e+03 | 1.13e+03 | +4.88 K |  |  |
| pareto | knee S/C | 4 | 4 | +0.0 % |  |  |
| pareto | max_H2 H2 (kmol/h) | 19.3 | 19.3 | +0.0 % |  |  |
| pareto | max_H2 fuel (MJ/kmol) | 154 | 154 | +0.0 % |  |  |
| pareto | max_H2 life rel. base | 3.55 | 4.33 | +21.9 % | **changed** | the regime itself is unchanged (same hot spot); the ratio rises because the v2 base at 21.6 % excess air is 3.3 K cooler and lives ~22 % longer |
| pareto | max_H2 T_wo,max | 1.18e+03 | 1.18e+03 | +0.00 K |  |  |
| pareto | max_H2 S/C | 4 | 4 | +0.0 % |  |  |
| pareto | min_fuel_overall H2 (kmol/h) | 9.87 | 9.87 | -0.0 % |  |  |
| pareto | min_fuel_overall fuel (MJ/kmol) | 121 | 121 | -0.0 % |  |  |
| pareto | min_fuel_overall life rel. base | 0.387 | 0.472 | +21.9 % | **changed** | the regime itself is unchanged (same hot spot); the ratio rises because the v2 base at 21.6 % excess air is 3.3 K cooler and lives ~22 % longer |
| pareto | min_fuel_overall T_wo,max | 1.13e+03 | 1.13e+03 | -0.00 K |  |  |
| pareto | min_fuel_overall S/C | 4 | 4 | +0.0 % |  |  |
| pareto | members dominating the base | 50 | 37 | -26.0 % | **changed** | the v2 base point (21.6 % excess air, cooler hot spot) is a better reference, so fewer Pareto members dominate it |
| uq | S1_base T_wo,max median | 1.15e+03 | 1.15e+03 | -1.40 K |  |  |
| uq | S1_base T_wo,max 90 % width | 42.2 | 41.9 | -0.34 K |  |  |
| uq | S1_base log10 t_r median | 7.87 | 7.91 | +0.04 dec |  |  |
| uq | S1_base life per kmol rel. base median | 1 | 1 | +0.0 % |  |  |
| uq | RQ3_knee T_wo,max median | 1.13e+03 | 1.13e+03 | +4.87 K |  |  |
| uq | RQ3_knee T_wo,max 90 % width | 56.5 | 57.8 | +1.26 K |  |  |
| uq | RQ3_knee log10 t_r median | 8.43 | 8.3 | -0.13 dec | **changed** | the v2 knee regime sits at a 5 K hotter hot spot |
| uq | RQ3_knee life per kmol rel. base median | 0.234 | 0.338 | +44.3 % | **changed** | the knee point is selected by distance to the utopia point and moved along the front; its own hot spot is 5 K hotter and the base denominator lengthened by ~22 % |
| uq | RQ3_min_life_iso_H2 T_wo,max median | 1.12e+03 | 1.11e+03 | -12.57 K | **changed** | the minimum-life regime now uses excess air above 20 %, which dilutes the flame and lowers the hot-spot temperature |
| uq | RQ3_min_life_iso_H2 T_wo,max 90 % width | 52 | 50.4 | -1.59 K |  |  |
| uq | RQ3_min_life_iso_H2 log10 t_r median | 8.54 | 8.9 | +0.36 dec | **changed** | the v2 minimum-life regime uses excess air above 20 % and runs 13 K cooler at the hot spot |
| uq | RQ3_min_life_iso_H2 life per kmol rel. base median | 0.209 | 0.1 | -52.0 % | **changed** | the minimum-life regime now uses the excess-air headroom above 20 % (leaner, cooler flame), lowering its hot spot by 13 K; the base denominator also lengthened by ~22 % |
| uq | S5_y4_holdslip_end T_wo,max median | 1.17e+03 | 1.17e+03 | -1.56 K |  |  |
| uq | S5_y4_holdslip_end T_wo,max 90 % width | 65.6 | 64.7 | -0.95 K |  |  |
| uq | S5_y4_holdslip_end log10 t_r median | 7.32 | 7.35 | +0.04 dec |  |  |
| uq | S5_y4_holdslip_end life per kmol rel. base median | 3.49 | 3.48 | -0.3 % |  |  |
| uq | variance share log10_t_r kinetics | 0.235 | 0.233 | -0.00  |  |  |
| uq | variance share log10_t_r heat_transfer | 0.0494 | 0.0504 | +0.00  |  |  |
| uq | variance share log10_t_r creep | 0.621 | 0.622 | +0.00  |  |  |
| uq | variance share log10_t_r measurement | 0.0935 | 0.0938 | +0.00  |  |  |
| uq | variance share T_wo_max_K kinetics | 0.505 | 0.513 | +0.01  |  |  |
| uq | variance share T_wo_max_K heat_transfer | 0.106 | 0.111 | +0.00  |  |  |
| uq | variance share T_wo_max_K creep | 0.209 | 0.211 | +0.00  |  |  |
| uq | variance share T_wo_max_K measurement | 0.197 | 0.202 | +0.01  |  |  |
| robustness | nominal_L_q S5y4_life_gt_2x_S1 | 0.908 | 0.906 | -0.00  |  |  |
| robustness | nominal_L_q knee_life_lt_0.5_base | 0.988 | 0.806 | -0.18  | **changed** | the v2 knee regime consumes 0.34 of base life instead of 0.23 (median), so fewer paired samples fall below the 0.5 threshold |
| robustness | nominal_L_q SC3.5_less_life_per_kmol_than_SC2.5 | 0.827 | 0.836 | +0.01  |  |  |
| robustness | nominal_L_q S2_life_per_kmol_le_S1 | 0.608 | 0.568 | -0.04  |  |  |
| robustness | nominal_L_q S3_life_per_kmol_le_S1 | 0.661 | 0.632 | -0.03  |  |  |
| robustness | load_dependent_L_q S5y4_life_gt_2x_S1 | 0.908 | 0.906 | -0.00  |  |  |
| robustness | load_dependent_L_q knee_life_lt_0.5_base | 0.998 | 0.936 | -0.06  | **changed** | the v2 knee regime consumes 0.34 of base life instead of 0.23 (median), so fewer paired samples fall below the 0.5 threshold |
| robustness | load_dependent_L_q SC3.5_less_life_per_kmol_than_SC2.5 | 0.827 | 0.836 | +0.01  |  |  |
| robustness | load_dependent_L_q S2_life_per_kmol_le_S1 | 0 | 0 | +0.00  |  |  |
| robustness | load_dependent_L_q S3_life_per_kmol_le_S1 | 0.000977 | 0.000977 | +0.00  |  |  |
| steam | alpha 0.0 min_life S/C | 4 | 3.99 | -0.1 % |  |  |
| steam | alpha 0.0 min_life life rel. base | 0.0724 | 0.0587 | -19.0 % | **changed** | the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 % |
| steam | alpha 0.0 min_fuel S/C | 2.76 | 2.66 | -3.5 % |  |  |
| steam | alpha 0.0 min_fuel life rel. base | 3.5 | 4.54 | +29.6 % | **changed** | the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 % |
| steam | alpha 0.25 min_life S/C | 4 | 4 | +0.0 % |  |  |
| steam | alpha 0.25 min_life life rel. base | 0.0901 | 0.0572 | -36.4 % | **changed** | the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 % |
| steam | alpha 0.25 min_fuel S/C | 3.29 | 2.92 | -11.1 % | **changed** | the minimum-fuel end of the iso-production front moved because the base fuel changed (21.6 % excess air) and the front is flat in S/C near its fuel optimum |
| steam | alpha 0.25 min_fuel life rel. base | 0.687 | 2.02 | +193.3 % | **changed** | the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 % |
| steam | alpha 0.5 min_life S/C | 4 | 4 | -0.0 % |  |  |
| steam | alpha 0.5 min_life life rel. base | 0.0663 | 0.0573 | -13.5 % | **changed** | the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 % |
| steam | alpha 0.5 min_fuel S/C | 3.34 | 3.63 | +8.7 % | **changed** | the minimum-fuel end of the iso-production front moved because the base fuel changed (21.6 % excess air) and the front is flat in S/C near its fuel optimum |
| steam | alpha 0.5 min_fuel life rel. base | 0.628 | 0.48 | -23.5 % | **changed** | the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 % |
| steam | alpha 0.75 min_life S/C | 4 | 4 | +0.0 % |  |  |
| steam | alpha 0.75 min_life life rel. base | 0.067 | 0.0569 | -15.2 % | **changed** | the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 % |
| steam | alpha 0.75 min_fuel S/C | 4 | 4 | +0.0 % |  |  |
| steam | alpha 0.75 min_fuel life rel. base | 0.225 | 0.276 | +22.7 % | **changed** | the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 % |
| steam | alpha 1.0 min_life S/C | 4 | 4 | +0.0 % |  |  |
| steam | alpha 1.0 min_life life rel. base | 0.0664 | 0.0572 | -13.9 % | **changed** | the v2 minimum-life points use excess air above 20 %, cooling the hot spot; the base denominator lengthened by ~22 % |
| steam | alpha 1.0 min_fuel S/C | 4 | 4 | +0.0 % |  |  |
| steam | alpha 1.0 min_fuel life rel. base | 0.224 | 0.276 | +22.9 % | **changed** | the minimum-fuel end of the iso-production front is a different member (flat fuel objective in S/C), and the base life denominator lengthened by ~22 % |
| surrogate | best test RMSE T_wo_max_K | 0.03 | 0.113 | +276.3 % | **changed** | same 2000 training points spread over a 25 % wider excess-air range (lower point density); the surrogate remains far inside the test thresholds (3 K, 0.05) |
| surrogate | best test RMSE CH4_slip_dry_pct | 0.0012 | 0.00125 | +3.7 % |  |  |
| surrogate | best test RMSE log10_t_r_hot | 0.00272 | 0.00566 | +108.0 % | **changed** | same 2000 training points spread over a 25 % wider excess-air range (lower point density); the surrogate remains far inside the test thresholds (3 K, 0.05) |

29 of 108 headline numbers changed beyond the threshold.
