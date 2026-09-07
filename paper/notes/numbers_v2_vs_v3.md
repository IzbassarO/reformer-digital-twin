# Numbers in the manuscript: v2 (in the text now) vs v3 (Centralloy G 4852)

Keyed to the LaTeX `\label{}` of the unit the number sits in, **not** to a printed table number:
the archived draft and the current manuscript number their tables differently, so labels are the only
stable handle. `paper/notes/reconcile_overleaf_vs_repo.md` maps labels between the two versions.

- **v2** = the numbers currently written into `paper/main.tex` and `paper/supplementary.tex`. This
  manuscript carries them inline; there is no `\input{numbers}` to regenerate.
- **v3** = the same pipeline on the Centralloy G 4852 data-sheet curve (`RunConfig` tag `v3`,
  *C* = 18.6, scatter measured from the printed lower scatter band).
- Operating space, sample sizes and every model parameter are identical between the two.

Because the numbers are inline, **every change below is a manual edit** -- there is no generator to
re-run against this manuscript. `rdt.paper_tables` still regenerates the *archived* draft's tables
under `paper/archive/draft_step21/`, and those carry the same values, so it can be used as a source
for the replacement cells.


## Summary of what moves

| | count |
|---|---|
| pipeline quantities bit-identical between v2 and v3 | 16 |
| moved but round to the same printed digits | 5 |
| changed in the printed digits | 38 |

Verified directly on the raw result files: across 2000 LHS runs, 4096 Saltelli runs, 2048 Monte Carlo
samples and the 20-row scenario summary, **every thermal, hydraulic and hydrogen column is
bit-identical** (`T_wo_max_K`, `T_out_K`, `P_out_bar`, `H2_net_kmol_h`, `CH4_slip_dry_pct`,
`sigma_hot_MPa`, `duty_W`, `T_fg_out_K`, `z_frac_T_wo_max`, `CH4_conversion`, `Q_comb_W`). Only
creep-derived columns move. The alloy enters the model only where `t_r(T_wo, sigma)` is evaluated.


## Frontmatter (abstract and highlights -- no `\label{}`)

Both blocks are inside `\begin{frontmatter}` of `paper/main.tex`.

| Quantity | v2 in the text | v3 | change | action |
|---|---|---|---|---|
| calibration residuals "4~K ... 12~K" | 4, 12 | unchanged | -- | none |
| "82\,\% of the hot-spot temperature variance" | 82 | unchanged (Sobol $S_T$ on $T_{wo,max}$ bit-identical) | -- | none |
| "load-following and renewable-following cost **no** tube life per kmol" | 0.991 / 0.983 | 1.010 / 1.013 | sign reverses | **rewrite** |
| ageing multiplier "2.3 within four years" | 2.255 | 1.959 | -13.1 % | 2.3 -> 2.0 |
| knee "+16\,\% hydrogen" | +16.24 | +15.18 | | +16 -> +15 |
| knee "$-18$\,\% fuel" | -17.59 | -17.31 | | -18 -> -17 |
| knee "0.28 of the base life cost" | 0.2801 | 0.2810 | +0.31 % | none, rounds the same |
| "creep-curve scatter contributes 62\,\%" | 62.2 | 85.9 | +38.0 % | 62 -> 86, and "contributes" -> "dominates" |

**Highlights block.** Bullet 3 ("Catalyst ageing costs 2.3 times more tube life than any flexible-load
scenario") becomes 2.0; bullets 1, 2 and 5 are unaffected. Bullet 4 mentions no number. Elsevier's
85-character limit still applies.

**Abstract wording to assert.** Load-following must lose its direction: on the G 4852 curve S2 and S3
cost 1.010 and 1.013 of the steady life per kmol, i.e. marginally *more* rather than nothing. Both
sit inside 2 % of unity and the robustness fractions support the claim in fewer than half the
samples, so the honest statement is that flexible operation changes the life cost per kmol by under
2 % in either direction and the sign is not resolved by the model.


## `sec:creep` — Tube wall and creep life

| Text in the manuscript | status |
|---|---|
| "$C = 22.96$" and the cubic master curve fitted to the Manaurite XM minimum curve of Yeh | **must change**: *C* = 18.6, Centralloy G 4852 lower scatter band |
| "This curve is a *placeholder*" / "carries no scatter information" | **must change**: the sheet prints an Average and a Lower Scatter Band curve at the 95 % confidence level, and the scatter is now measured from their separation |
| "the module is written so that the curve can be replaced by ... the NIMS creep data sheets" | **must change**: the replacement is done, from manufacturer data sheets |
| "the placeholder rupture curve halves the life every 11~K" (in `sec:intro`) | **must change**: recompute on G 4852 |
| "ratios of rupture times are far less sensitive to the curve than absolute values" | keep -- now demonstrated rather than asserted, see `data/creep_derived/alloy_comparison_v3.txt` |

**Wording to assert.** The master curve is fitted separately to the Average and Lower Scatter Band
curves of the Centralloy G 4852 data sheet (September 2009, Rev. 02), each with the constant printed
on that sheet, and the scatter model is the measured horizontal separation of the two,
`delta_log10_tr = (LMP_avg - LMP_min)/(scale*T)`, not an assumed 0.3 decades. Worth adding: the three
sheets do not share a constant (G 4852 18.6, G 4852 Micro 22.9, ET 45 Micro 19.3), so the 22.9 often
quoted for "the 4852 family" belongs to the Micro sheet alone.


## `sec:data` — Data, `\paragraph{Creep data}`

The paragraph describing the Yeh transcription and the NIMS template must be rewritten around the
three Schmidt + Clemens data sheets digitised by `rdt.creep_ingest_datasheet`, including the
non-redistribution note and the Wayback provenance. The synthetic-recovery claim ("recovered within
10 %") survives and was re-verified: *C* to 0.17 %, scatter to 6.4 % (`data/creep_derived/validation_v3.txt`).
The NIMS path stays as the remaining step, but for heat-to-heat statistics rather than for the master
curve itself.


## `sec:campaigns` — `\paragraph{RQ4}`, the creep uncertainty group

| Quantity | v2 | v3 | action |
|---|---|---|---|
| scatter on $\log_{10} t_r$ | 0.3 decades (assumed, "factor of two between heats") | 0.181 decades (measured band / 1.645) | **rewrite** |
| $C$ prior | uniform $\pm 1$ of 22.96 | uniform $\pm 1$ of 18.6 | **rewrite** |
| "with the master curve refitted for each $C$" | -- | no longer true | **rewrite** |

**Wording to assert.** The scatter is the median 0.298-decade separation of the data-sheet Average and
Lower Scatter Band curves at 1147 K, divided by 1.645 because the sheet calls the lower band the 95 %
confidence level. A digitised curve is already a function of LMP, so *C* now shifts only the
LMP-to-time conversion and the polynomial is carried over unchanged rather than refitted.


## `sec:verification` — Verification and validation

Calibration and Xu--Froment numbers come from `data/literature_validation/`, not the pipeline, and are unchanged. One pipeline number:

| Text | v2 | v3 | change | action |
|---|---|---|---|---|
| "the rupture time from the placeholder minimum curve is $10^8$~h at this point" | 7.64e7 h ($10^{7.88}$) | 1.14e7 h ($10^{7.06}$) | -85 % | $10^8$ -> $10^7$ |
| hot-spot temperature 1147~K, hoop stress 12.9~MPa | 1147 | 1147 | 0 | none |

The sentence that follows ("evaluating the life at the tube-average temperature instead of the hot
spot overstates it by a factor of about 110") is a ratio on the same curve and needs re-checking, but
it is a within-curve comparison so it moves far less than the absolute value.


## `tab:xf` / `tab:latham` / `tab:params` — verification and calibration tables

**Unchanged.** All three are built from `data/literature_validation/`, which the alloy change does not touch. Cell-for-cell identical between v2 and v3.


## `tab:sobol` — Sobol indices

The $T_{wo,max}$ and CH$_4$-slip columns are bit-identical; only the $\log_{10}t_r$ columns move.

| v2 row | v3 row |
|---|---|
| `steam-to-carbon & 0.29 $\pm$ 0.09 & 0.28 $\pm$ 0.05 & 0.25 $\pm$ 0.08 & 0.25 $\pm$ 0.04 & 0.04 $\pm$ 0.03 & 0.04 $\pm$ 0.01` | `steam-to-carbon & 0.29 $\pm$ 0.09 & 0.28 $\pm$ 0.05 & 0.25 $\pm$ 0.09 & 0.25 $\pm$ 0.05 & 0.04 $\pm$ 0.03 & 0.04 $\pm$ 0.01` |
| `load & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.01 & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.01 & 0.19 $\pm$ 0.07 & 0.19 $\pm$ 0.03` | `load & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.01 & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.01 & 0.19 $\pm$ 0.07 & 0.19 $\pm$ 0.04` |
| `firing & 0.34 $\pm$ 0.09 & 0.35 $\pm$ 0.06 & 0.29 $\pm$ 0.09 & 0.29 $\pm$ 0.05 & 0.54 $\pm$ 0.12 & 0.54 $\pm$ 0.10` | `firing & 0.34 $\pm$ 0.09 & 0.35 $\pm$ 0.06 & 0.28 $\pm$ 0.08 & 0.29 $\pm$ 0.06 & 0.54 $\pm$ 0.12 & 0.54 $\pm$ 0.08` |
| `inlet $T$ & 0.06 $\pm$ 0.05 & 0.06 $\pm$ 0.01 & 0.05 $\pm$ 0.04 & 0.05 $\pm$ 0.01 & 0.08 $\pm$ 0.06 & 0.09 $\pm$ 0.02` | `inlet $T$ & 0.06 $\pm$ 0.04 & 0.06 $\pm$ 0.01 & 0.05 $\pm$ 0.04 & 0.05 $\pm$ 0.01 & 0.08 $\pm$ 0.05 & 0.09 $\pm$ 0.02` |
| `inlet $p$ & 0.02 $\pm$ 0.02 & 0.02 $\pm$ 0.00 & 0.16 $\pm$ 0.06 & 0.16 $\pm$ 0.03 & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.01` | `inlet $p$ & 0.02 $\pm$ 0.02 & 0.02 $\pm$ 0.00 & 0.17 $\pm$ 0.07 & 0.17 $\pm$ 0.03 & 0.03 $\pm$ 0.03 & 0.03 $\pm$ 0.00` |
| `activity & 0.07 $\pm$ 0.04 & 0.08 $\pm$ 0.02 & 0.06 $\pm$ 0.05 & 0.06 $\pm$ 0.01 & 0.00 $\pm$ 0.01 & 0.00 $\pm$ 0.00` | `activity & 0.07 $\pm$ 0.05 & 0.08 $\pm$ 0.02 & 0.06 $\pm$ 0.04 & 0.06 $\pm$ 0.01 & 0.00 $\pm$ 0.01 & 0.00 $\pm$ 0.00` |
| `excess air & 0.19 $\pm$ 0.06 & 0.19 $\pm$ 0.03 & 0.16 $\pm$ 0.07 & 0.16 $\pm$ 0.03 & 0.11 $\pm$ 0.04 & 0.11 $\pm$ 0.02` | `excess air & 0.19 $\pm$ 0.07 & 0.19 $\pm$ 0.04 & 0.16 $\pm$ 0.06 & 0.16 $\pm$ 0.03 & 0.11 $\pm$ 0.06 & 0.11 $\pm$ 0.02` |

## `sec:rq1` — RQ1

| Text | v2 | v3 | change | action |
|---|---|---|---|---|
| surrogate RMSE on $\log_{10} t_r$, "0.0057 decades" | 0.00566 | 0.00494 | -12.7 % | 0.0057 -> 0.0049 |
| "speed-up of about 4700" | 4658 | 4712 | +1.2 % | none, still "about 4700" |
| $S_T$ of inlet pressure on $\log_{10} t_r$ | 0.164 | 0.173 | +5.6 % | 0.16 -> 0.17 |
| surrogate RMSE on $T_{wo,max}$, PICP, $z/L$ $R^2 = 0.64$ | | bit-identical | 0 | none |

## `sec:rq2` — RQ2

| Text | v2 | v3 | change | action |
|---|---|---|---|---|
| S2 daily, 0.991 | 0.991 | 1.010 | +1.9 % | **rewrite** |
| S3 renewable, 0.983 | 0.983 | 1.013 | +3.0 % | **rewrite** |
| S/C 2.5, 1.16 | 1.161 | 1.147 | -1.2 % | substitute |
| S/C 3.5, 0.83 | 0.834 | 0.847 | +1.6 % | substitute |
| ageing year 1, 0.71 | 0.715 | 0.756 | +5.8 % | substitute |
| ageing year 4, 2.25 | 2.255 | 1.959 | -13.1 % | substitute |
| hold $T_{out}$ year 4, 1.88 | 1.875 | 1.686 | -10.1 % | substitute |
| mild event, 2.4x | 2.378 | 2.056 | -13.5 % | substitute |
| severe event, 6.9x | 6.899 | 4.932 | -28.5 % | substitute |
| hot band, 11.1x | 11.115 | 7.415 | -33.3 % | substitute |
| annual-mean ratio, 0.982 | 0.982 | 0.987 | +0.5 % | substitute |
| "hot spot climbs from 1138 to 1164~K" | 1138 / 1164 | unchanged (temperatures bit-identical) | 0 | none |
| "share of the annual damage is 0.4--2.7\,\%" | 0.40--2.70 | 0.27--1.73 | | 0.3--1.7 |
| "the life consumed in year 4 is three times that in year 1" | 3.15 | 2.59 | | three -> two and a half |
| "annual life cost per kmol rises by at most 3\,\%" | 2.8 % | 1.8 % | | may tighten to 2 % |

**Wording to assert.** The sentence "the daily profile (S2) and the renewable-following profile (S3)
consume 0.991 and 0.983 of the S1 life per kmol H$_2$, i.e. slightly *less*, because ..." reverses:
on the G 4852 curve both come out marginally *dearer* (1.010 and 1.013). The mechanism given -- the
hot spot falls with load at constant outlet temperature -- is unchanged and still correct; what
changes is that it no longer quite outweighs the hours at overload on the flatter master curve. The
paragraph should say load-following is life-neutral to within about 1 % and drop the italicised
*less*. `tab:robust` supports this: S2 $\le$ S1 holds in 0.425 of samples.


## `tab:scenarios` — RQ2 scenario summary

Hydrogen, $\bar T_{wo,max}$, max $T_{wo,max}$ and firing are bit-identical for every row. Damage and relative-life columns:

| Scenario | $D$ v2 | $D$ v3 | rel. life v2 | rel. life v3 | avg-cond. v2 | avg-cond. v3 |
|---|---|---|---|---|---|---|
| S1_steady | 9.08e-05 | 0.00063 | 1.000 | 1.000 | 1.000 | 1.000 |
| S2_daily | 8.1e-05 | 0.000572 | 0.991 | 1.010 | 0.982 | 0.987 |
| S3_renewable | 7.61e-05 | 0.000544 | 0.983 | 1.013 | 0.982 | 0.987 |
| S4a_mild_overfire | 9.14e-05 | 0.000634 | 1.007 | 1.006 | 0.993 | 0.994 |
| S4b_severe_overfire | 9.11e-05 | 0.000632 | 1.004 | 1.003 | 0.996 | 0.997 |
| S4c_hot_band | 9.33e-05 | 0.000641 | 1.028 | 1.018 | 0.980 | 0.988 |
| S6_SC2.5 | 9.83e-05 | 0.000674 | 1.161 | 1.147 | 1.000 | 1.000 |
| S6_SC3.5 | 8.23e-05 | 0.00058 | 0.834 | 0.847 | 1.000 | 1.000 |
| S5_hold_CH4_slip_ageing_y1 | 6.49e-05 | 0.000477 | 0.715 | 0.756 | 0.997 | 0.998 |
| S5_hold_CH4_slip_ageing_y4 | 0.000205 | 0.00123 | 2.255 | 1.959 | 0.949 | 0.961 |
| S5_hold_T_out_ageing_y1 | 7e-05 | 0.000508 | 0.768 | 0.803 | 0.998 | 0.999 |
| S5_hold_T_out_ageing_y4 | 0.000169 | 0.00105 | 1.875 | 1.686 | 0.967 | 0.975 |

Absolute damage rises by a factor of about 7 in every row, because G 4852 gives a shorter absolute
rupture time at this operating point than the placeholder (1.14e7 h against 7.64e7 h at 1147 K and
12.9 MPa). The ratios move far less. **The caption must stop calling the curve the Yeh Manaurite XM
placeholder** and name the Centralloy G 4852 lower scatter band with *C* = 18.6.


## `sec:rq3` — RQ3

| Text | v2 | v3 | change | action |
|---|---|---|---|---|
| feasible members, 199 | 199.000 | 198.000 | -0.5 % | substitute |
| knee +16 % | 16.240 | 15.177 | -6.5 % | substitute |
| knee -18 % | -17.591 | -17.314 | -1.6 % | substitute |
| knee 0.28 | 0.280 | 0.281 | +0.3 % | substitute |
| iso-H2 min-life 0.083 | 0.083 | 0.123 | +48.5 % | **largest move in the paper** |
| iso-H2 fuel -16 % | -15.591 | -15.584 | -0.0 % | substitute |
| iso-H2 hot spot 1109 K | 1108.513 | 1108.355 | -0.0 % | substitute |
| max-H2 life 4.3x | 4.329 | 3.220 | -25.6 % | substitute |
| min-fuel S/C at alpha=0, 2.66 | 2.663 | 2.639 | -0.9 % | substitute |
| min-life life 0.059 | 0.059 | 0.094 | +60.7 % | substitute |
| min-life fuel +3.2 % | 3.212 | 4.084 | +27.2 % | substitute |
| "the knee to 3.3" ($\alpha = 0$) | 3.29 | 3.36 | | 3.3 -> 3.4 |
| "maximum discrepancy ... 0.03 in relative life" | 0.0286 | 0.0287 | | none |
| "upper bound of the steam-to-carbon ratio (3.96--4.00)" | 3.96--4.00 | 3.96--4.00 | 0 | none |

**Wording to assert.** `\isoLife` is the largest single move in the manuscript, 0.083 -> 0.123
(+48 %): at constant production the life cost falls to about 0.12 of the base, an eightfold rather
than a twelvefold reduction. The maximum-hydrogen corner costs 3.2 times the base life per kmol, not
4.3. The knee is re-optimised under G 4852 and no longer sits exactly at the inlet-temperature bound
(893 K against 900 K, excess air 5.4 against 5.0 %), so "the maximum inlet temperature" becomes "a
high inlet temperature".

**Pre-existing inaccuracy, unrelated to the alloy.** "For $\alpha \ge 0.5$ ... the ratio of 4.0 is
optimal for all three regimes" is wrong in both versions: at $\alpha = 0.5$ the minimum-fuel regime
sits at S/C 3.63 (v2) and 3.69 (v3). It first holds at $\alpha = 0.75$.


## `tab:regimes` — RQ3 representative regimes

The Pareto set is re-optimised under G 4852, so the regimes sit at slightly different operating points.

| Regime | rel. life v2 | rel. life v3 | change | $T_{wo,max}$ v2 | v3 | H$_2$ v2 | v3 |
|---|---|---|---|---|---|---|---|
| base | 1.000 | 1.000 | +0.0 % | 1147 | 1147 | 14.82 | 14.82 |
| min_life_iso_H2 | 0.083 | 0.123 | +48.5 % | 1109 | 1108 | 14.72 | 14.87 |
| knee | 0.280 | 0.281 | +0.3 % | 1130 | 1127 | 17.22 | 17.06 |
| max_H2 | 4.329 | 3.220 | -25.6 % | 1177 | 1177 | 19.26 | 19.26 |
| min_fuel_overall | 0.472 | 0.575 | +21.8 % | 1128 | 1128 | 9.87 | 9.87 |
| min_life_overall | 0.046 | 0.084 | +82.1 % | 1091 | 1092 | 9.18 | 9.56 |

Feasible members fall from 199 to 198 of 200. The caption's "(Yeh placeholder curve)" must be replaced.


## `sec:rq4` — RQ4

| Text | v2 | v3 | change | action |
|---|---|---|---|---|
| $T_{wo,max}$ 90 % width, 42 K | 41.913 | 41.913 | +0.0 % | none |
| $\log_{10}t_r$ p05, 7.04 | 7.036 | 5.938 | -15.6 % | substitute |
| $\log_{10}t_r$ p95, 8.73 | 8.725 | 8.200 | -6.0 % | substitute |
| creep share, 62 % | 62.197 | 85.853 | +38.0 % | substitute |
| heat-transfer share, 5 % | 5.042 | 1.846 | -63.4 % | substitute |
| kinetics share of $T_{wo}$, 51 % | 51.276 | 51.276 | +0.0 % | none |
| ageing fraction, 0.91 | 0.906 | 0.819 | -9.6 % | substitute |
| knee fraction, 0.81 | 0.806 | 0.885 | +9.8 % | substitute |
| S/C fraction, 0.84 | 0.836 | 0.856 | +2.5 % | substitute |
| load 0.70 under the $L_q$ variant, 1.71 | 1.711 | 1.664 | -2.7 % | substitute |
| "a factor of 50 in absolute rupture time" | 49 | 183 | | 50 -> about 180 |
| "the knee regime costs 0.20--0.64" | 0.201--0.643 | 0.213--0.561 | | 0.21--0.56 |
| "the fourth ageing year 1.8--9.0 times" | 1.841--8.974 | 1.657--6.242 | | 1.7--6.2 |
| "23\,\% to the kinetics" | 23.3 | 8.5 | | 23 -> 9 |
| "9\,\% to the wall-temperature measurement" | 9.4 | 3.4 | | 9 -> 3 |
| "it holds in 0.57--0.63 of the samples" | 0.568--0.632 | 0.425--0.484 | | 0.43--0.48, and "close to a coin toss" -> "in fewer than half the samples" |

**Wording to assert.** The 90 % interval of $\log_{10} t_r$ widens from 1.69 to 2.26 decades, a factor
of about 180 rather than 50, because the measured G 4852 band is wider in relation to its own
constant. The creep group carries 86 % of that variance, so "contributes" becomes "dominates". The
ordering of robustness reverses: the knee conclusion (0.885) overtakes the ageing conclusion (0.819)
as the most robust of the study.


## `tab:uq` — RQ4 Monte Carlo intervals

The $T_{wo,max}$ column is bit-identical at every point. $\log_{10}t_r$ drops by about 0.85 decades everywhere and the paired life ratios move a little.

| v2 row | v3 row |
|---|---|
| `S1 base & 1145 [1127, 1169] & 7.91 [7.04, 8.73] & 1.000 [1.000, 1.000]` | `S1 base & 1145 [1127, 1169] & 7.06 [5.94, 8.20] & 1.000 [1.000, 1.000]` |
| `RQ3 knee & 1132 [1107, 1165] & 8.30 [7.33, 9.20] & 0.338 [0.201, 0.643]` | `RQ3 knee & 1128 [1103, 1161] & 7.47 [6.26, 8.63] & 0.329 [0.213, 0.561]` |
| `RQ3 min-life iso-H$_2$ & 1110 [1088, 1138] & 8.90 [7.94, 9.80] & 0.100 [0.069, 0.159]` | `RQ3 min-life iso-H$_2$ & 1110 [1088, 1138] & 7.90 [6.69, 9.06] & 0.143 [0.107, 0.212]` |
| `S5 year 4, hold slip & 1165 [1138, 1203] & 7.35 [6.28, 8.30] & 3.478 [1.841, 8.974]` | `S5 year 4, hold slip & 1165 [1138, 1203] & 6.60 [5.34, 7.79] & 2.823 [1.657, 6.242]` |
| `S/C 2.5 & 1146 [1129, 1166] & 7.88 [7.04, 8.69] & 1.153 [0.916, 1.350]` | `S/C 2.5 & 1146 [1129, 1166] & 7.03 [5.93, 8.18] & 1.140 [0.941, 1.315]` |
| `S/C 3.5 & 1144 [1127, 1172] & 7.94 [7.04, 8.83] & 0.846 [0.587, 1.119]` | `S/C 3.5 & 1144 [1127, 1172] & 7.11 [5.94, 8.25] & 0.857 [0.623, 1.083]` |
| `load 0.70 & 1138 [1123, 1156] & 8.08 [7.22, 8.90] & 0.966 [0.638, 1.350]` | `load 0.70 & 1138 [1123, 1156] & 7.20 [6.09, 8.34] & 1.034 [0.732, 1.384]` |
| `load 0.85 & 1141 [1124, 1162] & 8.00 [7.14, 8.82] & 0.951 [0.787, 1.143]` | `load 0.85 & 1141 [1124, 1162] & 7.13 [6.02, 8.28] & 0.989 [0.841, 1.170]` |

## `tab:robust` — robustness fractions

| v2 row | v3 row |
|---|---|
| `ageing year 4 costs $>2\times$ S1 life per kmol & 0.906 & 0.906` | `ageing year 4 costs $>2\times$ S1 life per kmol & 0.819 & 0.819` |
| `knee regime $<0.5\times$ base & 0.806 & 0.936` | `knee regime $<0.5\times$ base & 0.885 & 0.974` |
| `S/C 3.5 cheaper in life than S/C 2.5 & 0.836 & 0.836` | `S/C 3.5 cheaper in life than S/C 2.5 & 0.856 & 0.856` |
| `S2 daily $\le$ S1 & 0.568 & 0.000` | `S2 daily $\le$ S1 & 0.425 & 0.000` |
| `S3 renewable $\le$ S1 & 0.632 & 0.001` | `S3 renewable $\le$ S1 & 0.484 & 0.000` |

## `sec:discussion` — Industrial implications and Limitations

| Text | v2 | v3 | action |
|---|---|---|---|
| "the twin predicts the fourth year to cost 2.3 times the first per kmol" | 2.25 vs S1 | 1.96 vs S1 | 2.3 -> 2.0 |
| "(a factor of 2.3 against 2\,\%)" | 2.25 vs 1.7 % | 1.96 vs 1.3 % | a factor of 2.0 against 1 % |
| "this ageing conclusion is the most robust of the study, and the conclusion on load-following the least" | ageing 0.91 highest | ageing 0.82, knee 0.88 | **swap the ranking** |
| "the recommended steam-to-carbon ratio moves from 4.0 to about 2.7" | 2.66 | 2.64 | none |

The Limitations subsection is the one that changes most in substance. "The rupture curve is a
placeholder", "no public curve with scatter information was available", "the uncertainty analysis
compensates with a factor-of-two heat-to-heat scatter" and "replacing the curve by a fit to the NIMS
data sheets is the first step for any absolute statement" are all now false.

**Wording to assert.** The limitation is no longer that the curve is a placeholder but that it is a
manufacturer *design band*, not a heat-resolved data set: the printed band is a 95 % confidence level
on the manufacturer's own population, not the heat-to-heat spread of the tubes actually installed,
and the alloy of the Latham plant is still unknown. The absolute rupture time is still not credible
(now $10^7$ h rather than $10^8$ h). This is also where the invariance evidence belongs: the material
conclusions keep their direction on all four curves tested
(`data/creep_derived/alloy_comparison_v3.txt`).


## `sec:conclusions` — Conclusions

| Bullet | v2 | v3 | action |
|---|---|---|---|
| 1, "within 14~K in gas temperature" | 14.0 | unchanged | none |
| 1, "outlet within 4 K, wall within 12 K" | 4, 12 | unchanged | none |
| 2, "explain 82\,\% of the variance" | 82 | unchanged | none |
| 3, "load-following ... cost less than 2\,\% of the steady life per kmol" | 0.991 / 0.983 | 1.010 / 1.013 | **sign reverses** |
| 3, "short over-firing events at most 3\,\%" | 2.8 % | 1.8 % | may tighten to 2 % |
| 3, "steam-to-carbon ratio $\pm 17$\,\%" | +16.1 / -16.6 | +14.7 / -15.3 | $\pm 15$\,\% |
| 3, "catalyst ageing multiplies the life cost by 2.3" | 2.25 | 1.96 | 2.3 -> 2.0 |
| 4, knee "+16\,\% / $-18$\,\% / 0.28" and iso "0.083" | | see `sec:rq3` | +15 / -17 / 0.28 / 0.12 |
| 4, "moving from 4.0 to 2.7 when steam raising is charged in full" | 2.66 | 2.64 | none |
| 5, "spans a factor of 50, dominated by creep-curve scatter (62\,\%)" | 49 / 62 | 183 / 86 | about 180 / 86 |
| 5, "hold in 81\,\% to 91\,\% of the samples" | 0.81--0.91 | 0.82--0.88 | 82 % to 88 % |

The closing paragraph ("The next steps are the transcription of the NIMS rupture data ... into the
prepared ingestion path") must change: the manufacturer route has been taken, and the remaining NIMS
motivation is heat-to-heat scatter rather than the master curve.


## `tab:ranges` — supplement, operating space

**Unchanged.** Both versions sample `data/design/parameter_ranges_v2.csv`; v3 uses the same file.


## `tab:surrogate` — supplement, surrogate performance

Only the $\log_{10}t_r$ row moves; every thermal target is bit-identical.

| v2 row | v3 row |
|---|---|
| `log10\_t\_r\_hot & -- & gp\_\_base & 0.9999 & 0.00566 & 0.00597 & 0.959 & 0.0111` | `log10\_t\_r\_hot & -- & gp\_\_base & 0.9999 & 0.00494 & 0.00516 & 0.958 & 0.00956` |

The surrogate is marginally *more* accurate on the flatter G 4852 curve. No wording change needed.


## `tab:steam` — supplement, steam-credit sweep

Every `Rel. life cost` cell moves; S/C, load, firing, excess air and $T_{wo}$ move only through the re-optimisation.

| v2 row | v3 row |
|---|---|
| `0.00 & min life & 3.99 & 0.96 & 0.909 & 24.6 & 0.059 & +3.2 & 1103` | `0.00 & min life & 4.00 & 0.96 & 0.930 & 24.6 & 0.094 & +4.1 & 1103` |
| `0.00 & knee & 3.29 & 0.97 & 0.850 & 5.1 & 0.828 & -6.6 & 1144` | `0.00 & knee & 3.36 & 0.97 & 0.851 & 5.0 & 0.783 & -6.3 & 1142` |
| `0.00 & min fuel & 2.66 & 1.00 & 0.888 & 5.0 & 4.539 & -8.6 & 1171` | `0.00 & min fuel & 2.64 & 1.00 & 0.887 & 5.0 & 3.647 & -8.7 & 1172` |
| `0.25 & min life & 4.00 & 0.96 & 0.943 & 24.9 & 0.057 & +2.3 & 1103` | `0.25 & min life & 4.00 & 0.96 & 0.908 & 24.7 & 0.093 & +0.5 & 1103` |
| `0.25 & knee & 3.80 & 0.94 & 0.850 & 5.0 & 0.365 & -7.1 & 1131` | `0.25 & knee & 3.84 & 0.94 & 0.850 & 5.0 & 0.416 & -7.0 & 1130` |
| `0.25 & min fuel & 2.92 & 0.99 & 0.865 & 5.0 & 2.016 & -8.9 & 1158` | `0.25 & min fuel & 2.96 & 0.99 & 0.866 & 5.0 & 1.706 & -8.9 & 1157` |
| `0.50 & min life & 4.00 & 0.96 & 1.039 & 25.0 & 0.057 & +5.5 & 1103` | `0.50 & min life & 4.00 & 0.96 & 0.908 & 24.8 & 0.093 & -2.8 & 1103` |
| `0.50 & knee & 4.00 & 0.95 & 0.850 & 13.8 & 0.103 & -7.6 & 1112` | `0.50 & knee & 4.00 & 0.95 & 0.850 & 11.0 & 0.195 & -8.4 & 1116` |
| `0.50 & min fuel & 3.63 & 0.95 & 0.850 & 5.0 & 0.480 & -10.5 & 1135` | `0.50 & min fuel & 3.69 & 0.94 & 0.850 & 5.0 & 0.504 & -10.4 & 1134` |
| `0.75 & min life & 4.00 & 0.96 & 0.911 & 25.0 & 0.057 & -6.6 & 1103` | `0.75 & min life & 4.00 & 0.96 & 0.909 & 24.9 & 0.092 & -6.8 & 1103` |
| `0.75 & knee & 4.00 & 0.96 & 0.850 & 13.5 & 0.105 & -12.2 & 1112` | `0.75 & knee & 4.00 & 0.96 & 0.850 & 14.3 & 0.143 & -12.0 & 1111` |
| `0.75 & min fuel & 4.00 & 0.93 & 0.850 & 5.0 & 0.276 & -14.6 & 1127` | `0.75 & min fuel & 4.00 & 0.93 & 0.850 & 5.0 & 0.342 & -14.6 & 1127` |
| `1.00 & min life & 4.00 & 0.96 & 0.958 & 25.0 & 0.057 & -7.3 & 1103` | `1.00 & min life & 4.00 & 0.96 & 0.902 & 23.7 & 0.096 & -12.7 & 1104` |
| `1.00 & knee & 4.00 & 0.96 & 0.850 & 15.3 & 0.087 & -17.7 & 1109` | `1.00 & knee & 4.00 & 0.96 & 0.850 & 14.3 & 0.144 & -18.0 & 1111` |
| `1.00 & min fuel & 4.00 & 0.93 & 0.850 & 5.0 & 0.276 & -20.4 & 1127` | `1.00 & min fuel & 4.00 & 0.93 & 0.850 & 5.0 & 0.341 & -20.4 & 1127` |

## `tab:variance` — supplement, variance decomposition

The $T_{wo,max}$ columns are bit-identical; the $\log_{10}t_r$ columns change substantially.

| v2 row | v3 row |
|---|---|
| `kinetics & 0.233 & 0.225 & 0.513 & 0.620` | `kinetics & 0.085 & 0.157 & 0.513 & 0.620` |
| `heat transfer & 0.050 & 0.114 & 0.111 & 0.302` | `heat transfer & 0.018 & 0.080 & 0.111 & 0.302` |
| `creep & 0.622 & 0.655 & 0.211 & 0.093` | `creep & 0.859 & 0.779 & 0.211 & 0.093` |
| `measurement & 0.094 & 0.036 & 0.202 & 0.092` | `measurement & 0.034 & 0.025 & 0.202 & 0.092` |
| `sum / interaction residual & 0.999 / +0.001 &  & 1.037 / -0.037 &` | `sum / interaction residual & 0.997 / +0.003 &  & 1.037 / -0.037 &` |

Note that this manuscript's `tab:variance` prints only the row sums, whereas the archived draft
printed sum **and** interaction residual. The surrounding text still says the shares "sum to within
0.04 of one", which remains true in v3 (0.997 and 1.037).


## v3 UQ run: per-curve variance shares, interval and robustness fractions

The base column comes from `data/uq/uq_summary_v3.json`. The other three curves are added by
re-scoring the stored Monte Carlo draws by their quantiles rather than re-running the sampling, which
reproduces the base-alloy values to 1e-4 and exactly for the shares and fractions.


**90 % interval of $\log_{10}t_r$ at the S1 base point (`tab:uq`, first row)**

| Curve | median | p05 | p95 | width (decades) | factor in $t_r$ |
|---|---|---|---|---|---|
| G 4852 (base, C=18.6) | 7.062 | 5.938 | 8.200 | 2.263 | 183 |
| G 4852 Micro (C=22.9) | 7.817 | 6.567 | 9.074 | 2.508 | 322 |
| ET 45 Micro (C=19.3) | 7.164 | 5.955 | 8.376 | 2.421 | 264 |
| Yeh legacy (C=22.96) | 7.899 | 6.659 | 9.143 | 2.484 | 305 |

**Variance shares, $\log_{10}t_r$ (`tab:variance`)**

| Curve | kinetics | heat transfer | creep | measurement | sum | interaction residual |
|---|---|---|---|---|---|---|
| G 4852 (base, C=18.6) | 0.0854 | 0.0185 | 0.8586 | 0.0344 | 0.997 | +0.0032 |
| G 4852 Micro (C=22.9) | 0.1013 | 0.0219 | 0.8324 | 0.0408 | 0.996 | +0.0036 |
| ET 45 Micro (C=19.3) | 0.0805 | 0.0173 | 0.8664 | 0.0324 | 0.997 | +0.0034 |
| Yeh legacy (C=22.96) | 0.1024 | 0.0222 | 0.8309 | 0.0412 | 0.997 | +0.0033 |

**Variance shares, $T_{wo,max}$ (`tab:variance`)**

| Curve | kinetics | heat transfer | creep | measurement | sum | interaction residual |
|---|---|---|---|---|---|---|
| G 4852 (base, C=18.6) | 0.5128 | 0.1106 | 0.2112 | 0.2025 | 1.037 | -0.0371 |
| G 4852 Micro (C=22.9) | 0.5128 | 0.1106 | 0.2112 | 0.2025 | 1.037 | -0.0371 |
| ET 45 Micro (C=19.3) | 0.5128 | 0.1106 | 0.2112 | 0.2025 | 1.037 | -0.0371 |
| Yeh legacy (C=22.96) | 0.5128 | 0.1106 | 0.2112 | 0.2025 | 1.037 | -0.0371 |

**The $T_{wo,max}$ shares are identical to every printed digit on all four curves.** The creep group
touches the thermal problem only through wall thickness and tube conductivity; neither the
Larson--Miller constant nor the scatter enters it. Worth stating in `sec:rq4`: that decomposition is a
property of the furnace and tube model alone and carries over to any alloy.

For $\log_{10}t_r$ the creep group carries 83--87 % on every *measured* curve, so the v2 figure of
62 % was the outlier caused by the assumed 0.3-decade scatter.


**Robustness fractions, nominal $L_q$ (`tab:robust`)**

| Statement | G 4852 (base, C=18.6) | G 4852 Micro (C=22.9) | ET 45 Micro (C=19.3) | Yeh legacy (C=22.96) |
|---|---|---|---|---|
| ageing year 4 costs >2x S1 | 0.819 | 0.905 | 0.840 | 0.906 |
| knee regime <0.5x base | 0.885 | 0.938 | 0.910 | 0.936 |
| S/C 3.5 cheaper than S/C 2.5 | 0.856 | 0.839 | 0.858 | 0.833 |
| S2 daily <= S1 | 0.425 | 0.553 | 0.429 | 0.565 |
| S3 renewable <= S1 | 0.484 | 0.619 | 0.486 | 0.632 |

**Robustness fractions, $L_q \propto$ load$^{0.5}$ (`tab:robust`)**

| Statement | G 4852 (base, C=18.6) | G 4852 Micro (C=22.9) | ET 45 Micro (C=19.3) | Yeh legacy (C=22.96) |
|---|---|---|---|---|
| ageing year 4 costs >2x S1 | 0.819 | 0.905 | 0.840 | 0.906 |
| knee regime <0.5x base | 0.974 | 0.993 | 0.986 | 0.993 |
| S/C 3.5 cheaper than S/C 2.5 | 0.856 | 0.839 | 0.858 | 0.833 |
| S2 daily <= S1 | 0.000 | 0.000 | 0.000 | 0.000 |
| S3 renewable <= S1 | 0.000 | 0.000 | 0.000 | 0.001 |

The three material conclusions hold in 82--94 % of samples on every curve. Only the two flexibility
statements cross the 0.5 line, and they do so with the alloy: 0.43/0.48 on G 4852 and ET 45 Micro
against 0.55/0.62 on the two high-$C$ curves. That is the quantitative backing for reporting the
flexibility result as "no material difference" rather than as a direction.

