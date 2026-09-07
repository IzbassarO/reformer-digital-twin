# v3 rewrite: verification record

State after `09_v3_rewrite.md` sections 1--14, the move of the surrogate-parity and variance-share
figures to the supplement, and the v3 re-rendering of those two figures.

## Build

| | pages | undefined references | undefined citations | duplicate labels |
|---|---|---|---|---|
| `paper/main.tex` | 20 | 0 | 0 | 0 (45 labels) |
| `paper/supplementary.tex` | 9 | 0 | 0 | 0 (23 labels) |

Figures: 10 in the main text, 7 in the supplement; no orphaned and no missing figure files.

## texcount

| | body | headers | captions | tables |
|---|---|---|---|---|
| `main.tex` | **8158** | **128** | **475** | **248** |
| `supplementary.tex` | **1410** | **82** | **281** | **246** |

For reference, before this pass: main 7956 / 128 / 530 / 241, supplement 1410 / 82 / 225 / 246.
The main body grew by 202 words rather than shrinking by the ~400 the rewrite file projected,
because several replacement paragraphs are longer than the text they replace (`sec:rq2` and
`sec:rq4` in particular add the two-competing-effects explanation and the four-curve invariance
result). The projected saving assumed the §4.3 paragraph merge, which the rewrite file lists as a
separate step and which was not part of sections 1--14.

## Forbidden strings

Grepped in the compiled PDFs, not the sources.

| term | `main.pdf` | `supplementary.pdf` |
|---|---|---|
| placeholder | absent | absent |
| Manaurite | absent | absent |
| 22.96 | absent | absent |
| 0.3 decades | absent | absent |

All four are absent everywhere, including the bibliographies.

## Computed [CHECK] values

Both from `rdt.creep` at the base operating point (1147 K, 12.9 MPa), with the v3 configuration
active (`creep.use_config(config.V3)` -> `centralloy_g_4852`).

**(a) life-halving temperature interval** --- `creep.LarsonMillerCurve.life_halving_dT(1147.0, 12.9)`

| curve | log10 t_r | dT_half | check t_r(T+dT)/t_r |
|---|---|---|---|
| legacy Yeh, C = 22.96 | 7.883 | **11.19 K** | 0.5034 |
| G 4852, C = 18.6 | 7.059 | **13.46 K** | 0.5040 |

The legacy value reproduces the published "11 K", so the code path is the right one. Substituted in
`sec:intro` as 13 K.

**(b) average-versus-hot-spot overstatement** --- `creep.life_along_tube` on the calibrated Plant A
tube, hot-spot rupture time against the rupture time at the length-averaged outer-wall temperature
and mean hoop stress. **Not substituted.**

| definition of "tube-average" | legacy | G 4852 |
|---|---|---|
| outer wall, length-average (trapezoid) | 127.7 | 56.6 |
| mid-wall, length-average | 123.4 | 54.4 |
| outer wall, arithmetic mean of samples | 130.7 | 57.7 |
| mid-wall, arithmetic mean of samples | 125.8 | 55.3 |
| hot spot on outer wall vs mid-wall average | 408.9 | 149.1 |

No definition reproduces the published 110 on the legacy curve, so the quantity the sentence refers
to is not identifiable from the code, and substituting 56.6 would be a guess. The sentence is left
at "about 110" pending a decision on which quantity it should state.

## Open items

1. **The (b) factor above.**
2. **"life-consumption ratios of up to two orders of magnitude above the base"** in `sec:rq2`. This
   is pre-existing and is carried over unchanged by the rewrite, but the LHS campaign does not
   support it: the maximum ratio is 89.5 on the legacy curve and 43.8 on G 4852, with no run above
   100 on either.
3. **Campaign figures still at v2.** `fig07_operating_maps`, `fig08_sobol`, `fig11_scenarios_summary`,
   `fig12_pareto`, `fig13_steam_credit` and `fig14_uq_distributions` in the main text are still the
   v2 renderings and now disagree with the v3 text around them. v3 versions of all of them are in
   `paper/figures/v3/`.
4. **Length.** The main body is 8158 words. If the IJHE target is the ~7550 the rewrite file
   projects, the §4.3 paragraph merge and the §1.2 / §4.1 trims it mentions are still outstanding.
