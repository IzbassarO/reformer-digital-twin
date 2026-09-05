# Bed-side heat-transfer coefficient for the Xu & Froment (1989) reformer tube: a discrepancy and a decision

**Context.** The 1-D tube model (`rdt.reactor1d`) was compared with the digitised Fig. 3 of Xu & Froment
(1989, Part II, *AIChE J.* 35:97-103) using the digitised outer-wall temperature as boundary condition
(0-11.12 m heated, adiabatic to 12 m), the Table 2 feed, a bed voidage of 0.526 and an equivalent ring
diameter of 9.44 mm fitted to the pressure curve, and effectiveness factors of 0.1 (Latham et al. 2011).
With the outer wall prescribed, the model's heat flux is set by the wall-to-gas temperature difference and
the bed-side (inner film) coefficient `alpha_i`; the results are controlled almost entirely by `alpha_i`.

**Four estimates of `alpha_i` at the conditions of this tube** (Re_p ~ 4000-5000 on the particle diameter
and superficial mass velocity, Pr ~ 0.71):

| Estimate | alpha_i [W/(m2 K)] | Source |
|---|---|---|
| Leva & Grummer (1948) as written in Latham et al. (2011) Eq. 20, multiplier 1 | 940 (inlet) to 1580 (11 m) | `reactor1d.leva_grummer_alpha` |
| Xu & Froment Part II Eqs. 11-12 as printed (De Wasch & Froment alpha_w0, Yagi-Kunii dynamic terms, Kunii-Smith lambda_er0) | 3100 (inlet) to 3700 (11 m) | `reactor1d.xu_froment_alpha_i` |
| Back-calculated from the digitised Fig. 3: wall conduction from T_wall_outer - T_wall_inner with lambda_tube = 29.6 W/(m K) (Latham 2008, Table 17; Xu & Froment give no value), then q''_i/(T_wall_inner - T_gas) | median 348, interquartile range 313-363, full range 299-378 over 0.5-11 m | `validate_xf1989.alpha_i_profile_from_digitised` |
| Leva-Grummer multiplied by f_htg = 0.4, the value that reproduces Fig. 3 within 11 K and 0.023 in conversion | about 380 to 630 | `xf1989_fit.yaml`, `leva_grummer.diagnostics_f_htg_sweep` |

The printed correlations therefore disagree with the paper's own figure by roughly a factor of ten, and
with the two independent empirical estimates (back-calculation and fitted Leva-Grummer) that agree with
each other. Used as printed, they give a gas about 115 K hotter and a conversion 0.22 higher than the
figure (RMSE 114 K, 0.24); Leva-Grummer with multiplier 1 gives +73 K and +0.14.

**Why the printed correlations extrapolate.** Both dynamic terms are linear in Re Pr:
`alpha_w = alpha_w0 + 0.444 Re Pr lambda_g/d_p` and `lambda_er = lambda_er0 + 0.14 lambda_g Re Pr`. They
were fitted to laboratory packed-bed data at particle Reynolds numbers of order 10^2-10^3 (De Wasch &
Froment 1972; Yagi & Kunii 1957). At Re Pr ~ 3000-3500 the second term drives the effective radial
conductivity to about 50 W/(m K), against a static Kunii-Smith value of 0.9-1.9 W/(m K), and the
2-D-to-1-D lumping `1/alpha_i = 1/alpha_w + d_ti/(8 lambda_er)` then returns 3000-3700 W/(m2 K). An
effective conductivity of 50 W/(m K) is not supported by any measurement in a 0.1 m reformer tube.

**Misprint hypothesis (to be verified against De Wasch & Froment 1972 and Froment & Bischoff 1979).** A
single displaced decimal would resolve the order of magnitude. If the coefficient 0.444 in `alpha_w` were
0.0444, `alpha_i` at the inlet becomes about 1200-1400 W/(m2 K), still 3-4 times the back-calculated
value. If instead the 0.14 in `lambda_er` were 0.014, `lambda_er` becomes 5.7-6.7 W/(m K) and `alpha_i`
about 435-515 W/(m2 K), within 30 % of the back-calculation. The static-only limit (both dynamic terms
zero) gives 50-80 W/(m2 K), too low. The `lambda_er` coefficient is therefore the more likely location
of the error, but this is a hypothesis, not a finding, and neither the original notation (which lists
`lambda_g` in W and `lambda_st` in kJ per hour) nor the Reynolds definition can be excluded as the cause.

**Decision.**
1. The Xu & Froment Part II case is treated as a *verification* of the 1-D model with a matched film
   coefficient: `alpha_i = 348 W/(m2 K)` (the back-calculated median) with Latham's inlet effectiveness
   profile (0.05 in the first 10 % of the length, 0.1 elsewhere). This configuration reproduces Fig. 3 with
   RMSE 14.0 K in gas temperature, 16.7 K in inner-wall temperature and 0.015 in the conversion increment
   x(z) - x(0.5 m) (`xf1989_fit.yaml`, `verification_matched_alpha`). It verifies kinetics, thermodynamics,
   the energy balance and the wall conduction chain, not the heat-transfer correlation.
2. The bed-side coefficient itself is *validated* against the Latham (2008) plant data (Appendix H, four
   operating points with measured tube-wall temperatures at two elevations), where the furnace model
   supplies the outer-wall condition and `f_htg` is an estimable parameter, as in Latham et al. (2011).
3. Until then, `alpha_i` in the range **[300, 3000] W/(m2 K)** (back-calculated value to the printed
   correlation) is carried into the uncertainty analysis as an epistemic parameter, with Leva-Grummer
   (about 1000) as the nominal correlation.

**Resolution (plant calibration, 2026-09-04).** Calibrating the coupled furnace-tube model on the Latham
(2008) plant cases A, B and C2 (Appendix H; outlet temperature, pressure, wet composition, flue-gas outlet
temperature and the two absolute tube-wall temperatures per case; `latham_fit.yaml`, `calibration`) gives
`f_htg = 1.74` (95 % interval 1.43-2.05) on the Leva-Grummer correlation, i.e. a plant-constrained bed-side
coefficient of **alpha_i = 1410 W/(m2 K) at mid-height of the Plant A tube** (926-1540 W/(m2 K) along the
tube; overall inner-area coefficient U about 860 W/(m2 K) including the 15 mm wall). This is within 4 % of
Latham et al. (2011)'s own f_htg = 1.68 and reproduces the two measured tube-wall temperatures to within
11 K and 1 K, respectively, with leave-one-case-out errors of the same size. The value is four times the
348 W/(m2 K) back-calculated from Xu & Froment's Fig. 3 for a ring catalyst at Re_p of about 4500, and two
to three times below the printed Xu & Froment chain. The two plants differ in catalyst shape (quadralobe
5.4 mm equivalent diameter vs 17 mm rings), tube diameter and mass velocity, so a single universal value is
not expected; for this project the Latham plant value (with its interval) is the validated coefficient for
the industrial top-fired reformer, the Xu & Froment value remains a verification-only setting, and the
range [300, 3000] W/(m2 K) stays in the uncertainty analysis as the epistemic envelope.
