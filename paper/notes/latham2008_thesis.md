# Notes: Latham (2008), M.Sc.E. thesis, Queen's University

**Source.** D. Latham, *Mathematical Modelling of an Industrial Steam Methane Reformer*, M.Sc.(Eng.) thesis,
Dept. of Chemical Engineering, Queen's University, Kingston, December 2008. Supervisors K. B. McAuley and
B. A. Peppley; industrial contact T. Raybold (Praxair). Open access: <http://hdl.handle.net/1974/1650>.
Local copy: `literature/latham2008_msce_thesis.pdf` (279 PDF pages). BibTeX key `latham2008`.

**Page convention.** "p." = printed page number as it appears in the thesis. PDF page = printed page + 23
(front matter is roman-numbered; printed p. 1 is PDF page 24). The journal version is `latham2011`
(*Fuel Process. Technol.* 92:1574-1586), which is paywalled; the thesis contains more detail and all the data.

---

## 1. Furnace and tube geometry (Chapter 1.3, pp. 6-8; Figures 2 and 3, pp. 7-8)

| Item | Value | Where |
|---|---|---|
| Furnace type | Top-fired, co-current (process gas and flue gas both flow downward), rectangular radiant box designed by Selas Fluid Processing Corp. | p. 4, p. 6 |
| Tubes | 7 rows x 48 tubes = **336 tubes** (the number 336 is also used explicitly as N_tubes in the obstacle-zone balance, p. 80) | p. 6, p. 80 |
| Burners | 8 rows x 12 burners = **96 burners**, in the roof between tube rows; the two rows adjacent to the walls are fired at a lower rate because they see only one row of tubes | p. 6 |
| Tube outer diameter | **14.6 cm** | p. 6 |
| Tube heated (exposed) length | **12.5 m** | p. 6 |
| Tube inner diameter / wall thickness | **Not stated anywhere in the thesis.** The symbols r_in and r_out appear in the List of Symbols (PDF p. 21) and in the wall-conduction terms (eqs. 30, 32, 48, 50) but no numerical values are given. Must be assumed (typical for 14.6 cm OD HP-modified tubes: ~1.0-1.5 cm wall) or taken from the 2011 paper. | -- |
| Tube material | **Not stated.** Only the tube thermal conductivity is given: k_tube = 106 500 J/(m h K) = 29.6 W/(m K), from Davis (2000), *Alloy Digest Sourcebook: Stainless Steel* (Table 17, p. 93). Tube emissivity 0.85 (company experience, Table 17, p. 94). | p. 93-94 |
| Flame length | 4.5-6 m (design); the fitted heat-release length L_q is bounded 3.05-6.10 m and the best fit is 6.10 m | p. 6, p. 96, p. 121 |
| Flue-gas tunnels ("coffin boxes") | Rectangular intrusions on the furnace floor between tube rows, 2.86 m high, openings 0.6 m above the floor; run front-to-back | p. 6 |
| Refractory | Thickness 0.305 m; ceramic insulation k = 1153 J/(m h K); ambient 22 degC; design heat loss ~2% of fired duty, implemented as an overall U_refrac = 5483 J/(m2 h K) | pp. 96-98 |
| Refractory emissivity | 0.60 (company experience) | Table 17, p. 94 |
| Plant capacity | 2.83 million std m3/day H2 (120 Mmol/day) at 2413 kPa(g); 71 200 kg/h export steam at 390 degC, 4580 kPa(g) | p. 6 |
| Peep-hole elevations | Upper 3.66 m and lower 8.53 m below the top of the tubes (these are the only TWT measurement locations) | p. 9 |

### Catalyst (p. 10; Table 17, pp. 93-94; Table 18, p. 96)

| Item | Value | Where |
|---|---|---|
| Type | Nickel on alumina, **Johnson Matthey Katalco 23-4Q**, quadralobe (four-hole) pellets (Figure 4, p. 10, cross-section only, no dimensions) | p. 10, p. 94, p. 126 |
| Pellet dimensions | **Not stated.** Only an *equivalent particle diameter* D_p = **5.40 mm** "calculated according to Twigg (1989, p. 101)" is given (Table 17, p. 93). | p. 93 |
| Packing density | rho_cat = 1100 kg_cat / m3_tube (Katalco 23-4Q product bulletin) | Table 17, p. 94 |
| Bed voidage | Adjustable parameter phi, base value 0.7 +/- 0.1 (bounds 0.5-1.0); **best-fit value 0.607** (Table 22, p. 121); other adequate fits 0.609 and 0.614 (Tables 32-33, pp. 251-253). Note: the model is insensitive to phi (changing it by 0.1 gives no visible change in profiles, p. 105), so this value is weakly identified. | pp. 96, 105, 121 |
| Catalyst loading by zone | The thesis mentions operators load different catalyst types at different bed depths but deliberately does not model it (p. 13). | p. 13 |

---

## 2. Base-case operating conditions

Four reconciled data sets from third-party reformer surveys of three geometrically identical plants
(Plants A, B, C1, C2) are given **in tables** in Appendix H (pp. 239-242). These are directly usable.

### Process-side feed (Table 30, p. 241)

| | Plant A | Plant B | Plant C1 | Plant C2 |
|---|---|---|---|---|
| Total feed to all tubes, gmol/h | 7.886E+06 | 7.073E+06 | 9.654E+06 | 5.355E+06 |
| Feed per tube (/336), kmol/h (derived) | 23.5 | 21.1 | 28.7 | 15.9 |
| x_H2 | 0.0004 | 0.0018 | 0.0027 | 0.0025 |
| x_CH4 | 0.2421 | 0.2487 | 0.2458 | 0.2401 |
| x_C2H6 | 0.0042 | 0.0039 | 0.0041 | 0.0058 |
| x_C3H8 | 0.0009 | 0.0009 | 0.0008 | 0.0013 |
| x_nC4 / x_nC5 / x_nC6 | 0.0005 / 0.0002 / 0.0003 | 0.0005 / 0.0002 / 0.0003 | 0.0004 / 0.0002 / 0.0001 | 0.0006 / 0.0002 / 0.0002 |
| x_N2 | 0.0006 | 0.0007 | 0.0011 | 0.0014 |
| x_CO2 | 0.0047 | 0.0053 | 0.0052 | 0.0042 |
| x_H2O | 0.7462 | 0.7377 | 0.7395 | 0.7437 |
| Inlet temperature, degC | 611.4 | 613.9 | 617.8 | 606.7 |
| Inlet pressure, kPa (abs) | 3006.0 | 2944.0 | 3026.7 | 2808.5 |
| Steam-to-carbon (all hydrocarbon C), derived | 2.89 | 2.79 | 2.85 | 2.86 |
| Steam / CH4, derived | 3.08 | 2.97 | 3.01 | 3.10 |
| Plant rate, % of capacity (p. 125) | 99 | 92 | 96 | 71 |

### Furnace-side feed (Table 29, pp. 239-241)

Three furnace inlet streams per data set. Stream 1 is natural-gas fuel (~95% CH4 with C2-C6, N2, CO2) at ~16-23 degC.
Stream 2 is PSA off-gas (approx. 24-27% H2, 8-10% CO, 18-20% CH4, 46-48% CO2) at ~16-31 degC.
Stream 3 is combustion air (78% N2, 1-3% H2O) preheated to 300-338 degC. Furnace pressure 101.3-103.4 kPa.

| | Plant A | Plant B | Plant C1 | Plant C2 |
|---|---|---|---|---|
| Fuel gas, gmol/h | 1.253E+05 | 1.341E+05 | 1.784E+05 | 1.042E+05 |
| PSA off-gas, gmol/h | 2.835E+06 | 2.523E+06 | 2.540E+06 | 1.878E+06 |
| Combustion air, gmol/h | 1.079E+07 | 1.058E+07 | 9.654E+06 | 7.412E+06 |
| Air temperature, degC | 331.1 | 315.6 | 337.8 | 299.8 |

Note that the off-gas and air flow rates are considered unreliable by the plant: the model carries
adjustable multipliers f_nOffGas (best fit 0.963) and f_nCombAir (best fit 1.050) on them (Table 22, p. 121).

---

## 3. Model structure

### Kinetics and effectiveness factors (Section 2.6, pp. 32-33; Section 2.9, p. 42; Appendix E, pp. 215-218; Table 17, pp. 94-95)

- **Xu and Froment (1989a)** Langmuir-Hinshelwood rate expressions for the three reactions (SMR to CO, SMR to CO2, WGS),
  written in full in Appendix E, eqs. 224-228 (pp. 215-218), with the standard DEN term. Rate constants and adsorption
  constants use Xu-Froment activation energies (240.1, 243.9, 67.13 kJ/mol) and adsorption enthalpies
  (-82.9, -70.65, -38.2, +88.68 kJ/mol for H2, CO, CH4, H2O); pre-exponential factors listed in Table 17 (pp. 94-95).
- A single **adjustable multiplier f_prx on the reforming pre-exponential factors** was introduced (Table 18, p. 96);
  best-fit value **0.200** (Table 22, p. 121), i.e. the fitted intrinsic activity is one fifth of Xu-Froment. f_prx is
  one of the two least estimable parameters (Table 21, p. 116).
- Process side is a **1-D pseudo-homogeneous plug-flow model** (catalyst and gas at the same temperature, assumption 14,
  p. 62) represented as a series of CSTR zones (p. 42). Diffusion limitation is handled with **constant effectiveness
  factors** taken from Wesenberg and Svendsen (2007): eta_refrm = **0.05 in tube segment 1** and **0.1 in all other
  segments**, eta_WGS = 0.1 everywhere (Table 17, p. 94). In the zone-count study all were set to 0.1 (p. 98).
- Equilibrium constants from Gibbs energy with Reid et al. (1977) thermochemistry (p. 217).
- Six species only (H2, CO, CO2, CH4, N2, H2O). Higher alkanes (C2-C6) in the feed are converted instantaneously at the
  inlet by an overall "water-cracking" reaction (eq. 45, p. 78) to CH4, CO2 and H2, and the heat effect is added to the
  top process zone (pp. 76-79).
- Pressure drop: **Ergun** friction factor (eq. 229, p. 218); Hicks was tried and under-predicted the plant pressure drop (pp. 42-43).
- Tube-to-process-gas heat transfer: **Leva and Grummer (1948)** correlation multiplied by adjustable f_htg (eq. 17, p. 35;
  best fit **f_htg = 1.680**, Table 22, p. 121). f_htg is the most estimable parameter (Table 21, p. 116).
- Furnace-gas-to-tube convection: Dittus-Boelter (eq. 12) times adjustable f_ctube (best fit 0.381, p. 121).

### Furnace model (Section 2.9, pp. 37-42; Section 3.3, pp. 72-86; Appendix B, pp. 155-167)

- **Hottel zone method** with a weighted-sum-of-gray-gases furnace atmosphere: one clear gas plus **three gray gases**
  (Taylor and Foster 1974 coefficients; absorption coefficients K1 = 0.300, K2 = 3.10, K3 = 42.9 1/m; Table 17, pp. 93-94).
  Classified as "long furnace model type 2" (Rhine and Tucker 1991), p. 41.
- Total-exchange areas computed offline by Monte Carlo ray tracing in **RADEX** (Lawson and Ziesler 1996) with 200 000 rays/m2
  (p. 129); RADEX validated against Hottel-Sarofim charts within 4-8% and against an analytical view factor within 1%
  (Appendix G, Tables 27-28, pp. 235-237).
- Zoning: the furnace is sliced into **horizontal vertical sections**; each section has refractory surface zones, one
  furnace-gas volume zone, and obstacle zones for the tubes (all 336 tubes at one elevation lumped as one "average tube"
  zone) and for the coffin boxes (Figures 7-9, pp. 52-53). With 10 sections there are 38 surface zones, 13 obstacle zones,
  10 volume zones -> 151 equations (Tables 13-14, pp. 73-75). The **final model uses 15 non-uniform sections: 10 in the top
  half and 5 in the bottom half**, giving **226 equations / 226 unknowns** (p. 130; Figure 24, p. 105). 20 and 40 uniform
  sections gave 299 and 594 equations with no significant benefit beyond 20 (p. 130).
- Furnace gas in perfect plug flow downward (assumption 5, p. 61); uniform composition after assumed complete
  isothermal pre-combustion of the fuel (assumption 2). Heat of combustion is distributed over the volume zones by a
  **discrete downward-opening parabola** over the heat-release length L_q with fraction alpha_top in the top zone
  (eqs. 36-42, pp. 66-68). L_q (3.05-6.10 m) and alpha_top (best fit 0.182) are fitted (Tables 18 and 22).
- Solved as a set of nonlinear algebraic equations by Newton-Raphson with numerical Jacobian (Appendix D, pp. 203-208).
  Runtime ~134 s for 10 sections, <4 min for 15 sections, on a 2.8 GHz Pentium 4 (p. 73, p. 128).

### How tube-wall temperatures are computed (eqs. 48 and 50, pp. 80-81)

- The **outer-wall temperature** of each vertical section is the temperature of the tube *obstacle zone* in the Hottel
  balance: radiation in from all zones + convection from furnace gas = conduction through the wall (eq. 48).
- The **inner-wall temperature** comes from the inner-tube-surface balance (eq. 50): cylindrical conduction
  2*pi*k_tube*N_tubes*dy*(T_out - T_in)/ln(r_out/r_in) equals convection to the process gas
  2*pi*r_in*N_tubes*dy*f_htg*h_tg*(T_in - T_gas).
- Thus TWT is an **average over all 336 tubes at a given elevation** (one radiative environment per elevation). Wall tubes
  and "gap tubes" (four per row) are known to run hotter but are not modelled; this is the first recommendation for
  future work (Section 5.2, Figures 43-44, pp. 135-137).

---

## 4. Where the MEASURED plant data are

All plant measurements in the thesis come from the quarterly third-party surveys (reconciled by plant mass/energy
balances, p. 89); the nightly hand-held pyrometer readings and hourly historian data are described (p. 9) but **not
tabulated**.

| Data | Location | Form | Directly usable? |
|---|---|---|---|
| **Process-side inlet conditions** (composition, flow, T, P) for Plants A, B, C1, C2 | Table 30, p. 241 (PDF 264) | Table | Yes |
| **Furnace-side inlet streams** (fuel, off-gas, air: composition, flow, T, P) | Table 29, pp. 239-241 (PDF 262-264) | Table | Yes |
| **Measured outputs**: process outlet T, P, molar flow, wet outlet composition (6 species), **upper and lower peep-hole TWT**, flue-gas outlet T | Table 31, p. 242 (PDF 265) | Table | Yes. Values: T_proc,out 832.4 / 834.0 / 848.9 / 835.8 degC; T_upper 822.8 / 809.2 / 838.0 / 802.5 degC; T_lower 858.9 / 857.7 / 878.4 / 853.1 degC; T_fur,out 1012.3 / 998.1 / 1005.9 / 940.9 degC; x_CH4,out 0.0505 / 0.0526 / 0.0454 / 0.0482 (wet) for A / B / C1 / C2 |
| Measurement uncertainties used as weights | Table 15, pp. 89-90 | Table | Yes (T_proc,out +/-2, peep-hole +/-3, T_fur,out +/-8 degC; P +/-36 kPa; x_H2 +/-0.01, x_CH4 +/-0.003) |
| **Model vs. plant, outputs not in the figures** (P_out, n_out, all six outlet fractions) | Table 23, pp. 124-125 | Table | Yes |
| **Full simulated profiles** (furnace gas, outer wall, inner wall, process gas vs. elevation, 0-12.5 m) with the 3 measured points overlaid, best-fit parameters, Plants A, B, C1, C2 | Figures 39-42, pp. 122-124 | Plots only | Needs digitizing. Some data labels are printed on the plots (e.g. Plant A: 832, 823, 859, 1012 / 1008, 867, 811, 835) |
| Same, unfitted (base-case) parameters | Figures 63-70, pp. 246-250 | Plots | Digitize |
| Same, alternative adequate parameter sets | Figures 71-78, pp. 251-255; parameter values Tables 32-33 | Plots + tables | Digitize |
| Base-case profile with Plant B inputs | Figure 25, p. 106 | Plot | Digitize |
| Sensitivity of profiles to each parameter (L_q, alpha_top, f_htg, f_prx, f_nOffGas, f_nCombAir, f_ctube) | Figures 26-33, pp. 107-114 | Plots | Digitize if needed |
| Zone-count convergence (10/20/40/15 sections) | Figures 18-24, pp. 99-105; runtimes Table 19, p. 104 | Plots + table | -- |

**There are no measured TWT profiles along the tube**, only two elevations (3.66 and 8.53 m). There are no measured
flue-gas temperature profiles inside the box, only the flue-gas outlet temperature. Simulated furnace-gas profiles
peak at roughly 1300-1400 degC near 1-3 m depth in Figures 39-42.

---

## 5. Stated model-vs-plant errors (Section 4.6, pp. 116-125; Conclusions, pp. 130-131; Abstract, p. ii)

- Best fit: heat-release length **6.10 m with 7 parameters** (Table 22, p. 121): f_htg 1.680, alpha_top 0.182,
  f_nOffGas 0.963, phi 0.607, f_nCombAir 1.050, f_prx 0.200, f_ctube 0.381. Objective J = 59.3 (Figure 38, p. 120).
- The absolute minimum J = 39.2 at L_q = 4.88 m was **rejected** because it drove f_ctube to 0.0004 (no furnace convection), p. 120.
- With the accepted best fit: **process outlet temperature within 4 degC, flue-gas outlet within 4 degC, upper and
  lower peep-hole TWT within 12 degC**; outlet pressure, composition and flow "accurately matched" (p. 130, Abstract).
- Any parameter set with J < 80 was deemed adequate; J < 80 corresponds to **7 degC on process outlet T, 27 degC on
  flue-gas outlet T, 12 degC on peep-hole TWT, 0.0044 on wet outlet CH4 fraction** (p. 121). Adequate fits exist at
  4.88 m (>= 4 parameters), 5.49 m (7 parameters) and 6.10 m (>= 5 parameters) (p. 125).
- Table 23 (pp. 124-125): worst outlet-pressure miss 89 kPa (Plant B), worst outlet CH4 miss 0.0106 (Plant B),
  worst molar flow miss ~3% (Plant B).
- Author's own caveat (p. 125): TWT data at more elevations, furnace-gas temperatures and wall temperatures would be
  needed to pin down the parameters, and the data span only 71-99% of capacity.

---

## 6. Relevance to our digital twin

1. **Complete, tabulated validation case.** Appendix H gives four fully specified inlet states and the matching measured
   outlet states of a 336-tube industrial top-fired reformer. This is the only open, tabulated industrial SMR data set
   in our bibliography and should become `data/literature_validation/latham2008_*.csv`.
2. **Missing geometry** (tube ID / wall thickness, pellet dimensions, tube alloy) must be assumed and documented; check
   `latham2011` when obtained through the library.
3. Latham's fitted **effective kinetics** (0.2 x Xu-Froment with eta = 0.05-0.1) and **Leva-Grummer x 1.68** give a
   reference point for our own parameter estimation; note the strong parameter correlations he reports.
4. His model produces only one average TWT per elevation and no uncertainty; both are gaps we address.
