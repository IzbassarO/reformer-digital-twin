# Notes: Yeh (2021), Applied Sciences 11(1):231

**Source.** C.-L. Yeh, "Effect of Burner Operation on the Catalyst Tube Lifetime of a Steam Methane Reformer:
A Numerical Study", *Applied Sciences* 11(1), 231 (published 29 Dec 2020), doi:10.3390/app11010231, CC BY 4.0.
Local copy: `literature/yeh2021_applsci_11_231.pdf` (20 pages). BibTeX key `yeh2021`.
Data and plant: Formosa Petrochemical Corporation, Taiwan (Acknowledgments, p. 19). Single author,
Dept. of Aeronautical Engineering, National Formosa University.

**Page convention.** Page numbers are the journal page numbers (1-20), which coincide with PDF page numbers.

---

## 1. What they did

### Simulation tool and physical models (Section 2.1, pp. 3-5)

- **ANSYS FLUENT v17**, 3-D steady RANS: standard k-epsilon turbulence, **discrete-ordinates (DO)** radiation,
  **finite-rate / eddy-dissipation (FRED)** combustion and reforming chemistry, standard wall functions, catalyst bed as a
  **porous zone** (p. 3).
- Reforming kinetics: Arrhenius constants of **Odegard, Johnsen and Karoliussen (1995)**, an SOFC-anode (Ni/zirconia)
  kinetic set, chosen because its temperature range matches the reformer (p. 5). Not Xu-Froment; no effectiveness factor
  or intraparticle diffusion is mentioned.

### Geometry and boundary conditions (Section 2.2, pp. 5-8; Figure 1, pp. 6-7)

- Reformer with **276 catalyst tubes and 432 burners**; half the box simulated by symmetry (138 tubes, 216 burners).
  Tube **OD 136 mm, wall thickness 13.4 mm**; burner diameter 197 mm (p. 5). Box length 37.68 m in the flow direction,
  divided into 6 groups of 23 tubes and 36 burners each (Figure 1d, p. 7; group ranges p. 13). Tube heated length and
  inner catalyst geometry are not stated. Burner flow is described as radial and tubes as axial (p. 7-8), and the
  furnace has a divergent roof with the flue outlet at one end, so hot gas convects toward the "downstream" end
  (pp. 11-12). This is a **side/terrace-fired-type layout, not a top-fired box like Latham's**.
- Burner inlet (premixed fuel + flue gas + air): 139 710 m3/h, 673.15 K, 10.45 kPa(g); mole fractions
  H2 0.0816, CH4 0.0474, N2 0.4906, O2 0.1282, CO2 0.2523 (p. 7-8).
- Tube inlet: 24 740 m3/h total, **912.75 K (639.6 degC)**, **2.1658 MPa(g)**; CH4 0.2029, H2O 0.60, H2 0.1286,
  CO2 0.0657, CO 0.00145, N2 0.00145 (p. 8). Implied steam-to-carbon = 0.60/0.2029 = **2.96**; note the unusually high
  H2 (12.9%) and CO2 (6.6%) in the feed, suggesting a pre-reformed or recycle-containing feed.
- Mesh: ~4 million cells (full prototype), ~1 million (one-group periodic model), GAMBIT; y* 20-60 (p. 8).

### Validation (Section 3.1, pp. 8-10; Figure 2, p. 9; Table 1, p. 10)

- Plant data: infrared-thermography **average outer-surface tube temperatures** and **outlet hydrogen yields**
  (H2 mole fraction at tube outlet) per tube, plus average front-wall and back-wall temperatures (923 and 934 degC).
- Full ("prototype") model: TWT error **within 4%**, H2 yield 0.67 vs measured ~0.69 (**within 3%**), walls 891/892 degC
  (3.5% / 4.5% error). Periodic single-group model: TWT far too low, H2 yield 0.13, wall errors ~23% (Table 1, p. 10).
  Conclusion: periodic boundaries are unusable for this furnace; all results use the full model.
- The measured per-tube TWT and H2 yield exist only as **plots** (Figure 2a,b, p. 9); no table.

### Tube life computation (Section 3.2, pp. 10-12)

- **Larson-Miller parameter**: LMP = T (C + log10 t_r) x 1e-3 with T in K, t_r in h, **C = 22.96** "as suggested by the
  manufacturer" (eq. 14-15, p. 10). The LMP-vs-stress master curve is **Figure 3 (p. 11), reproduced from the Manoir
  Industries data sheet for Manaurite XM** (ref. [34], p. 20), i.e. the alloy is **Manaurite XM (HP-Nb micro-alloyed,
  25Cr-35Ni-Nb-Ti)**. The curve itself is only a plot; no equation or table of LMP(stress) is given.
- **Stress**: thin-wall hoop formula S = 0.5 P D_o / t (eq. 16, p. 10) -> **11.9 MPa** at the design pressure
  24 kg/cm2 (2.35 MPa) with D_o 136 mm, t 13.4 mm. Mean-diameter or thick-wall (Lame) stress is not used; no thermal
  stress; no stress redistribution or creep-relaxation.
- **Design condition**: 925 degC and 24 kg/cm2 -> **design life ~5e6 h (~570 years)**, described as including a
  safety factor (p. 11). (For comparison, the industry design basis is 100 000 h; Latham 2008 p. 1.)
- **Simulated life**: for each tube, LMP evaluated at (a) the average and (b) the **maximum** outer-surface temperature
  and inner pressure from CFD (Figure 4, p. 12), giving per-tube lifetimes in **Figure 5 (p. 12)**. Life from average
  T/P is "the same order of magnitude as the design lifetime"; life from **maximum T/P shows upstream tubes rupturing
  at ~4 years** (p. 11). Periodic model would give ~1e12 years (p. 11), demonstrating the sensitivity. The T and P
  values behind Figures 4-5 are not tabulated.
- No time-integration of damage (life is a single steady-state number per tube), no life-fraction rule, no
  temperature history or start-up/shutdown cycling.

### TWT ranges obtained

- Only plotted (Figures 4a, 6, 7; pp. 12-16). Text gives: measured and simulated average outer TWTs agree within 4% of
  values comparable to the 891-934 degC wall temperatures (p. 9-10); design TWT 925 degC (p. 10). Turning off
  **one, two or three burner groups lowers average outer TWT by about 106, 163 and 215 degC** respectively (p. 13).
  TWTs are higher at the lower part of the tubes (weaker convection near the floor) and in the upstream part of the box (p. 14).

### Scenarios compared (Section 3.3, pp. 12-18; Table 2, p. 17; Figures 6-10)

Twelve steady-state burner on/off combinations: all on; one group off (1-6); two groups off (1&2, 3&4, 5&6);
three groups off (1&2&3, 4&5&6).

| Scenario | Average H2 yield (Table 2, p. 17) |
|---|---|
| All burners on | 0.670 |
| Group 1 / 2 / 3 / 4 / 5 / 6 off | 0.641 / 0.647 / 0.655 / 0.654 / 0.652 / 0.644 |
| Groups 1&2 / 3&4 / 5&6 off | 0.616 / 0.637 / 0.629 |
| Groups 1,2,3 / 4,5,6 off | 0.602 / 0.604 |

- H2 yield drops about **3.2%, 6.4%, 10%** for one, two, three groups off (p. 14).
- Lifetimes per tube per scenario are in **Figure 10 (p. 18), plot only**; no numbers are quoted except the ~4 year
  upstream rupture for the all-on case. Qualitative result: switching off **upstream** burners helps life most,
  switching off **central** groups helps least; the price is H2 yield (pp. 17-19).
- The life-vs-yield trade-off is presented scenario by scenario, **not as a Pareto front or optimisation**.

### Stated limitations (by the author)

- Periodic boundary conditions badly under-predict temperatures and pressures and would over-estimate life (Abstract; pp. 9-11, 19).
- Design-basis life is unrealistically long because "operational, human and environmental factors" shorten it (p. 11).
- Mesh size was capped by workstation memory (~4 M cells on a Core i7-8700, 64 GB) (p. 8).
- Life estimate is "preliminary" (word used for the stress estimate, p. 10) and considers creep only, not other failure
  mechanisms listed on p. 10.
- Data "available on request" only (Data Availability Statement, p. 19).

---

## 2. What they did NOT do (confirmed from the text)

| Item in our plan | Yeh (2021) | Evidence |
|---|---|---|
| Uncertainty quantification | **None.** Deterministic CFD, single LMP curve, single C constant; no error bars, no sensitivity study on C, alloy scatter, emissivity or kinetics. The only "uncertainty" discussion is the periodic-vs-full model comparison. | pp. 8-12 |
| ML surrogate / reduced model | **None.** Every scenario is a full 4 M-cell CFD run; no surrogate, no regression, no fast model for operations. | Sections 2-3 |
| Multi-objective / Pareto optimisation | **None.** Twelve hand-picked on/off scenarios; trade-off stated qualitatively. No optimiser, no Pareto set, no continuous decision variables (firing rate, S/C, throughput). | Section 3.3 |
| Open data | **None.** Plant data are proprietary (Formosa Petrochemical); results only in plots; "available on request". Only the inlet BCs (p. 7-8), Table 1 and Table 2 are numeric. | p. 19 |
| Time-dependent damage / life-fraction accumulation | None. Steady-state LMP life only. | Section 3.2 |
| Thick-wall or thermal stress, alloy-specific LMP scatter bands | None. Thin-wall hoop stress, single manufacturer curve. | p. 10 |
| Coupling to a 1-D process model for on-line use | None. CFD only. | -- |
| Physics-informed kinetics for industrial Ni catalyst | Used SOFC-anode kinetics (Odegard et al. 1995), no Xu-Froment, no effectiveness factor. | p. 5 |
| Validation against axial TWT profiles | Only tube-averaged outer TWT and outlet H2 per tube (IR thermography). | Figure 2, p. 9 |

---

## 3. Useful take-aways for our project

1. **Reusable creep-life recipe and numbers**: Manaurite XM, LMP with C = 22.96, S = 0.5 P D_o/t = 11.9 MPa at
   925 degC / 24 kg/cm2, design life ~5e6 h; and the demonstration that using maximum instead of average TWT changes
   life from ~centuries to ~4 years. We should reproduce this with our own LMP curve (NIMS HP40 data, `nims1991`) and
   report the spread as uncertainty.
2. **Burner-group trade-off data**: Table 2 gives a clean 12-point yield table that our surrogate/Pareto study can be
   compared against qualitatively (yield loss of ~3% per group turned off, upstream groups matter most).
3. **A concrete validation benchmark for CFD-level accuracy**: 4% on TWT, 3% on H2 yield, 3.5-4.5% on wall
   temperatures. Our 1-D model plus surrogate should aim for the same order while running in seconds instead of
   million-cell CFD.
4. Full inlet boundary conditions (p. 7-8) and tube dimensions (136 x 13.4 mm) are numeric and can be used as a
   second, side-fired-style test case, although the tube heated length is missing.
