# Literature validation data

Curated, hand-transcribed data from the open and library-accessed literature used to validate the
`rdt` 1-D reformer model. Datasets in this directory are CC BY 4.0 (see `../LICENSE.md`); the
underlying papers remain the copyright of their publishers and are not redistributed here.

Rules for every file: give the source table/figure/page for every block, never invent a number,
write `null` with a comment when the source does not give a value, and mark derived quantities as
derived.

## Files

| File | Source | Content | Status |
|---|---|---|---|
| `xu_froment_1989_industrial.yaml` | Xu & Froment (1989) Part II, *AIChE J.* 35:97-103, Table 2 (p. 100) and Figs. 3-7 (pp. 101-102). Key `xu1989b`. | Single-tube industrial **simulation** case: inlet 793.15 K, 29 bar, feed per tube, natural-gas analysis, tube ID/OD 0.1016/0.1322 m, heated length 11.12 m, ring catalyst dimensions, tortuosity 3.54. Result profiles exist only as plots. | Table 2 verified from page image. Figures need digitising. |
| `latham2011_cases.yaml` | Latham et al. (2011), *Fuel Process. Technol.* 92:1574-1586. Key `latham2011`. | Table 1 specs, Table 8 inputs (99 % and 75 % cases), Tables 5/9 tube-wall statistics (**shifted**), Tables 7/10 sim-minus-measured, Tables 4/6 adjustable parameters, model assumptions, and a cross-check against the thesis. | Transcribed from text layer. |
| (planned) `latham2008_appendixH.csv` | Latham (2008) thesis Appendix H, Tables 29-31, pp. 239-242 | Four absolute inlet/outlet states (Plants A, B, C1, C2). | Not yet transcribed; see `paper/notes/latham2008_thesis.md`. |

Kinetic constants (Xu & Froment Part I) live in `../kinetics/xu_froment_1989.yaml`, checked by
`tests/test_kinetics_data.py`.

## Provenance notes

- **Xu & Froment scans.** Both AIChE papers are scanned PDFs with an OCR text layer that garbles
  tables and equations. All numbers were read from the page images (`pypdf` image extraction,
  cropped and inspected), not from the OCR text.
- **K_H2O lower confidence limit.** Table 5 (Part I, p. 94) prints LL = 0.0317 for K_H2O,823.
  Every other 95 % interval in the table is symmetric about the estimate, which would give 0.3272.
  The printed value is recorded; treat it as suspect until checked against another edition.
- **Equivalent CH4 feed.** Part II Table 2 gives 5.168 kmol/h and H2O/CH4 = 3.358, but 135.00
  Nm3/h at 0 degC / 1 atm with the listed composition gives about 5.4 kmol/h of carbon-equivalent
  methane and 399.17 Nm3/h steam gives about 17.8 kmol/h, ratio about 3.45. The reference
  conditions for "Nm3" are not stated. Recorded as printed.
- **Latham 2011 temperatures are shifted.** Page 1579: "All temperatures in the publication have
  been vertically shifted to show the shape of the temperature profiles without indicating the
  true value of the temperatures." Only temperature *differences* from the paper are usable.
- **Thesis temperatures are absolute.** The thesis never mentions shifting; its Table 31 values are
  physically plausible; the one absolute composition the paper reveals (46.3 mol% wet H2 at 91 %)
  matches thesis Plant B exactly; and the paper-vs-thesis tube-wall differences at both elevations
  and both rates are reconciled by a single constant of about 525 degC. Absolute tube-wall
  validation must therefore use the thesis data, not the paper.
- **Paper vs thesis inlet states.** Paper Table 8 (99 %, 75 %) are new validation sets and match
  none of the thesis Appendix H states. The paper quotes fitting rates 99/94/91/75 %, the thesis
  99/92/96/71 %.
- **Tube wall thickness** is not given in the thesis or the 2011 paper. The Xu & Froment Part II
  tube (ID 0.1016 m, OD 0.1322 m, wall about 15 mm) is recorded as a documented reference value only.

## `latham2008_plant_cases.csv` (added 2026-09-04)

One row per plant case of Latham (2008), M.Sc.E. thesis, Appendix H: Table 29 furnace-side inputs
(pp. 239-241), Table 30 process-side inputs (p. 241), Table 31 measured outputs (p. 242); plant rates
99/92/96/71 % from p. 125. Built by a script from the thesis text; the Plant B anchor (wet H2 outlet
0.463134 = 46.3 mol%, matching Latham et al. 2011 Table 7) is reproduced.

Conventions and conversions:
- Temperatures: thesis degC + 273.15 -> K. **All temperatures are absolute** (unlike the 2011 paper).
- Pressures: kPa x 1000 -> Pa. Flows: gmol/h / 1000 -> kmol/h; per tube = total / 336.
- Compositions: mole fraction x 100 -> mol%. Feed C6plus = thesis n-C6; iso-C4, iso-C5, neo-C5 are zero in all cases.
- Tube-wall temperatures: peep holes 3.66 m (upper) and 8.53 m (lower) from the top of the tubes (p. 9);
  fractions 0.293 and 0.682 of the **exposed (heated) tube length 12.5 m** (p. 6), which is what "height" means here.
- Uncertainties (`sigma_*`) are the "uncertainty in value" of thesis Table 15 (pp. 89-90), used as one standard
  deviation: +/-2 degC outlet T, +/-3 degC tube wall, +/-8 degC flue-gas outlet, +/-36 kPa, +/-227 kgmol/h,
  +/-0.01 H2, +/-0.003 CH4 mole fraction. The thesis gives no per-tube spread; the 2011 paper's standard
  deviations refer to shifted data of different cases.
- Columns `*_given` for heat of combustion, flue-gas flow and excess air are **empty**: the thesis does not state
  them. Columns `*_derived` are computed from the Table 29 streams: complete combustion of fuel gas + PSA off-gas
  in the air stream (O2 in air by difference from the listed N2, CO2, H2O), LHV at 298.15 K with gri30
  enthalpies (CH4, C2H6, C3H8, H2, CO) and tabulated values for n-C4/C5/C6 (2657.3, 3244.9, 3855.1 kJ/mol);
  the mixed inlet temperature is the enthalpy-weighted mixing temperature of the three unburnt streams
  (C4-C6 lumped into C3H8 for cp).
- Known oddities, not corrected: Plant C1 process feed (9.654E+06 gmol/h) equals its air flow to four digits,
  suspicious; thesis Table 23 (p. 124) has the Plant B / C1 columns for P_out and n_out swapped relative to
  Table 31; Table 31 (Appendix H) is used.
