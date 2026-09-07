# Creep-derived data

Files here feed `rdt.creep` (hoop stress, Larson-Miller master curve, Robinson life fraction).

## `yeh2021_manaurite_xm.yaml`
Larson-Miller constant and rupture curves for Manaurite XM (Manoir Industries HP-Nb micro-alloy) as used by
Yeh (2021), *Appl. Sci.* 11:231 (`literature/yeh2021_applsci_11_231.pdf`). The master curve exists in the paper
only as Figure 3 (p. 11, reproduced from the Manoir data sheet, ref. [34]); the points stored here were extracted
by pixel analysis of the figure image (axis calibration checked against the 900/1000/1100 degC guide lines and
the 10 and 30 MPa gridlines) and are **approximate placeholders (TODO: replace by WebPlotDigitizer export or by
the Manoir data sheet)**. Yeh's design-point check (925 degC, 11.9 MPa -> about 5e6 h) is reproduced by the
*minimum* curve, which is therefore the default. This is a **placeholder alloy curve** until the NIMS sheets are
processed.

## NIMS Creep Data Sheets 16B (HK40-type) and 38A (HP40-type)
`nims_creep_data_schema.yaml` is an empty schema. NIMS terms: the sheets are CC BY-NC 4.0, obtained after free
registration on the MatNavi portal; **raw NIMS rupture data are not redistributed in this repository**. The
workflow is: download locally, fill the schema (heat_id, T_K, sigma_MPa, t_r_h, elongation, source), fit the
Larson-Miller master curve with `rdt.creep.LarsonMillerCurve`, and store **only the derived fit** (constant,
polynomial coefficients, fit range, number of points) in a new YAML file here, citing the DOI.

## Schmidt + Clemens Centralloy data sheets (`sources.yaml`, `datasheets/`)
`sources.yaml` describes the three manufacturer stress-rupture sources that replace the Yeh placeholder:
Centralloy **G 4852** (GX40NiCrSiNb35-25, HP-Nb), **G 4852 Micro** (GX45NiCrSiNbTi35-25, HP micro-alloy)
and **ET 45 Micro** (GX45NiCrSiNb45-35). All three are "September 2009, Rev. 02" sheets; the digitised
chart is page 6, *Parametric stress rupture strength*, which prints an *Average* and a *Lower Scatter
Band* curve (the latter stated on the sheet to represent the 95 % confidence level).

The Larson-Miller constant **differs between the sheets** and is read off each one individually:
G 4852 **C = 18.6**, G 4852 Micro **C = 22.9**, ET 45 Micro **C = 19.3**.

Two notes on the alloy identities:

* ET 45 Micro is **not** an HK40-type alloy. It is a 45Ni-35Cr micro-alloyed grade, far more highly
  alloyed than ACI HK40 (GX40CrNi25-20). It is kept as the high-alloy end member of the comparison
  only; for a genuine HK40 curve use NIMS Creep Data Sheet 16B.
* Only G 4852 Micro carries the C = 22.9 that is often quoted for "the 4852 family".

**The PDFs are copyright Schmidt + Clemens and are not redistributed here.**
`data/creep_derived/datasheets/` is git-ignored. To reproduce the digitisation, place the three files
named in `sources.yaml` there and run `python -m rdt.creep_ingest_datasheet --all`.

### Derived files
| File | Written by | Contents |
|---|---|---|
| `centralloy_g_4852_rupture_curve.csv` | `rdt.creep_ingest_datasheet` | Average and lower-scatter-band curves, LMP 25.32-32.39 |
| `centralloy_g_4852_micro_rupture_curve.csv` | `rdt.creep_ingest_datasheet` | LMP 31.33-39.29 |
| `centralloy_et_45_micro_rupture_curve.csv` | `rdt.creep_ingest_datasheet` | LMP 26.07-34.59 |
| `validation_v3.txt` | `rdt.creep_validate` | The three validation checks (see below) |

Columns: `alloy, curve (average|minimum), lmp, stress_mpa, source_file, page, method (vector|raster),
extracted_on`. The extractor has two independent paths - analytic Bezier flattening of the vector chart,
and a 600 dpi raster trace - which agree to 0.06 % median and 1.6 % worst case in stress on all six
curves. `rdt.creep.alloy_curves_from_csv` fits `log10(sigma) = P3(LMP)` to each curve with the alloy's
own constant, and derives the scatter from the **horizontal** separation of the two curves at constant
stress, `delta_log10_tr = (LMP_avg - LMP_min)/(scale T)`. At 1147 K that gives

| Alloy | Scatter over the digitised range | Median |
|---|---|---|
| G 4852 | 0.243-0.364 decades | 0.298 |
| G 4852 Micro | 0.427-0.558 decades | 0.461 |
| ET 45 Micro | 0.337-0.491 decades | 0.414 |

replacing the assumed constant 0.3 decades of v1/v2 (still selectable as `legacy_yeh_manaurite_xm`).

`python -m rdt.creep_validate` checks that the fitted lower-scatter-band curves reproduce the quoted
100 000 h strengths (18.3 MPa for G 4852 and 21.2 MPa for G 4852 Micro, both matched at about 930 C
rather than at 900 C), that `ingest_nims` still recovers C and scatter within 10 % on synthetic data,
and that the base operating point (1147 K, 12.9 MPa) is interpolation on every digitised curve.
`python -m rdt.alloy_comparison` re-scores the headline result ratios on each alloy.
