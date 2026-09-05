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
