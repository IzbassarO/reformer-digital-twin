# Reconciliation: repository sectioned draft vs consolidated Overleaf pair

Compared without modifying either version.

| | repository draft | Overleaf pair |
|---|---|---|
| Entry point | `paper/main.tex` + `paper/sections/*.tex` + `paper/tables/*.tex` + `paper/numbers.tex` | `chem/main.tex` (single file) |
| Class options | `[preprint,12pt]`, `lineno` active | `[final,3p,times]`, `lineno` commented out |
| Bibliography | `\bibliography{references}` + `main.bbl` | `thebibliography` embedded, no `.bib` |
| Numbers | macros from `paper/numbers.tex` (auto-generated) | hard-coded in the text |
| Length | 37 pages | ~22 pages |
| Tables in main text | 11 | 8 |
| Tables in supplement | 1 | 4 |
| Figures in main text | 12 | 12 |
| Figures in supplement | 5 (unlabelled) | 5 (all labelled) |
| TODOs | 44 (40 in the main document, 4 in the supplement) | 3 (all in `main.tex`) |

Both versions are built on the **v2** result files. The Overleaf pair does **not** contain the
Centralloy G 4852 (v3) numbers.

---

## 1. Section structure

The numbered section tree is **identical** in both: 1 Introduction (1.1 Previous work,
1.2 Contribution and research questions), 2 Methods (2.1--2.7), 3 Results (3.1--3.5),
4 Discussion (4.1--4.3), 5 Conclusions, followed by the unnumbered back matter.

The only structural difference in the back matter is that the repository draft has a **Funding**
section that the Overleaf pair does not — that section is one of my earlier frontmatter edits and
is re-applied to the Overleaf file in step 3.

The supplements are organised differently:

| repository `supplementary.tex` | Overleaf `supplementary.tex` |
|---|---|
| 1 Digitised verification data | 1 Digitisation of the verification curves |
| 2 Additional figures | 2 Operating space and campaign design |
| 3 Operating space | 3 Surrogate performance |
| 4 Version note | 4 Steam-credit sweep |
| | 5 Variance decomposition |
| | 6 Supplementary figures |

The repository supplement's "Version note" (the v1-vs-v2 comparison paragraph) has **no counterpart**
in the Overleaf supplement.

---

## 2. Where each table lives

Keyed to the LaTeX label, since the printed numbers differ between versions.

| Content | repository draft | Overleaf | moved? |
|---|---|---|---|
| Xu--Froment verification | main, `tab:xf` | main, `tab:xf` | -- |
| Latham plant cases | main, `tab:latham` | main, `tab:latham` | -- |
| Calibrated parameters | main, `tab:cal` | main, **`tab:params`** | renamed |
| Sobol indices | main, `tab:sobol` | main, `tab:sobol` | -- |
| Best surrogate per target | main, `tab:surr` | **supplement, `tab:surrogate`** | **moved to supplement** |
| RQ2 scenario summary | main, `tab:scen` | main, **`tab:scenarios`** | renamed |
| RQ3 regimes | main, `tab:regimes` | main, `tab:regimes` | -- |
| Steam-credit sweep | main, `tab:steam` | **supplement, `tab:steam`** | **moved to supplement** |
| RQ4 Monte Carlo intervals | main, `tab:uq` | main, `tab:uq` | -- |
| Variance decomposition | main, `tab:var` | **supplement, `tab:variance`** | **moved to supplement** |
| Robustness fractions | main, `tab:robust` | main, `tab:robust` | -- |
| Operating ranges | supplement, `tab:ranges` | supplement, `tab:ranges` | -- |

Eight tables in the Overleaf main text, three moved to the supplement, exactly as described.

## 3. Where each figure lives

| File | repository draft | Overleaf |
|---|---|---|
| `fig01_schematic` | main, `fig:schematic` | main, `fig:schematic` |
| `fig03_validation_xf1989` | main, `fig:xf` | main, `fig:xf` |
| `fig05_latham_calibrated` | main, `fig:latham` | main, `fig:latham` |
| `fig06_creep_life_profile` | main, `fig:creepprofile` | main, **`fig:life`** |
| `fig07_operating_maps` | main, `fig:maps` | main, `fig:maps` |
| `fig08_sobol` | main, `fig:sobol` | main, `fig:sobol` |
| `fig09_surrogate_parity` | main, `fig:parity` | main, **`fig:surrogate`** |
| `fig11_scenarios_summary` | main, `fig:scenarios` | main, `fig:scenarios` |
| `fig12_pareto` | main, `fig:pareto` | main, `fig:pareto` |
| `fig13_steam_credit` | main, `fig:steam` | main, `fig:steam` |
| `fig14_uq_distributions` | main, `fig:uq` | main, `fig:uq` |
| `fig15_variance_shares` | main, `fig:variance` | main, **`fig:shares`** |
| `fig00_equilibrium_conversion` | supplement, no label | supplement, `fig:equilibrium` |
| `fig01_kinetics` | supplement, no label | supplement, `fig:kinetics` |
| `fig02_reactor_profiles` | supplement, no label | supplement, `fig:profiles` |
| `fig04_latham_baseline` | supplement, no label | supplement, `fig:baseline` |
| `fig10_scenarios_timeseries` | supplement, no label | supplement, `fig:timeseries` |

No figure moved between main text and supplement; five labels were renamed and the five supplement
figures gained labels.

### Figure files in `chem/figures/` and their use

17 PDF files are present.

- **Referenced by `chem/main.tex` (12):** `fig01_schematic`, `fig03_validation_xf1989`,
  `fig05_latham_calibrated`, `fig06_creep_life_profile`, `fig07_operating_maps`, `fig08_sobol`,
  `fig09_surrogate_parity`, `fig11_scenarios_summary`, `fig12_pareto`, `fig13_steam_credit`,
  `fig14_uq_distributions`, `fig15_variance_shares`.
- **Referenced by `chem/supplementary.tex` (5):** `fig00_equilibrium_conversion`, `fig01_kinetics`,
  `fig02_reactor_profiles`, `fig04_latham_baseline`, `fig10_scenarios_timeseries`.
- **Referenced by neither: none.** Every file is used exactly once.
- **Referenced but missing from `figures/`: none.**
- No file is used by both documents, so the two can be uploaded to separate Overleaf projects
  without pruning.

Note that these are the **v2** figure PDFs. The v3 re-run wrote its figures to
`paper/figures/v3/`, which is a separate set and is not what `chem/figures/` holds.

---

## 4. TODOs unique to each version

**Overleaf (3, all in `main.tex`):**

1. `university e-mail of the second author` — in the second author's `\ead`.
2. `Zenodo DOI` — Data and code availability.
3. `Add funding and acknowledgements, if any.` — Acknowledgements.

**Repository draft (44).** All three Overleaf TODOs have a counterpart in the repository draft
except the first: the repository has `\todo{Second author name}` (the *name* is missing there,
whereas the Overleaf already names Madina Sissenbay and is missing only her e-mail). The
Acknowledgements TODO was already resolved in the repository by my earlier frontmatter edit.

The other **41 repository TODOs have no counterpart in the Overleaf pair** — they were resolved
during consolidation. By file: `methods.tex` 14, `results.tex` 11, `discussion.tex` 6,
`supplementary.tex` 4, `introduction.tex` 3, `backmatter.tex` 3 remaining
(`GitHub URL`, `Second author`, CRediT roles, journal AI wording), `main.tex` 1.

The GitHub URL is a concrete example: the repository draft still has `\todo{GitHub URL}` while the
Overleaf gives `https://github.com/IzbassarO/reformer-digital-twin`.

---

## 5. Content in the repository draft but absent from the Overleaf pair

27 sentences have no close counterpart. Excluding the three that are my own earlier frontmatter
edits (Funding section, Acknowledgements text, competing-interest sentence), the substantive
omissions are:

| Location | Content dropped in consolidation |
|---|---|
| 1 Introduction | the retube-cost sentence (was a TODO) |
| 1.1 Previous work | Latham's model described as "three-dimensional"; the Yeh sentence is rewritten; the neural-network/optimisation survey sentence is compressed |
| 2.1 Reactor tube model | the higher-alkane inlet rule (was a TODO) |
| 2.3 Tube wall and creep life | "wall thickness not given by Latham, taken as 15 mm" (was a TODO); the NIMS-replacement sentence |
| 2.4 Data | Plant C1 exclusion rationale; the NIMS template/schema and synthetic-validation sentences |
| 3.1 Verification and validation | the displaced-decimal hypothesis; the pressure-RMSE explanation; the "fixing $\alpha_{top}$ raises $\chi^2$ by 14 %, above our 10 % threshold" model-selection sentence; the expected-failure/4.5 K test note |
| 3.2 RQ1 | the mechanism sentence for why the hot spot rises with load at constant outlet temperature |
| 4.1 Industrial implications | "the plant's steam balance, not the reformer, determines whether running steam-rich is a life-saving measure or a fuel penalty" |
| 4.2 Flame length | "The three published operating points cannot separate the two hypotheses." |
| 4.3 Limitations | the 20 K tube-to-tube spread figure; the thermal-stress caveat; the NIMS-replacement sentence; the ramp-rate justification |
| 5 Conclusions | the closing "next steps" sentence (NIMS transcription, per-tube extension, pyrometer survey) |
| Data and code availability | "The two Monte Carlo result files above 5 MB are provided in the archive only." |
| Supplement | the entire "Version note" section (v1 vs v2 differences) |

Most of these were TODO-bearing or hedging sentences, which is consistent with the Overleaf pair
being a deliberate tightening. Two are worth a second look before submission because they carry
reasoning rather than hedging: the **model-selection justification** in 3.1 (why the three-parameter
fit is adopted) and the **steam-balance conclusion** in 4.1.

Conversely, the Overleaf adds material the repository draft does not have, including a
tube-geometry paragraph (146 mm OD, 12.2/12.4/15.3 mm wall variants), fuller digitisation detail in
2.4, the $z/L$ surrogate $R^2 = 0.64$ and the "average-temperature evaluation overstates life by a
factor of about 110" comparison with Yeh.

---

## 6. Numbers that differ between the two texts

Both texts quote v2 results, and the **tables agree cell for cell** for `tab:xf`, `tab:latham`,
`tab:cal`/`tab:params`, `tab:sobol`, `tab:scen`/`tab:scenarios`, `tab:regimes`, `tab:uq`,
`tab:robust`, `tab:surr`/`tab:surrogate` and `tab:steam` — the apparent differences are notation
only (`8.1e-05` vs `$8.10\times10^{-5}$`, `0.4` vs `0.400`).

Genuine differences:

| Quantity | repository draft | Overleaf | comment |
|---|---|---|---|
| Hot-spot temperature at the start of the ageing campaign | **1137 K** (`\TwoAgeStart`) | **1138 K** | Real inconsistency. The macro reads `T_wo_start` of the first campaign year; the Overleaf quotes the year-1 *mean* from `tab:scenarios` (1138). One of the two should be chosen. |
| Surrogate speed-up | 4658 (`\speedup`) | "about 4700" | rounding |
| Bed-side coefficient at mid-height | 1408 W m$^{-2}$ K$^{-1}$ (`\alphaIplant`) | "about 1400" | rounding |
| Flexible-operation robustness fraction | 0.63 (`\robFlexNom`, S3 only) | "0.57--0.63" | the Overleaf quotes the range over S2 and S3 and is the more complete statement |
| `tab:ranges` base inlet temperature | 884.6 K | 884.5 K | rounding of the same calibrated value |
| `tab:ranges` base inlet pressure | 30.06 bar | 30.1 bar | rounding |
| `tab:var` / `tab:variance` | sum **and** interaction residual (`0.999 / +0.001`, `1.037 / -0.037`) | sum only (`0.999`, `1.037`) | the interaction residual is dropped in the Overleaf table, although the text still refers to the shares summing to within 0.04 of one |

Everything else that the two texts share — 82 %, 4 K, 12 K, 2.25, +16 %, $-18$ %, 0.28, 0.083,
62 %, 42 K, 7.04--8.73, the factor of 50, 0.4--2.7 %, 14.0 K, 348 W m$^{-2}$ K$^{-1}$, $10^8$ h,
all robustness fractions — is identical.

---

## 7. Name spellings

The Overleaf `main.tex` carries **both** spellings: `\author[su]{Izbassar Orynbassar}` in the
frontmatter but `\textbf{Izbassar Orynbassarov:}` in the CRediT statement. Step 3 resolves this to
`Orynbassar` throughout. The second author's name is left exactly as the Overleaf has it,
**Madina Sissenbay**, and is not normalised against the e-mail address.

---

## 8. Word counts (`texcount -inc`)

Captions and tables are separated by extracting the `figure` and `table` environments into standalone
files and counting them apart; texcount lumps both into "words outside text".

### Current manuscript, `paper/main.tex`

| Category | Words |
|---|---|
| Body text | **8469** |
| Section headers | **128** |
| Figure captions (12 environments) | **530** |
| Tables (8 environments, captions + cells) | **317** |
| *outside-text total, as texcount reports it* | *847* |

Also 226 inline maths, 5 displayed, 37 headers, 20 floats. Compiled length 22 pages.

### Current supplement, `paper/supplementary.tex`

| Category | Words |
|---|---|
| Body text | 715 |
| Section headers | 35 |
| Figure captions (5 environments) | 225 |
| Tables (4 environments) | 143 |

Compiled length 6 pages.

### Archived draft, `paper/archive/draft_step21/main.tex`

| Category | Words |
|---|---|
| Body text | 7495 |
| Section headers | 128 |
| Figure captions | 589 |
| Tables (9 files, 11 environments) | 432 |

Compiled length 37 pages.

The consolidated manuscript has **more** body text than the draft (8469 against 7495) while being 15
pages shorter. The page count fell because of the layout — `final,3p,times` in place of
`preprint,12pt` with line numbers — not because prose was cut. What was removed was TODO-bearing and
hedging material (Section 5 above); what was added was the extra methodological detail also listed
there. Caption and table words fell (530 against 589, 317 against 432) mainly because three tables
moved to the supplement.
