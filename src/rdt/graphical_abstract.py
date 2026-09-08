"""Three-panel graphical abstract for the IJHE submission, built from the v3 result files.

``python -m rdt.graphical_abstract [--version v3]``

Every number is read from the same files the manuscript figures use; none is written into this
module. The three panels answer, left to right, what is varied, what the twin computes, and what it
costs:

* **Panel 1** -- simplified cross-section of a top-fired box reformer with the five operating inputs
  as labelled arrows. Drawn here rather than reused: ``paper/figures/fig01_schematic.pdf`` was
  committed as a rendered artefact and its source is not in the repository, so panel (a) of that
  figure cannot be imported. The geometry below follows the caption of Fig.~1 and the tube geometry
  of :mod:`rdt.latham_cases`, and carries no number of its own.
* **Panel 2** -- outer-wall temperature against ``z/L`` for the calibrated Plant A case, the same
  calculation that produces Fig.~5 (:func:`rdt.latham_cases.run_case` with
  :func:`rdt.operate.calibrated_params`). The hot spot is located by ``argmax`` of the computed
  profile, not asserted. The shaded band is the 5--95\\% Monte Carlo interval of the hot-spot
  temperature from the RQ4 campaign, applied as a uniform offset to the profile: the campaign stores
  scalars per sample, not profiles, so the band states the spread of the peak and must be read that
  way rather than as a per-elevation envelope.
* **Panel 3** -- life consumption per kmol H2 relative to the base, with 5--95\\% whiskers, for five
  operating points of the RQ4 table (``uq_summary_<tag>.json``, field
  ``life_rate_per_kmol_H2_rel_base``).

Design constraints of the journal: single row at 13 x 5 cm, sans-serif, nothing below 8 pt at that
size, no title, no author names, no logos. The colours are those of Figs.~5 and 8 --- the Okabe-Ito
blue and vermillion of :mod:`rdt.plotting` --- plus grey.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from rdt import config as C
from rdt import plotting as P

OUT_DIR = C.ROOT / "paper" / "figures"
STEM = "graphical_abstract"

#: 13 x 5 cm, the single-row aspect ratio the journal asks for
WIDTH_CM, HEIGHT_CM = 13.0, 5.0
#: panel widths, left to right
PANEL_WIDTHS = (0.25, 0.35, 0.40)
#: 600 dpi raster for the TIFF; 13 cm at 600 dpi is 3071 px, well above the 1328 px floor
RASTER_DPI = 600
MIN_PX = (1328, 531)

#: two colours plus grey, taken from the palette of Figs. 5 and 8
BLUE, VERMILLION = P.PALETTE[0], P.PALETTE[1]
GREY = "0.55"

#: the five operating points of panel 3, keyed as they appear in uq_summary_<tag>.json
BARS: Tuple[Tuple[str, str], ...] = (
    ("RQ3_min_life_iso_H2", "min life, iso-H$_2$"),
    ("RQ3_knee", "knee"),
    ("S1_base", "base"),
    ("load_0.70", "load-following"),
    ("S5_y4_holdslip_end", "aging year 4"),
)
CAPTION = "Aging costs more tube life than flexibility"


# ---------------------------------------------------------------------------
# Data, all of it read from the v3 result files
# ---------------------------------------------------------------------------
def wall_profile() -> Dict[str, np.ndarray]:
    """Calibrated Plant A outer-wall profile -- the calculation behind Fig.~5."""
    from rdt import latham_cases as lc
    from rdt import operate as op

    df = lc.load_cases()
    row = df[df.case == "Plant_A"].iloc[0]
    res = lc.run_case(row, op.calibrated_params())
    z, T = np.asarray(res.z_frac, float), np.asarray(res.T_wo, float)
    k = int(np.argmax(T))
    return {"z_frac": z, "T_wo": T, "hot_z_frac": float(z[k]), "hot_T": float(T[k])}


def uq_intervals(cfg: C.RunConfig) -> Dict[str, Dict[str, float]]:
    """The RQ4 table: median and 5--95\\% interval of the relative life cost, per operating point."""
    d = json.loads(Path(cfg.uq_summary).read_text())["intervals"]
    return {k: v["life_rate_per_kmol_H2_rel_base"] for k, v in d.items()}


def hotspot_band(cfg: C.RunConfig) -> Dict[str, float]:
    """5--95\\% Monte Carlo interval of the hot-spot temperature at the base point."""
    t = json.loads(Path(cfg.uq_summary).read_text())["intervals"]["S1_base"]["T_wo_max_K"]
    return {"median": float(t["median"]), "p05": float(t["p05"]), "p95": float(t["p95"])}


def collect(cfg: C.RunConfig) -> Dict[str, object]:
    prof = wall_profile()
    band = hotspot_band(cfg)
    iv = uq_intervals(cfg)
    missing = [k for k, _ in BARS if k not in iv]
    if missing:
        raise KeyError(f"{cfg.uq_summary} has no interval for {missing}")
    return {"profile": prof, "band": band, "intervals": iv,
            "bars": [(lab, iv[k]) for k, lab in BARS]}


# ---------------------------------------------------------------------------
# Panels. Nothing is drawn below MIN_PT, the journal's type floor at 13 x 5 cm.
# ---------------------------------------------------------------------------
#: every piece of type in the figure is at least this size
MIN_PT = 8.0

#: axes rectangles in figure coordinates. The panel *regions* are 25 / 35 / 40 % of the width; the
#: boxes below sit inside them, leaving room for the tick and category labels that draw outside.
#: rdt.plotting sets 7 pt ticks for the manuscript figures; the abstract lifts everything to the floor
RC_8PT = {"font.size": MIN_PT, "axes.labelsize": MIN_PT, "xtick.labelsize": MIN_PT,
          "ytick.labelsize": MIN_PT, "legend.fontsize": MIN_PT}

AXES_RECTS = {"furnace": (0.005, 0.10, 0.215, 0.86),
              "profile": (0.330, 0.30, 0.205, 0.58),
              "bars":    (0.740, 0.30, 0.245, 0.58)}


def _furnace(ax) -> None:
    """Top-fired box reformer, after the caption of Fig.~1. Schematic: carries no data."""
    from matplotlib.patches import FancyArrow, Polygon, Rectangle

    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    x0, x1 = 5.9, 9.9
    ax.add_patch(Rectangle((x0, 1.1), x1 - x0, 7.2, fc="none", ec=GREY, lw=0.9))        # firebox
    ax.add_patch(Rectangle((x0, 8.3), x1 - x0, 0.6, fc=GREY, ec=GREY, alpha=0.35))      # roof
    for x in (7.0, 8.8):                                                                # catalyst tubes
        ax.add_patch(Rectangle((x - 0.28, 1.7), 0.56, 6.4, fc=BLUE, ec=BLUE, alpha=0.35, lw=0.6))
    for x in (6.4, 7.9, 9.4):                                                           # flames
        ax.add_patch(Polygon([[x - 0.32, 8.3], [x + 0.32, 8.3], [x, 6.2]],
                             closed=True, fc=VERMILLION, ec="none", alpha=0.55))
    for x in (6.4, 9.4):                                                                # flue gas down
        ax.annotate("", xy=(x, 2.1), xytext=(x, 5.6),
                    arrowprops=dict(arrowstyle="-|>", color=GREY, lw=0.7, shrinkA=0, shrinkB=0))
    ax.add_patch(Rectangle((x0, 1.1), x1 - x0, 0.55, fc=GREY, ec=GREY, alpha=0.30))     # tunnels

    for i, lab in enumerate(("load", "S/C", "excess air", "firing", "activity")):
        y = 8.05 - 1.62 * i
        ax.text(0.0, y, lab, ha="left", va="center", fontsize=MIN_PT)
        ax.add_patch(FancyArrow(4.6, y, 1.0, 0.0, width=0.05, head_width=0.36, head_length=0.36,
                                length_includes_head=True, fc="k", ec="k", lw=0))


def _profile(ax, prof, band) -> None:
    z, T = prof["z_frac"], prof["T_wo"]
    lo = T + (band["p05"] - band["median"])
    hi = T + (band["p95"] - band["median"])
    ax.fill_between(z, lo, hi, color=BLUE, alpha=0.20, lw=0)
    ax.plot(z, T, color=BLUE, lw=1.4)
    ax.plot([prof["hot_z_frac"]], [prof["hot_T"]], "o", color=VERMILLION, ms=4.5, zorder=5)
    ax.set_ylim(float(lo.min()) - 15, float(hi.max()) + 70)
    ax.annotate(f"hot spot, $z/L$ = {prof['hot_z_frac']:.2f}",
                xy=(prof["hot_z_frac"], prof["hot_T"] + 6), xytext=(0.0, float(hi.max()) + 14),
                fontsize=MIN_PT, color=VERMILLION, ha="left", va="bottom",
                arrowprops=dict(arrowstyle="-", color=VERMILLION, lw=0.6, shrinkA=1, shrinkB=3))
    ax.text(0.97, 0.04, "5–95 % MC", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=MIN_PT, color=BLUE)
    ax.set_xticks([0.0, 0.5, 1.0])
    P.tidy(ax, "$z/L$ (–)", "Outer-wall $T$ (K)")
    #: the arrow that carries the profile into the life calculation, kept inside this panel
    ax.annotate("", xy=(1.02, 0.40), xytext=(0.78, 0.40), xycoords="axes fraction",
                textcoords="axes fraction", annotation_clip=False,
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.0))
    ax.text(0.90, 0.45, "Larson–Miller", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=MIN_PT, clip_on=False)


def _bars(ax, bars) -> None:
    labels = [b[0] for b in bars]
    med = np.array([b[1]["median"] for b in bars], float)
    lo = np.array([b[1]["p05"] for b in bars], float)
    hi = np.array([b[1]["p95"] for b in bars], float)
    y = np.arange(len(bars))[::-1]
    colours = [GREY if lab == "base" else (VERMILLION if m > 1.0 else BLUE)
               for lab, m in zip(labels, med)]
    ax.barh(y, med, height=0.60, color=colours, alpha=0.9, lw=0)
    paired = np.isclose(lo, hi)                       # the base is paired: no interval to draw
    ax.errorbar(med[~paired], y[~paired],
                xerr=[med[~paired] - lo[~paired], hi[~paired] - med[~paired]],
                fmt="none", ecolor="k", elinewidth=0.8, capsize=1.8)
    ax.axvline(1.0, color="k", lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=MIN_PT)
    ax.set_ylim(-0.7, len(bars) - 0.3)
    ax.set_xlim(0.72 * lo.min(), 3.1 * hi.max())
    ax.set_xticks([0.1, 1.0]); ax.set_xticklabels(["0.1", "1"])
    for yy, m, h in zip(y, med, hi):
        ax.text(h * 1.26, yy, f"{m:.2f}", ha="left", va="center", fontsize=MIN_PT)
    P.tidy(ax, "Life per kmol H$_2$, rel. base (–)", None)
    # the label is wider than this panel at 8 pt; shift it left so it stays on the canvas
    ax.xaxis.set_label_coords(0.26, -0.17)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------
def build(cfg: C.RunConfig, out_dir: Optional[Path] = None) -> Dict[str, object]:
    import matplotlib.pyplot as plt

    data = collect(cfg)
    out_dir = Path(out_dir or OUT_DIR); out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / STEM
    import matplotlib as mpl
    with P.style(), mpl.rc_context(RC_8PT):
        fig = plt.figure(figsize=(WIDTH_CM * P.CM, HEIGHT_CM * P.CM))
        _furnace(fig.add_axes(AXES_RECTS["furnace"]))
        _profile(fig.add_axes(AXES_RECTS["profile"]), data["profile"], data["band"])
        _bars(fig.add_axes(AXES_RECTS["bars"]), data["bars"])
        # the caption belongs under panel 3; at 8 pt it is wider than that panel alone, so it is
        # centred on the right-hand part of the figure rather than clipped at the edge
        fig.text(0.68, 0.035, CAPTION, ha="center", va="bottom", fontsize=MIN_PT, style="italic")
        fig.savefig(stem.with_suffix(".pdf"), format="pdf")
        png = stem.with_suffix(".png")
        fig.savefig(png, format="png", dpi=RASTER_DPI)
        checks = {"overflow": _overflow(fig), "min_pt": _min_font_pt(fig)}
        plt.close(fig)
    tif = _to_tiff(png, stem.with_suffix(".tif"))
    data["files"] = {"pdf": stem.with_suffix(".pdf"), "tif": tif, "png": png}
    data["checks"] = checks
    return data


def _texts(fig):
    return [o for o in fig.findobj(match=lambda x: hasattr(x, "get_text") and hasattr(x, "get_fontsize"))
            if (o.get_text() or "").strip()]


def _overflow(fig) -> List[str]:
    """Text whose rendered box leaves the canvas -- the defect a clipped caption would cause."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    out = []
    for o in _texts(fig):
        b = o.get_window_extent(renderer=r)
        if b.x0 < -1 or b.y0 < -1 or b.x1 > fig.bbox.width + 1 or b.y1 > fig.bbox.height + 1:
            out.append(o.get_text().strip())
    return out


def _min_font_pt(fig) -> float:
    """Smallest type actually rendered, in points at the figure's own size."""
    return min(float(o.get_fontsize()) for o in _texts(fig))


def _to_tiff(png: Path, tif: Path) -> Path:
    """LZW-compressed TIFF at :data:`RASTER_DPI`, flattened onto white (TIFF has no alpha here)."""
    from PIL import Image

    im = Image.open(png).convert("RGBA")
    flat = Image.new("RGB", im.size, (255, 255, 255))
    flat.paste(im, mask=im.split()[3])
    flat.save(tif, format="TIFF", compression="tiff_lzw", dpi=(RASTER_DPI, RASTER_DPI))
    return tif


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def verify(cfg: C.RunConfig, data: Dict[str, object]) -> str:
    from PIL import Image

    L: List[str] = []
    p = L.append
    prof, band = data["profile"], data["band"]
    p(f"Graphical abstract, {cfg.tag} result files. No number is written into the module.")
    p("=" * 96)
    p("PANEL 2 -- calibrated Plant A outer-wall profile (the calculation behind Fig. 5)")
    p(f"  source            rdt.latham_cases.run_case(Plant_A, rdt.operate.calibrated_params())")
    p(f"  nodes             {len(prof['z_frac'])}")
    p(f"  hot spot          z/L = {prof['hot_z_frac']:.4f}  (argmax of the computed profile)")
    p(f"  hot-spot T        {prof['hot_T']:.2f} K")
    p(f"  profile range     {prof['T_wo'].min():.2f} to {prof['T_wo'].max():.2f} K")
    p(f"  band source       {Path(cfg.uq_summary).name}, intervals.S1_base.T_wo_max_K")
    p(f"  band              median {band['median']:.3f} K, 5 % {band['p05']:.3f} K, "
      f"95 % {band['p95']:.3f} K  (offsets {band['p05']-band['median']:+.2f}/"
      f"{band['p95']-band['median']:+.2f} K)")
    p("  the band is the spread of the PEAK applied as a uniform offset: the RQ4 campaign stores")
    p("  scalars per sample, not profiles, so it is not a per-elevation envelope")
    p("")
    p("PANEL 3 -- relative life cost, 5-95 % (the RQ4 table)")
    p(f"  source            {Path(cfg.uq_summary).name}, intervals.*.life_rate_per_kmol_H2_rel_base")
    p(f"  {'bar':<22} {'median':>8} {'5 %':>8} {'95 %':>8}   key")
    for (key, lab), (_, iv) in zip(BARS, data["bars"]):
        p(f"  {lab:<22} {iv['median']:8.4f} {iv['p05']:8.4f} {iv['p95']:8.4f}   {key}")
    p("")
    p("PANEL 1 -- schematic only, carries no number")
    p("  fig01_schematic.pdf was committed as a rendered artefact with no source in the repository,")
    p("  so panel (a) could not be reused; the cross-section is redrawn from the Fig. 1 caption.")
    p("")
    p("OUTPUT")
    p("-" * 96)
    f = data["files"]
    im = Image.open(f["tif"])
    px_ok = im.size[0] >= MIN_PX[0] and im.size[1] >= MIN_PX[1]
    p(f"  {f['pdf'].name:<28} {f['pdf'].stat().st_size:>9} bytes   vector")
    p(f"  {f['tif'].name:<28} {f['tif'].stat().st_size:>9} bytes   {im.size[0]} x {im.size[1]} px, "
      f"{im.info.get('compression', 'n/a')}, {tuple(int(round(d)) for d in im.info.get('dpi', (0, 0)))} dpi")
    p(f"  size              {WIDTH_CM:g} x {HEIGHT_CM:g} cm, aspect "
      f"{WIDTH_CM/HEIGHT_CM:.2f}:1 (13:5 = {13/5:.2f}:1)")
    p(f"  pixel floor       {MIN_PX[0]} x {MIN_PX[1]} required -> {'OK' if px_ok else 'FAILS'}")
    p(f"  panel widths      {', '.join(f'{100*w:.0f} %' for w in PANEL_WIDTHS)}")
    ch = data["checks"]
    p(f"  smallest type     {ch['min_pt']:.1f} pt at {WIDTH_CM:g} cm "
      f"(floor {MIN_PT:g} pt) -> {'OK' if ch['min_pt'] >= MIN_PT else 'FAILS'}")
    p(f"  text off canvas   {', '.join(ch['overflow']) if ch['overflow'] else 'none'}")
    p(f"  colours           {BLUE} (blue), {VERMILLION} (vermillion), {GREY} (grey) -- "
      "the palette of Figs. 5 and 8")
    return "\n".join(L)


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--version", default="v3", choices=sorted(C.BY_TAG))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    a = ap.parse_args(argv)
    cfg = C.BY_TAG[a.version]
    from rdt import pipeline as pl
    pl.apply_config(cfg)
    data = build(cfg, a.out_dir)
    print(verify(cfg, data))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
