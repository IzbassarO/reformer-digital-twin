"""Journal-style plotting helpers shared by all figures.

Single column 8.5 cm, double column 17.8 cm; 8 pt fonts; Okabe-Ito colour-blind-safe palette; no in-figure
titles (captions carry the description); SI units in axis labels; every figure saved as vector PDF and PNG.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Optional, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt

CM = 1.0 / 2.54
SINGLE = 8.5 * CM
DOUBLE = 17.8 * CM
#: Okabe & Ito (2008) palette
PALETTE = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
RC = {
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"], "mathtext.fontset": "dejavusans",
    "axes.prop_cycle": mpl.cycler(color=PALETTE), "axes.linewidth": 0.6, "lines.linewidth": 1.0, "lines.markersize": 3.5,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.direction": "out", "ytick.direction": "out",
    "legend.frameon": False, "legend.handlelength": 1.6, "axes.grid": False, "figure.dpi": 150, "savefig.dpi": 300,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
}


@contextmanager
def style():
    with mpl.rc_context(RC):
        yield


def figure(width: float = SINGLE, height: Optional[float] = None, nrows: int = 1, ncols: int = 1, **kw):
    """Create a figure of journal width; height defaults to 0.75 x width for a single panel row."""
    height = height or (0.75 * width if nrows == 1 else 0.6 * width * nrows / max(ncols, 1) + 1.0 * CM)
    return plt.subplots(nrows, ncols, figsize=(width, height), **kw)


def save(fig, stem: Path | str, formats: Sequence[str] = ("pdf", "png")) -> list[Path]:
    """Save ``fig`` as ``<stem>.pdf`` and ``<stem>.png`` (tight bounding box) and close it."""
    stem = Path(stem); stem.parent.mkdir(parents=True, exist_ok=True); out = []
    for fmt in formats:
        p = stem.with_suffix(f".{fmt}"); fig.savefig(p, bbox_inches="tight", pad_inches=0.02); out.append(p)
    plt.close(fig)
    return out


def panel_label(ax, text: str, x: float = -0.18, y: float = 1.02):
    ax.text(x, y, text, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom", ha="left")


def tidy(ax, xlabel: Optional[str] = None, ylabel: Optional[str] = None, legend: bool = False, **legend_kw):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if legend: ax.legend(**legend_kw)
    ax.tick_params(length=2.5)
    return ax
