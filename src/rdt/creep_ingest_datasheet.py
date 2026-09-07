"""Digitise the *Parametric stress rupture strength* chart of a manufacturer material data sheet.

The Schmidt + Clemens Centralloy sheets (see ``data/creep_derived/sources.yaml``) print, on page 6, one
chart carrying two curves: **Average** and **Lower Scatter Band** (the latter stated on the sheet to be
the 95 % confidence level; it is stored here as ``minimum``). This module turns those two curves into
``data/creep_derived/<alloy>_rupture_curve.csv`` with columns

    alloy, curve, lmp, stress_mpa, source_file, page, method, extracted_on

Two extraction paths are provided, selected by ``method``:

``vector``
    The chart is line art. :mod:`pdfplumber` exposes every stroked path in ``page.curves`` and
    ``page.lines``; the two rupture curves are the only thick strokes inside the plot frame and are
    identified by their stroke colour (see ``curve_stroke_cmyk`` in ``sources.yaml``). Cubic Bezier
    segments (``c``/``v``/``y`` operators) are flattened analytically, so the result is exact up to the
    flattening tolerance -- there is no pixel quantisation at all.
``raster``
    Fallback for sheets whose chart is a flattened image. The plot region is rendered at 600 dpi with
    ``pdftoppm`` (or :mod:`pypdfium2` when poppler is not installed), the two curves are separated by
    colour, and one value is traced per pixel column.

Device-to-data conversion is driven by **two calibration ticks per axis**, given in PDF user-space points
with the origin at the top-left of the page (pdfplumber's ``x`` / ``top``). ``sources.yaml`` stores the
ticks read off the printed gridlines. Note that these charts are rotated 90 degrees on the page: the LMP
axis (linear) runs vertically and the stress axis (log10) runs horizontally. The calibration is therefore
expressed per axis with an explicit page direction (``x_pt`` or ``top_pt``) rather than assuming a layout.

CLI::

    python -m rdt.creep_ingest_datasheet --all
    python -m rdt.creep_ingest_datasheet --key G4852 --method raster
    python -m rdt.creep_ingest_datasheet --key G4852 --page 6 \
        --lmp-ticks 25,762.4,32,244.1 --stress-ticks 100,181.7,1,555.6
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data" / "creep_derived"
SOURCES_YAML = DATA_DIR / "sources.yaml"
CURVE_KINDS = ("average", "minimum")
RASTER_DPI = 600


# ---------------------------------------------------------------------------
# Source specification
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Axis:
    """A chart axis calibrated by two ticks on one page direction.

    ``direction`` is ``"x"`` (page abscissa) or ``"top"`` (page ordinate, increasing downwards);
    ``transform`` is ``"linear"`` or ``"log10"``. ``ticks`` are two ``(value, coordinate_pt)`` pairs.
    """

    name: str
    direction: str
    transform: str
    ticks: Tuple[Tuple[float, float], Tuple[float, float]]

    def to_data(self, coord):
        """Map a page coordinate [pt] to the axis value."""
        (v0, c0), (v1, c1) = self.ticks
        if c1 == c0:
            raise ValueError(f"axis {self.name}: calibration ticks share a coordinate ({c0})")
        y0, y1 = (np.log10(v0), np.log10(v1)) if self.transform == "log10" else (v0, v1)
        y = y0 + (np.asarray(coord, float) - c0) * (y1 - y0) / (c1 - c0)
        return 10.0**y if self.transform == "log10" else y

    def to_page(self, value):
        """Inverse of :meth:`to_data`: axis value to page coordinate [pt]."""
        (v0, c0), (v1, c1) = self.ticks
        y0, y1 = (np.log10(v0), np.log10(v1)) if self.transform == "log10" else (v0, v1)
        y = np.log10(np.asarray(value, float)) if self.transform == "log10" else np.asarray(value, float)
        return c0 + (y - y0) * (c1 - c0) / (y1 - y0)


@dataclass(frozen=True)
class SourceSpec:
    """One data-sheet chart: file, page, calibration and curve colours (from ``sources.yaml``)."""

    key: str
    alloy: str
    designation: str
    datasheet_file: str
    page: int
    C: float
    chart_title: str
    lmp: Axis
    stress: Axis
    frame: Dict[str, float]
    stroke_cmyk: Dict[str, Tuple[float, ...]]
    datasheet_dir: Path

    @property
    def pdf_path(self) -> Path:
        return self.datasheet_dir / self.datasheet_file

    @property
    def slug(self) -> str:
        return self.alloy.lower().replace(" ", "_").replace("(r)", "").replace("+", "")

    @property
    def out_csv(self) -> Path:
        return DATA_DIR / f"{self.slug}_rupture_curve.csv"


def _axis_from_yaml(name: str, d: dict) -> Axis:
    ticks = d["ticks"]
    if len(ticks) != 2:
        raise ValueError(f"axis {name}: expected exactly two calibration ticks, got {len(ticks)}")
    direction = "x" if "x_pt" in ticks[0] else "top"
    key = f"{direction}_pt"
    return Axis(name=name, direction=direction, transform=d.get("transform", "linear"),
                ticks=tuple((float(t["value"]), float(t[key])) for t in ticks))  # type: ignore[arg-type]


def load_sources(path: Path = SOURCES_YAML) -> Dict[str, SourceSpec]:
    """Read ``sources.yaml`` into :class:`SourceSpec` objects keyed by their ``key`` field."""
    d = yaml.safe_load(Path(path).read_text())
    base = ROOT / d.get("datasheet_dir", "data/creep_derived/datasheets")
    out = {}
    for s in d["sources"]:
        out[s["key"]] = SourceSpec(
            key=s["key"], alloy=s["alloy"], designation=s.get("designation", ""), datasheet_file=s["datasheet_file"],
            page=int(s["page"]), C=float(s["larson_miller_C"]), chart_title=s.get("chart_title", ""),
            lmp=_axis_from_yaml("lmp", s["axis"]["lmp"]), stress=_axis_from_yaml("stress_MPa", s["axis"]["stress_MPa"]),
            frame={k: float(v) for k, v in s["plot_frame_pt"].items()},
            stroke_cmyk={k: tuple(float(c) for c in v) for k, v in s["curve_stroke_cmyk"].items()},
            datasheet_dir=base)
    return out


# ---------------------------------------------------------------------------
# Path flattening
# ---------------------------------------------------------------------------
def flatten_path(path: Sequence[tuple], n_per_segment: int = 240) -> np.ndarray:
    """Flatten a pdfplumber path (``m``/``l``/``c``/``v``/``y``/``h`` operators) to an ``(n, 2)`` polyline.

    PDF's three cubic operators differ only in which control points are implied: ``c`` gives both,
    ``v`` reuses the current point as the first control point, ``y`` reuses the endpoint as the second.
    """
    pts: List[np.ndarray] = []
    cur: Optional[np.ndarray] = None
    start: Optional[np.ndarray] = None
    t = np.linspace(0.0, 1.0, n_per_segment)[1:, None]
    for op in path:
        kind, args = op[0], [np.array(a, float) for a in op[1:]]
        if kind == "m":
            cur = args[0]; start = cur.copy(); pts.append(cur.copy()); continue
        if cur is None:
            raise ValueError("path does not start with a moveto")
        if kind == "l":
            pts.append(args[0]); cur = args[0]; continue
        if kind == "h":
            if start is not None:
                pts.append(start.copy()); cur = start.copy()
            continue
        if kind == "c":
            p1, p2, p3 = args
        elif kind == "v":
            p1, (p2, p3) = cur, args
        elif kind == "y":
            (p1, p3), p2 = args, args[1]
        else:
            continue
        pts.extend((1 - t) ** 3 * cur + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3)
        cur = p3
    return np.asarray(pts, float).reshape(-1, 2)


# ---------------------------------------------------------------------------
# Colour helpers
# ---------------------------------------------------------------------------
def _cmyk_to_rgb(cmyk: Sequence[float]) -> np.ndarray:
    c, m, y, k = (list(cmyk) + [0.0, 0.0, 0.0, 0.0])[:4]
    return np.array([255 * (1 - c) * (1 - k), 255 * (1 - m) * (1 - k), 255 * (1 - y) * (1 - k)], float)


def _colour_distance(a: Optional[Sequence[float]], b: Sequence[float]) -> float:
    if a is None or not isinstance(a, (tuple, list)) or len(a) != len(b):
        return np.inf
    return float(np.max(np.abs(np.asarray(a, float) - np.asarray(b, float))))


# ---------------------------------------------------------------------------
# Vector extraction
# ---------------------------------------------------------------------------
def _inside_frame(obj, frame: Dict[str, float], pad: float = 6.0) -> bool:
    return (obj["x0"] >= frame["x0"] - pad and obj["x1"] <= frame["x1"] + pad
            and obj["top"] >= frame["top"] - pad and obj["bottom"] <= frame["bottom"] + pad)


def extract_vector(spec: SourceSpec, min_linewidth: float = 1.0, colour_tol: float = 0.12,
                   n_per_segment: int = 240, min_span_frac: float = 0.30) -> Dict[str, np.ndarray]:
    """Read the two curves from the vector content of the page.

    Returns ``{curve_kind: (n, 2) array of (lmp, stress_MPa)}``. Candidates are the stroked paths inside
    the plot frame that are thicker than the gridlines and that span at least ``min_span_frac`` of the
    frame diagonal -- the span test is what rejects the two **legend swatches**, which are short strokes
    drawn inside the frame in exactly the two curve colours and at a similar width. Each survivor is then
    matched to a curve kind by stroke colour, the longest one winning any tie.
    """
    import pdfplumber

    fw = spec.frame["x1"] - spec.frame["x0"]
    fh = spec.frame["bottom"] - spec.frame["top"]
    min_span = min_span_frac * float(np.hypot(fw, fh))

    with pdfplumber.open(spec.pdf_path) as pdf:
        page = pdf.pages[spec.page - 1]
        candidates = []
        for obj in list(page.curves) + list(page.lines):
            if obj.get("fill") or (obj.get("linewidth") or 0.0) < min_linewidth:
                continue
            if not _inside_frame(obj, spec.frame):
                continue
            if np.hypot(obj["x1"] - obj["x0"], obj["bottom"] - obj["top"]) < min_span:
                continue
            path = obj.get("path") or [("m", (obj["x0"], obj["top"])), ("l", (obj["x1"], obj["bottom"]))]
            candidates.append((obj, flatten_path(path, n_per_segment)))

    if not candidates:
        raise RuntimeError(f"{spec.pdf_path.name} page {spec.page}: no stroked path spanning the plot frame")

    out: Dict[str, np.ndarray] = {}
    used = set()
    for kind in CURVE_KINDS:
        target = spec.stroke_cmyk[kind]
        scored = [(i, _colour_distance(o.get("stroking_color"), target), np.hypot(o["x1"] - o["x0"], o["bottom"] - o["top"]))
                  for i, (o, _) in enumerate(candidates)]
        scored = [(i, d, span) for i, d, span in scored if i not in used and d <= colour_tol]
        if not scored:
            raise RuntimeError(f"{spec.pdf_path.name} page {spec.page}: no stroke matching the "
                               f"{kind!r} colour {target} within {colour_tol}")
        i = min(scored, key=lambda s: (s[1], -s[2]))[0]
        used.add(i)
        out[kind] = _to_data(candidates[i][1], spec)
    return _order_and_check(out, spec)


def _to_data(poly: np.ndarray, spec: SourceSpec) -> np.ndarray:
    """Convert an ``(n, 2)`` polyline of page ``(x, top)`` points to ``(lmp, stress_MPa)``."""
    page = {"x": poly[:, 0], "top": poly[:, 1]}
    lmp = spec.lmp.to_data(page[spec.lmp.direction])
    stress = spec.stress.to_data(page[spec.stress.direction])
    arr = np.column_stack([lmp, stress])
    arr = arr[np.argsort(arr[:, 0])]
    keep = np.concatenate([[True], np.diff(arr[:, 0]) > 1e-9])
    return arr[keep]


def _order_and_check(curves: Dict[str, np.ndarray], spec: SourceSpec) -> Dict[str, np.ndarray]:
    """Verify that the *average* curve lies above the *minimum* curve; swap them if the colours lied."""
    a, m = curves["average"], curves["minimum"]
    lo = max(a[:, 0].min(), m[:, 0].min()); hi = min(a[:, 0].max(), m[:, 0].max())
    if not np.isfinite(lo) or hi <= lo:
        raise RuntimeError(f"{spec.key}: the two extracted curves do not overlap in LMP")
    grid = np.linspace(lo, hi, 64)
    da = np.interp(grid, a[:, 0], np.log10(a[:, 1])) - np.interp(grid, m[:, 0], np.log10(m[:, 1]))
    if np.median(da) < 0:
        curves = {"average": m, "minimum": a}
    elif np.min(da) < -1e-3:
        raise RuntimeError(f"{spec.key}: the two curves cross; check the calibration and stroke colours")
    return curves


# ---------------------------------------------------------------------------
# Raster extraction (fallback for sheets whose chart is a flattened image)
# ---------------------------------------------------------------------------
def _crop_origin(spec: SourceSpec) -> Tuple[float, float]:
    """Page ``(x, top)`` [pt] of the top-left corner of the rendered image.

    Renderers (both poppler and pdfium) rasterise the **CropBox**, whereas pdfplumber reports object
    coordinates against the **MediaBox**. These sheets carry a 15 mm bleed, so the two differ by 42.52 pt
    in each direction; ignoring it shifts the whole raster trace. Returned as an offset to subtract.
    """
    import pdfplumber

    with pdfplumber.open(spec.pdf_path) as pdf:
        page = pdf.pages[spec.page - 1]
        media, crop = page.mediabox, page.cropbox
    if crop is None or media is None:
        return 0.0, 0.0
    return float(crop[0] - media[0]), float(media[3] - crop[3])


def _render_page(spec: SourceSpec, dpi: int, out_dir: Path) -> np.ndarray:
    """Render the page to an RGB array at ``dpi``, preferring ``pdftoppm`` and falling back to pypdfium2."""
    if shutil.which("pdftoppm"):
        stem = out_dir / "page"
        subprocess.run(["pdftoppm", "-r", str(dpi), "-f", str(spec.page), "-l", str(spec.page),
                        "-png", str(spec.pdf_path), str(stem)], check=True, capture_output=True)
        png = sorted(out_dir.glob("page*.png"))
        if not png:
            raise RuntimeError("pdftoppm produced no output")
        from PIL import Image
        return np.asarray(Image.open(png[0]).convert("RGB"), dtype=np.uint8)
    import pypdfium2 as pdfium

    pdf = pdfium.PdfDocument(str(spec.pdf_path))
    try:
        img = pdf[spec.page - 1].render(scale=dpi / 72.0).to_pil().convert("RGB")
        return np.asarray(img, dtype=np.uint8)
    finally:
        pdf.close()


def _find_ink(crop: np.ndarray, seed: np.ndarray, max_thickness: float, min_pixels: int,
              max_candidates: int = 64) -> np.ndarray:
    """Rendered ink colour of a curve: the *thin* exact colour nearest the declared stroke colour.

    Two things make this harder than nearest-colour matching. First, the CMYK stroke colour declared in
    the PDF is only a seed -- the rasteriser applies its own colour management, so the ink lands a fair
    distance away (DeviceCMYK 1.00/0.51/0/0 renders as RGB 0/113/187, not the naive 0/125/255). Second,
    the plot area is not white: it carries a pale blue gradient whose shades sit *closer* to the declared
    scatter-band cyan than the real light-blue ink does, so nearest-colour alone locks onto the
    background and the trace leaks along the frame.

    The discriminator used here is geometric rather than chromatic: a stroke covers only about its own
    line width in each column it occupies, whereas a background fills hundreds of rows per column.
    Candidates are therefore examined nearest-seed-first and the first one thinner than
    ``max_thickness`` pixels per occupied column is accepted. Exact (unquantised) colours are used, so
    the accepted mask is the stroke interior and a tight tolerance can be applied afterwards.
    """
    key = (crop[:, :, 0].astype(np.int64) << 16) | (crop[:, :, 1].astype(np.int64) << 8) | crop[:, :, 2].astype(np.int64)
    vals, counts = np.unique(key, return_counts=True)
    vals = vals[counts >= min_pixels]
    if not vals.size:
        raise RuntimeError("no colour is common enough to be a curve")
    rgb = np.stack([(vals >> 16) & 255, (vals >> 8) & 255, vals & 255], axis=1).astype(float)
    for i in np.argsort(np.linalg.norm(rgb - seed[None, :], axis=1))[:max_candidates]:
        mask = key == vals[i]
        n_cols = int(mask.any(axis=0).sum())
        if n_cols and mask.sum() / n_cols <= max_thickness:
            return rgb[i]
    raise RuntimeError(f"no sufficiently thin colour near seed {seed.round(0)}")


def _trace_runs(mask: np.ndarray, max_run: int, min_run: int, max_jump: float) -> Tuple[np.ndarray, np.ndarray]:
    """Trace one row per column through ``mask`` by run continuity.

    For every column the matching rows are grouped into consecutive runs, keeping only those between
    ``min_run`` and ``max_run`` pixels thick: the floor drops anti-aliasing speckle, and the ceiling
    drops the **legend swatches**, which are tall vertical bars drawn inside the plot frame in exactly
    the two curve colours. The trace then walks outwards from the cleanest column, always taking the run
    whose centre is nearest the previously accepted one and refusing steps larger than ``max_jump`` per
    column, so neither a gridline nor a legend fragment can capture it.
    """
    per_col: List[List[float]] = []
    for j in range(mask.shape[1]):
        idx = np.flatnonzero(mask[:, j])
        runs: List[float] = []
        if idx.size:
            for grp in np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1):
                if min_run <= len(grp) <= max_run:
                    runs.append(float(grp.mean()))
        per_col.append(runs)

    singles = [j for j, r in enumerate(per_col) if len(r) == 1]
    if not singles:
        raise RuntimeError("no column with an unambiguous run to seed the trace")
    seed_col = singles[len(singles) // 2]
    centres: Dict[int, float] = {seed_col: per_col[seed_col][0]}
    for rng in (range(seed_col + 1, len(per_col)), range(seed_col - 1, -1, -1)):
        last, last_j = centres[seed_col], seed_col
        for j in rng:
            if not per_col[j]:
                continue
            c = min(per_col[j], key=lambda r: abs(r - last))
            if abs(c - last) > max_jump * max(1, abs(j - last_j)):
                continue
            centres[j], last, last_j = c, c, j
    cols = np.array(sorted(centres), float)
    return cols, np.array([centres[int(j)] for j in cols], float)


def extract_raster(spec: SourceSpec, dpi: int = RASTER_DPI, colour_tol: float = 25.0,
                   linewidth_pt: float = 2.72) -> Dict[str, np.ndarray]:
    """Trace the two curves from a 600 dpi rendering of the plot region, separating them by colour.

    One value is traced per pixel column of the plot region (the page-abscissa direction).
    """
    with tempfile.TemporaryDirectory() as td:
        rgb = _render_page(spec, dpi, Path(td))
    s = dpi / 72.0
    ox, ot = _crop_origin(spec)
    f = spec.frame
    x0, x1 = max(int(round((f["x0"] - ox) * s)), 0), min(int(round((f["x1"] - ox) * s)), rgb.shape[1])
    t0, t1 = max(int(round((f["top"] - ot) * s)), 0), min(int(round((f["bottom"] - ot) * s)), rgb.shape[0])
    crop = rgb[t0:t1, x0:x1, :].astype(float)
    if crop.size == 0:
        raise RuntimeError(f"{spec.key}: empty plot region at {dpi} dpi")
    lw_px = linewidth_pt * s
    max_run, min_run = max(4, int(round(4.0 * lw_px))), max(2, int(round(0.4 * lw_px)))
    max_jump = max(4.0, 1.5 * lw_px)

    out: Dict[str, np.ndarray] = {}
    for kind in CURVE_KINDS:
        seed = _cmyk_to_rgb(spec.stroke_cmyk[kind])
        ink = _find_ink(crop, seed, max_thickness=max_run, min_pixels=max(200, crop.shape[1] // 4))
        mask = np.linalg.norm(crop - ink[None, None, :], axis=2) <= colour_tol
        try:
            cols, rows = _trace_runs(mask, max_run, min_run, max_jump)
        except RuntimeError as e:
            raise RuntimeError(f"{spec.key}: raster trace failed for {kind!r} (ink {ink.round(0)}): {e}") from e
        if len(cols) < 20:
            raise RuntimeError(f"{spec.key}: raster trace found only {len(cols)} columns for {kind!r}")
        poly = np.column_stack([(cols + x0) / s + ox, (rows + t0) / s + ot])
        out[kind] = _to_data(poly, spec)
    return _order_and_check(out, spec)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def digitise(spec: SourceSpec, method: str = "auto", dpi: int = RASTER_DPI, verbose: bool = True) -> pd.DataFrame:
    """Digitise one data sheet and return the tidy frame written to ``spec.out_csv``."""
    if not spec.pdf_path.exists():
        raise FileNotFoundError(
            f"{spec.pdf_path} not found. The Schmidt + Clemens PDFs are copyright and are not "
            f"redistributed with this repository; place {spec.datasheet_file} in {spec.datasheet_dir}.")
    if method not in ("auto", "vector", "raster"):
        raise ValueError("method must be 'auto', 'vector' or 'raster'")

    used = method
    if method in ("auto", "vector"):
        try:
            curves = extract_vector(spec); used = "vector"
        except Exception as e:  # noqa: BLE001 - fall back to the raster path
            if method == "vector":
                raise
            if verbose:
                print(f"  vector extraction failed ({type(e).__name__}: {e}); falling back to raster")
            curves = extract_raster(spec, dpi=dpi); used = "raster"
    else:
        curves = extract_raster(spec, dpi=dpi)

    stamp = _dt.datetime.now(_dt.timezone.utc).date().isoformat()
    rows = []
    for kind in CURVE_KINDS:
        arr = curves[kind]
        rows.append(pd.DataFrame({"alloy": spec.alloy, "curve": kind, "lmp": arr[:, 0], "stress_mpa": arr[:, 1],
                                  "source_file": spec.datasheet_file, "page": spec.page, "method": used,
                                  "extracted_on": stamp}))
    df = pd.concat(rows, ignore_index=True)
    spec.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(spec.out_csv, index=False)

    if verbose:
        print(f"{spec.alloy}  ({spec.designation}, C = {spec.C})")
        print(f"  {spec.datasheet_file} page {spec.page} -> {spec.out_csv.relative_to(ROOT)}  [{used}]")
        for kind in CURVE_KINDS:
            a = curves[kind]
            print(f"    {kind:8s}: {len(a):4d} points   LMP {a[:, 0].min():6.3f} .. {a[:, 0].max():6.3f}   "
                  f"stress {a[:, 1].min():7.3f} .. {a[:, 1].max():7.3f} MPa")
        print(f"    sha256({spec.datasheet_file}) = {_sha256(spec.pdf_path)[:16]}...")
    return df


def digitise_all(method: str = "auto", dpi: int = RASTER_DPI, sources: Path = SOURCES_YAML,
                 verbose: bool = True) -> Dict[str, pd.DataFrame]:
    return {k: digitise(s, method=method, dpi=dpi, verbose=verbose) for k, s in load_sources(sources).items()}


def _parse_ticks(text: str) -> List[dict]:
    v = [float(t) for t in text.split(",")]
    if len(v) != 4:
        raise argparse.ArgumentTypeError("ticks must be 'value1,coord1,value2,coord2'")
    return [{"value": v[0], "_c": v[1]}, {"value": v[2], "_c": v[3]}]


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--all", action="store_true", help="digitise every source in sources.yaml")
    p.add_argument("--key", help="digitise a single source by its key (e.g. G4852)")
    p.add_argument("--method", default="auto", choices=("auto", "vector", "raster"))
    p.add_argument("--dpi", type=int, default=RASTER_DPI, help="raster fallback resolution (default 600)")
    p.add_argument("--page", type=int, help="override the chart page number")
    p.add_argument("--lmp-ticks", type=_parse_ticks, metavar="v1,c1,v2,c2",
                   help="override the LMP calibration ticks (value,pt,value,pt)")
    p.add_argument("--stress-ticks", type=_parse_ticks, metavar="v1,c1,v2,c2",
                   help="override the stress calibration ticks (value,pt,value,pt)")
    p.add_argument("--sources", type=Path, default=SOURCES_YAML)
    a = p.parse_args(argv)

    specs = load_sources(a.sources)
    if a.all:
        targets = list(specs.values())
    elif a.key:
        if a.key not in specs:
            p.error(f"unknown key {a.key!r}; available: {', '.join(specs)}")
        targets = [specs[a.key]]
    else:
        p.error("give --all or --key")

    import dataclasses

    for spec in targets:
        over = {}
        if a.page:
            over["page"] = a.page
        for axis_name, ticks in (("lmp", a.lmp_ticks), ("stress", a.stress_ticks)):
            if ticks:
                axis: Axis = getattr(spec, axis_name)
                over[axis_name] = dataclasses.replace(
                    axis, ticks=tuple((t["value"], t["_c"]) for t in ticks))  # type: ignore[arg-type]
        if over:
            spec = dataclasses.replace(spec, **over)
        digitise(spec, method=a.method, dpi=a.dpi)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
