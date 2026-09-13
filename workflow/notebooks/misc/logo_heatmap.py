"""Render the AnoExpress logo (docs/logo.png) as an expression-style heatmap from docs/logo_img2.png.

Usage (from repo root):
    uv run --with numpy --with pillow --with matplotlib python workflow/notebooks/misc/logo_heatmap.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

DOCS = Path(__file__).resolve().parents[3] / "docs"
SEED = 7
NROWS = 42  # heatmap rows spanning the logo artwork
PAD = 4  # cells around the artwork, including the border ring
GAP_KEEP = 0.45  # fraction of the blank gap between mosquito and text to keep


def load_ink(path):
    """Return a 0-1 ink array (1 = black) cropped to the artwork, frame removed."""
    im = Image.open(path).convert("RGBA")
    im = Image.alpha_composite(Image.new("RGBA", im.size, "white"), im).convert("L")
    w, h = im.size
    im = im.crop((30, 30, w - 30, h - 30))  # drop the 23px frame
    ink = 1 - np.asarray(im, dtype=float) / 255
    ys, xs = np.where(ink > 0.5)
    ink = ink[ys.min() : ys.max() + 1, xs.min() : xs.max() + 1]

    # shrink the empty gap between the mosquito and the wordmark
    empty = np.where(ink.max(0) < 0.1)[0]
    runs = np.split(empty, np.where(np.diff(empty) != 1)[0] + 1)
    gap = max(runs, key=len)
    drop = gap[int(len(gap) * GAP_KEEP) :]
    return np.delete(ink, drop, axis=1)


def to_grid(ink, nrows):
    ncols = round(nrows * ink.shape[1] / ink.shape[0])
    im = Image.fromarray((ink * 255).astype(np.uint8))
    return np.asarray(im.resize((ncols, nrows), Image.BOX), dtype=float) / 255


def logo_values(grid, rng, split_col, lo=0.2, hi=0.45, sd=0.16):
    """Background ~ N(0, sd) 'fold changes'; mosquito strongly down, text strongly up, border ring down."""
    t = np.clip((grid - lo) / (hi - lo), 0, 1)
    bg = rng.normal(0, sd, grid.shape).clip(-0.55, 0.55)
    fg = rng.uniform(0.7, 1.0, grid.shape)
    fg[:, :split_col] *= -1
    vals = (1 - t) * bg + t * fg

    ring = -rng.uniform(0.7, 1.0, grid.shape)
    vals[[0, -1], :] = ring[[0, -1], :]
    vals[:, [0, -1]] = ring[:, [0, -1]]
    return vals


def render(vals, out, cmap="RdBu_r", facecolor="white", gap=0.6):
    nr, nc = vals.shape
    fig, ax = plt.subplots(figsize=(nc / 8, nr / 8))
    ax.pcolormesh(vals[::-1], cmap=cmap, vmin=-1, vmax=1, edgecolors=facecolor, linewidth=gap)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.03, facecolor=facecolor)
    plt.close(fig)
    print("wrote", out)


def main():
    ink = load_ink(DOCS / "logo_img2.png")
    grid = np.pad(to_grid(ink, NROWS), PAD)
    # first all-blank column after the mosquito marks the mosquito/text boundary
    split = int(np.where(grid[:, PAD:].max(0) < 0.1)[0][0] + PAD)
    render(logo_values(grid, np.random.default_rng(SEED), split), DOCS / "logo.png")


if __name__ == "__main__":
    main()
