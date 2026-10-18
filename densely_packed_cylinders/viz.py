#!/usr/bin/env python3

"""Figures for the densely packed cylinders pipeline.

Reproduces:
    - Fig. 1: inner diameter histogram
    - Fig. 2: initial configuration (before packing)
    - Fig. 3: final axon packing (periodic tiling 3x3)
    - Fig. 4: lookup table
    - Fig. 5: axon binary mask
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# ------------------------------------------------------------------
# Fig. 1 — inner diameter histogram
# ------------------------------------------------------------------
def plot_diameter_histogram(ri, roi_name="CC", shrink_factor=1.0,
                            edges=None, outfile=None):
    """Bar plot of the inner diameter histogram."""
    if edges is None:
        edges = np.arange(0, 9.01, 0.2) * shrink_factor

    counts, _ = np.histogram(2 * ri, bins=edges)
    pct = counts / counts.sum() * 100

    fig, ax = plt.subplots(figsize=(6, 3), dpi=150)
    ax.bar(edges[1:], pct, width=edges[1] - edges[0], align="edge",
           color="steelblue", edgecolor="k", linewidth=0.4)
    ax.set_xlim(0, 9)
    ax.set_ylim(0, 30)
    ax.set_xticks(range(0, 10))
    ax.set_yticks(range(0, 31, 5))
    ax.set_xlabel("Inner Diameter (µm)")
    ax.set_ylabel("Frequency (%)")
    ax.set_title(f"{roi_name}, shrinkage factor = {shrink_factor:.2f}")
    ax.set_box_aspect(0.5)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    if outfile is not None:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, ax


# ------------------------------------------------------------------
# Fig. 2 — initial configuration
# ------------------------------------------------------------------
def plot_initial(x, y, rs, outfile=None):
    """Initial configuration with rescaled radii."""
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    for xi, yi, ri in zip(x, y, rs):
        ax.add_patch(Circle((xi, yi), ri, fill=False,
                            edgecolor="steelblue", linewidth=0.8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Initial configuration")
    fig.tight_layout()
    if outfile is not None:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, ax


# ------------------------------------------------------------------
# Fig. 3 — final axon packing (periodic 3x3 tiling)
# ------------------------------------------------------------------
def plot_packing(xc, yc, rc, outfile=None, title="Axon Packing"):
    """Plot the final packing with 3x3 periodic tiling."""
    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for xi, yi, ri in zip(xc + i, yc + j, rc):
                ax.add_patch(Circle((xi, yi), ri, fill=False,
                                    edgecolor="k", linewidth=1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title)
    fig.tight_layout()
    if outfile is not None:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, ax


# ------------------------------------------------------------------
# Fig. 4 — lookup table
# ------------------------------------------------------------------
def plot_lookup_table(A, N_axons, outfile=None, title="Lookup Table"):
    """
    Plot the lookup table: each pixel colored by axon label.
    Background pixels black, multi-axon pixels white.
    """
    Nmax = 10 ** int(np.ceil(np.log10(N_axons + 1)))

    # Start from the raw A (labels, not normalized)
    A2 = A.astype(float).copy()

    # Map single-axon labels to their axon index (1..N_axons)
    # Multi-axon pixels (> Nmax) get the axon of the outermost circle
    single = (A > 0) & (A < Nmax)
    multi = A >= Nmax

    A2[single] = A[single]            # 1..N_axons
    A2[multi] = A[multi] % Nmax       # keep the outer axon label

    # Normalize to [0, 1] for imshow
    A2 = A2 / N_axons

    # Background black, multi-axon white
    A2[A == 0] = np.nan
    A2[multi] = 1.0

    fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
    cmap = plt.cm.viridis.copy()
    cmap.set_bad("black")
    ax.imshow(np.rot90(A2), cmap=cmap, interpolation="nearest",
              vmin=0, vmax=1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title)
    fig.tight_layout()
    if outfile is not None:
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, ax

# ------------------------------------------------------------------
# Fig. 5 — Axon binary mask
# ------------------------------------------------------------------
def axon_binary_mask(xc, yc, rc, n=1024, outfile=None):
    """
    Binary mask of the packing: 1 inside an axon, 0 in the EAS.
    Periodic, wrapped into [0, 1).
    
    Returns a (n, n) uint8 array; 1 = axon, 0 = EAS. 
    Save as .npy with np.save(out/"axon_mask.npy", mask) if you want it for downstream MC.
    """
    mask = np.zeros((n, n), dtype=np.uint8)
    xs = np.arange(n) / n
    Y, X = np.meshgrid(xs, xs, indexing="ij")

    for xi, yi, ri in zip(xc, yc, rc):
        dx = X - xi
        dy = Y - yi
        dx -= np.round(dx)
        dy -= np.round(dy)
        mask[(dx**2 + dy**2) <= ri**2] = 1

    if outfile is not None:
        fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
        ax.imshow(mask, origin="lower", cmap="gray_r",
                  extent=[0, 1, 0, 1])
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("Axon binary mask")
        fig.tight_layout()
        fig.savefig(outfile, dpi=200, bbox_inches="tight")
        plt.close(fig)

    return mask