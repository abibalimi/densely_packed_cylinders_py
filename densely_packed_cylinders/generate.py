#!/usr/bin/env python3

"""Orchestrate the Donev C++ packer: full pipeline from histogram to packing."""

import os
import subprocess
import numpy as np

from .initposition import init_position
from .histogram import load_cc_histogram
from .lookuptable import create_lookup_table
from .viz import (
    axon_binary_mask, plot_diameter_histogram, plot_initial, plot_packing, plot_lookup_table
)

INPUT_TEMPLATE = """int eventspercycle = 20;
int N = {N};
double initialpf = 0.01;
double maxpf = {maxpf};
double temp = 0.2;
double growthrate = 0.001;
double maxpressure = 10000000;
double maxcollisionrate = 100000000;
double bidispersityratio = 0.5;
double bidispersityfraction = 0.5;
double massratio = 1.;
int hardwallBC = {hardwallBC};
char* readfile = read.dat
char* writefile = write.dat
char* datafile = stats.dat
"""


def generate_packing(root, histo_mat, roi=1, shrink_factor=1.0,
                     maxdensity=0.76, gratio=0.585, hardwallBC=0,
                     seed=9, naxons=500):
    """
    Full pipeline: histogram → initial positions → Donev packer → packing.

    Parameters
    ----------
    root : project root containing spheres_poly/ and packing/ dirs
    histo_mat : path to CC_diameter_histogram.mat
    roi : ROI index (1 = genu)
    shrink_factor : histological shrink factor
    maxdensity : target outer packing fraction (e.g. 0.76)
    gratio : g-ratio (inner / outer)
    hardwallBC : 0 periodic, 1 hard wall
    seed : RNG seed
    naxons : number of axons to generate
    Returns
    -------
    dict with keys: xc, yc, rc, res, gratio, ri, roi_name
    """
    root = os.path.abspath(root)
    spheres_dir = os.path.join(root, "spheres_poly")
    packing_dir = os.path.join(root, "packing")
    os.makedirs(packing_dir, exist_ok=True)

    # ---- 1. Load histogram, build per-axon inner radii ----
    ri, roi_name = load_cc_histogram(histo_mat, roi, naxons, shrink_factor)
    N = len(ri)

    # ---- 2. Init positions and write read.dat ----
    read_path = os.path.join(spheres_dir, "read.dat")
    rinit = ri / gratio  # initial OUTER radii (before rescale)
    xinit, yinit, rs = init_position(rinit, read_path, seed=seed)

    # ---- 3. Write inputN ----
    inputN = os.path.join(spheres_dir, "inputN")
    with open(inputN, "w") as f:
        f.write(INPUT_TEMPLATE.format(
            N=N, maxpf=maxdensity, hardwallBC=hardwallBC))

    # ---- 4. Compile and run the C++ packer ----
    binary = os.path.join(spheres_dir, "spheres")
    if not os.path.exists(binary):
        subprocess.run(
            ["g++", "-O3", "-ffast-math", "-o", "spheres",
             "neighbor.C", "spheres.C", "box.C", "sphere.C",
             "event.C", "heap.C", "read_input.C"],
            cwd=spheres_dir, check=True)

    subprocess.run(["./spheres", "inputN"], cwd=spheres_dir, check=True)

    # ---- 5. Parse write.dat ----
    write_path = os.path.join(spheres_dir, "write.dat")
    M = np.loadtxt(write_path, skiprows=6)
    xc = M[:, 0]
    yc = M[:, 1]
    rc = M[:, 2] / 2.0

    # ---- 6. Compute physical FOV size ----
    res = _physical_fov(histo_mat, roi, gratio, rc.mean(), shrink_factor)

    # ---- 7. Save ----
    np.savez(os.path.join(packing_dir, "packing.npz"),
             xc=xc, yc=yc, rc=rc, gratio=gratio, res=res)
    
    r_min_axon_physical = (rc * res).min()
    n_grid = compute_n_grid(res, r_min_axon_physical)
    print(f"r_min physical = {(rc * res).min():.3f} µm,  res = {res:.2f} µm,  n_grid = {n_grid}")
    A,B,_ = create_lookup_table(xc, yc, rc, n_grid)
    np.save(os.path.join(packing_dir, "B.npy"), B)
    np.save(os.path.join(packing_dir, "A.npy"), A)

    print(f"Packing complete:")
    print(f"  ROI: {roi_name}")
    print(f"  N cylinders: {N}")
    print(f"  Cylindrical volume fraction: {np.sum(np.pi * rc**2):.4f}")
    print(f"  Intra-axonal volume fraction: "
          f"{np.sum(np.pi * (rc * gratio)**2):.4f}")
    print(f"  FOV size (physical): {res:.4f} µm")

    # ---- 8. Visualize ----
    out = os.path.join(packing_dir, "figures")
    os.makedirs(out, exist_ok=True)
    plot_diameter_histogram(ri, roi_name=roi_name,
                            shrink_factor=shrink_factor,
                            outfile=out+"/histogram.png")

    plot_initial(xinit, yinit, rs, outfile=out+"/initial.png")

    plot_packing(xc, yc, rc, outfile=out+"/packing.png")
    
    axon_binary_mask(xc, yc, rc, outfile=out+"/axon_binary_mask.png")
    
    plot_lookup_table(A, N_axons=len(rc), outfile=out+"/lookup_table.png")

    return {"xc": xc, "yc": yc, "rc": rc, "res": res,
                "gratio": gratio, "ri": ri, "roi_name": roi_name}

    
def compute_n_grid(L, r_min, safety=2.0):
    """Next power of 2 >= safety * L / r_min."""
    n = int(np.ceil(safety * L / r_min))
    return 1 << (n - 1).bit_length()


def _physical_fov(histo_mat, roi, gratio, rc_mean, shrink_factor):
    """Compute physical FOV size (µm) from histogram and packing."""
    from scipy.io import loadmat
    m = loadmat(histo_mat, squeeze_me=True, struct_as_record=False)
    entry = m["cc"][roi - 1]
    diam = np.asarray(entry.diameter, dtype=float).ravel() * shrink_factor
    pct = np.asarray(entry.frequency, dtype=float).ravel()
    ri_per_bin = diam / 2.0
    mean_r_outer = np.sum((ri_per_bin / gratio) * pct) / np.sum(pct)
    return mean_r_outer / rc_mean