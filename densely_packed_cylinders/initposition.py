#!/usr/bin/env python3

"""Initial position + radius generation for the Donev C++ packer.

Port of initposition.m.
"""

import numpy as np


def init_position(rinit, read_filename="spheres_poly/read.dat",
                  dens0=0.01, seed=None):
    """
    Initialize positions and radii for the Donev packer and write
    read.dat.

    Parameters
    ----------
    rinit : (N,) array of inner radii (µm)
    read_filename : path to write the read.dat file
    dens0 : initial density parameter (dimensionless)
    seed : RNG seed (None = random)

    Returns
    -------
    x, y : (N,) initial positions in [0, 1)
    rs   : (N,) rescaled radii (dimensionless)
    """
    rng = np.random.default_rng(seed)
    rinit = np.asarray(rinit, dtype=float).ravel()
    N = len(rinit)

    # Rescale radii
    rmax = rinit.max()
    Rscale = rmax * np.sqrt(N / dens0)
    rinit_sorted = np.sort(rinit)[::-1]
    rs = rinit_sorted / Rscale

    # Random positions with no overlap (periodic unit box)
    x = np.zeros(N)
    y = np.zeros(N)

    def dist2(x1, y1, x2, y2):
        dx = min(abs(x1 - x2), 1 - abs(x1 - x2))
        dy = min(abs(y1 - y2), 1 - abs(y1 - y2))
        return dx**2 + dy**2

    for n in range(N):
        while True:
            xc = rng.random()
            yc = rng.random()
            ok = True
            for k in range(n):
                if dist2(xc, yc, x[n_prev := k], y[k]) <= (rs[n] + rs[k])**2:
                    ok = False
                    break
            if ok:
                x[n] = xc
                y[n] = yc
                break

    # Write read.dat
    with open(read_filename, "w") as f:
        f.write(f"{2}\n")                        # dimension
        f.write(f"{N}\n")                        # N
        f.write(f"{N}\n")                        # N
        f.write(f"{2 * rs.max():e}\n")           # max diameter
        f.write("10 01\n")                       # box unit vectors
        f.write("T T\n")                         # flags
        for n in range(N):
            f.write(f"{rs[n]:e} {rs[n]:e} {1.0:f} {x[n]:f} {y[n]:f}\n")

    return x, y, rs