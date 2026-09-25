#!/usr/bin/env python3

"""Create the axon lookup table used by the MC simulator."""

import numpy as np


def create_lookup_table(xc, yc, r, n=1024):
    """
    Port of createlookuptable.m.

    Parameters
    ----------
    n : grid resolution (n x n)
    xc, yc : cylinder centers in [0, 1)
    r : cylinder radii in [0, 1)

    Returns
    -------
    A : (n, n) uint32, axon labels (encoding up to 2 overlapping axons)
    B : (n, n) uint16, number of axons per pixel
    Nmax : int, base-10 scaling constant
    """
    A = np.zeros((n, n), dtype=np.int64)
    B = np.zeros((n, n), dtype=np.int64)
    N_ax = len(r)
    Nmax = 10 ** int(np.ceil(np.log10(N_ax + 1)))

    for i in range(N_ax):
        x0, y0, ri = xc[i], yc[i], r[i]

        iimin = int(np.ceil((x0 - ri) * n))
        iimax = int(np.ceil((x0 + ri) * n))
        jjmin = int(np.ceil((y0 - ri) * n))
        jjmax = int(np.ceil((y0 + ri) * n))

        for ii in range(iimin, iimax + 1):
            ti = ii % n
            for jj in range(jjmin, jjmax + 1):
                tj = jj % n
                if _inside_circle(ii, jj, x0, y0, ri, n):
                    if A[ti, tj] == 0:
                        A[ti, tj] = i + 1
                        B[ti, tj] = 1
                    elif A[ti, tj] < Nmax:
                        A[ti, tj] = A[ti, tj] * Nmax + (i + 1)
                        B[ti, tj] = 2
                    else:
                        print("Warning: >2 circles in one pixel.")

    return A.astype(np.uint32), B.astype(np.uint16), Nmax


def _inside_circle(ii, jj, xc, yc, r, n):
    x = max(ii - 1, min(xc * n, ii))
    y = max(jj - 1, min(yc * n, jj))
    d2 = (x - xc * n) ** 2 + (y - yc * n) ** 2
    return d2 < (r * n) ** 2