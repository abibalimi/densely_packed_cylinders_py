#!/usr/bin/env python3

"""Load corpus callosum diameter histogram from MATLAB .mat file."""

import numpy as np
from scipy.io import loadmat


def load_cc_histogram(mat_path, roi=1, naxons=1000, shrink_factor=1.0):
    """
    Load inner axonal diameter histogram from the CC .mat file.

    Parameters
    ----------
    mat_path : path to CC_diameter_histogram.mat
    roi : ROI index (1 = genu of corpus callosum)
    naxons : number of axons to generate
    shrink_factor : multiply histological diameters by this factor

    Returns
    -------
    ri : (M,) inner radii (µm), one entry per axon
    roi_name : name of the ROI
    """
    m = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    cc = m["cc"]
    entry = cc[roi - 1]

    pct = np.asarray(entry.frequency, dtype=float).ravel()
    diam = np.asarray(entry.diameter, dtype=float).ravel()

    # Number of axons per bin (integer counts, floored)
    N_target = naxons
    ct = np.floor(N_target * pct).astype(int)

    # Build the per-axon inner radius list
    ri = np.zeros(ct.sum())
    idx = 0
    for k in range(len(ct)):
        ri[idx:idx + ct[k]] = diam[k] / 2.0 * shrink_factor
        idx += ct[k]

    return ri, str(entry.name)