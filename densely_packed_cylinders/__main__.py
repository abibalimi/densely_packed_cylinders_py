#!/usr/bin/env python3

"""Entry point:

    python -m densely_packed_cylinders --root . 
    --histo diameter_histogram/CC_diameter_histogram.mat 
    --roi 1 --phi 0.76 --gratio 0.585 --seed 7 
"""

import argparse
from .generate import generate_packing


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", required=True,
                   help="Project root containing spheres_poly/")
    p.add_argument("--histo", required=True,
                   help="Path to CC_diameter_histogram.mat")
    p.add_argument("--roi", type=int, default=1)
    p.add_argument("--shrink", type=float, default=1.0)
    p.add_argument("--phi", type=float, default=0.76)
    p.add_argument("--gratio", type=float, default=0.585)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--res", type=int, default=None)
    args = p.parse_args()

    generate_packing(
        root=args.root,
        histo_mat=args.histo,
        roi=args.roi,
        shrink_factor=args.shrink,
        maxdensity=args.phi,
        gratio=args.gratio,
        seed=args.seed,
        L=args.res)


if __name__ == "__main__":
    main()