# densely_packed_cylinders_py

Python wrapper for the [Donev](https://doi.org/10.1016/j.jcp.2004.08.014) C++ packing algorithm, generating randomly 
and densely packed cylinders with the freedom of tuning axonal water fraction, inner diameter distribution, and g-ratio. Axonal diameter histogram comes from [Aboitiz, et al., Brain Research, 1992](https://doi.org/10.1016/0006-8993(92)90178-C).
Part of the MATLAB demo from the NYU-DiffusionMRI
[monte-carlo-simulation-recipes](https://github.com/NYU-DiffusionMRI/monte-carlo-simulation-recipes.git) repository.

## What it does

1. Loads the axon diameter histogram (`CC_diameter_histogram.mat`)
2. Writes `read.dat` and `inputN` for the Donev C++ packer
3. Compiles and runs `spheres`
4. Parses `write.dat` (axon centers and outer radii)
5. Builds the lookup table `A`, `B` for MC simulations
6. Saves `packing.npz`, `A.npy`, `B.npy`


## Structure

```text
densely_packed_cylinders_py/
|-- README.md
|-- requirements.txt
|-- pyproject.toml
|-- .gitignore
|-- densely_packed_cylinders/
|   |-- __init__.py
|   |-- __main__.py
|   |-- histogram.py
|   |-- initposition.py
|   |-- generate.py
|   |-- lookuptable.py
|   `-- viz.py
|-- spheres_poly/
|   |-- input
|   |-- neighbor.C
|   |-- neighbor.h
|   |-- spheres.C
|   |-- spheres.h
|   |-- box.C
|   |-- box.h
|   |-- sphere.C
|   |-- sphere.h
|   |-- event.C
|   |-- event.h
|   |-- heap.C
|   |-- heap.h
|   |-- read_input.C
|   `-- read_input.h
|-- diameter_histogram/
|    `-- CC_diameter_histogram.mat
|-- packing/
|   |-- figures/
|    `-- histogram.png
|   |-- A.npy
|   |-- B.npy
`   |-- packing.npy

```

## Install

```bash
bash setup.sh
```

## Run

```bash
python -m densely_packed_cylinders \
    --root . \
    --histo diameter_histogram/CC_diameter_histogram.mat \
    --roi 1 \
    --phi 0.5 \
    --gratio 0.585 \
    --naxons 4000 \
    --seed 9
```

## Arguments

| Flag | Meaning | Default |
| :--- | :--- | :--- |
| `--root` | project root containing spheres_poly/ | *required* |
| `--histo` | path to CC_diameter_histogram.mat | *required* |
| `--roi` | ROI index (1 = genu) | `1` |
| `--shrink` | histological shrink factor | `1.0` |
| `--phi` | target packing fraction | `0.5` |
| `--gratio` | inner / outer radius ratio | `0.585` |
| `--naxons` | number of axons | `5000` |
| `--seed` | RNG seed | `9` |

## Outputs

Saved to `packing/`:

*   **packing.npz** — xc, yc, rc, gratio, L, n_grid
*   **A.npy, B.npy** — lookup tables (A = labels, B = occupancy)
*   **histogram.png, packing.png, lookup_table.png, axon_mask.png**

## Notes

*   `xc`, `yc`, `rc` are normalized to `[0, 1)`; physical radii are `rc * L`.
*   `B` is the combined 0/1 lookup table for MC walkers (`0` = EAS, `>0` = obstacle).
*   `L` is the physical box size (FOV), inferred from the histogram and the packing. 
*   `n_grid` is derived automatically from `L` and the smallest axon radius (`n_grid = next_pow2(2 * L / r_min)`).
*   The package has no downstream dependencies: it produces a packing, not a simulation.

## References

* **Random packing generation**
  - [Donev, et al., J Comput Phys, 2005](https://doi.org/10.1016/j.jcp.2004.08.014) and [code](https://cims.nyu.edu/~donev/Packing/C++/)

* **Axonal diameter histogram in corpus callosum**
  - [Aboitiz, et al., Brain Research, 1992](https://doi.org/10.1016/0006-8993(92)90178-C)

