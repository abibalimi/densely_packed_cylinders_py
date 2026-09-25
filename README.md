# densely_packed_cylinders_py

Python wrapper for the Donev C++ packing algorithm, generating densely
packed aligned cylinders with radii drawn from the corpus callosum axon
diameter histogram. Port of the MATLAB demo from the NYU-DiffusionMRI
`monte-carlo-simulation-recipes` repository.

## What it does

1. Loads the axon diameter histogram (`CC_diameter_histogram.mat`)
2. Writes `read.dat` and `inputN` for the Donev C++ packer
3. Compiles and runs `spheres`
4. Parses `write.dat` (axon centers and outer radii)
5. Builds the lookup table `A`, `B` for MC simulations
6. Saves `packing.npz`, `A.npy`, `B.npy`
