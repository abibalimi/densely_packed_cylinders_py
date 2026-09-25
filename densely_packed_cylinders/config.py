#!/usr/bin/env python3

"""Configuration dataclass for the densely packed cylinders pipeline."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


@dataclass
class Config:
    # --- Paths ---
    root: Path = Path(".")
    spheres_dir: Path = Path("spheres_poly")
    histogram_mat: Path = Path("diameter_histogram/CC_diameter_histogram.mat")
    output_dir: Path = Path("packing")

    # --- Histogram / axon radii ---
    roi: int = 1                    # 1 = genu, 2 = body, 3 = splenium
    shrink_factor: float = 1.0      # histological shrinkage
    n_target: int = 500             # target axon count before flooring

    # --- Packing parameters ---
    maxdensity: float = 0.76        # target outer packing fraction
    gratio: float = 0.585           # inner / outer radius ratio
    hardwallBC: int = 0             # 0 = periodic, 1 = hard wall
    dens0: float = 0.01             # initial density for placement

    # --- Lookup table ---
    lut_size: int = 1000            # n x n

    # --- RNG ---
    seed: Optional[int] = None

    # --- C++ build ---
    compile_flags: list = field(default_factory=lambda: [
        "g++", "-O3", "-ffast-math", "-o", "spheres",
        "neighbor.C", "spheres.C", "box.C", "sphere.C",
        "event.C", "heap.C", "read_input.C",
    ])
    binary_name: str = "spheres"

    # --- Derived paths (filled by pipeline) ---
    @property
    def spheres_path(self) -> Path:
        return self.root / self.spheres_dir

    @property
    def histogram_path(self) -> Path:
        return self.root / self.histogram_mat

    @property
    def output_path(self) -> Path:
        return self.root / self.output_dir

    @classmethod
    def from_yaml(cls, path: Path) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        # Convert path-like strings to Path
        for key in ("root", "spheres_dir", "histogram_mat", "output_dir"):
            if key in data and data[key] is not None:
                data[key] = Path(data[key])
        return cls(**data)

    def to_yaml(self, path: Path) -> None:
        data = {k: (str(v) if isinstance(v, Path) else v)
                for k, v in self.__dict__.items()}
        with open(path, "w") as f:
            yaml.safe_dump(data, f, sort_keys=False)