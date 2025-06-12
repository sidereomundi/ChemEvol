"""Input routines translated from ``src/io.f90``.

Each ``leggi*`` function reads a data file from the ``YIELDSBA``
directory using :func:`numpy.loadtxt`. The Fortran COMMON blocks
containing the arrays are stored inside an :class:`IORoutines` object.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent


@dataclass
class IORoutines:
    data: Dict[str, np.ndarray] = field(default_factory=dict)
    basepath: Path = field(init=False)

    def __post_init__(self) -> None:
        self.basepath = BASE_DIR / "YIELDSBA"

    def _load(self, filename: str, usecols: slice | None = None) -> np.ndarray:
        path = self.basepath / filename
        arr = np.loadtxt(path, skiprows=1)
        if usecols is not None:
            arr = arr[:, usecols]
        self.data[filename] = arr
        return arr

    def leggi(self) -> np.ndarray:
        """Load Ba yields from ``CristalloBa2.dat``."""
        return self._load("CristalloBa2.dat", slice(0, None))

    def leggiSr(self) -> np.ndarray:
        """Load Sr yields from ``CristalloSr.dat``."""
        return self._load("CristalloSr.dat", slice(0, None))

    def leggiY(self) -> np.ndarray:
        """Load Y yields from ``CristalloY.dat``."""
        return self._load("CristalloY.dat", slice(0, None))

    def leggiEu(self) -> np.ndarray:
        """Load Eu yields from ``CristalloEu.dat``."""
        return self._load("CristalloEu.dat", slice(0, None))

    def leggiZr(self) -> np.ndarray:
        """Load Zr yields from ``CristalloZr.dat``."""
        return self._load("CristalloZr.dat", slice(0, None))

    def leggiLa(self) -> np.ndarray:
        """Load La yields from ``CristalloLa.dat``."""
        return self._load("CristalloLa.dat", slice(0, None))

    def leggiRb(self) -> np.ndarray:
        """Load Rb yields from ``CristalloRb.dat``."""
        return self._load("CristalloRb.dat", slice(0, None))

    def leggiLi(self) -> np.ndarray:
        """Load Li yields from ``KarakasLi.dat``."""
        return self._load("KarakasLi.dat", slice(0, None))

    # ------------------------------------------------------------------
    # Additional helpers ------------------------------------------------

    def load_yield_grid(self, filenames: list[str], basepath: str = "DATI") -> "InterpolationData":
        """Read a set of yield tables forming a mass/metallicity grid.

        Parameters
        ----------
        filenames : list of str
            Sequence of files containing the yields for different
            metallicities. Each file should start with a line of the form
            ``"Z=<value>"`` which is used to build the metallicity grid.
        basepath : str, optional
            Directory containing the input files.
        """

        from .interpolation import InterpolationData

        mass_grid = None
        tables = []
        zetas = []

        for fname in filenames:
            path = BASE_DIR / basepath / fname
            with open(path, "r") as fh:
                first = fh.readline()
                if "=" in first:
                    zetas.append(float(first.split("=")[1]))
                else:
                    raise ValueError(f"Could not parse metallicity from {fname}")
            arr = np.loadtxt(path, skiprows=2)
            masses = arr[:, 0]
            yields = arr[:, 1:].T  # elements x mass
            if mass_grid is None:
                mass_grid = masses
            else:
                if not np.allclose(mass_grid, masses):
                    raise ValueError("Inconsistent mass grids in input files")
            tables.append(yields)

        W = np.stack(tables, axis=2)
        return InterpolationData(massa=mass_grid, zeta=np.array(zetas), W=W)

