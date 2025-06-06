"""Input routines translated from ``src/io.f90``.

Each ``leggi*`` function reads a data file from the ``YIELDSBA``
directory using :func:`numpy.loadtxt`. The Fortran COMMON blocks
containing the arrays are stored inside an :class:`IORoutines` object.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Dict

import numpy as np


@dataclass
class IORoutines:
    basepath: str = "YIELDSBA"
    data: Dict[str, np.ndarray] = field(default_factory=dict)

    def _load(self, filename: str, usecols: slice | None = None) -> np.ndarray:
        path = os.path.join(self.basepath, filename)
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

