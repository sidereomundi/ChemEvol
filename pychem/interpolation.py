"""Interpolation utilities translated from ``src/interpolation.f90``.

Only a minimal subset is ported. The Fortran COMMON blocks storing
mass grids and yield tables are represented by attributes of the
:class:`InterpolationData` container. The heavy use of global state in
Fortran is replaced by object-oriented design in Python.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple


@dataclass
class InterpolationData:
    """Container for the stellar yield grid.

    Parameters mirror the arrays used in ``src/interpolation.f90`` but are
    represented with NumPy arrays instead of COMMON blocks.
    """

    massa: np.ndarray = field(default_factory=lambda: np.empty(0))
    zeta: np.ndarray = field(default_factory=lambda: np.empty(0))
    W: np.ndarray = field(default_factory=lambda: np.empty((0, 0, 0)))


@dataclass
class Interpolator:
    data: InterpolationData

    def polint(
        self, xa: np.ndarray, ya: np.ndarray, x: float
    ) -> Tuple[float, float]:
        """Polynomial interpolation of order len(xa)-1.

        This is a direct translation of the `polint` subroutine. NumPy is used
        to perform vectorized operations.
        """
        n = len(xa)
        c = ya.astype(float).copy()
        d = ya.astype(float).copy()
        ns = np.argmin(np.abs(x - xa))
        y = ya[ns]
        ns -= 1
        dy = 0.0
        for m in range(1, n):
            for i in range(n - m):
                ho = xa[i] - x
                hp = xa[i + m] - x
                w = c[i + 1] - d[i]
                den = ho - hp
                if den == 0:
                    raise ZeroDivisionError("polint failure")
                den = w / den
                d[i] = hp * den
                c[i] = ho * den
            if 2 * (ns + 1) < n - m:
                dy = c[ns + 1]
            else:
                dy = d[ns]
                ns -= 1
            y += dy
        return y, dy

    def interp(
        self, mass: float, zeta: float, binmax: float
    ) -> Tuple[np.ndarray, float]:
        """Bilinear interpolation of the yield table.

        This is a greatly simplified version of the Fortran ``interp`` routine.
        The yield grid :data:`InterpolationData.W` is assumed to have the shape
        ``(n_elements, n_mass, n_zeta)`` where the second and third axes
        correspond to the mass and metallicity grids stored in
        :data:`InterpolationData.massa` and :data:`InterpolationData.zeta`.

        Parameters
        ----------
        mass : float
            Stellar mass ``H`` in the original code.
        zeta : float
            Metallicity value.
        binmax : float
            Unused in this minimal implementation. It controls some special
            behaviour in the Fortran version related to binary stars.
        """

        masses = self.data.massa
        zetas = self.data.zeta
        W = self.data.W

        if masses.size < 2 or zetas.size < 2:
            raise ValueError("Interpolation grid is not initialised")

        # Locate the bracketing indices for ``mass`` and ``zeta``
        i = np.searchsorted(masses, mass) - 1
        j = np.searchsorted(zetas, zeta) - 1
        i = np.clip(i, 0, len(masses) - 2)
        j = np.clip(j, 0, len(zetas) - 2)

        m1, m2 = masses[i], masses[i + 1]
        z1, z2 = zetas[j], zetas[j + 1]
        fm = 0.0 if m2 == m1 else (mass - m1) / (m2 - m1)
        fz = 0.0 if z2 == z1 else (zeta - z1) / (z2 - z1)

        # Bilinear interpolation for each element
        q = (
            W[:, i, j] * (1 - fm) * (1 - fz)
            + W[:, i + 1, j] * fm * (1 - fz)
            + W[:, i, j + 1] * (1 - fm) * fz
            + W[:, i + 1, j + 1] * fm * fz
        )

        # ``Hecore`` in the Fortran code is related to the He core mass.
        # We use a very rough approximation here just to populate the
        # variable.
        hecore = 0.1 * mass

        return q, float(hecore)
