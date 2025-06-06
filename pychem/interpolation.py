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
    """Container for yield tables used in interpolation."""

    massa: np.ndarray = field(default_factory=lambda: np.empty(0))
    W: np.ndarray = field(default_factory=lambda: np.empty((0, 0, 0)))


@dataclass
class Interpolator:
    data: InterpolationData

    def polint(self, xa: np.ndarray, ya: np.ndarray, x: float) -> Tuple[float, float]:
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

    def interp(self, mass: float, zeta: float, binmax: float) -> Tuple[np.ndarray, float]:
        """Placeholder for the complex `interp` Fortran routine.

        Parameters map directly to the Fortran variables ``H``, ``zeta`` and
        ``BINMAX``. The actual interpolation of stellar yields is highly
        domain specific. Here we simply return zeros with the correct
        dimensionality.
        """
        q = np.zeros(33)
        hecore = 0.0
        # Real implementation would select the proper grid points from
        # ``self.data.W`` and call :func:`polint` several times.
        return q, hecore

