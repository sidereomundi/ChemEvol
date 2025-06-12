"""Main module providing the ``MinGCE`` routine in Python.

The Fortran implementation heavily relies on COMMON blocks and global
state. Here we encapsulate the state inside the :class:`GCEModel` class.
Only a very small subset of the original logic is implemented to
illustrate the structure.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

from .interpolation import Interpolator, InterpolationData
from .io_routines import IORoutines
from .tau import tau


@dataclass
class GCEModel:
    io: IORoutines = field(default_factory=IORoutines)
    interpolator: Interpolator = field(init=False)

    def __post_init__(self) -> None:
        """Initialise the interpolation tables from the data files."""

        # Use a few of the ``Cris`` tables shipped with the repository to build
        # an example mass/metallicity grid.  In the original Fortran code many
        # more files are read at this stage.
        files = ["Cris0.dat", "Cris004.dat", "Cris02.dat"]
        interp_data = self.io.load_yield_grid(files, basepath="DATI")
        self.interpolator = Interpolator(interp_data)

    def MinGCE(
        self,
        endoftime: int,
        sigmat: float,
        sigmah: float,
        psfr: float,
        pwind: float,
        delay: int,
        time_wind: int,
    ) -> None:
        """Very small Python reimplementation of ``MinGCE``.

        Only a tiny fraction of the original algorithm is reproduced here: a
        few yield tables are loaded, the lifetime function :func:`tau` is
        evaluated and the interpolation routine is exercised.  The parameters
        are kept for API compatibility but are not used.
        """

        # Load one of the heavy element yield tables
        ba = self.io.leggi()

        # Compute a stellar lifetime for a 5 Msun star
        lifetime = tau(5.0, tautype=1, binmax=0.0)

        # Interpolate the lighter element yields at the same mass and a
        # metallicity of 0.004
        q, hecore = self.interpolator.interp(5.0, 0.004, 0.0)

        print("Loaded Ba yields:", ba.shape)
        print("Lifetime for 5 Msun:", lifetime)
        print("Interpolated Q vector (first 5):", q[:5])
        print("He core mass:", hecore)
