import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
import pytest
from pychem.interpolation import Interpolator, InterpolationData


def test_simple_interpolation():
    masses = np.array([1.0, 3.0])
    zetas = np.array([0.0, 2.0])
    W = np.empty((1, 2, 2))
    for i, m in enumerate(masses):
        for j, z in enumerate(zetas):
            W[0, i, j] = m + z
    interp = Interpolator(InterpolationData(massa=masses, zeta=zetas, W=W))
    q, hecore = interp.interp(2.0, 1.0, 0.0)
    assert q[0] == pytest.approx(3.0)
    assert hecore == pytest.approx(0.2)
