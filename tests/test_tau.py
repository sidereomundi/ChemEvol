import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from pychem import tau


def test_tau_typical_values():
    assert tau(1.0, 1, 0.0) == 10000.0
    assert tau(0.5, 2, 0.0) == 50000.0


def test_tau_binmax_effect():
    assert tau(1.0, 1, -10.0) == 14000.0
