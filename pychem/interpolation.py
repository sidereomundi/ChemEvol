"""Interpolation routines translated from ``src/interpolation.f90``."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np

from .io_routines import FortranState


@dataclass
class InterpolationData:
    """Compatibility container for tests and simple interpolation mode."""

    massa: np.ndarray = field(default_factory=lambda: np.empty(0))
    zeta: np.ndarray = field(default_factory=lambda: np.empty(0))
    W: np.ndarray = field(default_factory=lambda: np.empty((0, 0, 0)))


@dataclass
class Interpolator:
    state: FortranState | None = None
    data: InterpolationData | None = None

    def __post_init__(self) -> None:
        if isinstance(self.state, InterpolationData):
            # Backward-compatible positional construction used by tests:
            # Interpolator(InterpolationData(...))
            self.data = self.state
            self.state = None
        if self.state is None and self.data is None:
            self.data = InterpolationData()

    def polint(self, xa: np.ndarray, ya: np.ndarray, x: float) -> Tuple[float, float]:
        n = len(xa)
        c = ya.astype(float).copy()
        d = ya.astype(float).copy()
        ns = int(np.argmin(np.abs(x - xa)))
        y = float(ya[ns])
        ns -= 1
        dy = 0.0
        for m in range(1, n):
            for i in range(n - m):
                ho = xa[i] - x
                hp = xa[i + m] - x
                w = c[i + 1] - d[i]
                den = ho - hp
                if den == 0.0:
                    return y, dy
                den = w / den
                d[i] = hp * den
                c[i] = ho * den
            if 2 * (ns + 1) < n - m:
                dy = float(c[ns + 1])
            else:
                dy = float(d[ns])
                ns -= 1
            y += dy
        return float(y), float(dy)

    def _simple_interp(self, mass: float, zeta: float) -> Tuple[np.ndarray, float]:
        masses = self.data.massa
        zetas = self.data.zeta
        W = self.data.W

        i = np.searchsorted(masses, mass) - 1
        j = np.searchsorted(zetas, zeta) - 1
        i = np.clip(i, 0, len(masses) - 2)
        j = np.clip(j, 0, len(zetas) - 2)

        m1, m2 = masses[i], masses[i + 1]
        z1, z2 = zetas[j], zetas[j + 1]
        fm = 0.0 if m2 == m1 else (mass - m1) / (m2 - m1)
        fz = 0.0 if z2 == z1 else (zeta - z1) / (z2 - z1)

        q = (
            W[:, i, j] * (1 - fm) * (1 - fz)
            + W[:, i + 1, j] * fm * (1 - fz)
            + W[:, i, j + 1] * (1 - fm) * fz
            + W[:, i + 1, j + 1] * fm * fz
        )
        return q, 0.1 * mass

    def litio(self, zcerc2: float, massa: float) -> float:
        state = self.state
        indice_m = 1
        for i in range(1, 16):
            if massa >= state.massaLi[i]:
                indice_m = i
        indice_m = min(indice_m, 14)

        massav = np.array([state.massaLi[indice_m], state.massaLi[indice_m + 1]], dtype=float)
        yd = np.array([state.YLi[indice_m, 1], state.YLi[indice_m + 1, 1]], dtype=float)
        qli, _ = self.polint(massav, yd, massa)
        return float(qli)

    def _bario_component(self, grid: np.ndarray, indice_m: int, indice_z: int, massa: float, zv: np.ndarray, massav: np.ndarray, zcerc: float) -> float:
        yd = np.array([grid[indice_m, indice_z], grid[indice_m, indice_z + 1]], dtype=float)
        yd2 = np.array([grid[indice_m + 1, indice_z], grid[indice_m + 1, indice_z + 1]], dtype=float)
        q1, _ = self.polint(zv, yd, zcerc)
        q2, _ = self.polint(zv, yd2, zcerc)
        y, _ = self.polint(massav, np.array([q1, q2], dtype=float), massa)
        return float(y)

    def bario(self, zcerc2: float, massa: float) -> Tuple[float, float, float, float, float, float, float]:
        state = self.state
        zcerc = float(zcerc2)

        indice_m = 1
        indice_z = 1
        for i in range(1, 6):
            if massa >= state.massaba[i]:
                indice_m = i
        indice_m = min(indice_m, 4)

        if zcerc < state.zbario[1]:
            zcerc = state.zbario[1]
        if zcerc > state.zbario[9]:
            zcerc = state.zbario[9]

        for i in range(1, 10):
            if zcerc >= state.zbario[i]:
                indice_z = i
        indice_z = min(indice_z, 8)

        zv = np.array([state.zbario[indice_z], state.zbario[indice_z + 1]], dtype=float)
        massav = np.array([state.massaba[indice_m], state.massaba[indice_m + 1]], dtype=float)

        qba = self._bario_component(state.ba, indice_m, indice_z, massa, zv, massav, zcerc)
        qy = self._bario_component(state.yt, indice_m, indice_z, massa, zv, massav, zcerc)
        qsr = self._bario_component(state.sr, indice_m, indice_z, massa, zv, massav, zcerc)
        qeu = self._bario_component(state.eu, indice_m, indice_z, massa, zv, massav, zcerc)
        qzr = self._bario_component(state.zr, indice_m, indice_z, massa, zv, massav, zcerc)
        qla = self._bario_component(state.la, indice_m, indice_z, massa, zv, massav, zcerc)
        qrb = self._bario_component(state.rb, indice_m, indice_z, massa, zv, massav, zcerc)
        return qba, qsr, qy, qeu, qzr, qla, qrb

    def interp(self, mass: float, zeta: float, binmax: float) -> Tuple[np.ndarray, float]:
        if self.state is None:
            return self._simple_interp(mass, zeta)

        s = self.state
        elem = 33
        nmax = 23

        H = float(mass)
        q = np.zeros(elem + 1, dtype=float)
        q1 = np.zeros(elem + 1, dtype=float)
        q2 = np.zeros(elem + 1, dtype=float)
        cosn2 = np.ones(nmax + 1, dtype=float)

        qia = np.zeros(elem + 1, dtype=float)
        qia[1] = 0.0
        qia[2] = 0.048
        qia[3] = 0.143
        qia[4] = 1.16e-6
        qia[5] = 1.40e-6
        qia[6] = 0.00202
        qia[7] = 0.0425 * 1.2
        qia[8] = 0.154
        qia[9] = 0.6
        qia[10] = 0.0
        qia[11] = 0.0
        qia[12] = 0.0846
        qia[13] = 0.0119
        qia[14] = 0.0
        qia[15] = 6.345e-4
        qia[16] = 1.19e-3
        qia[17] = 1.28e-5
        qia[18] = 4.16e-4
        qia[19] = 7.49e-5
        qia[20] = 5.67e-3
        qia[21] = 6.21e-3
        qia[22] = 1.56e-3
        qia[23] = 1.83e-2

        cosn2[3] = 1.0
        cosn2[7] = 6.0
        cosn2[8] = 1.0
        cosn2[13] = 1.0
        cosn2[23] = 1.0
        cosn2[16] = 2.0
        cosn2[17] = 1.15
        cosn2[18] = 15.0
        cosn2[20] = 3.0
        cosn2[22] = 0.3
        cosn2[19] = 4.5

        if binmax >= 0.0:
            if binmax > 0.0:
                ratio = (H + binmax) / binmax
                H = binmax
            else:
                ratio = 0.0

            aa = np.zeros(16, dtype=float)
            if H > 8.0:
                aa[1:6] = [0.0, 0.02e-4, 0.02e-2, 0.02e-1, 0.02]
            else:
                aa[1:6] = [0.0, 0.004, 0.008, 0.02, 0.04]

            z = 1
            for j in range(1, 6):
                if aa[j] <= zeta:
                    z = j
            zz = z if z == 5 else z + 1
            met = np.array([aa[z], aa[zz]], dtype=float)

            k = 1
            for j in range(1, s.ninputyield + 1):
                if s.massa[j] < H:
                    k = j
            k = max(1, min(k, s.ninputyield - 1))
            kk = k + 1
            mm = np.array([s.massa[k], s.massa[kk]], dtype=float)

            for j in range(1, 24):
                dd = np.array([s.W[j, k, z], s.W[j, kk, z]], dtype=float)
                q1[j], _ = self.polint(mm, dd, H)
                q1[j] *= cosn2[j]

            if z < 5:
                for j in range(1, 24):
                    dd = np.array([s.W[j, k, zz], s.W[j, kk, zz]], dtype=float)
                    q2[j], _ = self.polint(mm, dd, H)
                    q2[j] *= cosn2[j]
                for j in range(1, 24):
                    dd = np.array([q1[j], q2[j]], dtype=float)
                    q[j], _ = self.polint(met, dd, zeta)
            else:
                q[1:24] = q1[1:24]

            aa[6:11] = [0.0, 1.0e-8, 1.0e-5, 4.0e-3, 2.0e-2]
            aa[11:16] = [0.0, 0.001, 0.004, 0.008, 2.0e-2]

            k = 1
            for j in range(1, 14):
                if s.massac[j] < H:
                    k = j
            k = max(1, min(k, 12))
            kk = k + 1
            mm = np.array([s.massac[k], s.massac[kk]], dtype=float)

            z = 6
            for j in range(6, 11):
                if aa[j] <= zeta:
                    z = j
            zz = min(z + 1, 10)

            for j in range(1, 15):
                dd = np.array([s.W[j, k, z], s.W[j, kk, z]], dtype=float)
                q1[j], _ = self.polint(mm, dd, H)

            if z < 10:
                for j in range(1, 15):
                    dd = np.array([s.W[j, k, zz], s.W[j, kk, zz]], dtype=float)
                    q2[j], _ = self.polint(mm, dd, H)

            aa[1:6] = [0.0, 0.001, 0.004, 0.02, 0.05]
            z = 1
            for j in range(1, 6):
                if aa[j] <= zeta:
                    z = j
            zz = z if z == 5 else z + 1
            met = np.array([aa[z], aa[zz]], dtype=float)

            k = 1
            for j in range(1, s.ninputyield + 1):
                if s.massa[j] < H:
                    k = j
            k = max(1, min(k, s.ninputyield - 1))
            kk = k + 1
            mm = np.array([s.massa[k], s.massa[kk]], dtype=float)

            for idx in (9, 21, 13):
                dd = np.array([s.W[idx, k, z], s.W[idx, kk, z]], dtype=float)
                q1[idx], _ = self.polint(mm, dd, H)
            if z < 5:
                for idx in (9, 21, 13):
                    dd = np.array([s.W[idx, k, zz], s.W[idx, kk, zz]], dtype=float)
                    q2[idx], _ = self.polint(mm, dd, H)
                dd = np.array([q1[9], q2[9]], dtype=float)
                q[9], _ = self.polint(met, dd, zeta)
            else:
                q[9] = q1[9]

            qbar = 1.0e-30
            qeu2 = 1.0e-30
            qla = 1.0e-30
            qsrr = 1.0e-30
            qy = 1.0e-30
            qzr = 1.0e-30
            qrb = 1.0e-30

            mlow1 = 10.0
            mup1 = 30.0
            value1 = 0.8e-6

            cost_ba = 1.0
            cost_sr = 3.16 * 88.0 / 138.0
            cost_la = 0.136
            cost_zr = 2.53 * 90.0 / 138.0
            cost_eu = 0.117 * 151.0 / 138.0
            cost_y = 1.625 * 89.0 / 138.0 / 3.0
            cost_rb = 3.16 * 86.0 / 138.0

            if mlow1 <= H <= mup1:
                qbar = value1 * cost_ba
                qeu2 = value1 * cost_eu
                qla = value1 * cost_la
                qsrr = value1 * cost_sr
                qy = value1 * cost_y
                qzr = value1 * cost_zr
                qrb = value1 * cost_rb

            q[24] = qla
            q[25] = qbar
            q[26] = qeu2
            q[27] = qsrr
            q[28] = qy
            q[29] = qzr
            q[30] = qrb

            aa[1:4] = [1.4e-2, 1.0e-3, 1.0e-5]
            if 15.0 < H < 80.0 and zeta > 1.0e-30:
                k = 4
                for j in range(1, 4):
                    if s.MBa[j] <= H <= s.MBa[j + 1]:
                        k = j
                if k == 4:
                    if H >= 40.0:
                        k = kk = 4
                    else:
                        k = kk = 1
                        
                else:
                    kk = k + 1
                mm = np.array([s.MBa[k], s.MBa[kk]], dtype=float)

                if zeta < 1.0e-5:
                    for elem_idx, grid in ((25, s.WBa), (27, s.WSr), (28, s.WY), (24, s.WLa), (29, s.WZr), (30, s.WRb), (26, s.WEu)):
                        dd = np.array([grid[k, 3], grid[kk, 3]], dtype=float)
                        qbar_i, _ = self.polint(mm, dd, H)
                        q[elem_idx] += qbar_i
                elif 1.0e-5 <= zeta < 1.4e-2:
                    z = 1
                    for j in range(1, 3):
                        if aa[j + 1] < zeta <= aa[j]:
                            z = j
                    zz = z + 1
                    met = np.array([aa[zz], aa[z]], dtype=float)
                    for elem_idx, grid in ((25, s.WBa), (27, s.WSr), (28, s.WY), (24, s.WLa), (29, s.WZr), (30, s.WRb), (26, s.WEu)):
                        dd = np.array([grid[k, zz], grid[kk, zz]], dtype=float)
                        q1v, _ = self.polint(mm, dd, H)
                        dd = np.array([grid[k, z], grid[kk, z]], dtype=float)
                        q2v, _ = self.polint(mm, dd, H)
                        qbar_i, _ = self.polint(met, np.array([q1v, q2v], dtype=float), zeta)
                        q[elem_idx] += qbar_i
                elif zeta > 1.4e-2:
                    for elem_idx, grid in ((25, s.WBa), (27, s.WSr), (28, s.WY), (24, s.WLa), (29, s.WZr), (30, s.WRb), (26, s.WEu)):
                        dd = np.array([grid[k, 1], grid[kk, 1]], dtype=float)
                        qbar_i, _ = self.polint(mm, dd, H)
                        q[elem_idx] += qbar_i

            if 1.3 <= H <= 3.0:
                qbar, qsrr, qy, qeu2, qzr, qla, qrb = self.bario(zeta, H)

            if 1.0 <= H <= 6.0:
                qli = self.litio(zeta, H)
                if qli > 1.0e-20:
                    q[31] = qli

            if 1.3 <= H <= 3.0:
                q[24] = qla / 2.0
                q[25] = qbar / 2.0
                q[26] = qeu2 / 2.0
                q[27] = qsrr / 2.0
                q[28] = qy / 2.0
                q[29] = qzr / 2.0
                q[30] = qrb / 2.0

            if 12.0 <= H <= 50.0:
                q[9] = 0.07
            else:
                q[9] = 1.0e-20

            hecore = 0.0
            for i in range(1, 32):
                if H < 0.5:
                    q[i] = 0.0
                    if i == 14:
                        q[i] = H
                if binmax > 0.0:
                    q[i] = qia[i] * ratio + q[i]
                hecore += q[i]
        else:
            q = np.zeros(elem + 1, dtype=float)
            if binmax >= -8.0:
                q[31] = 2.0e-6 * 4.0
            else:
                value1 = 0.8e-6 * 20.0
                q[24] = value1 * 0.136
                q[25] = value1 * 1.0
                q[26] = value1 * (0.117 * 151.0 / 138.0)
                q[27] = value1 * (3.16 * 88.0 / 138.0)
                q[28] = value1 * (1.625 * 89.0 / 138.0 / 3.0)
                q[29] = value1 * (2.53 * 90.0 / 138.0)
                q[30] = value1 * (3.16 * 86.0 / 138.0)
            hecore = float(np.sum(q[1:32]))

        self.state.Q[:] = q
        return q[1:34].copy(), float(hecore)
