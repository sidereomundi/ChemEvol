"""Main MinGCE routine translated from ``src/main.f90``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys

import numpy as np

from .interpolation import Interpolator
from .io_routines import BASE_DIR, FortranState, IORoutines
from .tau import tau

try:
    from mpi4py import MPI  # type: ignore
except Exception:  # pragma: no cover - optional dependency
    MPI = None


@dataclass
class GCEModel:
    io: IORoutines = field(default_factory=IORoutines)
    state: FortranState = field(default_factory=FortranState)
    interpolator: Interpolator = field(init=False)

    UM1: float = 0.0
    UM2: float = 0.0
    UM3: float = 0.0
    UM4: float = 0.0
    A: float = 0.0
    B: float = 0.0
    C: float = 0.0
    D: float = 0.0
    M1: float = 0.0
    M2: float = 0.0
    M3: float = 0.0
    M4: float = 0.0

    def __post_init__(self) -> None:
        self.interpolator = Interpolator(state=self.state)

    def _initialize_from_fortran_tables(self, lowmassive: int = 1, mm: int = 0) -> None:
        self.io.load_main_tables(self.state, lowmassive=lowmassive, mm=mm)

    def _imf_scalo(self) -> None:
        self.UM1 = -1.35
        self.UM2 = -1.7
        mi = 0.1
        self.M1 = 2.0
        ms = 100.0
        self.M3 = 1.0
        zita = 0.3
        azita = ms ** (1.0 + self.UM2)
        bzita = self.M1 ** (1.0 + self.UM2)
        czita = self.M1 ** (1.0 + self.UM1)
        dzita = self.M3 ** (1.0 + self.UM1)
        self.B = zita / (
            self.M1 ** (self.UM2 - self.UM1) * (czita - dzita) / (1.0 + self.UM1)
            + (azita - bzita) / (1.0 + self.UM2)
        )
        self.A = self.B * self.M1 ** (self.UM2 - self.UM1)

    def _imf_salpeter(self) -> None:
        self.UM1 = -1.35
        mi = 0.1
        ms = 80.0
        azita = ms ** (1.0 + self.UM1)
        bzita = mi ** (1.0 + self.UM1)
        self.A = (1.0 + self.UM1) / (azita - bzita)

    def _multi_scalo(self, mmax: float, mmin: float) -> float:
        if mmax <= self.M1:
            return ((mmin**self.UM1 - mmax**self.UM1) / self.UM1) * self.A
        return ((mmin**self.UM2 - mmax**self.UM2) / self.UM2) * self.B

    def _multi_salpeter(self, mmax: float, mmin: float) -> float:
        return ((mmin**self.UM1 - mmax**self.UM1) / self.UM1) * self.A

    def _build_mass_bins(self, imf: int, tautype: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, int]:
        amu_vals = np.loadtxt(BASE_DIR / "DATI" / "amu.dat").flatten()
        amu = np.zeros(116, dtype=float)
        amu[1:116] = amu_vals[:115]

        mstars = np.zeros(2000, dtype=float)
        mstars1 = np.zeros(2000, dtype=float)
        multi1 = np.zeros(2000, dtype=float)
        binmax = np.zeros(2000, dtype=float)
        binmax1 = np.zeros(2000, dtype=float)
        multi2 = np.zeros(2000, dtype=float)

        norm = 0.0
        ss = 0
        binmass = 115
        ss2 = binmass - 1

        for jj in range(1, binmass):
            mstars1[jj] = amu[binmass - jj]
            mstars1[jj + 1] = amu[binmass - jj + 1]

            if imf == 1:
                multi1[jj] = self._multi_scalo(mstars1[jj], mstars1[jj + 1])
            elif imf == 3:
                multi1[jj] = self._multi_salpeter(mstars1[jj], mstars1[jj + 1])
            else:
                raise NotImplementedError("IMFKroupa setup is not present in src/main.f90")

            mstars[jj] = 0.5 * (mstars1[jj] + mstars1[jj + 1])
            binmax[jj] = 0.0

            if 3.0 <= mstars[jj] <= 16.0:
                mmu = mstars[jj]
                mumin = 0.8 / mmu
                mumin2 = 1.0 - 8.0 / mmu
                if mumin2 > mumin:
                    mumin = mumin2

                mux = 0.5 - mumin
                xmu = np.zeros(12, dtype=float)
                xmu[1] = mumin
                xmu[2] = mumin + 0.01 * mux
                xmu[3] = mumin + 0.02 * mux
                xmu[4] = mumin + 0.05 * mux
                xmu[5] = mumin + 0.1 * mux
                xmu[6] = mumin + 0.2 * mux
                xmu[7] = mumin + 0.3 * mux
                xmu[8] = mumin + 0.4 * mux
                xmu[9] = mumin + 0.6 * mux
                xmu[10] = mumin + 0.8 * mux
                xmu[11] = 0.5

                for j3 in range(1, 11):
                    ss += 1
                    ss2 = ss + binmass - 1
                    mu1 = xmu[j3]
                    mu2 = xmu[j3 + 1]
                    mstars[ss2] = mstars[jj] * (xmu[j3] + xmu[j3]) / 2.0
                    multi1[ss2] = (8.0 * (mu2**3 - mu1**3)) * multi1[jj] * 0.09
                    binmax[ss2] = mstars[jj] - mstars[ss2]
                    norm += multi1[ss2] * mstars[jj]

                multi1[jj] = 0.95 * multi1[jj]

            if 2.0 <= mstars[jj] <= 8.0:
                ss += 1
                ss2 = ss + binmass - 1
                mstars[ss2] = mstars[jj]
                binmax[ss2] = -1.0
                multi1[ss2] = 0.015 * multi1[jj]

            norm += multi1[jj] * mstars[jj]

        tdead_raw = np.zeros(ss2 + 1, dtype=float)
        for j in range(1, ss2 + 1):
            tdead_raw[j] = tau(max(mstars[j], 1.0e-8), tautype, binmax[j])

        order = np.argsort(tdead_raw[1 : ss2 + 1])
        for idx, o in enumerate(order, start=1):
            src = o + 1
            mstars1[idx] = mstars[src]
            binmax1[idx] = binmax[src]
            multi2[idx] = multi1[src]

        for jj in range(1, ss2 + 1):
            mstars[jj] = mstars1[jj]
            binmax[jj] = binmax1[jj]
            multi1[jj] = multi2[jj]

        tdead = np.zeros(2001, dtype=float)
        tdead[ss2 + 1] = 1.0e30
        for jj in range(1, ss2 + 1):
            tdead[jj] = tau(max(mstars[jj], 1.0e-8), tautype, binmax[jj])

        return mstars, binmax, multi1, tdead, norm, ss2

    def _mpi_ctx(self) -> tuple[object | None, int, int]:
        if MPI is None:
            return None, 0, 1
        comm = MPI.COMM_WORLD
        return comm, int(comm.Get_rank()), int(comm.Get_size())

    def MinGCE(
        self,
        endoftime: int,
        sigmat: float,
        sigmah: float,
        psfr: float,
        pwind: float,
        delay: int,
        time_wind: int,
        use_mpi: bool = True,
        show_progress: bool = True,
    ) -> None:
        nmax = 15000
        elem = 33

        self._initialize_from_fortran_tables(lowmassive=1, mm=0)

        imf = 1
        tautype = 1
        if imf == 1:
            self._imf_scalo()
        elif imf == 3:
            self._imf_salpeter()

        mstars, binmax, multi1, tdead, norm, ss2 = self._build_mass_bins(imf, tautype)

        allv = np.zeros(nmax + 2, dtype=float)
        gas = np.zeros(nmax + 2, dtype=float)
        stars = np.zeros(nmax + 2, dtype=float)
        remn = np.zeros(nmax + 2, dtype=float)
        hot = np.zeros(nmax + 2, dtype=float)
        wind = np.zeros(nmax + 2, dtype=float)
        oldstars = np.zeros(nmax + 2, dtype=float)
        zeta = np.zeros(nmax + 2, dtype=float)
        snianum = np.zeros(nmax + 2, dtype=float)
        spalla = np.zeros(nmax + 2, dtype=float)
        sfr_hist = np.zeros(nmax + 2, dtype=float)

        qqn = np.full((elem + 1, nmax + 2), 1.0e-20, dtype=float)
        ini = np.zeros(elem + 1, dtype=float)
        ini[31] = 1.0e-9

        winds = np.ones(32, dtype=float)
        winds[9] = 1.0

        superf = 20000.0
        threshold = 0.1
        sigmasun = 50.0
        kappa = 1.5
        rm = 8.0

        comm, rank, size = self._mpi_ctx()
        mpi_active = bool(use_mpi and comm is not None and size > 1)

        out_dir = BASE_DIR / "RISULTATI2"
        f_fis = None
        f_mod = None
        if rank == 0:
            out_dir.mkdir(parents=True, exist_ok=True)
            f_fis = open(out_dir / "fis.encesmin.dat", "w", encoding="ascii")
            f_mod = open(out_dir / "modencesmin.dat", "w", encoding="ascii")
            f_fis.write("time all gas stars remn hot zeta SFR nume SFR2 SIaN SIar\n")
            f_mod.write(
                "time all gas star SFR run Hen C12 O16 N14 C13 Ne Mg Si Fe S14 C13S S32 Ca Remn Zn K Sc Ti V Cr Mn Co Ni La Ba Eu Sr Y Zr Rb Li H He4\n"
            )

        elem_idx_no14 = np.array([i for i in range(1, elem) if i != 14], dtype=int)

        t = 0
        last_progress_step = 0
        progress_stride = max(1, endoftime // 200) if endoftime > 0 else 1

        def _print_progress(step: int) -> None:
            if not (show_progress and rank == 0 and endoftime > 0):
                return
            pct = 100.0 * step / endoftime
            bar_w = 30
            fill = int(bar_w * step / endoftime)
            bar = "#" * fill + "-" * (bar_w - fill)
            sys.stdout.write(f"\rProgress [{bar}] {pct:6.2f}% ({step}/{endoftime})")
            sys.stdout.flush()

        while True:
            t += 1

            t3 = t
            while True:
                if sigmat != 0.0:
                    allv[t3] = allv[t3 - 1] + sigmah * superf / (2.5 * sigmat) * np.exp(-((t3 - delay) ** 2) / (2.0 * sigmat**2))
                else:
                    allv[t3] = allv[t3 - 1]

                gas[t3] = allv[t3] - stars[t3] - remn[t3] - hot[t3] - wind[t3]
                if t3 >= endoftime:
                    break
                t3 += 1

            if gas[t] / superf > threshold and sigmah != 0.0:
                sfr = (
                    psfr
                    * (gas[t] / (superf * sigmah)) ** kappa
                    * (sigmah / sigmasun) ** (kappa - 1.0)
                    * (8.0 / rm)
                    * (superf / 1000.0)
                    * sigmah
                )
            else:
                sfr = 0.0
            sfr_hist[t] = sfr

            if gas[t] / superf >= threshold:
                hecores = np.zeros(ss2 + 2, dtype=float)
                mstars1_eff = np.zeros(ss2 + 2, dtype=float)
                qispecial = np.zeros((elem + 1, ss2 + 2), dtype=float)
                oldstars_contrib = 0.0

                if mpi_active:
                    local_bins = range(1 + rank, ss2 + 1, size)
                else:
                    local_bins = range(1, ss2 + 1)

                for jj in local_bins:
                    if tdead[jj] + t > 13500.0:
                        oldstars_contrib += multi1[jj] * sfr

                    q, hecore = self.interpolator.interp(mstars[jj], zeta[t], binmax[jj])
                    mstars1_eff[jj] = (binmax[jj] + mstars[jj]) if (binmax[jj] > 0.0) else mstars[jj]
                    hecores[jj] = hecore
                    qispecial[1:elem, jj] = q[: elem - 1]

                if mpi_active:
                    oldstars[t] += comm.allreduce(oldstars_contrib, op=MPI.SUM)
                    comm.Allreduce(MPI.IN_PLACE, hecores, op=MPI.SUM)
                    comm.Allreduce(MPI.IN_PLACE, mstars1_eff, op=MPI.SUM)
                    comm.Allreduce(MPI.IN_PLACE, qispecial, op=MPI.SUM)
                else:
                    oldstars[t] += oldstars_contrib

                t3 = t
                starstot = sfr * norm
                difftot = sfr * norm
                hecoretot = 0.0
                snian = 0.0
                jj = 1
                qacc = np.zeros(elem + 1, dtype=float)

                while True:
                    next_dead = tdead[jj + 1] if (jj + 1) <= ss2 else 1.0e30
                    died_now = t3 >= (t + tdead[jj])
                    next_died = t3 >= (t + next_dead)

                    if died_now:
                        dm = multi1[jj] * sfr
                        qacc[1:elem] += qispecial[1:elem, jj] * dm
                        hecoretot += hecores[jj] * dm
                        difftot -= (mstars1_eff[jj] - hecores[jj]) * dm
                        starstot -= mstars1_eff[jj] * dm
                        if binmax[jj] > 0.0:
                            snian += dm

                    if (not next_died) and gas[t] > 0.0:
                        qqn[elem_idx_no14, t3] = (
                            qqn[elem_idx_no14, t3] + qacc[elem_idx_no14] - qqn[elem_idx_no14, t] * difftot / gas[t]
                        )
                        qqn[31, t3] = (
                            qqn[31, t3] - qqn[31, t] * (mstars1_eff[jj] + hecores[jj]) * multi1[jj] * sfr / gas[t]
                        )

                    if died_now:
                        if next_died:
                            jj += 1
                            if jj > ss2:
                                break
                            continue
                        else:
                            stars[t3] += starstot
                            q14 = qacc[14]
                            remn[t3] += q14
                            if q14 < 0.0:
                                break
                            snianum[t3] += snian
                            if t3 >= endoftime:
                                break
                            t3 += 1
                            jj += 1
                            if jj > ss2:
                                break
                    else:
                        stars[t3] += starstot
                        q14 = qacc[14]
                        remn[t3] += q14
                        if q14 < 0.0:
                            break
                        snianum[t3] += snian
                        if t3 >= endoftime:
                            break
                        t3 += 1

            t3 = t
            if t > time_wind and gas[t] > 0.0:
                windist = pwind * sfr
            else:
                windist = 0.0

            while True:
                wind[t3] += windist

                for ii in range(1, 32):
                    if ii != 14 and sfr > 0.0 and gas[t] > 0.0:
                        qqn[ii, t3] = qqn[ii, t3] - qqn[ii, t] * winds[ii] * windist / gas[t]
                        if qqn[ii, t3] <= 1.0e-20:
                            qqn[ii, t3] = 1.0e-20

                zeta[t3] = 0.0
                for i in range(2, elem - 1):
                    if i != 14:
                        zeta[t3] += qqn[i, t3]

                if sfr > 0.0 and gas[t3] > 0.0:
                    zeta[t3] = zeta[t3] / gas[t3]
                else:
                    zeta[t3] = 1.0e-20

                qqn[elem, t3] = gas[t3] * 0.241 + qqn[1, t3]
                qqn[elem - 1, t3] = gas[t3] * (0.759 - zeta[t3]) - qqn[1, t3]

                if t > 1:
                    qqn[31, t3] = qqn[31, t3] + (allv[t] - allv[t - 1]) * ini[31]
                    denom = max(qqn[elem - 1, t3], 1.0e-30)
                    feh = np.log10(max(qqn[9, t3] / denom, 1.0e-30))
                    spalla[t3] = 10 ** (-9.50 + 1.24 * (feh - (-2.75)) + np.log10(denom))
                    qqn[31, t3] = qqn[31, t3] + spalla[t3] - spalla[t3 - 1]
                else:
                    qqn[31, t3] = qqn[31, t3] + allv[t] * ini[31]

                if t3 >= endoftime:
                    break
                t3 += 1

            if rank == 0:
                vals = [
                    float(t),
                    allv[t],
                    gas[t],
                    stars[t],
                    sfr,
                    oldstars[t],
                ] + [qqn[i, t] for i in range(1, 34)]
                f_mod.write(" ".join(f"{v: .5e}" for v in vals) + "\n")

            if t >= endoftime:
                break

            if t - last_progress_step >= progress_stride:
                _print_progress(t)
                last_progress_step = t

        if rank == 0:
            if show_progress and endoftime > 0:
                _print_progress(endoftime)
                sys.stdout.write("\n")
                sys.stdout.flush()

            for t in range(1, endoftime + 1):
                sfr2 = 20.0 * (stars[t] - stars[t - 1])
                siar = snianum[t] - snianum[t - 1]
                vals = [
                    float(t),
                    allv[t],
                    gas[t],
                    stars[t],
                    remn[t],
                    hot[t],
                    zeta[t],
                    sfr_hist[t],
                    1.0,
                    sfr2,
                    snianum[t],
                    siar,
                ]
                f_fis.write(" ".join(f"{v: .5e}" for v in vals) + "\n")

            f_mod.close()
            f_fis.close()

            print("MinGCE full translation run complete")
            if mpi_active:
                print("MPI ranks:", size)
            print("ninputyield:", self.state.ninputyield)
            print("final gas:", gas[endoftime] if endoftime > 0 else gas[0])
            print("final zeta:", zeta[endoftime] if endoftime > 0 else zeta[0])
