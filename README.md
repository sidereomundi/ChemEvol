# ChemEvol

This repository contains a Python translation of the `MinGCE` chemical
evolution model alongside the original Fortran sources.
The Python implementation lives in `pychem/` and includes:

- Fortran-style table loading and shared state
- translated interpolation/yield routines
- full MinGCE time-evolution loop
- optional MPI parallelism for the heavy per-bin interpolation stage
- optional console progress bar

The repository ships with a small set of yield tables in the `DATI/` and
`YIELDSBA/` directories so that the example driver can run without external
resources.

## Requirements

* Python 3.10+
* `numpy`
* Optional: `mpi4py` for MPI parallel runs
* `gfortran` to build the original Fortran code
* `meson` (required by `f2py` on Python 3.12+)

Install the Python dependencies using `pip`:

```bash
pip install -r requirements.txt
```

For MPI support:

```bash
pip install mpi4py
```

## Running

1. Basic run (driver defaults):

   ```bash
   python -m pychem.driver
   ```

2. Custom run from Python:

   ```bash
   python - <<'PY'
   from pychem.main import GCEModel
   m = GCEModel()
   m.MinGCE(
       endoftime=2000,
       sigmat=2000.0,
       sigmah=54.0,
       psfr=0.2,
       pwind=0.0,
       delay=7000,
       time_wind=1000000,
       use_mpi=False,
       show_progress=True,
   )
   PY
   ```

3. MPI run:

   ```bash
   mpiexec -n 8 python -c "from pychem.main import GCEModel; m=GCEModel(); m.MinGCE(2000,2000.0,54.0,0.2,0.0,7000,1000000,use_mpi=True,show_progress=True)"
   ```

Output files are written to `RISULTATI2/modencesmin.dat` and
`RISULTATI2/fis.encesmin.dat`.

See `docs/RUNNING.md` for a detailed run guide.
See `docs/PERFORMANCE.md` for performance tuning and MPI scaling guidance.

## Running the tests

The test suite uses `pytest`. Install the test dependencies and execute
`pytest` from the repository root:

```bash
pip install numpy pytest
pytest
```

## Fortran version

The original Fortran implementation is still available under the `src/`
folder. To compile it you need a Fortran compiler such as `gfortran`.
The provided `Makefile` builds the executable `GCE_min.x`:

```bash
make
./GCE_min.x
```

To compile the same sources into a Python module using `f2py` you can run:

```bash
make python
```
This produces `gce.so` which exposes the Fortran routines to Python.
Building with Python 3.12 requires the external `meson` build
system in addition to `gfortran`. Install the prerequisites with

```bash
sudo apt-get install gfortran
pip install meson
```

## Repository layout

```
pychem/       # Python port of the main routines
src/          # Original Fortran code
DATI/         # Example yield tables used by the Python demo
YIELDSBA/     # Additional tables for the heavy elements
docs/         # Usage documentation
```
