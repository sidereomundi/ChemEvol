# ChemEvol

This repository contains a minimal Python reimplementation of the `MinGCE`
chemical evolution model alongside the original Fortran sources.
The Python code lives in the `pychem/` package and mirrors several of the
Fortran modules.

The repository ships with a small set of yield tables in the `DATI/` and
`YIELDSBA/` directories so that the example driver can run without external
resources.

## Requirements

* Python 3.10+
* `numpy`
* (optional) `gfortran` to build the original Fortran code

Install the Python dependencies using `pip`:

```bash
pip install -r requirements.txt
```

## Running the Python demo

1. Install the Python dependencies:

   ```bash
   pip install -r requirements.txt
   ```

2. Execute the driver which loads a few tables, performs an interpolation
   and prints some diagnostic output:

   ```bash
   python -m pychem.driver
   ```

   You can also invoke the package directly:

   ```bash
   python -m pychem
   ```

This command initialises a `GCEModel` instance, reads several files from
`DATI/`, evaluates the stellar lifetime function and runs the interpolation
routine. The printed values should match those obtained from the Fortran
version on the same input data.

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

## Repository layout

```
pychem/       # Python port of the main routines
src/          # Original Fortran code
DATI/         # Example yield tables used by the Python demo
YIELDSBA/     # Additional tables for the heavy elements
```

The Python code is intentionally lightweight and only implements a subset of
the full Fortran functionality, but the structure should make it clear how
additional routines can be translated.
