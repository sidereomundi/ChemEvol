# Running ChemEvol (Python)

## Requirements

- Python 3.10+
- `numpy`
- Optional for MPI runs: `mpi4py` and an MPI launcher (`mpiexec`/`mpirun`)

Install:

```bash
pip install -r requirements.txt
pip install mpi4py
```

## Quick Start

Run the default driver:

```bash
python -m pychem.driver
```

This uses the driver defaults and writes results under `RISULTATI2/`.

## Custom Serial Run

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

## MPI Parallel Run

`MinGCE(..., use_mpi=True)` parallelizes the expensive per-mass-bin
interpolation stage across MPI ranks.

```bash
mpiexec -n 8 python -c "from pychem.main import GCEModel; m=GCEModel(); m.MinGCE(2000,2000.0,54.0,0.2,0.0,7000,1000000,use_mpi=True,show_progress=True)"
```

Notes:

- Rank 0 writes output files and prints progress.
- If `mpi4py` is unavailable, the code falls back to serial mode.

## Progress Bar

Use `show_progress=True` (default) to print progress in the main time loop.
Set `show_progress=False` to disable.

## Outputs

Runs produce:

- `RISULTATI2/modencesmin.dat`: abundances/evolution table
- `RISULTATI2/fis.encesmin.dat`: physical summary table

## Main Parameters

- `endoftime`: total time steps (Myr in original model logic)
- `sigmat`: width parameter for infall law
- `sigmah`: normalization for infall/SFR terms
- `psfr`: star formation efficiency factor
- `pwind`: wind strength multiplier
- `delay`: delay term used in infall profile
- `time_wind`: time threshold after which wind can start
- `use_mpi`: enable MPI work sharing
- `show_progress`: enable console progress bar
