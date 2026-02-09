# Performance Guide

This document covers practical performance tuning for the Python MinGCE run.

## Where Time Is Spent

The dominant cost is the per-time-step interpolation/ejecta stage over stellar
mass bins (`jj` loop). This is the part parallelized with MPI.

## Serial vs MPI

- Use serial (`use_mpi=False`) for quick tests and very small `endoftime`.
- Use MPI (`use_mpi=True`) for production runs with larger `endoftime`.
- Rank 0 handles file output and progress display; all ranks contribute to the
  heavy interpolation work.

## Recommended MPI Rank Counts

Start with:

- `endoftime <= 500`: 1-2 ranks
- `500 < endoftime <= 5000`: 4-8 ranks
- `endoftime > 5000`: 8-16 ranks

Then benchmark on your machine. Past a certain rank count, communication
(`Allreduce`) and filesystem overhead dominate.

## Benchmarking

Use wall-clock timing for realistic runs:

```bash
time python -c "from pychem.main import GCEModel; m=GCEModel(); m.MinGCE(2000,2000.0,54.0,0.2,0.0,7000,1000000,use_mpi=False,show_progress=False)"
```

MPI benchmark:

```bash
time mpiexec -n 8 python -c "from pychem.main import GCEModel; m=GCEModel(); m.MinGCE(2000,2000.0,54.0,0.2,0.0,7000,1000000,use_mpi=True,show_progress=False)"
```

Run each command 2-3 times and compare median runtime.

## Practical Tuning Tips

- Disable progress while benchmarking:
  - `show_progress=False`
- Keep output writing on default rank 0 only (already implemented).
- Use larger `endoftime` for stable speedup measurements.
- Avoid launching more MPI ranks than physical CPU cores.
- If your cluster/node is NUMA-heavy, test rank pinning options from your MPI
  launcher for better cache locality.

## Expected Scaling Behavior

- Good speedup at low-to-moderate ranks.
- Diminishing returns as rank count increases due to:
  - collective communication (`Allreduce`)
  - synchronized per-step progress in global simulation state
  - output and memory bandwidth limits

## Correctness While Optimizing

After any performance change, validate:

```bash
pytest -q
```

And run one short serial + MPI smoke test to ensure outputs are generated:

```bash
python -m pychem.driver
mpiexec -n 2 python -c "from pychem.main import GCEModel; m=GCEModel(); m.MinGCE(10,2000.0,54.0,0.2,0.0,7000,1000000,use_mpi=True,show_progress=False)"
```
