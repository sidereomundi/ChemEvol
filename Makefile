FC = gfortran
FFLAGS = -O2 -ffixed-form -I.

# Run ``make`` to build the standalone executable ``GCE_min.x``.
# ``make python`` compiles the Fortran sources into ``gce.so`` using f2py.

# Fortran source files used by both the executable and the Python module
# ``SRC`` intentionally excludes ``src/driver.f90`` which is only required
# for the standalone executable.  The same list is reused by the ``python``
# target defined below.
SRC = \
src/io.f90 \
src/interpolation.f90 \
src/tau.f90 \
src/main.f90

# Additional source only for the standalone program
DRIVER = src/driver.f90

# Objects compiled in the same order as the source files
OBJ = $(SRC:.f90=.o) $(DRIVER:.f90=.o)

all: GCE_min.x

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Ensure modules are built in correct order
src/main.o: src/io.o src/interpolation.o
src/driver.o: src/main.o src/interpolation.o src/io.o

GCE_min.x: $(OBJ) Makefile
	$(FC) $(FFLAGS) -o $@ $(OBJ)

PY_SRC = src/interpolation.f90 src/io.f90 src/tau.f90
PY_MOD = gce

# Build the Fortran sources into a Python extension using ``f2py``. This
# produces ``gce.so`` which exposes the legacy routines to Python.

$(PY_MOD).so: $(PY_SRC)
	f2py -c --fcompiler=gfortran --f90flags="-ffixed-form -ffixed-line-length-none -fdollar-ok -I$(CURDIR)" -m $(PY_MOD) $(PY_SRC)

# Alias so ``make python`` will produce the module
python: $(PY_MOD).so

clean:
	rm -f $(OBJ) GCE_min.x $(PY_MOD).so
