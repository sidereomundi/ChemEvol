FC = gfortran
FFLAGS = -O2 -ffixed-form -I.

# Source files ordered to satisfy module dependencies
SRC = \
src/io.f90 \
src/interpolation.f90 \
src/tau.f90 \
src/main.f90 \
src/driver.f90

# Objects compiled in the same order as SRC
OBJ = $(SRC:.f90=.o)

all: GCE_min.x

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Ensure modules are built in correct order
src/main.o: src/io.o src/interpolation.o
src/driver.o: src/main.o src/interpolation.o src/io.o

GCE_min.x: $(OBJ) Makefile
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm -f $(OBJ) GCE_min.x
