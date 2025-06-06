FC = gfortran
FFLAGS = -O2 -ffixed-form -I.

SRC = src/main.f90 src/interpolation.f90 src/io.f90 src/driver.f90 src/tau.f90
OBJ = $(SRC:.f90=.o)

all: GCE_min.x

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

GCE_min.x: $(OBJ) Makefile
	$(FC) $(FFLAGS) -o $@ $(OBJ)

clean:
	rm -f $(OBJ) GCE_min.x
