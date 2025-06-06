F77 = gfortran

FLAGS  =  -O2 


.f.o:
	$(F77) $(FLAGS) -c $<



GCE_min: GCE_min.o Makefile
	$(F77) $(FLAGS) -o GCE_min.x  GCE_min.o

