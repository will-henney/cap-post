F90C=gfortran
F90FLAGS=-O3
LFITS=/usr/local/lib/libcfitsio.a
-include local.mk # optionally override compiler and flags with this file

EXECS= \
     cubet cubevr

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

.SUFFIXES: .f90
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<

cubet.o:  wfitsutils.o
cubevr.o:  wfitsutils.o

cubet: cubet.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubevr: cubevr.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)
