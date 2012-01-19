F90C=gfortran
F90FLAGS=-O3
LFITS=/usr/local/lib/libcfitsio.a
-include local.mk # optionally override compiler and flags with this file

EXECS= \
     cubet cubevr cubeemiss makemap makerotmap makerotmapsmall

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

.SUFFIXES: .f90
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<

cubet.o:  wfitsutils.o
cubevr.o:  wfitsutils.o
cubeemiss.o:  wfitsutils.o emissmod.o
emissmod.o: em2levmod.o
makemap.o: wfitsutils.o emissmod.o
makerotmap.o: cuberotate.o wfitsutils.o emissmod.o
makerotmapsmall.o: cuberotate.o wfitsutils.o emissmod.o

cubet: cubet.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubevr: cubevr.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubeemiss: cubeemiss.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makemap: makemap.o emissmod.o em2levmod.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makerotmapsmall: makerotmapsmall.o cuberotate.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makerotmap: makerotmap.o cuberotate.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)


