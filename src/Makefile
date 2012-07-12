F90C=gfortran
F90FLAGS=-O3
LFITS=/usr/local/lib/libcfitsio.a
-include local.mk # optionally override compiler and flags with this file

EXECS= \
     cubeav cubet cubevr cubewind cubeemiss makemap makerotmap makerotmapsmall makevcubes

all: $(EXECS)

clean: 
	rm -fr $(EXECS) *.o *.mod

.SUFFIXES: .f90
.f90.o:  
	$(F90C) $(F90FLAGS) -c $<

cubeav.o:  wfitsutils.o findcolumn_standalone.o
cubet.o:  wfitsutils.o
cubevr.o:  wfitsutils.o
cubewind.o:  wfitsutils.o
cubeemiss.o:  wfitsutils.o emissmod.o
emissmod.o: em2levmod.o
makemap.o: wfitsutils.o emissmod.o
makerotmap.o: cuberotate.o wfitsutils.o emissmod.o
makerotmapsmall.o: cuberotate.o wfitsutils.o emissmod.o
makevcubes.o: wfitsutils.o psfutils.o fft.o
psfutils.o: fft.o


makevcubes: makevcubes.o wfitsutils.o psfutils.o fft.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubeav: cubeav.o wfitsutils.o findcolumn_standalone.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubet: cubet.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubevr: cubevr.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubewind: cubewind.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

cubeemiss: cubeemiss.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makemap: makemap.o emissmod.o em2levmod.o wfitsutils.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makerotmapsmall: makerotmapsmall.o cuberotate.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)

makerotmap: makerotmap.o cuberotate.o wfitsutils.o emissmod.o em2levmod.o
	$(F90C) $(F90FLAGS) -o $@ $^ $(LFITS)


