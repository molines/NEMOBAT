## -------------------------------------------------------------------------
##  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 228 $
##  $Id: Makefile 228 2009-04-29 14:19:25Z forge $
## -------------------------------------------------------------------------

FC=ifort
NCDF=/usr/local

FFLAGS=-O2 -I$(NCDF)/include -L$(NCDF)/lib -lnetcdf

EXEC=grd2nc mk_coor smith_sandwell2nc revertetopo2 find_zoom

all: $(EXEC)

grd2nc: grd2nc.f90
	$(FC)  grd2nc.f90 -o grd2nc $(FFLAGS)

mk_coor: mk_coor.f90
	$(FC)  mk_coor.f90 -o mk_coor $(FFLAGS)

revertetopo2: revertetopo2.f90
	$(FC)  revertetopo2.f90 -o revertetopo2 $(FFLAGS)

find_zoom: find_zoom.f90
	$(FC)  find_zoom.f90 -o find_zoom $(FFLAGS)


smith_sandwell2nc: smith_sandwell2nc.f90
	$(FC)  smith_sandwell2nc.f90 -o smith_sandwell2nc $(FFLAGS) -assume byterecl -convert big_endian

clean:
	\rm -f *.o *~ *.mod

cleanall : clean
	 \rm -f $(EXEC)
