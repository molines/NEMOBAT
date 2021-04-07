## -------------------------------------------------------------------------
##  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 228 $
##  $Id: Makefile 228 2009-04-29 14:19:25Z forge $
## -------------------------------------------------------------------------
include make.macro

EXEC=grd2nc.x mk_coor.x smith_sandwell2nc.x revertetopo2.x find_zoom.x NSIDC_map_trf.x \
    test_xy.x mergebat.x apply_history.x apply_patch.x

all: $(EXEC)

grd2nc.x: grd2nc.f90
	$(FC)  grd2nc.f90 -o grd2nc.x $(FFLAGS)

mk_coor.x: mk_coor.f90
	$(FC)  mk_coor.f90 -o mk_coor.x $(FFLAGS)

revertetopo2.x: revertetopo2.f90
	$(FC)  revertetopo2.f90 -o revertetopo2.x $(FFLAGS)

find_zoom.x: find_zoom.f90
	$(FC)  find_zoom.f90 -o find_zoom.x $(FFLAGS)


smith_sandwell2nc.x: smith_sandwell2nc.f90
	$(FC)  smith_sandwell2nc.f90 -o smith_sandwell2nc.x $(FFLAGS) -assume byterecl -convert big_endian

NSIDC_map_trf.x: NSIDC_map_trf.f90
	$(FC) NSIDC_map_trf.f90 -o NSIDC_map_trf.x $(FFLAGS)

mergebat.x: mergebat.f90
	$(FC) mergebat.f90 -o mergebat.x $(FFLAGS)

apply_history.x: apply_history.f90
	$(FC) apply_history.f90 -o apply_history.x $(FFLAGS)

apply_patch.x: apply_patch.f90
	$(FC) apply_patch.f90 -o apply_patch.x $(FFLAGS)

test_xy.x: test_xy.f90
	$(FC) test_xy.f90 -o test_xy.x $(FFLAGS)

clean:
	\rm -f *.o *~ *.mod 

cleanall : clean
	 \rm -f $(EXEC)
