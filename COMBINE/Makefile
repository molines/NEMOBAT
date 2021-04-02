## -------------------------------------------------------------------------
##  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 225 $
##  $Id: Makefile 225 2009-04-29 13:52:46Z molines $
## -------------------------------------------------------------------------
#
#
CPP     = /lib/cpp
FC      = ifort
LD      = ifort
dirsub = ../SUBROUTINES
#
# Regles de construction
#
.SUFFIXES:	.F90
.DEFAULTS:
.F90.o:
	@echo default rule .F90.o used 
	$(FC) -c $(FFLAGS)  $*.F90 -o $*.o
#
NCDF_PATH       = /usr/local/
#
FFLAGS		=  -I$(NCDF_PATH)/include 
#
LDFLAGS	 	= -L$(NCDF_PATH)/lib -lnetcdf
# resultats --> les executables
#
combine:	combine.f90   $(dirsub)/rewrite_bathy_nc.o $(dirsub)/lbc.o
	$(LD) $(FFLAGS) -o $@ combine.f90 $(dirsub)/rewrite_bathy_nc.o $(dirsub)/lbc.o   $(LDFLAGS)

$(dirsub)/rewrite_bathy_nc.o: $(dirsub)/rewrite_bathy_nc.f90
	$(FC) -c $(FFLAGS) $(dirsub)/rewrite_bathy_nc.f90 -o $(dirsub)/rewrite_bathy_nc.o

$(dirsub)/lbc.o: $(dirsub)/lbc.f90
	$(FC) -c $(FFLAGS) $(dirsub)/lbc.f90 -o $(dirsub)/lbc.o

clean: 
	/bin/rm -f *.o *~

#
