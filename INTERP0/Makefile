#
## $Date: 2009-04-29 15:53:40 +0200 (mer. 29 avril 2009) $
## $Rev: 226 $
## $Id: Makefile 226 2009-04-29 13:53:40Z molines $

# change compiler from pgf90 to ifort (compile on kaaba)
#
#****************************************************************
#
#  cpp, compilation, chargement
#
include make.macro

dirsub = ../SUBROUTINES
#
# Regles de construction
#
.SUFFIXES:	.F90
.DEFAULTS:
.F90.o:
	@echo default rule .F90.o used 
	$(CPP)   $*.F90 -o $*.f90
	$(FC) -c $(FFLAGS)  $*.f90 -o $*.o
#
# resultats --> les executables
#
#  etopo2 option: needs to define "etopo2" 

all : batinterp.exe #batinterp1.exe batinterp2.exe
batinterp2.exe:	batinterpnew.F90  $(dirsub)/scilib.o  $(dirsub)/lbc.o
	$(LD) -Detopo2 $(FFLAGS) -o $@ batinterpnew.F90 $(dirsub)/scilib.o $(dirsub)/lbc.o $(LDFLAGS)
#
#   gebco1 option (default)
batinterp1.exe:	batinterpnew.F90  $(dirsub)/scilib.o $(dirsub)/lbc.o
	$(LD) $(FFLAGS) -o $@ batinterpnew.F90 $(dirsub)/scilib.o $(dirsub)/lbc.o $(LDFLAGS)

#
batinterp.exe: batinterp.F90 $(dirsub)/scilib.o $(dirsub)/lbc.o
	$(LD) $(FFLAGS) -o $@ batinterp.F90 $(dirsub)/scilib.o $(dirsub)/lbc.o $(LDFLAGS)

batinterp_grd.exe:	batinterpnew_grd.F90  $(dirsub)/scilib.o $(dirsub)/lbc.o
	$(LD) $(FFLAGS) -o $@ batinterpnew_grd.F90 $(dirsub)/scilib.o $(dirsub)/lbc.o $(LDFLAGS)

$(dirsub)/scilib.o: $(dirsub)/scilib.F
	$(FC) -c $(FFLAGS) $(dirsub)/scilib.F -o $(dirsub)/scilib.o

$(dirsub)/lbc.o: $(dirsub)/lbc.f90
	$(FC) -c $(FFLAGS) $(dirsub)/lbc.f90 -o $(dirsub)/lbc.o

clean: 
	/bin/rm -f *.o *~ batinterp.exe batinterp1.exe batinterp2.exe
#
