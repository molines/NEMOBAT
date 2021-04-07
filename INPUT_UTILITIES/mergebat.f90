PROGRAM mergebat
  !!======================================================================
  !!                     ***  PROGRAM  mergebat  ***
  !!=====================================================================
  !!  ** Purpose : Produce a merged bathymetry between 2 files, using
  !!               a distance as weight
  !!
  !!  ** Method  :
  !!
  !! History :  4.0  : 03/2017  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
   USE netcdf
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE 
  INTEGER(KIND=4) :: ji,jj
  INTEGER(KIND=4) :: npiglo, npjglo

  INTEGER(KIND=4) :: ierr, ncid, id

  REAL(KIND=4)                              :: rlimit, alfa
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  Vorig, Vpatch, Vfinal, tdist

  CHARACTER(LEN=80)  :: cf_orig="eGREENLAND025.L75_bathymetry.nc"
  CHARACTER(LEN=80)  :: cf_patch="eGREENLAND025.L75_bathy_meter_002.nc2"
  CHARACTER(LEN=80)  :: cf_final="eGREENLAND025.L75_bathymetry_patched.nc"
  CHARACTER(LEN=80)  :: cf_dist="dist.coast"
  CHARACTER(LEN=80)  :: cmd

  rlimit=80 !km
  cmd="cp "//TRIM(cf_orig)//" "//TRIM(cf_final)
  CALL SYSTEM(cmd)

  ierr=NF90_OPEN(cf_orig,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_DIMID(ncid,"x",id ) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr= NF90_INQ_DIMID(ncid,"y",id ) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)
  PRINT *, ' NPIGLO : ', npiglo
  PRINT *, ' NPJGLO : ', npjglo
  ALLOCATE( Vorig(npiglo,npjglo), Vpatch(npiglo,npjglo), Vfinal(npiglo,npjglo), tdist(npiglo,npjglo))

  ierr= NF90_INQ_VARID(ncid,"Bathymetry",id) ; ierr = NF90_GET_VAR(ncid,id,Vorig)
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cf_patch,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,"Bathymetry",id) ; ierr = NF90_GET_VAR(ncid,id,Vpatch)
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cf_dist,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,"Tcoast",id) ; ierr = NF90_GET_VAR(ncid,id,tdist)
  tdist=tdist/1000.  ! from meters to km
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cf_final,NF90_WRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,"Bathymetry",id) 
  Vfinal=Vorig
  DO jj=2,npjglo
     DO ji = 2, npiglo
     IF ( Vpatch(ji,jj) > 0 ) THEN
     IF ( tdist(ji,jj) <= rlimit)  THEN
       alfa =  tdist(ji,jj)/rlimit
     ELSE
       alfa = 1.
     ENDIF
     Vfinal(ji,jj) = Vpatch(ji,jj) * alfa + Vorig(ji,jj) *(1 - alfa)
     ELSE IF (Vorig(ji,jj) > 0 ) THEN
     Vfinal(ji,jj) = -50.
  
     ENDIF
     ENDDO
  ENDDO
 ierr = NF90_PUT_VAR(ncid,id,Vfinal)

  ierr = NF90_CLOSE(ncid)

END PROGRAM mergebat
