PROGRAM findzoom
  !----------------------------------------------------------------------------
  !           ***  PROGRAM findzoom  ***
  !
  !  Purpose : Give a geographical zone and a file name on input
  !            Output the zone in grid point indexes.
  !
  !  Method : Assume longitude between -180 +180; assume
  !           lon lat dimensions and variables.
  !           Work with regular grid.
  !
  ! history : Original J.M. Molines 26/04/2009
  !----------------------------------------------------------------------------
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
  !!  $Rev: 228 $
  !!  $Id: find_zoom.f90 228 2009-04-29 14:19:25Z forge $
  !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER :: ncid,idvar,istatus
  INTEGER :: narg, iargc
  INTEGER :: ji,jj
  INTEGER :: imin, imax, jmin,jmax, ni,nj
  REAL(KIND=4) :: xmin, xmax, ymin, ymax
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rlon, rlat
  CHARACTER(LEN=80) :: cfilin, cdum

  narg=iargc()
  IF ( narg /= 5 ) THEN
   PRINT *, 'USAGE: findzoom  nc_file xmin xmax ymin ymax'
   PRINT *, '      nc_file : name of a regular nc file with lon/lat'
   PRINT *, '      xmin xmax : limits in longitude '
   PRINT *, '      ymin ymax : limits in latitude '
  END IF
  CALL getarg(1, cfilin)
  CALL getarg(2, cdum) ; READ(cdum,*) xmin
  CALL getarg(3, cdum) ; READ(cdum,*) xmax
  CALL getarg(4, cdum) ; READ(cdum,*) ymin
  CALL getarg(5, cdum) ; READ(cdum,*) ymax

  ! open file and get longitude and latitude
  istatus = NF90_OPEN(cfilin,NF90_NOWRITE,ncid)
  istatus = NF90_INQ_DIMID(ncid,'lon',idvar)
  istatus = NF90_INQUIRE_DIMENSION(ncid,idvar,len=ni)
  istatus = NF90_INQ_DIMID(ncid,'lat',idvar)
  istatus = NF90_INQUIRE_DIMENSION(ncid,idvar,len=nj)

  ALLOCATE ( rlon(ni), rlat(nj) )

  istatus = NF90_INQ_VARID(ncid,'lon',idvar)
  istatus = NF90_GET_VAR(ncid,idvar,rlon)
  istatus = NF90_INQ_VARID(ncid,'lat',idvar)
  istatus = NF90_GET_VAR(ncid,idvar,rlat)
  istatus = NF90_CLOSE(ncid)

  ! search for longitude range
  ji=1
  DO WHILE ( rlon(ji) < xmin )
   ji=ji+1
  ENDDO
  imin=ji-1
   
  ji=imin
  DO WHILE ( rlon(ji) < xmax )
   ji=ji+1
  ENDDO
  imax=ji-1

  ! search for latitude range
  jj=1
  DO WHILE ( rlat(jj) < ymin )
   jj=jj+1
  ENDDO
  jmin=jj-1
   
  jj=jmin
  DO WHILE ( rlat(jj) < ymax )
   jj=jj+1
  ENDDO
  jmax=jj-1
  PRINT *, imin, imax, jmin, jmax
  
END PROGRAM findzoom
