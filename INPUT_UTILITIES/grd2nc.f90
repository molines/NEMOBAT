PROGRAM grd2nc
  !----------------------------------------------------------------------------
  !                     ***  PROGRAM grd2nc  ***
  !
  !  Purpose : convert grd type file to standard ncdf map with lon,lat
  !
  !  Method : use a namelist for defining the name of the dimension an variables
  !           on both input (grd) and output files (nc)
  !
  ! history : Original Jean-Marc Molines, April 2009
  !----------------------------------------------------------------------------
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
  !!  $Rev: 228 $
  !!  $Id: grd2nc.f90 228 2009-04-29 14:19:25Z forge $
  !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: idepth
  INTEGER        , DIMENSION(:)  , ALLOCATABLE :: ndim
  INTEGER                :: istatus, ncid,ivarid !: net cdf stuff
  INTEGER                :: id_lon, id_lat, id_vlon, id_vlat   !: net cdf stuff
  INTEGER                :: id_data              !: net cdf stuff
  INTEGER                :: ji,jj , ii           !: dummy loop index
  INTEGER                :: i1, i2, j1, j2
  INTEGER                :: numnam=10            !: namlist logical unit
  INTEGER                :: narg, jarg, iargc
  INTEGER                :: iside, ixy 
  INTEGER                :: niout, njout, ni, nj

  REAL(KIND=8)                        :: dx,dy, x1,x2, y1, y2
  REAL(KIND=4), DIMENSION(:),ALLOCATABLE       :: z, rlon, rlat
  REAL(KIND=8), DIMENSION(:),ALLOCATABLE       :: x_range, y_range, spacing

  CHARACTER(LEN=80)      :: cdum                  !: for reading line arguments
  CHARACTER(LEN=80)      :: cnamlist = 'namelist' !: name of namelist file (ARG)

  ! Input Grid                       default name
  CHARACTER(LEN=80)      :: cnamegrd = 'dummy.grd' !: name of grd_file
  CHARACTER(LEN=80)      :: cside    = 'side'     !: name side dimension
  CHARACTER(LEN=80)      :: cxysize  = 'xysize'   !: name xysize dimension
  !
  CHARACTER(LEN=80)      :: cx_range = 'x_range'  !: name of x_range variable
  CHARACTER(LEN=80)      :: cy_range = 'y_range'  !: name of y_range variable
  CHARACTER(LEN=80)      :: cz_range = 'z_range'  !: name of z_range variable
  CHARACTER(LEN=80)      :: cspacing = 'spacing'  !: name of spacing variable
  CHARACTER(LEN=80)      :: cdim     = 'dimension' !: name of dimension variable
  CHARACTER(LEN=80)      :: cdatain  = 'z'        !: name of data  variable
  !
  ! Output grid                       default name
  CHARACTER(LEN=80)      :: cnamenc  = 'dummy.nc' !: name of nc file
  CHARACTER(LEN=80)      :: clon     = 'lon'      !: name of longitude dim/var
  CHARACTER(LEN=80)      :: clat     = 'lat'      !: name of latitude dim/var
  CHARACTER(LEN=80)      :: cdataout = 'Bathy'    !: name of data  variable


  NAMELIST/input_grd/ cnamegrd, cside, cxysize, cx_range, cy_range, cz_range, &
     &                cspacing, cdim, cdatain

  NAMELIST/output_grd/cnamenc,  clon, clat, cdataout

  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *, ' Usage : grd2nc namelist_name [-p] [-zoom x1 x2 y1 y2 ] '
     PRINT *, ' All required informations are given in the namelist file'
     PRINT *, ' specified in the argument of the command '
     PRINT *, ' Using grd2nc -p gives you a template namefile that you '
     PRINT *, ' use as a starting point : namelist_template'
     STOP
  ENDIF

  CALL getarg(1,cdum)
  IF ( cdum == '-p' ) THEN
     CALL print_template 
     STOP
  ENDIF
  ! open and read namelist
  cnamlist=cdum
  OPEN(numnam, FILE=cnamlist) ; REWIND (numnam)
  READ(numnam,input_grd)
  READ(numnam,output_grd)
  CLOSE(numnam)

  istatus=NF90_OPEN(cnamegrd,NF90_NOWRITE,ncid)
  ! dimension and allocation
  istatus=NF90_INQ_DIMID (ncid,cxysize,ivarid)
  istatus=NF90_INQUIRE_DIMENSION (ncid,ivarid,LEN=ixy)
  ALLOCATE (z(ixy))

  istatus=NF90_INQ_DIMID (ncid,cside,ivarid)
  istatus=NF90_INQUIRE_DIMENSION (ncid,ivarid,LEN=iside)
  ALLOCATE (ndim(iside), x_range(iside), y_range(iside), spacing(iside) )
  !
  istatus=NF90_INQ_VARID (ncid,cspacing,ivarid)
  istatus=NF90_GET_VAR (ncid,ivarid,spacing)
  dx=spacing(1) ; dy=spacing(2)
  !
  istatus=NF90_INQ_VARID (ncid,cx_range,ivarid)
  istatus=NF90_GET_VAR (ncid,ivarid,x_range)
  !
  istatus=NF90_INQ_VARID (ncid,cy_range,ivarid)
  istatus=NF90_GET_VAR (ncid,ivarid,y_range)
  !
  istatus=NF90_INQ_VARID (ncid,cdim,ivarid)
  istatus=NF90_GET_VAR (ncid,ivarid,ndim)
  !
  ni= ndim(1)
  nj= ndim(2)
  IF ( ni * nj /= ixy ) THEN
    PRINT *,' This grd file is weird: number of grid point incoherent with size of domain'
    STOP
  ENDIF

  x1=x_range(1) ; x2=x_range(2)
  y1=y_range(1) ; y2=y_range(2)

  ! Look for extra arguments on the command line
  IF ( narg > 1 ) THEN
   jarg=2
   DO WHILE ( jarg <= narg )
     CALL getarg( jarg, cdum )
     IF ( cdum == '-zoom' ) THEN ! get x1 x2 y1 y2
        jarg=jarg+1
        CALL getarg( jarg, cdum) ; READ(cdum,*) x1 ; jarg=jarg+1
        CALL getarg( jarg, cdum) ; READ(cdum,*) x2 ; jarg=jarg+1
        CALL getarg( jarg, cdum) ; READ(cdum,*) y1 ; jarg=jarg+1
        CALL getarg( jarg, cdum) ; READ(cdum,*) y2 ; jarg=jarg+1
     ENDIF
   END DO
  END IF

  ALLOCATE(idepth(ni,nj), rlon(ni), rlat(nj) )
  print *, ni,nj, x_range, y_range
  DO ji=1, ni
    rlon(ji) = x_range(1) + (ji-1)*dx
  END DO
  DO jj=1, nj
    rlat(jj) = y_range(1) + (jj-1)*dy
  END DO
  ! verification :
  PRINT *,' Should be 0 for x ', x_range(2) - rlon(ni)
  PRINT *,' Should be 0 for y ', y_range(2) - rlat(nj)
  !
  i1 =     (x1 - x_range(1))/dx +1
  i2 =ni + (x2 - x_range(2))/dx 
  j1 =     -(y2 - y_range(2))/dy +1
  j2 =nj - (y1 - y_range(1))/dy 
! j2 =     (y2 - y_range(2))/dy +1
! j1 =     (y1 - y_range(1))/dy 
  niout=i2-i1 +1
  njout=j2-j1 +1
  !
  istatus=NF90_INQ_VARID(ncid,cdatain,ivarid)
  istatus=NF90_GET_VAR(ncid,ivarid,z)
  istatus=NF90_CLOSE( ncid )
  !
  ii=1
  DO jj=1,nj
     DO ji=1,ni
        idepth(ji,jj) = z(ii)
        ii = ii+ 1
     END DO
  END DO
  ! output results on netcdf file
  istatus=NF90_CREATE(cnamenc, NF90_CLOBBER, ncid )
  ! define dimension
  istatus=NF90_DEF_DIM(ncid,clon,niout, id_lon)
  istatus=NF90_DEF_DIM(ncid,clat,njout, id_lat)
  ! define variables
  istatus=NF90_DEF_VAR(ncid,clon,NF90_FLOAT,(/id_lon/), id_vlon )
  istatus=NF90_DEF_VAR(ncid,clat,NF90_FLOAT,(/id_lat/), id_vlat )
  istatus=NF90_DEF_VAR(ncid,cdataout,NF90_SHORT,(/id_lon,id_lat/), id_data )
  ! eventually attributes there ...
  
  istatus=NF90_ENDDEF(ncid)
  ! put variable
  istatus=NF90_PUT_VAR(ncid, id_vlon, rlon(i1:i2) )
  istatus=NF90_PUT_VAR(ncid, id_vlat, rlat(j1:j2) )
  istatus=NF90_PUT_VAR(ncid, id_data, idepth(i1:i2,j2:j1:-1) )
  istatus=NF90_CLOSE (ncid)
CONTAINS
  SUBROUTINE print_template
    ! print template namelist
    OPEN(numnam,FILE='namelist_template')
    WRITE(numnam,'("!! Template namelist for grd2nc ")' )
    WRITE(numnam,'("!!edit to fit your own requirement ")' )
    WRITE(numnam,'("!!template produced by grd2nc -p ")' )
    WRITE(numnam,'("!")' )
    WRITE(numnam,'("! input grid grd file")')
    WRITE(numnam,'("&input_grd")')
    WRITE(numnam,'("  cnamegrd  =  ",a10,"  ! name of grd input file")') TRIM(cnamegrd)
    WRITE(numnam,'("  cside     =  ",a10,"  ! name of side dimension")') TRIM(cside)
    WRITE(numnam,'("  cxysize   =  ",a10,"  ! name of xysize dimension")') TRIM(cxysize)
    WRITE(numnam,'("  cx_range  =  ",a10,"  ! name of x_range  variable")') TRIM(cx_range)
    WRITE(numnam,'("  cy_range  =  ",a10,"  ! name of y_range  variable")') TRIM(cy_range)
    WRITE(numnam,'("  cz_range  =  ",a10,"  ! name of z_range  variable")') TRIM(cz_range)
    WRITE(numnam,'("  cspacing  =  ",a10,"  ! name of spacing  variable")') TRIM(cspacing)
    WRITE(numnam,'("  cdim      =  ",a10,"  ! name of dimension  variable")') TRIM(cdim)
    WRITE(numnam,'("  cdatain   =  ",a10,"  ! name of input   variable")') TRIM(cdatain)
    WRITE(numnam,'("/")')
    WRITE(numnam,'("! output grid netcdf file")')
    WRITE(numnam,'("&output_grd")')
    WRITE(numnam,'("  cnamenc   =  ",a10,"  ! name of nc output file")') TRIM(cnamenc)
    WRITE(numnam,'("  clon      =  ",a10,"  ! name of longitude dimension")') TRIM(clon)
    WRITE(numnam,'("  clat      =  ",a10,"  ! name of latitude dimension")') TRIM(clat)
    WRITE(numnam,'("  cdataout  =  ",a10,"  ! name of output   variable")') TRIM(cdataout)
    WRITE(numnam,'("/")')
    CLOSE(numnam)
  END SUBROUTINE print_template


END PROGRAM grd2nc
