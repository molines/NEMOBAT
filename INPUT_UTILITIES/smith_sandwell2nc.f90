PROGRAM smith_sandwell2nc
  !!---------------------------------------------------------------------------
  !!               *** PROGRAM  smith_sandwell2nc  ***
  !!
  !!  Purpose: convert Smith Sandwell native img file into standard netcdf
  !!
  !!  Method: Read the native binary file ( coded in big_endian) 
  !!          - create longitude and latitude ( taking base information on README file)
  !!            (img data are on Mercator grid, 1' at equator)
  !!          - shift the data to put them in the frame -180 180 to be the same
  !!            as etopo1 etc...
  !!
  !! history: J.M. Molines, April 2009.
  !!---------------------------------------------------------------------------
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
  !!  $Rev: 228 $
  !!  $Id: smith_sandwell2nc.f90 228 2009-04-29 14:19:25Z forge $
  !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE

  INTEGER, PARAMETER :: jpi=21600, jpj=17280
  REAL(kind=8), PARAMETER :: rplatmin=-80.738 , rplatmax=80.738
  REAL(kind=8), PARAMETER :: rplonmin=-180. 
  REAL(KIND=8) :: dlon,rplonmax, y,ymin,ymax, pi, zphi
  REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: rlon, rlat
  INTEGER(kind=2), DIMENSION (jpi,jpj) :: idata
  INTEGER(kind=2), DIMENSION (jpi) ::     itmp
  INTEGER :: ji,jj,ji0,numdta=10
  INTEGER :: istatus, ncid, idvar,idx, idy , idlon, idlat

  ALLOCATE (rlon(jpi), rlat(jpj) )
  OPEN(numdta, FILE='topo_11.1.img', FORM='UNFORMATTED', ACCESS='DIRECT', RECL=jpi*2 )
  DO jj=jpj,1,-1
   READ(numdta,REC=jj) (idata(ji,jpj-jj+1),ji=1,jpi)
  ENDDO
  
  
  ! create rlon :
  dlon=1.d0/60.d0  ! 1' in degrees
  DO ji=1,jpi
   rlon(ji)=rplonmin+dlon/2. + (ji -1)*dlon
  END DO
  ! create rlat : Mercator projection:
  pi=acos(-1.d0)
  zphi=rplatmax*pi/180.d0          ! in radians
  ymax=log(tan(pi/4.d0+zphi/2.d0)) ! mercator
  zphi=rplatmin*pi/180.d0          ! in radians
  ymin=log(tan(pi/4.d0+zphi/2.d0)) ! mercator

  DO jj=1,jpj
    y=(jj-1.d0)/(jpj-1.d0)*(ymax - ymin ) + ymin ! y is evenly spaced from min to max
    zphi= 2.d0*atan( exp(y))-pi/2.d0 ! latitude correspondig to y (radians)
    rlat(jj)=zphi*180.d0/pi
  END DO
  ! which ji for rlon=0
  DO ji=1,jpi
   IF ( rlon(ji) >= 0 ) THEN
    EXIT
   ENDIF
  END DO
  ji0=ji ; print *, rlon(ji0)

  ! data are 0-360 put them -180 180
  DO jj=1,jpj
    itmp(:)=idata(:,jj)
    idata(ji0:jpi,jj)=itmp(1:jpi-ji0+1)
    idata(1:ji0-1,jj)=itmp(jpi-ji0:jpi)
  END DO
  
  ! create netcdf
  istatus = NF90_CREATE ('topo_11.1.nc',NF90_CLOBBER,ncid)
  istatus = NF90_DEF_DIM(ncid,'lon',jpi,idx)
  istatus = NF90_DEF_DIM(ncid,'lat',jpj,idy)
  istatus = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx/), idlon )
  istatus = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idy/), idlat )
  istatus = NF90_DEF_VAR(ncid,'z',NF90_SHORT,(/idx,idy/), idvar )
  istatus = NF90_ENDDEF(ncid)

  istatus = NF90_PUT_VAR(ncid,idlon,rlon)
  istatus = NF90_PUT_VAR(ncid,idlat,rlat)
  istatus = NF90_PUT_VAR(ncid,idvar,idata)
  istatus = NF90_CLOSE(ncid)

END PROGRAM smith_sandwell2nc
