PROGRAM mkcoor
 !-----------------------------------------------------------------------------
 !       *** PROGRAM mkcoor ***
 !
 ! Purpose: create a pseudo coordinate file with glamt, gphit,e1t,e2t
 !          from lon(ji) and lat(jj)
 !
 ! Method : read lon,lat from input file, build 2D glamt,gphit,
 !          use dist function to approximate e1t, e2t:
 !          e1t(ji) = 1/2*dist (T(ji+1), T(ji-1))
 !          e2t(jj) = 1/2*dist (T(jj+1), T(jj-1))
 !
 ! history:
 !       Original J.M. Molines 25/04/2009
 !-----------------------------------------------------------------------------
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
  !!  $Rev: 228 $
  !!  $Id: mk_coor.f90 228 2009-04-29 14:19:25Z forge $
  !! -------------------------------------------------------------------------
 USE netcdf
 IMPLICIT NONE
 INTEGER :: narg, iargc
 INTEGER :: ji,jj
 INTEGER :: ni,nj
 INTEGER :: ncid, idx, idy, idvar, istatus
 INTEGER :: idlam, idphi, ide1, ide2

 REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: rlon, rlat
 REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE :: gphit, glamt, e1t, e2t
 CHARACTER(LEN=80) :: cfilin, cfilout='pseudo_coordinates.nc'
 narg=iargc()
 IF ( narg == 0 ) THEN
   PRINT * ,' USAGE: mk_coor  input_file '
   PRINT *, '    input file is a file with only LON(ji) and LAT(jj) variables'
   PRINT *, '   produce pseudo_coordinates.nc with variables name '
   PRINT *, ' gphit, glamt, e1t, e2t + nav_lon, nav_lat '
   STOP
 ENDIF
 
 CALL getarg(1, cfilin )
 
 istatus=NF90_OPEN(cfilin,NF90_NOWRITE,ncid)
 istatus=NF90_INQ_DIMID(ncid,'lon',idx)
 istatus=NF90_INQUIRE_DIMENSION(ncid,idx,len=ni)
 istatus=NF90_INQ_DIMID(ncid,'lat',idy)
 istatus=NF90_INQUIRE_DIMENSION(ncid,idy,len=nj)

 ALLOCATE (rlon(ni), rlat(nj), glamt(ni,nj), gphit(ni,nj), e1t(ni,nj), e2t(ni,nj))

 istatus=NF90_INQ_VARID(ncid,'lon',idvar)
 istatus=NF90_GET_VAR(ncid,idvar,rlon)
 istatus=NF90_INQ_VARID(ncid,'lat',idvar)
 istatus=NF90_GET_VAR(ncid,idvar,rlat)
 istatus=NF90_CLOSE(ncid)
 
 ! build glamt
 DO jj=1,nj
   glamt(:,jj)=rlon(:)
 END DO
 ! build gphit
 DO ji=1,ni
   gphit(ji,:)=rlat(:)
 END DO

 !build e1
 DO jj=1,nj
  DO ji=2,ni-1
   e1t(ji,jj)=0.5* dist(rlon(ji-1),rlon(ji+1),rlat(jj),rlat(jj)  )*1000.
  END DO
 END DO
 ! fill first and last row
 e1t(1,:) =e1t(2,:)
 e1t(ni,:)=e1t(ni-1,:)
 
 !build e2
 DO jj=2,nj-1
  DO ji=1,ni
   e2t(ji,jj)=0.5* dist(rlon(ji),rlon(ji),rlat(jj-1),rlat(jj+1)  )*1000.
  END DO
 END DO
 ! fill first and last row
 e2t(:,1) =e2t(:,2)
 e2t(:,nj)=e2t(:,nj-1)

 ! create outputfile
 istatus=NF90_CREATE(cfilout,NF90_CLOBBER,ncid)
 ! dims ... as in nemo
 istatus=NF90_DEF_DIM(ncid,'x',ni,idx)
 istatus=NF90_DEF_DIM(ncid,'y',nj,idy)
 ! vars
 istatus=NF90_DEF_VAR(ncid,'glamt',NF90_FLOAT,(/idx,idy/),idlam)
 istatus=NF90_DEF_VAR(ncid,'gphit',NF90_FLOAT,(/idx,idy/),idphi)
 istatus=NF90_DEF_VAR(ncid,'e1t',NF90_FLOAT,(/idx,idy/),ide1)
 istatus=NF90_DEF_VAR(ncid,'e2t',NF90_FLOAT,(/idx,idy/),ide2)
 !
 istatus=NF90_ENDDEF(ncid)
 istatus=NF90_PUT_VAR(ncid,idlam,glamt)
 istatus=NF90_PUT_VAR(ncid,idphi,gphit)
 istatus=NF90_PUT_VAR(ncid,ide1,e1t)
 istatus=NF90_PUT_VAR(ncid,ide2,e2t)
 istatus=NF90_CLOSE(ncid)
 
 
CONTAINS

  FUNCTION dist(plona,plonb,plata,platb)
    !!----------------------------------------------------------
    !!           ***  FUNCTION  DIST  ***  
    !!
    !!  ** Purpose : Compute the distance (km) between
    !!          point A (lona, lata) and B(lonb,latb)
    !!
    !!  ** Method : Compute the distance along the orthodromy
    !!        
    !! * history : J.M. Molines in CHART, f90, may 2007
    !!----------------------------------------------------------
    IMPLICIT NONE
    ! Argument 
    REAL(KIND=8), INTENT(in) :: plata, plona, platb, plonb
    REAL(KIND=8) :: dist
    ! Local variables
    REAL(KIND=8),SAVE ::  zlatar, zlatbr, zlonar, zlonbr
    REAL(KIND=8) ::  zpds
    REAL(KIND=8),SAVE :: zux, zuy, zuz
    REAL(KIND=8) :: zvx, zvy, zvz

    REAL(KIND=8), SAVE :: prevlat=-1000., prevlon=-1000, zr, zpi, zconv
    LOGICAL :: lfirst=.TRUE.

    ! initialise some values at first call
    IF ( lfirst ) THEN
       lfirst=.FALSE.
       ! constants
       zpi=ACOS(-1.)
       zconv=zpi/180.  ! for degree to radian conversion
       ! Earth radius
       zr=(6378.137+6356.7523)/2.0 ! km
    ENDIF

    ! compute these term only if they differ from previous call
    IF ( plata /= prevlat .OR. plona /= prevlon) THEN
       zlatar=plata*zconv
       zlonar=plona*zconv
       zux=COS(zlonar)*COS(zlatar)
       zuy=SIN(zlonar)*COS(zlatar)
       zuz=SIN(zlatar)
       prevlat=plata
       prevlon=plona
    ENDIF

    zlatbr=platb*zconv
    zlonbr=plonb*zconv
    zvx=COS(zlonbr)*COS(zlatbr)
    zvy=SIN(zlonbr)*COS(zlatbr)
    zvz=SIN(zlatbr)

    zpds=zux*zvx+zuy*zvy+zuz*zvz

    IF (zpds >= 1.) THEN
       dist=0.
    ELSE
       dist=zr*ACOS(zpds)
    ENDIF
  END FUNCTION dist
 
END PROGRAM mkcoor

