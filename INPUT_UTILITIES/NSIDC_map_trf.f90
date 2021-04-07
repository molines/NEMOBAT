PROGRAM NSIDC_map_trf
  !!======================================================================
  !!                     ***  PROGRAM  NSIDC_map_trf >  ***
  !!=====================================================================
  !!  ** Purpose : Transform  Lon/lat into NSDIC sterographic projection
  !!               or NSDIC stereographic projection to lon/lat
  !!
  !!  ** Method  : Use routine from https://github.com/nsidc/polar_stereo_legacy.git
  !!               Use double precision 
  !!
  !! History :  1.0  : 04/2021  : J.M. Molines : Port routine to F90 Dr norm
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  !! Coddyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE


  INTEGER(KIND=4) :: ijarg, iargc, narg       ! line parser
  INTEGER(KIND=4) :: npiglo, npjglo, npx, npy ! dimensions
  INTEGER(KIND=4) :: ihem   ! hemisphere : 1= North 2=South
  INTEGER(KIND=4) :: isense ! sense of the transformation

  ! netcdf files
  INTEGER(KIND=4) :: ierr, ncid, id   ! input
  INTEGER(KIND=4) :: ncio, idxt, idxu, idxv, idxf, idxb
  INTEGER(KIND=4) ::       idyt, idyu, idyv, idyf, idyb
  INTEGER(KIND=4) ::       idlont, idlonu, idlonv, idlonf, idlonb
  INTEGER(KIND=4) ::       idlatt, idlatu, idlatv, idlatf, idlatb


  REAL(KIND=8) :: dslat, de, de2, dre, dpi, dsgn, delta       ! physical constant (global)
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dlont, dlatt   ! lon lat array 
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: dxt, dyt       ! X Y arrays

  CHARACTER(LEN=80) :: cf_lonlat   ! lon lat file
  CHARACTER(LEN=80) :: cf_xy       ! X Y file
  CHARACTER(LEN=80) :: cldum       ! dummy char variable
  CHARACTER(LEN=80) :: clhem       ! hemisphere
  !-----------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  NSIDC_map_trf -lonlat LONLAT-file -xy XY-file  -conv SENSE '
     PRINT *,'             [-hem N/S ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This program allow to switch between NSIDC stereographic coordinates'
     PRINT *,'        to Longitude/Latitude or reverse, according to SENSE.'
     PRINT *,'          NEMO mesh_hgr conventions are assumed for LONLAT-file, and NSIDC '
     PRINT *,'        conventions are assumed for XY-file'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -lonlat LONLAT-file : pass the name of the file with Longitude/latitude' 
     PRINT *,'             This can be either an input or an output file depending on SENSE.'
     PRINT *,'       -xy XY-file : pass the name of the file with x,y NSIDC coordinates.' 
     PRINT *,'             This can be either an input or an output file depending on SENSE.'
     PRINT *,'       -conv SENSE: Indicate the sense of the convertion :'
     PRINT *,'             SENSE = 1 : convert LONLAT to XY'
     PRINT *,'             SENSE = 2 : convert XY to LONLAT. In this case, hemisphere should'
     PRINT *,'                         be specified with option -hem'
     PRINT *,'             SENSE = 3 : convert BEDMACHINE file (with only x,y)  to lon lat'
     PRINT *,'                         In this case, hemisphere should be specified with '
     PRINT *,'                         option -hem'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'        -hem N/S : specify the hemisphere (N or S) in which the XY file will'
     PRINT *,'                   transformed.'
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'      ' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : either LONLAT or XY file'
     PRINT *,'         variables : for LONLAT : gphi?, glam? (? in {t u v f})'
     PRINT *,'                     for XY  : x?,y? (? in {t u v f} )'
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  clhem = 'N/A'
  isense = -1
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-lonlat'   ) ; CALL getarg(ijarg, cf_lonlat ) ; ijarg=ijarg+1
     CASE ( '-xy'       ) ; CALL getarg(ijarg, cf_xy     ) ; ijarg=ijarg+1
     CASE ( '-conv'     ) ; CALL getarg(ijarg, cldum     ) ; ijarg=ijarg+1 ; READ(cldum,*) isense
        ! option
     CASE ( '-hem'      ) ; CALL getarg(ijarg, clhem     ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! Some sanity check
  IF ( isense == -1 ) THEN 
     PRINT *,' ERROR : you must specify the conversion sense with -conv'
     PRINT *,'         1 : lonlat --> XY '
     PRINT *,'         2 : XY  --> lonlat'
     PRINT *,'         3 : XY  --> lonlat (BEDMACHINE file)'
     STOP 99
  ENDIF

  SELECT CASE (clhem)
  CASE ( 'N'   ) 
     dsgn  = +1.d0
     delta = 45.d0
  CASE ( 'S'   ) 
     dsgn  = -1.d0
     delta =  0.d0
  CASE ( 'N/A' )
     IF ( isense > 1  ) THEN
        PRINT *,' ERROR : for XY --> lonlat conversion, the hemisphere must be specified'
        PRINT *,'         use options -hem N or -hem S !'
        STOP 98
     ENDIF
  END SELECT

  ! Physical constant initialisation
  dslat = 70.d0         ! degree = latitude (N/S) of tangeant cone
  dre   = 6378.273d0    ! km : earth radius
  de2   = .006693883d0  ! Earth elllipticity
  dpi   = ACOS(-1.d0)   ! pi
  de    = SQRT(de2)     

  ! Open files according to sense
  SELECT CASE (isense)
  CASE ( 1 )   ! LONLAT to XY
     ierr = NF90_OPEN( cf_lonlat, NF90_NOWRITE, ncid )
     ierr = NF90_INQ_DIMID( ncid,'x',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo) 
     ierr = NF90_INQ_DIMID( ncid,'y',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo) 
     ALLOCATE ( dlont(npiglo,npjglo), dlatt(npiglo,npjglo) )   ! input arrays
     ALLOCATE ( dxt(npiglo,npjglo)  , dyt(npiglo,npjglo)   )   ! output arrays
     CALL CreateXY
     ! The program will convert all the NEMO glam. and gphi. into x. and y. 
     ! T points
     CALL toxy('T')
     ! U points
     CALL toxy('U')
     ! V points
     CALL toxy('V')
     ! F points
     CALL toxy('F')
     ierr = NF90_CLOSE( ncio) 
  CASE ( 2  )   ! XY to LONLAT
     ierr = NF90_OPEN( cf_xy, NF90_NOWRITE, ncid )
     ierr = NF90_INQ_DIMID( ncid,'x',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npx) 
     ierr = NF90_INQ_DIMID( ncid,'y',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npy) 
     ALLOCATE ( dlont(npx,npy), dlatt(npx,npy) )   ! output arrays
     ALLOCATE ( dxt(npx,npy)  , dyt(npx,npy)   )   ! input arrays
     CALL CreateLONLAT
     ! The program will convert all the x. and y. to glam. and gphi.
     ! T points
     CALL tolonlat('T')
     ! U points
     CALL tolonlat('U')
     ! V points
     CALL tolonlat('V')
     ! F points
     CALL tolonlat('F')
     ierr = NF90_CLOSE( ncio) 
  CASE ( 3  )   ! XY  BEDMACHINE to LONLAT
     ierr = NF90_OPEN( cf_xy, NF90_NOWRITE, ncid )
     ierr = NF90_INQ_DIMID( ncid,'x',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npx) 
     ierr = NF90_INQ_DIMID( ncid,'y',id ) ; ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npy) 
     ALLOCATE ( dlont(npx,npy), dlatt(npx,npy) )   ! output arrays
     ALLOCATE ( dxt(npx,1)  , dyt(1,npy)   )   ! input arrays
     CALL CreateLONLAT_BED
     CALL tolonlat('B')
     ierr = NF90_CLOSE( ncio) 
  CASE DEFAULT;   PRINT *,' ERROR : specify -conv 1 or -conv 2 '; STOP
  END SELECT
  !-----------------------------------------------------------

CONTAINS
  SUBROUTINE CreateXY
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateXY  ***
    !!
    !! ** Purpose :  Create output file for reciving sterographic coordinates
    !!               computed from longitude and latitudes
    !!
    !! ** Method  :  Netcdf primitives. Take care of using NO CLOBBER to
    !!               avoid erasing existing file in case of error in the
    !!               argument. 
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: idx, idy 

    ierr = NF90_CREATE(cf_xy,OR(NF90_NETCDF4,NF90_NOCLOBBER), ncio )
    IF ( ierr /= NF90_NOERR) THEN
       PRINT *,' ERROR ', NF90_STRERROR(ierr)
       PRINT *,' rm ',TRIM(cf_xy),' ???'
       STOP
    ENDIF
    ! define dimensions
    ierr = NF90_DEF_DIM(ncio, 'x', npiglo, idx )
    ierr = NF90_DEF_DIM(ncio, 'y', npjglo, idy )
    ! define variables
    ierr = NF90_DEF_VAR(ncio, 'xt', NF90_INT, (/idx,idy/), idxt )
    ierr = NF90_PUT_ATT (ncio, idxt,'long_name','Cartesian x-coordinate')
    ierr = NF90_PUT_ATT (ncio, idxt,'standard_name','projection_x_coordinate')
    ierr = NF90_PUT_ATT (ncio, idxt,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'xu', NF90_INT, (/idx,idy/), idxu )
    ierr = NF90_PUT_ATT (ncio, idxu,'long_name','Cartesian x-coordinate')
    ierr = NF90_PUT_ATT (ncio, idxu,'standard_name','projection_x_coordinate')
    ierr = NF90_PUT_ATT (ncio, idxu,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'xv', NF90_INT, (/idx,idy/), idxv )
    ierr = NF90_PUT_ATT (ncio, idxv,'long_name','Cartesian x-coordinate')
    ierr = NF90_PUT_ATT (ncio, idxv,'standard_name','projection_x_coordinate')
    ierr = NF90_PUT_ATT (ncio, idxv,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'xf', NF90_INT, (/idx,idy/), idxf )
    ierr = NF90_PUT_ATT (ncio, idxf,'long_name','Cartesian x-coordinate')
    ierr = NF90_PUT_ATT (ncio, idxf,'standard_name','projection_x_coordinate')
    ierr = NF90_PUT_ATT (ncio, idxf,'units','meter')

    ierr = NF90_DEF_VAR(ncio, 'yt', NF90_INT, (/idx,idy/), idyt )
    ierr = NF90_PUT_ATT (ncio, idyt,'long_name','Cartesian y-coordinate')
    ierr = NF90_PUT_ATT (ncio, idyt,'standard_name','projection_y_coordinate')
    ierr = NF90_PUT_ATT (ncio, idyt,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'yu', NF90_INT, (/idx,idy/), idyu )
    ierr = NF90_PUT_ATT (ncio, idyu,'long_name','Cartesian y-coordinate')
    ierr = NF90_PUT_ATT (ncio, idyu,'standard_name','projection_y_coordinate')
    ierr = NF90_PUT_ATT (ncio, idyu,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'yv', NF90_INT, (/idx,idy/), idyv )
    ierr = NF90_PUT_ATT (ncio, idyv,'long_name','Cartesian y-coordinate')
    ierr = NF90_PUT_ATT (ncio, idyv,'standard_name','projection_y_coordinate')
    ierr = NF90_PUT_ATT (ncio, idyv,'units','meter')
    ierr = NF90_DEF_VAR(ncio, 'yf', NF90_INT, (/idx,idy/), idyf )
    ierr = NF90_PUT_ATT (ncio, idyf,'long_name','Cartesian y-coordinate')
    ierr = NF90_PUT_ATT (ncio, idyf,'standard_name','projection_y_coordinate')
    ierr = NF90_PUT_ATT (ncio, idyf,'units','meter')
    ierr = NF90_ENDDEF(ncio)

  END SUBROUTINE CreateXY

  SUBROUTINE CreateLONLAT
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateLONLAT  ***
    !!
    !! ** Purpose :  Create output file for reciving longitude and latitude
    !!               computed from x, y 
    !!
    !! ** Method  :  Netcdf primitives. Take care of using NO CLOBBER to
    !!               avoid erasing existing file in case of error in the
    !!               argument. 
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: idx, idy 

    ierr = NF90_CREATE(cf_lonlat,OR(NF90_NETCDF4,NF90_NOCLOBBER), ncio )
    ! define dimensions
    ierr = NF90_DEF_DIM(ncio, 'x', npx, idx )
    ierr = NF90_DEF_DIM(ncio, 'y', npy, idx )
    ! define variables
    ierr = NF90_DEF_VAR(ncio, 'glamt', NF90_DOUBLE, (/idx,idy/), idlont )
    ierr = NF90_DEF_VAR(ncio, 'glamu', NF90_DOUBLE, (/idx,idy/), idlonu )
    ierr = NF90_DEF_VAR(ncio, 'glamv', NF90_DOUBLE, (/idx,idy/), idlonv )
    ierr = NF90_DEF_VAR(ncio, 'glamf', NF90_DOUBLE, (/idx,idy/), idlonf )

    ierr = NF90_DEF_VAR(ncio, 'gphit', NF90_DOUBLE, (/idx,idy/), idlatt )
    ierr = NF90_DEF_VAR(ncio, 'gphiu', NF90_DOUBLE, (/idx,idy/), idlatu )
    ierr = NF90_DEF_VAR(ncio, 'gphiv', NF90_DOUBLE, (/idx,idy/), idlatv )
    ierr = NF90_DEF_VAR(ncio, 'gphif', NF90_DOUBLE, (/idx,idy/), idlatf )
    ierr = NF90_ENDDEF(ncio)

  END SUBROUTINE CreateLONLAT

  SUBROUTINE CreateLONLAT_BED
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE CreateLONLAT_BED  ***
    !!
    !! ** Purpose :  Create output file for reciving longitude and latitude
    !!               computed from x, y  from a BEDMACHINE File
    !!
    !! ** Method  :  Netcdf primitives. Take care of using NO CLOBBER to
    !!               avoid erasing existing file in case of error in the
    !!               argument. 
    !!----------------------------------------------------------------------
    INTEGER(KIND=4) :: idx, idy

    ierr = NF90_CREATE(cf_lonlat,OR(NF90_NETCDF4,NF90_NOCLOBBER), ncio )
    IF ( ierr /= NF90_NOERR) THEN
       PRINT *,' ERROR :', NF90_STRERROR(ierr)
       PRINT *,' rm ',TRIM(cf_lonlat),' ???'
       STOP 97
    ENDIF
    ! define dimensions
    ierr = NF90_DEF_DIM(ncio, 'x', npx, idx )
    ierr = NF90_DEF_DIM(ncio, 'y', npy, idy )
    ! define variables
    ierr = NF90_DEF_VAR(ncio, 'nav_lon', NF90_DOUBLE, (/idx,idy/), idlont )

    ierr = NF90_DEF_VAR(ncio, 'nav_lat', NF90_DOUBLE, (/idx,idy/), idlatt )
    ierr = NF90_ENDDEF(ncio)

  END SUBROUTINE CreateLONLAT_BED

  SUBROUTINE toxy ( cd_typ)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE toxy  ***
    !!
    !! ** Purpose :  Perform transformation from lon lat to xy for cd_typ point 
    !!
    !! ** Method  :  wrapper   for calling mapll
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_typ

    INTEGER(KIND=4)   :: ji, jj
    INTEGER(KIND=4)   :: id, idx, idy
    CHARACTER(LEN=80) :: cl_lon, cl_lat
    !! ----------------------------------------------------------------------
    SELECT CASE (cd_typ)
    CASE ('T') 
       cl_lon = "glamt" ; cl_lat = "gphit"
       idx    =  idxt   ; idy    = idyt
    CASE ('U') 
       cl_lon = "glamu" ; cl_lat = "gphiu"
       idx    =  idxu   ; idy    = idyu
    CASE ('V') 
       cl_lon = "glamv" ; cl_lat = "gphiv"
       idx    =  idxv   ; idy    = idyv
    CASE ('F') 
       cl_lon = "glamf" ; cl_lat = "gphif"
       idx    =  idxf   ; idy    = idyf
    END SELECT

    ierr = NF90_INQ_VARID(ncid, cl_lon, id ) ; ierr = NF90_GET_VAR(ncid, id, dlont, start=(/1,1,1/), count=(/npiglo, npjglo,1/) )
    ierr = NF90_INQ_VARID(ncid, cl_lat, id ) ; ierr = NF90_GET_VAR(ncid, id, dlatt, start=(/1,1,1/), count=(/npiglo, npjglo,1/) )
    dlont=(dlont+delta)*dpi/180.d0
    dlatt=dlatt*dpi/180.d0

    IF ( dlatt(1,1 ) > 0 ) THEN
       dsgn  = 1.d0
       delta = 45.d0
    ELSE
       dsgn  = -1.d0
       delta =  0.d0
    ENDIF

    DO jj=1,npjglo
       DO ji = 1, npiglo
          CALL mapll(dxt(ji,jj), dyt(ji,jj), dlatt(ji,jj), dlont(ji,jj))
       ENDDO
    ENDDO
    !pass to meters
    dxt = dxt * 1000.d0
    dyt = dyt * 1000.d0
    ierr = NF90_PUT_VAR(ncio, idx, dxt(:,:),start=(/1,1/), count=(/npiglo, npjglo/) )
    ierr = NF90_PUT_VAR(ncio, idy, dyt(:,:),start=(/1,1/), count=(/npiglo, npjglo/) )

  END SUBROUTINE toxy

  SUBROUTINE tolonlat ( cd_typ)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE tolonlat  ***
    !!
    !! ** Purpose :  Perform transformation from xy  to lonlat for cd_typ point 
    !!
    !! ** Method  :  wrapper   for calling mapll
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_typ

    INTEGER(KIND=4)   :: ji, jj
    INTEGER(KIND=4)   :: id, idlon, idlat
    CHARACTER(LEN=80) :: cl_x, cl_y

    !! ----------------------------------------------------------------
    SELECT CASE (cd_typ)
    CASE ('T') 
       cl_x = "xt" ; cl_y = "yt"
       idlon    =  idlont   ; idlat    = idlatt
    CASE ('U') 
       cl_x = "xu" ; cl_y = "yu"
       idlon    =  idlonu   ; idlat    = idlatu
    CASE ('V') 
       cl_x = "xv" ; cl_y = "yv"
       idlon    =  idlonv   ; idlat    = idlatv
    CASE ('F') 
       cl_x = "xf" ; cl_y = "yf"
       idlon    =  idlonf   ; idlat    = idlatf
    CASE ('B')   ! bedmachine file only X
       cl_x = "x" ; cl_y = "y"
       idlon    =  idlont   ; idlat    = idlatt
       ierr = NF90_INQ_VARID(ncid, cl_x, id ) ; ierr = NF90_GET_VAR(ncid, id, dxt(:,1), start=(/1/), count=(/npx/) )
       ierr = NF90_INQ_VARID(ncid, cl_y, id ) ; ierr = NF90_GET_VAR(ncid, id, dyt(1,:), start=(/1/), count=(/npy/) )
       ! transform to km
       dxt = dxt/1000.d0
       dyt = dyt/1000.d0

       DO jj=1,npy
          DO ji = 1, npx
             CALL mapxy(dxt(ji,1), dyt(1,jj), dlatt(ji,jj), dlont(ji,jj) )
          ENDDO
       ENDDO
       ! transform to degree and offset
       dlont = dlont*180.d0/dpi - delta
       dlatt = dlatt*180.d0/dpi

       ierr = NF90_PUT_VAR(ncio, idlon, dlont(:,:),start=(/1,1/), count=(/npx, npy/) )
       ierr = NF90_PUT_VAR(ncio, idlat, dlatt(:,:),start=(/1,1/), count=(/npx, npy/) )
       RETURN

    END SELECT


    ierr = NF90_INQ_VARID(ncid, cl_x, id ) ; ierr = NF90_GET_VAR(ncid, id, dxt, start=(/1,1/), count=(/npx, npy/) )
    ierr = NF90_INQ_VARID(ncid, cl_y, id ) ; ierr = NF90_GET_VAR(ncid, id, dyt, start=(/1,1/), count=(/npx, npy/) )
    ! transform to km
    dxt = dxt/1000.d0
    dyt = dyt/1000.d0

    DO jj=1,npy
       DO ji = 1, npx
          CALL mapxy(dxt(ji,jj), dyt(ji,jj), dlatt(ji,jj), dlont(ji,jj) )
       ENDDO
    ENDDO
    ! transform to degree and offset
    dlont = dlont*180.d0/dpi - delta
    dlatt = dlatt*180.d0/dpi
    ierr = NF90_PUT_VAR(ncio, idlon, dlont(:,:),start=(/1,1/), count=(/npx, npy/) )
    ierr = NF90_PUT_VAR(ncio, idlat, dlatt(:,:),start=(/1,1/), count=(/npx, npy/) )

  END SUBROUTINE tolonlat

  SUBROUTINE mapxy (ddx,ddy,ddlat,ddlon)
    !!*****************************************************************************
    !!                                                                            *
    !!                                                                            *
    !!    DESCRIPTION:                                                            *
    !!                                                                            *
    !!    This subroutine converts from Polar Stereographic (ddx,ddy) coordinates *
    !!    to geodetic latitude and longitude for the polar regions. The equations *
    !!    are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *
    !!    Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *
    !!    Printing Office.  See JPL Technical Memorandum 3349-85-101 for further  *
    !!    details.                                                                *
    !!                                                                            *
    !!                                                                            *
    !!    ARGUMENTS:                                                              *
    !!                                                                            *
    !!    Variable    Type        I/O    Description                              *
    !!                                                                            *
    !!    ddx        REAL*4        I     Polar Stereographic ddx Coordinate (km)  *
    !!    ddy        REAL*4        I     Polar Stereographic ddy Coordinate (km)  *
    !!    ddlat      REAL*4        O     Geodetic Latitude (degrees, +90 to -90)  *
    !!    ddlon      REAL*4        O     Geodetic Longitude (degrees, 0 to 360)   *
    !!                                                                            *
    !!                                                                            *
    !!                  Written by C. S. Morris - April 29, 1985                  *
    !!                  Revised by C. S. Morris - December 11, 1985               *
    !!                                                                            *
    !!                  Revised by V. J. Troisi - January 1990
    !!                  dsgn - provide hemisphere dependency (+/- 1)
    !!
    !!*****************************************************************************
    REAL(KIND=8), INTENT(in ) :: ddx, ddy
    REAL(KIND=8), INTENT(out) :: ddlat, ddlon

    REAL(KIND=8) :: dl_sl
    REAL(KIND=8) :: dl_rho, dl_cm, dl_tan, dl_chi
    !!*****************************************************************************
    dl_sl = dslat/180.d0*dpi
    dl_rho=SQRT(ddx**2+ddy**2)
    IF (dl_rho  <= 0.1) THEN
       ddlat=dpi/2.d0*dsgn
       ddlon=0.0
       RETURN
    ENDIF
    dl_cm=COS(dl_sl)/SQRT(1.0-de2*(SIN(dl_sl)**2))
    dl_tan=TAN((dpi/4.0)-(dl_sl/(2.0)))/((1.0-de*SIN(dl_sl))/(1.0+de*SIN(dl_sl)))**(de/2.0)
    IF ( ABS(dslat-90.) < 1.d-5 ) THEN
       dl_tan=dl_rho*SQRT((1.+de)**(1.+de)*(1.-de)**(1.-de))/2./dre
    ELSE
       dl_tan=dl_rho*dl_tan/(dre*dl_cm)
    END IF
    dl_chi=(dpi/2.0)-2.0*ATAN(dl_tan)
    ddlat=dl_chi+((de2/2.0)+(5.0*de2**2.0/24.0)+(de2**3.0/12.0))*SIN(2*dl_chi)+& 
         &         ((7.0*de2**2.0/48.0)+(29.0*de2**3/240.0))*SIN(4.0*dl_chi)+&
         &          (7.0*de2**3.0/120.0)*SIN(6.0*dl_chi)
    ddlat=dsgn*ddlat
    ddlon=ATAN2(dsgn*ddx,-dsgn*ddy)
    ddlon=dsgn*ddlon

  END SUBROUTINE mapxy

  SUBROUTINE MAPLL (ddx,ddy,ddlat,ddlon)
    !!*****************************************************************************
    !!                                                                            *
    !!                                                                            *
    !!    DESCRIPTION:                                                            *
    !!                                                                            *
    !!    This subroutine converts from geodetic latitude and longitude to Polar  *
    !!    Stereographic (ddx,ddy) coordinates for the polar regions.  The equations   *
    !!    are from Snyder, J. P., 1982,  Map Projections Used by the U.S.         *
    !!    Geological Survey, Geological Survey Bulletin 1532, U.S. Government     *
    !!    Printing Office.  See JPL Technical Memorandum 3349-85-101 for further  *
    !!    details.                                                                *
    !!                                                                            *
    !!                                                                            *
    !!    ARGUMENTS:                                                              *
    !!                                                                            *
    !!    Variable    Type        I/O    Description                              *
    !!                                                                            *
    !!    ddlat       REAL*4        I     Geodetic Latitude (degrees, +90 to -90)  *
    !!    ddlon      REAL*4        I     Geodetic Longitude (degrees, 0 to 360)   *
    !!    ddx          REAL*4        O     Polar Stereographic ddx Coordinate (km)    *
    !!    ddy          REAL*4        O     Polar Stereographic ddy Coordinate (km)    *
    !!                                                                            *
    !!                                                                            *
    !!                  Written by C. S. Morris - April 29, 1985                  *
    !!                  Revised by C. S. Morris - December 11, 1985               *
    !!                                                                     	      *
    !!                  Revised by V. J. Troisi - January 1990                    *
    !!                  dsgn - provides hemisphere dependency (+/- 1)              *
    !!                  Revised by Xiaoming Li - October 1996                     *
    !!                  Corrected equation for dl_rho                                *
    !!*****************************************************************************
    REAL(KIND=8), INTENT(out) :: ddx, ddy
    REAL(KIND=8), INTENT(in ) :: ddlat, ddlon

    REAL(KIND=8) :: dl_sl
    REAL(KIND=4) :: dl_rho, dl_cm, dl_tan, dl_chi, dl_tc, dl_mc 
    !!*****************************************************************************
    !     Compute ddx and ddy in grid coordinates.
    IF ( ABS(ddlat) >= dpi/2.d0 ) THEN
       ddx=0.0
       ddy=0.0
       RETURN
    ENDIF
    dl_tan=TAN(dpi/4.-ddlat/2.)/((1.-de*SIN(ddlat))/(1.+de*SIN(ddlat)))**(de/2.)
    IF (ABS(90.-dslat)  < 1.d-5) THEN
       dl_rho=2.*dre*dl_tan/((1.+de)**(1.+de)*(1.-de)**(1.-de))**(1/2.)
    ELSE
       dl_sl=dslat*dpi/180.d0
       dl_tc=TAN(dpi/4.-dl_sl/2.)/((1.-de*SIN(dl_sl))/(1.+de*SIN(dl_sl)))**(de/2.)
       dl_mc=COS(dl_sl)/SQRT(1.0-de2*(SIN(dl_sl)**2))
       dl_rho=dre*dl_mc*dl_tan/dl_tc
    END IF
    ddy=-dl_rho*dsgn*COS(dsgn*ddlon)
    ddx= dl_rho*dsgn*SIN(dsgn*ddlon)
  END SUBROUTINE MAPLL
END PROGRAM NSIDC_map_trf
