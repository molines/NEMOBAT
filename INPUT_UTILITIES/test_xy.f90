PROGRAM NSIDC_map_trf
  !!======================================================================
  !!                     ***  PROGRAM  NSIDC_map_trf >  ***
  !!=====================================================================
  !!  ** Purpose : Transform  Lon/lat into NSDIC sterographic projection
  !!               or NSDIC stereographic projection to lon/lat
  !!
  !!  ** Method  : Use routine from https://github.com/nsidc/polar_stereo_legacy.git
  !!
  !! History :  1.0  : 04/2021  : J.M. Molines : Port routine to F90 Dr norm
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !! Coddyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  REAL(KIND=8) :: dslat, de, de2, dre, dpi, dsgn, delta
  REAL(KIND=8) :: dlat, dlon, dx, dy

  INTEGER(KIND=4) :: ihem   ! hemisphere : 1= North 2=South
  INTEGER(KIND=4) :: ijarg, iargc, narg
  INTEGER(KIND=4) :: isense
  INTEGER(KIND=4) :: npiglo, npjglo, npx, npy


  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=80) :: clhem
  !-----------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  test_xy   -x x  -y y       [-hem N/S ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This program allow to switch between NSIDC stereographic coordinates'
     PRINT *,'        to Longitude/Latitude or reverse, according to SENSE.'
     PRINT *,'          NEMO mesh_hgr conventions are assumed for LONLAT-file, and NSIDC '
     PRINT *,'        conventions are assumed for XY-file'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'         -x x '
     PRINT *,'         -y y '
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
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-x'   ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) dx ; print *, dx
     CASE ( '-y'   ) ; CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1 ; READ(cldum,*) dy ; print *, dy
        ! option
     CASE ( '-hem'      ) ; CALL getarg(ijarg, clhem     ) ; ijarg=ijarg+1 ; print *, TRIM(clhem)
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO


  SELECT CASE (clhem)
  CASE ( 'N'   ) 
     dsgn  = +1.d0
     delta = 45.d0
  CASE ( 'S'   ) 
     dsgn  = -1.d0
     delta =  0.d0
  END SELECT

  ! Physical constant initialisation
  dslat = 70.d0         ! degree = latitude (N/S) of tangeant cone
  dre   = 6378.273d0    ! km : earth radius
  de2   = 0.006693883d0  ! Earth elllipticity
  dpi   = ACOS(-1.d0)
  de    = SQRT(de2)

  CALL mapxy(dx/1000.d0,dy/1000.d0, dlat, dlon )
  dlon=dlon*180.d0/dpi - delta
  dlat=dlat *180.d0/dpi
  print *, ' X   Y    : ', dx, dy
  print *, ' LON LAT  : ', dlon, dlat
  dlon = (dlon+delta)*dpi/180.d0
  dlat = dlat*dpi/180.d0
  CALL mapll(dx,dy, dlat,dlon)
  print *, ' X   Y    : ', dx*1000.d0, dy*1000.d0

CONTAINS


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
    REAL(KIND=8) :: zrho, zcm, ztan, zchi
    !!*****************************************************************************
    dl_sl = dslat/180.d0*dpi
    zrho=SQRT(ddx**2+ddy**2)
    IF (zrho  <= 0.1) THEN
       ddlat=dpi/2.d0*dsgn
       ddlon=0.0d0
       RETURN
    ENDIF
    zcm=COS(dl_sl)/SQRT(1.0d0-de2*(SIN(dl_sl)**2))
    ztan=TAN((dpi/4.0d0)-(dl_sl/(2.0d0)))/((1.0-de*SIN(dl_sl))/(1.0+de*SIN(dl_sl)))**(de/2.0)
    IF ( ABS(dslat-90.) < 1.d-5 ) THEN
       ztan=zrho*SQRT((1.+de)**(1.+de)*(1.-de)**(1.-de))/2./dre
       print *, ' PIF'
    ELSE
       ztan=zrho*ztan/(dre*zcm)
       print *, ' PAF'
    END IF
    zchi=(dpi/2.0)-2.0*ATAN(ztan)
    ddlat=zchi+((de2/2.0)+(5.0*de2**2.0/24.0)+(de2**3.0/12.0))*SIN(2*zchi)+& 
        &         ((7.0*de2**2.0/48.0)+(29.0*de2**3/240.0))*SIN(4.0*zchi)+&
        &          (7.0*de2**3.0/120.0)*SIN(6.0*zchi)
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
    !!                  Corrected equation for zrho                                *
    !!*****************************************************************************
    REAL(KIND=8), INTENT(out) :: ddx, ddy
    REAL(KIND=8), INTENT(in ) :: ddlat, ddlon

    REAL(KIND=8) :: dl_sl
    REAL(KIND=4) :: zrho, zcm, ztan, zchi, ztc, zmc 
    !!*****************************************************************************
    !     Compute ddx and ddy in grid coordinates.
    IF ( ABS(ddlat) >= dpi/2.d0 ) THEN
       ddx=0.0
       ddy=0.0
       RETURN
    ENDIF
    ztan=TAN(dpi/4.-ddlat/2.)/((1.-de*SIN(ddlat))/(1.+de*SIN(ddlat)))**(de/2.)
    IF (ABS(90.-dslat)  < 1.d-5) THEN
       zrho=2.*dre*ztan/((1.+de)**(1.+de)*(1.-de)**(1.-de))**(1/2.)
    ELSE
       dl_sl=dslat*dpi/180.d0
       ztc=TAN(dpi/4.-dl_sl/2.)/((1.-de*SIN(dl_sl))/(1.+de*SIN(dl_sl)))**(de/2.)
       zmc=COS(dl_sl)/SQRT(1.0-de2*(SIN(dl_sl)**2))
       zrho=dre*zmc*ztan/ztc
    END IF
    ddy=-zrho*dsgn*COS(dsgn*ddlon)
    ddx= zrho*dsgn*SIN(dsgn*ddlon)
  END SUBROUTINE MAPLL
END PROGRAM NSIDC_map_trf
