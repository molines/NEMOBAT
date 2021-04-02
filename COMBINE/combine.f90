!-------------------------------------------------------------------
PROGRAM combine
  !-------------------------------------------------------------------
  ! read  two bathymetry files interpolated on an ORCA grid:
  !    - ETOPO2
  !    - gebco
  ! be careful! gebco goes only up to 88N. Where it is flagged
  ! (moce = -1) take etopo2. 
  !  
  !-------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: combine.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  REAL,   DIMENSION(:,:),ALLOCATABLE :: bathy1, bathy2, batnew, gphit
  CHARACTER(len=140)                 :: chfilin1,chfilin2,chanswer,chfilout
  CHARACTER(len=80 )                 :: chtreatment,chcomment

  INTEGER                            :: NCID_IN, NCID_OUT, XID, YID, TID
  INTEGER                            :: IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN, IDVAR3_IN
  INTEGER                            :: IDVAR4_IN
  INTEGER                            :: LONID, LATID
  INTEGER                            :: STATUS
  INTEGER ,DIMENSION(2)              :: DIM_msk  
  INTEGER ,DIMENSION(:,:) ,ALLOCATABLE :: moce 

  INTEGER                            :: ji,jj, jpi, jpj, nperio 
  REAL                               :: zmin,zmax,zeps
  INTEGER  :: narg, iargc

  narg=iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' USAGE: combine gebco.nc  etopo_filterd.nc'
     PRINT *,'  Both files have same grid (model grid)'
     PRINT *,'  output is done on combine_etopo1_gebco1.nc'
     STOP
  ENDIF
  CALL getarg(1,chfilin1)
  CALL getarg(2,chfilin2)

  !-------------------------------------------------------------------
  ! 0 - DEFINITIONS of variables 
  ! ------------------------------------------------------------------
  !
  zmin = 200
  zmax = 300
  chtreatment = 'Bathy after merging etopo2 and gebco1'
  WRITE(chcomment,'("zmin to keep gebco",i4," zmax:",i4)' ),INT(zmin), INT(zmax)

  !
  ! ------------------------------------------------------------------
  ! 1 -  open bathymetry files
  ! ------------------------------------------------------------------
  !                            default value for file name: 
  chfilout = 'combine_etopo1_gebco1.nc'

  !                             read actual bathymetry filename
  !  READ dimensions in netcdf file
  !
  STATUS = NF90_OPEN(TRIM(chfilin1), NF90_NOWRITE, NCID_IN)
  IF (STATUS /= NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin1)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  STATUS = NF90_INQ_DIMID(NCID_IN, 'x', LONID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpi)
  STATUS = NF90_INQ_DIMID(NCID_IN, 'y', LATID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpj)
  PRINT *,' reading netcdf file, dimension jpi=',jpi,' jpj=',jpj
  IF (jpi <=0 .OR.jpj <= 0) THEN
     PRINT *, ' dimensions read in .nc file is zero'
     STOP
  ENDIF
  STATUS=NF90_INQ_VARID(NCID_IN, 'Bathymetry'   , IDVAR3_IN)
  !
  ! Allocate arrays
  ! ------------------------------------
  ALLOCATE (bathy1(jpi,jpj)  ,  bathy2(jpi,jpj), batnew(jpi,jpj) )
  ALLOCATE (moce(jpi,jpj) ,gphit(jpi,jpj) )
  !
  ! ------------------------------------------------------------------
  ! 2 -  read bathymetry file
  ! ------------------------------------------------------------------
  !   Reading the bathymetry in  meters
  !   for Gebco file, we also read flag moce. 
  !
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR3_IN, bathy1  )
  STATUS=NF90_INQ_VARID(NCID_IN, 'moce'   , IDVAR1_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR1_IN, moce  )
  STATUS=NF90_CLOSE(NCID_IN)

  STATUS = NF90_OPEN(TRIM(chfilin2), NF90_NOWRITE, NCID_IN)
  IF (STATUS /= NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin2)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  STATUS=NF90_INQ_VARID(NCID_IN, 'Bathymetry'   , IDVAR3_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR3_IN, bathy2  )
  STATUS=NF90_INQ_VARID(NCID_IN, 'nav_lat'   , IDVAR4_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR4_IN, gphit  )
  !
  !          the bathymetry must be positive for OPA model.
  !          if it is not, convert it to positive values
  bathy1 = ABS(bathy1)
  bathy2 = ABS(bathy2)

  PRINT *, ' Bathy read ok, max =',MAXVAL(bathy1)
  PRINT *, ' Bathy read ok, min =',MINVAL(bathy2)

  ! ------------------------------------------------------------------
  !  3 -   combines bathymetries 
  ! -------------------------------------------------------------------
  DO jj=1,jpj
     DO ji=1,jpi
        IF (moce(ji,jj) .EQ. -1 .OR. gphit(ji,jj)  < -58. ) THEN
           !               there was no gebco1 information in that area
           batnew(ji,jj)=bathy2(ji,jj)
        ELSE
           IF      ( bathy1(ji,jj) <= zmin ) THEN
              batnew(ji,jj)=bathy1(ji,jj)
           ELSE IF ( bathy1(ji,jj) <= zmax ) THEN
              zeps = (bathy1(ji,jj)-zmin)/(zmax-zmin)
              batnew(ji,jj)=bathy1(ji,jj)*(1-zeps)+bathy2(ji,jj)*zeps
           ELSE
              batnew(ji,jj)= bathy2(ji,jj)
           ENDIF
        ENDIF
     ENDDO
  ENDDO

  ! -------------------------------------------------------------------
  !
  !   Apply lateral boundary conditions
  !   calling lbc ( ptab,kjpi,kjpj, kjpk, ktype,ksgn,nperio)
  !   The bathymetry file is at "T" points, no sign change: ktype = 1, ksgn = 1.
  !   it is a 2D array so kjpk = 1
  ! -------------------------------------------------------------------
  !
  IF (jpi== 182 .AND.jpj==149) THEN 
     nperio= 4
     PRINT *, 'This is an ORCA R2 configuration, nperio=4'
  ELSEIF (jpi== 722 .AND.jpj==511) THEN
     nperio= 6
     PRINT *, 'This is an ORCA R05 configuration, nperio=6'
  ELSEIF (jpi== 1442 .AND.jpj==1021) THEN
     nperio= 4
     PRINT *, 'This is an ORCA R025 configuration, nperio=4'
  ELSE
     PRINT *, 'your configuration is unknown: you need to choose nperio '
     PRINT *, '     jperio= 0, closed,      jperio= 1, cyclic east-west'
     PRINT *, '     jperio= 2, equatorial symmetric'
     PRINT *, '     jperio= 3, north fold with T-point pivot; jperio= 4, the same + cyclic east-west'   
     PRINT *, '     jperio= 5, north fold with F-point pivot; jperio= 6, the same + cyclic east-west'   
     READ (5,*) nperio
  ENDIF
  CALL lbc ( batnew,jpi,jpj, 1, 1 ,1, nperio)   

  ! -------------------------------------------------------------------
  !  Write new bathy file
  ! -------------------------------
  CALL rewrite_bathy_nc(batnew, jpi,jpj,chfilin2,chfilout, chtreatment,chcomment)

  STOP

  !*********************************************************


END PROGRAM combine

