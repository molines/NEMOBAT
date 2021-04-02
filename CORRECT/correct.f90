PROGRAM correct
  !-------------------------------------------------------------------
  !   Apply hand corrections from a "journal" file generated 
  !   by IDL or cdfbathy
  !  
  !-------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: correct.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  REAL,   DIMENSION(:,:),ALLOCATABLE :: bathy, batold
  CHARACTER(len=140)                 :: chfilin1, chfilout
  CHARACTER(len=80 )                 :: chtreatment,chcomment

  INTEGER                            :: NCID_IN, NCID_OUT, XID, YID, TID
  INTEGER                            :: IDVAR2, IDVAR3,  IDVAR2_IN, IDVAR3_IN
  INTEGER                            :: LONID, LATID
  INTEGER                            :: istatus
  INTEGER ,DIMENSION(2)              :: DIM_msk  
  INTEGER ,DIMENSION(:,:) ,ALLOCATABLE :: moce 
  INTEGER                            :: ji,jj, jpi, jpj, nperio , ipt

  INTEGER :: narg, iargc

  narg=iargc()
  IF ( narg /= 2 ) THEN
     PRINT *,' USAGE: correct bathyorig bathymodified '
     STOP
  ENDIF

  CALL getarg(1,chfilin1)
  CALL getarg(2,chfilout)
  !-------------------------------------------------------------------
  ! 0 - DEFINITIONS of variables 
  ! ------------------------------------------------------------------
  !
  chtreatment = 'correct bathy'
  !
  ! ------------------------------------------------------------------
  ! 1 -  open bathymetry files
  ! ------------------------------------------------------------------

  !  READ dimensions in netcdf file
  !
  istatus = NF90_OPEN(TRIM(chfilin1), NF90_NOWRITE, NCID_IN)
  IF (istatus /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(istatus)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin1)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF

  istatus = NF90_INQ_DIMID(NCID_IN, 'x', LONID)      
  istatus = NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpi)
  istatus = NF90_INQ_DIMID(NCID_IN, 'y', LATID)      
  istatus = NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpj)
  PRINT *,' reading netcdf file, dimension jpi=',jpi,' jpj=',jpj
  IF (jpi <= 0.OR.jpj <= 0) THEN
     PRINT *, ' dimensions read in .nc file is zero'
     STOP
  ENDIF
  istatus=NF90_INQ_VARID(NCID_IN, 'Bathymetry'   , IDVAR3_IN)
  !
  ! Allocate arrays
  ! ------------------------------------
  ALLOCATE (bathy(jpi,jpj), batold(jpi,jpj) )
  !
  ! ------------------------------------------------------------------
  ! 2 -  read bathymetry file
  ! ------------------------------------------------------------------
  !   Reading the bathymetry in  meters
  !
  istatus = NF90_GET_VAR(NCID_IN, IDVAR3_IN, bathy  )
  !
  !          the bathymetry must be positive for OPA model.
  !          if it is not, convert it to positive values
  bathy = ABS(bathy)

  PRINT *, ' Bathy read ok, max =',MAXVAL(bathy)
  PRINT *, ' Bathy read ok, min =',MINVAL(bathy)

  ! ------------------------------------------------------------------
  !  3 -   corrects bathymetries 
  !        applies same modifications as the old bathy, excepted for 
  !        indonesian throughflow where modifications are replaced by Ariane's.
  ! -------------------------------------------------------------------
  batold(:,:) = bathy(:,:)
  INCLUDE  'hand_modif.h90'
  !  now count how many points were changed
  !
  ipt = 0
  DO jj=1,jpj
     DO ji=1,jpi
        IF (batold(ji,jj) /= bathy(ji,jj)) ipt=ipt+1
     ENDDO
  ENDDO
  PRINT *, ' Number of points changed:', ipt

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
  CALL lbc ( bathy,jpi,jpj, 1, 1 ,1, nperio)   

  ! -------------------------------------------------------------------
  !  Write new bathy file
  ! -------------------------------
  CALL rewrite_bathy_nc(bathy, jpi,jpj,chfilin1,chfilout, chtreatment,chcomment)

  STOP

  !*********************************************************


END PROGRAM correct

