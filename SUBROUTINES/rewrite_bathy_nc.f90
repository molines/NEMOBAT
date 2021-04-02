SUBROUTINE REWRITE_BATHY_NC(bathy, jpi,jpj,chfilin,chfilout,chtreatment,chcomment)
  !----------------------------------------------------------------
  !     subroutine rewrite_bathy_nc
  !  rewrite a new bathy file by copying variables from a
  !  previous one (only bathy will change)
  !
  !-----------------------------------------------------------------
  ! 
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: rewrite_bathy_nc.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE      

  INTEGER                ,INTENT(in) :: jpi,jpj
  REAL,DIMENSION(jpi,jpj),INTENT(in) :: bathy
  CHARACTER(len=140) ,INTENT(in)     :: chfilin,chfilout
  CHARACTER(len=80),  INTENT(in)     :: chcomment,chtreatment
  !
  !  LOCAL VARIABLES
  !
  INTEGER                            :: jpiread,jpjread
  INTEGER                            :: NCID_IN, NCID_OUT, XID, YID
  INTEGER                            :: IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN, IDVAR3_IN
  INTEGER                            :: LONID, LATID
  INTEGER                            :: STATUS
  INTEGER ,DIMENSION(2)              :: DIM_msk
  !      
  REAL                               ::   bat_min,bat_max
  REAL,   DIMENSION(:,:),ALLOCATABLE ::   nav_lon,nav_lat
  !-----------------------------------------------------------------
  !          EXECUTION 
  !-----------------------------------------------------------------
  !  OPEN INPUT FILE AND CHECK 
  !
  PRINT *, ' -------------------------------------------'
  PRINT *, '             REWRITE_BATHY_NC routine'
  PRINT *, ' -------------------------------------------'
  PRINT *, 'Opening  the input file:' ,TRIM(chfilin)
  STATUS = NF90_OPEN(TRIM(chfilin), NF90_NOWRITE, NCID_IN)
  IF (STATUS .NE. NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  STATUS = NF90_INQ_DIMID(NCID_IN, 'x', LONID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpiread)
  STATUS = NF90_INQ_DIMID(NCID_IN, 'y', LATID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpjread)
  PRINT *,'--    Reading input file, dimension jpi=',jpiread,' jpj=',jpjread
  IF (jpi.NE.jpiread.OR.jpj.NE.jpjread) THEN
     PRINT *, ' dimensions read in .nc file do not agree with parameters'
     STOP
  ENDIF
  !
  !   ALLOCATE arrays for longitude and latitude 
  ! 
  ALLOCATE(nav_lon(jpi,jpj),nav_lat(jpi,jpj))
  !
  !            read longitude and latitude 
  !
  STATUS=NF90_INQ_VARID(NCID_IN, 'nav_lon', IDVAR1_IN)
  STATUS=NF90_INQ_VARID(NCID_IN, 'nav_lat', IDVAR2_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR1_IN, nav_lon)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR2_IN, nav_lat)
  !
  !   calculate max and min value of input array bathy 
  ! 
  bat_max= MAXVAL(bathy)
  bat_min= MINVAL(bathy)

  !--------------------------------------------------------------
  !  NOW deal with output file 
  !--------------------------------------------------------------
  PRINT *, 'Opening  the output file:' ,TRIM(chfilout)
  STATUS = NF90_CREATE(TRIM(chfilout), NF90_CLOBBER, NCID_OUT)
  IF (STATUS .NE. NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the output file ', TRIM(chfilout)
     PRINT *, ' '
     STOP 'We stop here'
  ELSE
     PRINT *,'--   Output file open OK'
  ENDIF
  !
  !  Definition of dimensions.
  !  note we need an unlimited time axis for today (2003) version of IOIPSL
  !  so that the bathymetry can be read by the OPA code ... (!!!)
  !
  STATUS = NF90_DEF_DIM(NCID_OUT, 'x', JPI, XID)
  STATUS = NF90_DEF_DIM(NCID_OUT, 'y', JPJ, YID)
  !
  !  Definition of variables
  !
  STATUS=NF90_DEF_VAR(NCID_OUT,'nav_lon',NF90_FLOAT,(/XID,YID/),IDVAR1) 
  STATUS=NF90_DEF_VAR(NCID_OUT,'nav_lat',NF90_FLOAT,(/XID,YID/),IDVAR2) 
  !                                     ---------------  
  !                                     Dim_msk has 2 dimensions  XID and YID
  STATUS=NF90_DEF_VAR(NCID_OUT,'Bathymetry' ,NF90_FLOAT,(/XID,YID/),IDVAR3) 

  !
  !  Attributs
  !
  !       Copy of attributes from the input file  

  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR1_IN,'units'    ,NCID_OUT,IDVAR1)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR1_IN,'long_name',NCID_OUT,IDVAR1)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR1_IN,'valid_min',NCID_OUT,IDVAR1)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR1_IN,'valid_max',NCID_OUT,IDVAR1)

  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR2_IN,'units'    ,NCID_OUT,IDVAR2)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR2_IN,'long_name',NCID_OUT,IDVAR2)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR2_IN,'valid_min',NCID_OUT,IDVAR2)
  STATUS=NF90_COPY_ATT(NCID_IN,IDVAR2_IN,'valid_max',NCID_OUT,IDVAR2)

  STATUS=NF90_PUT_ATT(NCID_OUT, IDVAR3,'units', 'meters') 
  STATUS=NF90_PUT_ATT(NCID_OUT, IDVAR3,'long_name', 'Median depth by area')
  STATUS=NF90_PUT_ATT(NCID_OUT,IDVAR3,'valid_min',bat_min) 
  STATUS=NF90_PUT_ATT(NCID_OUT,IDVAR3,'valid_max',bat_max)
  STATUS=NF90_PUT_ATT(NCID_OUT, NF90_GLOBAL,'Title','Bathymetry ORCA0.25')
  STATUS=NF90_PUT_ATT(NCID_OUT, NF90_GLOBAL,'Treatment',chtreatment)
  STATUS=NF90_PUT_ATT(NCID_OUT, NF90_GLOBAL,'Comment', chcomment)

  STATUS=NF90_ENDDEF(NCID_OUT)
  !
  !  Save variables
  !
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR1, nav_lon)
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR2, nav_lat)
  !
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR3, bathy  )

  !
  !  Close files
  !
  STATUS=NF90_CLOSE(NCID_OUT)
  STATUS=NF90_CLOSE(NCID_IN)

  PRINT * ,'--   Output file written:',TRIM(chfilout)

  DEALLOCATE(nav_lon,nav_lat)
  RETURN
END SUBROUTINE REWRITE_BATHY_NC
