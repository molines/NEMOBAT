!
!     subroutine rewrite_bathy_nc
!  rewrite a new bathy file by copying variables from a
!  previous one (only bathy will change)
!
!-----------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: WBATHGOOD.F90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
SUBROUTINE WBATHGOOD(bathy, jpi,jpj,chfilin,chfilout,chtreatment,chcomment)
!-----------------------------------------------------------------
! 
IMPLICIT NONE      

INTEGER                ,INTENT(in) :: jpi,jpj
REAL(KIND = 4),DIMENSION(jpi,jpj),INTENT(in) :: bathy
character(len=140) ,INTENT(in)     :: chfilin,chfilout
character(len=80),  INTENT(in)     :: chcomment,chtreatment
!
!  LOCAL VARIABLES
!
INTEGER                            :: jpiread,jpjread, istart1, icount1
INTEGER                            :: NCID_IN, NCID_OUT, XID, YID, ZID, TID, NT
INTEGER                            :: IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN, IDVAR3_IN
INTEGER                            :: IDVAR4, IDVAR5
INTEGER                            :: LONID, LATID
INTEGER                            :: STATUS
INTEGER ,dimension(2)              :: DIM_msk
INTEGER ,dimension(4)              :: DIM_xyzt, icount, istart
! 
REAL(KIND = 4) 				   ::   bat_min,bat_max, deptht, time_counter
REAL(KIND = 4) 				   ::   lon_min, lon_max, lat_min, lat_max
REAL(KIND = 4),   DIMENSION(:,:),ALLOCATABLE ::   nav_lon,nav_lat
REAL(KIND = 4),   DIMENSION(:,:,:,:),ALLOCATABLE ::   bathy4
include 'netcdf.inc'
!-----------------------------------------------------------------
!          EXECUTION 
!-----------------------------------------------------------------
!  OPEN INPUT FILE AND CHECK 
!
       PRINT *, ' -------------------------------------------'
       PRINT *, '             REWRITE_BATHY_NC routine'
       PRINT *, ' -------------------------------------------'
       PRINT *, 'Opening  the input file:' ,trim(chfilin)
       STATUS = NF_OPEN(trim(chfilin), NF_NOWRITE, NCID_IN)
       IF (STATUS .NE. NF_NOERR) THEN
           PRINT *, NF_STRERROR(STATUS)
           PRINT *, ' '
           PRINT *, 'Could not open the file:' ,trim(chfilin)
           PRINT *, ' '
           STOP 'We stop here'
       ENDIF

      STATUS = NF_INQ_DIMID(NCID_IN, 'x', LONID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LONID, jpiread)
      STATUS = NF_INQ_DIMID(NCID_IN, 'y', LATID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LATID, jpjread)
      print *,'--    Reading input file, dimension jpi=',jpiread,' jpj=',jpjread
      if (jpi.ne.jpiread.or.jpj.ne.jpjread) then
          print *, ' dimensions read in .nc file do not agree with parameters'
          stop
      endif
!
!   ALLOCATE arrays for longitude and latitude 
! 
         ALLOCATE(nav_lon(jpi,jpj),nav_lat(jpi,jpj))
!
!            read longitude and latitude 
!
      STATUS=NF_INQ_VARID(NCID_IN, 'nav_lon', IDVAR1_IN)
      STATUS=NF_INQ_VARID(NCID_IN, 'nav_lat', IDVAR2_IN)
      STATUS = NF_GET_VAR_REAL(NCID_IN, IDVAR1_IN, nav_lon)
      STATUS = NF_GET_VAR_REAL(NCID_IN, IDVAR2_IN, nav_lat)
!
!   calculate max and min value of input array bathy 
! 
       bat_max= MAXVAL(bathy)
       bat_min= MINVAL(bathy)
       lat_max= MAXVAL(nav_lat)
       lat_min= MINVAL(nav_lat)
       lon_min= MINVAL(nav_lon)
       lon_max= MAXVAL(nav_lon)
       ALLOCATE ( bathy4( JPI, JPJ, 1, 1))
       bathy4(:,:,1,1)=bathy
       NT = 1
       deptht=1
       time_counter=1
       PRINT *, bathy4(50,50,1,1)
!--------------------------------------------------------------
!  NOW deal with output file 
!--------------------------------------------------------------
       PRINT *, 'Opening  the output file:' ,trim(chfilout)
       STATUS = NF_CREATE(trim(chfilout), NF_CLOBBER, NCID_OUT)
       IF (STATUS .NE. NF_NOERR) THEN
           PRINT *, NF_STRERROR(STATUS)
           PRINT *, ' '
           PRINT *, 'Could not open the output file ', trim(chfilout)
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
       STATUS = NF_DEF_DIM(NCID_OUT, 'x', JPI, XID)
       STATUS = NF_DEF_DIM(NCID_OUT, 'y', JPJ, YID)
       STATUS = NF_DEF_DIM(NCID_OUT, 'deptht', 1, ZID)
       STATUS = NF_DEF_DIM(NCID_OUT, 'time_counter', NF_UNLIMITED, TID)
!       STATUS = NF_DEF_DIM(NCID_OUT, 'time_counter', 1, TID)
       DIM_msk(1)=XID
       DIM_msk(2)=YID
       DIM_xyzt= (/ XID, YID, ZID, TID /)
!
!  Definition of variables
!
      STATUS=NF_DEF_VAR(NCID_OUT,'nav_lon',NF_FLOAT,2,DIM_msk,IDVAR1) 
      STATUS=NF_DEF_VAR(NCID_OUT,'nav_lat',NF_FLOAT,2,DIM_msk,IDVAR2)
      STATUS=NF_DEF_VAR(NCID_OUT,'deptht',NF_FLOAT,1,ZID,IDVAR3)
      STATUS=NF_DEF_VAR(NCID_OUT,'time_counter',NF_FLOAT,1,TID,IDVAR4)
!                                     ---------------  
!                                     Dim_msk has 2 dimensions  XID and YID
      STATUS=NF_DEF_VAR(NCID_OUT,'Bathymetry',NF_FLOAT,4,DIM_xyzt,IDVAR5) 

!
!  Attributs
!
!       Copy of attributes from the input file  


      STATUS=NF_PUT_ATT_TEXT(NCID_OUT,IDVAR1,'units'    , 12,'degrees_east')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT,IDVAR1,'long_name',9,'Longitude')
      STATUS=NF_PUT_ATT_REAL(NCID_OUT,IDVAR1,'valid_min',NF_FLOAT,1,lon_min)
      STATUS=NF_PUT_ATT_REAL(NCID_OUT,IDVAR1,'valid_max',NF_FLOAT,1,lon_max)
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR1,'nav_model', 12,'Default grid')

      STATUS=NF_PUT_ATT_TEXT(NCID_OUT,IDVAR2,'units'    , 13,'degrees_north')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT,IDVAR2,'long_name',8,'Latitude')
      STATUS=NF_PUT_ATT_REAL(NCID_OUT,IDVAR2,'valid_min',NF_FLOAT,1,lat_min)
      STATUS=NF_PUT_ATT_REAL(NCID_OUT,IDVAR2,'valid_max',NF_FLOAT,1,lat_max)
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR2,'nav_model', 12,'Default grid')
      
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR3,'units', 6, 'meters') 
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR3,'long_name', 5, 'Depth')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR3,'nav_model', 12,'Default grid')

      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR4,'units', 34, 'seconds since 0000-01-01 00:00:00 ')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR4,'calendar', 6, 'noleap')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR4,'long_name', 9, 'Time axis')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR4,'Title', 4,'Time')
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR4,'time_origin',20,'0000-JAN-01 00:00:00')
      
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR5,'units', 6, 'meters') 
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR5,'long_name', 11, 'Bathymetry')
      STATUS=NF_PUT_ATT_REAL(NCID_OUT, IDVAR5,'valid_min',NF_FLOAT,1,bat_min) 
      STATUS=NF_PUT_ATT_REAL(NCID_OUT, IDVAR5,'valid_max',NF_FLOAT,1,bat_max)
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, IDVAR5,'Title', 15,'Bathymetry ROSS')

      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, NF_GLOBAL,'Treatment',len(trim(chtreatment)),chtreatment)
      STATUS=NF_PUT_ATT_TEXT(NCID_OUT, NF_GLOBAL,'Comment',  len(trim(chcomment)),chcomment)

      STATUS=NF_ENDDEF(NCID_OUT)
!
!  Save variables
!      PRINT *, bathy4(50,50,1,1)
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR1, nav_lon)
      IF (STATUS .NE. NF_NOERR) PRINT *,'7: ERROR, ierr=',STATUS

      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR2, nav_lat)
      IF (STATUS .NE. NF_NOERR) PRINT *,'7: ERROR, ierr=',STATUS

      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR3, deptht)
      IF (STATUS .NE. NF_NOERR) PRINT *,'7: ERROR, ierr=',STATUS

      istart1= 1
      icount1= 1
!      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR4,time_counter)
      STATUS=NF_PUT_VARA_REAL(NCID_OUT, IDVAR4, istart1, icount1, time_counter)
      IF (STATUS .NE. NF_NOERR) PRINT *,'7: ERROR, ierr=',STATUS

      
      istart= (/ 1 , 1 , 1 , 1 /)
      icount= (/ JPI , JPJ , 1 , 1 /)
!      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR5, bathy4 )
      STATUS=NF_PUT_VARA_REAL(NCID_OUT, IDVAR5, istart, icount, bathy4 )
      IF (STATUS .NE. NF_NOERR) PRINT *,'7: ERROR, ierr=',STATUS

!
!  Close files
!
      STATUS=NF_CLOSE(NCID_OUT)
      STATUS=NF_CLOSE(NCID_IN)

      print * ,'--   Output file written:',trim(chfilout)

      DEALLOCATE(nav_lon,nav_lat)
      RETURN
      END
