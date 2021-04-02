!-------------------------------------------------------------------
PROGRAM verify 
!-------------------------------------------------------------------
!  verify hand correction in a bathy file 
!  (history: this program written because apperently hand-corrections
!  were not correctly taken into account in OPABAT025 bathy)
!  
!-------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: verify.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  IMPLICIT NONE
  REAL,   DIMENSION(:,:),ALLOCATABLE :: bathy,  batnew
  CHARACTER(len=140)                 :: chfilin1,chanswer,chfilout
  CHARACTER(len=80 )                 :: chtreatment,chcomment

  INTEGER                            :: NCID_IN, NCID_OUT, XID, YID, TID
  INTEGER                            :: IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN, IDVAR3_IN
  INTEGER                            :: LONID, LATID
  INTEGER                            :: STATUS
  INTEGER ,dimension(2)              :: DIM_msk  
  INTEGER ,dimension(:,:) ,ALLOCATABLE :: moce 

  INTEGER                            :: ji,jj, jpi, jpj, nperio 


include "netcdf.inc"
!-------------------------------------------------------------------
! 0 - DEFINITIONS of variables 
! ------------------------------------------------------------------
!
  chtreatment = 'verify bathy'

!
! ------------------------------------------------------------------
! 1 -  open bathymetry files
! ------------------------------------------------------------------
!                            default value for file name: 
   chfilin1  = 'OPABAT025-OLD/bathy_meter_treated_ORCA_R025_checked_coast_nobug.nc'
   chfilin1  = 'Bathy05_AMT_ori+TP_modif.nc'
   chfilout = 'toto.nc'

!                             read actual bathymetry filename
   print *,' file name for bathy to verify, ''y'' for default:',trim(chfilin1)
     read(5,*) chanswer
     if (trim(chanswer).ne.'y') chfilin1=chanswer
   print * ,' file for output, ''y'' for default:',trim(chfilout)
     read(5,*) chanswer
     if (trim(chanswer).ne.'y') chfilout=chanswer

!  READ dimensions in netcdf file
!
       STATUS = NF_OPEN(trim(chfilin1), NF_NOWRITE, NCID_IN)
       IF (STATUS .NE. NF_NOERR) THEN
           PRINT *, NF_STRERROR(STATUS)
           PRINT *, ' '
           PRINT *, 'Could not open the file:' ,trim(chfilin1)
           PRINT *, ' '
           STOP 'We stop here'
       ENDIF
      STATUS = NF_INQ_DIMID(NCID_IN, 'x', LONID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LONID, jpi)
      STATUS = NF_INQ_DIMID(NCID_IN, 'y', LATID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LATID, jpj)
      print *,' reading netcdf file, dimension jpi=',jpi,' jpj=',jpj
      if (jpi.le.0.or.jpj.le.0) then
          print *, ' dimensions read in .nc file is zero'
          stop
      endif
      STATUS=NF_INQ_VARID(NCID_IN, 'Bathymetry'   , IDVAR3_IN)
!
! Allocate arrays
! ------------------------------------
 ALLOCATE (bathy(jpi,jpj), batnew(jpi,jpj) )
!
! ------------------------------------------------------------------
! 2 -  read bathymetry file
! ------------------------------------------------------------------
!   Reading the bathymetry in  meters
!   for Gebco file, we also read flag moce. 
!
       STATUS = NF_GET_VAR_REAL(NCID_IN, IDVAR3_IN, bathy  )
       STATUS=NF_INQ_VARID(NCID_IN, 'moce'   , IDVAR1_IN)

!
!          the bathymetry must be positive for OPA model.
!          if it is not, convert it to positive values
  bathy = abs(bathy)

       print *, ' Bathy read ok, max =',maxval(bathy)
       print *, ' Bathy read ok, min =',minval(bathy)


! ------------------------------------------------------------------
!  3 -   verifies bathymetries 
! -------------------------------------------------------------------
include 'orca05_hand_modif_2.h90'
! -------------------------------------------------------------------
!
!   Apply lateral boundary conditions
!   calling lbc ( ptab,kjpi,kjpj, kjpk, ktype,ksgn,nperio)
!   The bathymetry file is at "T" points, no sign change: ktype = 1, ksgn = 1.
!   it is a 2D array so kjpk = 1
! -------------------------------------------------------------------
!
    if (jpi== 182 .AND.jpj==149) then 
         nperio= 4
         print *, 'This is an ORCA R2 configuration, nperio=4'
    elseif (jpi== 722 .AND.jpj==511) then
         nperio= 6
         print *, 'This is an ORCA R05 configuration, nperio=6'
    elseif (jpi== 1442 .AND.jpj==1021) then
         nperio= 4
         print *, 'This is an ORCA R025 configuration, nperio=4'
    else
         print *, 'your configuration is unknown: you need to choose nperio '
         print *, '     jperio= 0, closed,      jperio= 1, cyclic east-west'
         print *, '     jperio= 2, equatorial symmetric'
         print *, '     jperio= 3, north fold with T-point pivot; jperio= 4, the same + cyclic east-west'   
         print *, '     jperio= 5, north fold with F-point pivot; jperio= 6, the same + cyclic east-west'   
         read (5,*) nperio
    endif
!    CALL lbc ( batnew,jpi,jpj, 1, 1 ,1, nperio)   

! -------------------------------------------------------------------
!  Write new bathy file
! -------------------------------
!  CALL rewrite_bathy_nc(batnew, jpi,jpj,chfilin2,chfilout, chtreatment,chcomment)

STOP

!*********************************************************


END PROGRAM verify

