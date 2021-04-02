PROGRAM batcoast
  !----------------------------------------------------------------------------
  ! Written by E. Remy from batmod.F (AmT) included in the OPABAT program,  
  ! (LODYC) and checkpog.f written by M. Bremond and L. Siefridt (?)
  !
  !
  ! This subroutines modifies a first guess bathymetric field in meters 
  ! and subsequently in levels.
  !
  ! In the present version the modifications are:
  !   - apply a coast mask;
  !   - force each ocean point to have at least 3 levels;
  !   - force periodicity or a line of zeroes east and west.
  !
  !----------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: batcoast.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  !-------------------------------------------------------------------                                       
  INTEGER                :: JPI,JPJ 
  !-------------------------------------------------------------------
  REAL,   DIMENSION(:,:),ALLOCATABLE    ::  batg, batmp
  INTEGER,   DIMENSION(:,:),ALLOCATABLE ::  mcoast
  INTEGER,   DIMENSION(:),ALLOCATABLE   ::  indi, indj
  INTEGER kiter,ji,jj,jk,kbat0,kmask1,kbattest,numcoast,nperio, ktodo
  REAL                                  :: spval,batg_min,batg_max

  INTEGER                               ::  NCID_IN, NCID_OUT, XID, YID, TID
  INTEGER                               ::  IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN, IDVAR3_IN
  INTEGER                               ::  LONID, LATID
  INTEGER                               ::  STATUS
  INTEGER                               ::  DIM_msk(2)
  CHARACTER(len=140)                    ::  chfilin
  CHARACTER(len=140)                    ::  chfilout, chfilmask,chanswer, chfilmsk
  CHARACTER(len=8)                      ::  chform
  CHARACTER(len=80 )                    ::  chtreatment,chcomment,chname

  INTEGER           :: narg, iargc, jarg
  CHARACTER(LEN=80) :: cdum
  LOGICAL           :: lstop=.false.
  

  narg=iargc()
  IF ( narg == 0 ) THEN 
   PRINT *, 'USAGE batcoast -if file_in -of file_out -msk maskfile.nc -asc ascii_coast file'
   PRINT *, '       file_in file_out and maskfile are netcdf file '
   PRINT *, '       if using an ascii coast file, mask is not required/used.'
   STOP
  ENDIF

  jarg=1
  chfilin='' ; chfilout='' ; chfilmsk='' ; chfilmask=''
  DO WHILE ( jarg <= narg )
   CALL getarg(jarg,cdum) ; jarg=jarg+1
    SELECT CASE ( cdum )
      CASE ( '-if' )
        CALL getarg(jarg,chfilin)   ; jarg=jarg+1
      CASE ( '-of' ) 
        CALL getarg(jarg,chfilout)  ; jarg=jarg+1
      CASE ( '-msk' )
        CALL getarg(jarg,chfilmsk)  ; jarg=jarg+1
      CASE ( '-asc' )
        CALL getarg(jarg,chfilmask) ; jarg=jarg+1
      CASE DEFAULT 
        PRINT *, 'Error: option ',TRIM(cdum),' not available in batcoast'
        STOP
      END SELECT
   END DO
  ! Check that all required file names are given on the command line :
  lstop=.false.
  IF ( chfilin == '' ) THEN ; PRINT * ,' You must supply an input file with -if option'   ; lstop=.true. ; ENDIF
  IF ( chfilout == '' ) THEN ; PRINT * ,' You must supply an output file with -of option' ; lstop=.true. ; ENDIF
  IF ( chfilmsk == '' .AND. chfilmask == ''  ) THEN 
                            PRINT * ,' You must supply a valid mask with either -msk or -asc option' ; lstop=.true. ; ENDIF
  IF ( lstop ) STOP
  !
  !  DEFINITION of the treatment ot be applied to the bathymetry 
  !
  chtreatment = 'Gebco1 (z<200m) + smoothed etopo1 (z>300m) with mask  '
  chcomment   = 'Original mask of MERCATOR+ modif amt'
  !
  chname = 'Bathymetry'
  !
  !      opens INPUT files 
  !      ------------------
  !
  STATUS = NF90_OPEN(TRIM(chfilin), NF90_NOWRITE, NCID_IN)
  IF (STATUS .NE. NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF

  STATUS = NF90_INQ_DIMID(NCID_IN, 'x', LONID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpi)
  STATUS = NF90_INQ_DIMID(NCID_IN, 'y', LATID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpj)
  PRINT *,' reading netcdf file, dimension jpi=',jpi,' jpj=',jpj
  IF (jpi.LE.0.OR.jpj.LE.0) THEN
     PRINT *, ' dimensions read in .nc file is zero'
     STOP
  ENDIF
  !
  !   ALLOCATE ARRAYS 
  !
  ALLOCATE(batg(jpi,jpj),mcoast(jpi,jpj),batmp(jpi,jpj))

  STATUS=NF90_INQ_VARID(NCID_IN, chname   , IDVAR3_IN)
  !
  ! Reading the bathymetry in  meters
  !
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR3_IN, batg   )
  STATUS = NF90_CLOSE(NCID_IN)
  !
  !  the bathymetry must be positive for OPA model.
  !   if it is not, convert it to positive values
  batg=ABS(batg) 
  PRINT *, ' Bathy read ok, max =',MAXVAL(batg)
  PRINT *, ' Bathy read ok, min =',MINVAL(batg)

  !
  ! Reading the coastal mask
  !
  IF ( chfilmask /= '' ) THEN
     numcoast=10
     OPEN(numcoast,form='formatted',file=chfilmask,status='old')
     WRITE(chform,102) jpi
102  FORMAT ( '(',i4,'i1)')
     PRINT *,'reading coast mask with format:',chform
     DO jj=jpj,1,-1
        READ (numcoast,chform)(mcoast(ji,jj),ji=1,jpi)
     ENDDO
     CLOSE(numcoast)      
     PRINT *, ' Coast mask read ok, max =',MAXVAL(mcoast)
     PRINT *, ' Coast mask read ok, min =',MINVAL(mcoast)
  ELSE IF ( chfilmsk /= '' ) THEN
     STATUS = NF90_OPEN(TRIM(chfilmsk), NF90_NOWRITE, NCID_IN)
     STATUS = NF90_INQ_VARID(NCID_IN, 'tmask', IDVAR3_IN)
     STATUS = NF90_GET_VAR(NCID_IN,IDVAR3_IN,mcoast,start=(/1,1,1/), count=(/jpi,jpj,1/) )
     STATUS = NF90_CLOSE(NCID_IN)
     PRINT *, ' Coast mask read ok, max =',MAXVAL(mcoast)
     PRINT *, ' Coast mask read ok, min =',MINVAL(mcoast)
  ELSE
     PRINT *,' no mask files. The program shold have stopped before ...'
  ENDIF
  !
  ! Fill the bathymetry where the coastal mask indicates land
  ! and compute a value where the coastal mask indicates ocean while the
  !         bathymetry indicates land. 
  !  we will need to find a suitable depth to put in that area, using extrapolation.
  !
  spval=99999.
  kbat0=0
  kmask1=0
  kbattest=0
  !
  DO jj=1,jpj
     DO ji=1,jpi
        !                case mask is land and batg is not
        IF((mcoast(ji,jj).EQ.0).AND.batg(ji,jj).NE.0)THEN
           kbat0=kbat0+1
           batg(ji,jj)=0.0
        ELSE IF((batg(ji,jj).EQ.0.).AND.mcoast(ji,jj).NE.0) THEN
           kmask1=kmask1+1
           batg(ji,jj)=spval
        ENDIF
     ENDDO
  ENDDO
  PRINT *, ' The coastal mask has been applied. '
  PRINT *, kbat0,' ocean points were put to 0'
  PRINT *, kmask1,' points were land in first guess bathy but ocean in the mask'
  !
  !   We will check further these later points. allocates arrays 
  ALLOCATE ( indi(kmask1) , indj (kmask1) )
  jk=1
  DO jj=1,jpj
     DO ji=1,jpi
        IF(  batg(ji,jj).EQ.spval  ) THEN
           indi(jk) = ji
           indj(jk)= jj
           jk=jk+1
        ENDIF
     ENDDO
  ENDDO
  !
  PRINT *, ' indi and indj have been constructed, ',jk-1, ' elements'
  !      DO jk=1,kmask1
  !       print 103, indi(jk), indj(jk), batg(indi(jk), indj(jk))
  !      enddo 
  ! 103  format ( 'i:',i4,'  j:',i4,' old bathy:',f8.1)

  ktodo = 1
  kiter=1
  !                 creates an array batmp that is zero over land and points to be filled.
  batmp = batg; 
  WHERE( batmp <= 0.  ) batmp = spval
  DO WHILE (ktodo .NE.0 .AND. kiter < 500)
     CALL ntrfil(batmp,jpi,jpj,spval)
     kbattest=0
     ktodo = 0
     DO jj=1,jpj
        DO ji=1,jpi
           IF(  batg(ji,jj).EQ.spval  ) THEN
              IF (batmp (ji,jj) > 0. .AND. batmp(ji,jj) .NE. spval ) THEN
                 batg(ji,jj)=batmp(ji,jj)
                 kbattest=kbattest+1
              ELSE
                 ktodo = ktodo+1
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     PRINT * , kbattest,' new ocean points extrapolated from neighbors'
     PRINT * , ktodo,   ' points remain to do at iteration:',kiter
     kiter=kiter+1
  ENDDO
  !
  !  final verification 
  !   
  kiter = 0
  DO jj=1,jpj
     DO ji=1,jpi
        !                case mask is land and batg is not
        IF((mcoast(ji,jj).EQ.0).AND.batg(ji,jj).NE.0)THEN
           PRINT *, ' mismatch at ji,jj',ji,jj,' mcoast = 0 but not batg '
           kiter = kiter+1
        ELSE IF((batg(ji,jj).EQ.0.).AND.mcoast(ji,jj).NE.0) THEN
           PRINT *, ' mismatch at ji,jj',ji,jj,' mcoast = 1  but batg = 0 '
           kiter = kiter+1
        ENDIF
     ENDDO
  ENDDO
  IF (kiter .NE. 0) STOP
  !
  !
  !   Apply lateral boundary conditions
  !   calling lbc ( ptab,kjpi,kjpj, kjpk, ktype,ksgn,nperio)
  !   The bathymetry file is at "T" points, no sign change: ktype = 1, ksgn = 1.
  !   it is a 2D array so kjpk = 1
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
  CALL lbc ( batg,jpi,jpj, 1, 1 ,1, nperio) 
  !
  !    printout new bathy where it has been extrapolated (points designated as ocean
  !    in the mask that were land in batg).
  !
  DO jk=1,kmask1
     PRINT 101, indi(jk), indj(jk), batg(indi(jk), indj(jk))
  ENDDO
101 FORMAT ( 'i:',i4,'  j:',i4,' new bathy:',f8.1)

  !  CREATE A NEW bathymetry .nc file using attributes and coodinates of the old one
  !
  CALL rewrite_bathy_nc (batg,jpi,jpj,chfilin,chfilout,chtreatment,chcomment)

  STOP

END PROGRAM batcoast
!!C---------------------------------------------------------------------

SUBROUTINE NTRFIL(PTIN,KIM,KJM,PVAL)
  !!C---------------------------------------------------------------------
  !!C
  !!C                       ROUTINE NTRFIL
  !!C                     *******************
  !!C
  !!C  PURPOSE :
  !!C  ---------
  !!C     FILL THE ARRAY PTIN WHERE THE VALUE ARE UNAVAILABLE.
  !!C                    i.e. WHERE, THERE IS THE VALUE PVAL.
  !!   METHOD :
  !!   -------
  !!     WE CALCULATE VALUES FOR ONE ROW MORE WHICH THE 8 VALUES AROUND.
  !!
  !!   INPUT :
  !!   -----
  !!      ARGUMENT
  !!              PTIN            : INPUT ARRAY
  !!              KIM,KJM         : DIMENSIONS
  !!              PVAL                : VALUE FOR UNAIVALABLE DATA
  !!
  !!   OUTPUT :  PTIN
  !!   ------
  !!
  !!   WORKSPACE :
  !!   ---------
  !!      LOCAL
  !!              ZWORK           : WORK ARRAY
  !!
  !
  IMPLICIT NONE
  INTEGER ji,jj,kim,kjm
  REAL PTIN(KIM,KJM)
  REAL, DIMENSION (:,:), ALLOCATABLE :: zwork, zmask
  REAL PVAL,znbpt
  !
  !  REMARK: dynamic allocation to avoid the stack limit on automatic arrays
  !  with ifort
  !  (stack limit 8 Mbytes by default on kaaba - here two arrays of 12 Mbytes)
  ALLOCATE( ZWORK(0:(kim+1),0:(kjm+1)),ZMASK(0:(kim+1),0:(kjm+1)))

  !
  !
  ! 1. INITIALIZATION
  ! -----------------
  !
  ZWORK=0.
  ZMASK=0.
  !
  DO JJ=1,KJM
     DO JI=1,KIM
        ZWORK(JI,JJ)=PTIN(JI,JJ)
        IF(ABS(PTIN(JI,JJ)-PVAL).GT.(1.E-5)) ZMASK(JI,JJ)=1.
     END DO
  END DO
  !
  !
  ! 2. FILLING
  ! ----------
  !
  DO JJ=1,KJM
     DO JI=1,KIM
        !
        ZNBPT=ZMASK(JI-1,JJ+1)+ZMASK(JI,JJ+1)+ZMASK(JI+1,JJ+1)+  &
             ZMASK(JI-1,JJ  )               +ZMASK(JI+1,JJ  )+  &
             ZMASK(JI-1,JJ-1)+ZMASK(JI,JJ-1)+ZMASK(JI+1,JJ-1)

        PTIN(JI,JJ)=ZMASK(JI,JJ)*ZWORK(JI,JJ)+(1.-ZMASK(JI,JJ))* &
             ((1.-MIN(1.,ZNBPT))*PVAL+MIN(1.,ZNBPT)*      &
             (ZWORK(JI-1,JJ+1)*ZMASK(JI-1,JJ+1)+          &
             ZWORK(JI  ,JJ+1)*ZMASK(JI  ,JJ+1)+          &
             ZWORK(JI+1,JJ+1)*ZMASK(JI+1,JJ+1)+          &
             ZWORK(JI-1,JJ  )*ZMASK(JI-1,JJ  )+          &
             ZWORK(JI+1,JJ  )*ZMASK(JI+1,JJ  )+          &
             ZWORK(JI-1,JJ-1)*ZMASK(JI-1,JJ-1)+          &
             ZWORK(JI  ,JJ-1)*ZMASK(JI  ,JJ-1)+          &
             ZWORK(JI+1,JJ-1)*ZMASK(JI+1,JJ-1))/         &
             MAX(1.,ZNBPT))
     END DO
  END DO
  !
  !
  !
  RETURN
END SUBROUTINE NTRFIL
