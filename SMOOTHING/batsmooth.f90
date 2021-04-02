PROGRAM batsmooth
  !-------------------------------------------------------------------
  ! read an OPA Bathymetry file and smoothes it 
  ! with a shapiro filter (METHLIS=1) or with a hanning filter 
  ! of order one in two directions (METHLIS=2)
  ! the number of times the smoothing is applied is determined by  NBRPASSE
  !
  !  history: 
  !  - Anne MArie Treguier, CLIPPER project, using shapiro filter routine 
  !           from the SPEM code (1997)
  !    Note this routine could filter better near the periodic boundaries 
  !    ( here does not filter the boundaries)
  !  - MERCATOR project, fortran 90, add Hanning filter
  !  - Anne marie Treguier, process bathy in meters only, 
  !    netcdf I/O.
  !  
  !-------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 18:28:10 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 230 $
   !!  $Id: batsmooth.f90 230 2009-04-29 16:28:10Z forge $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER                            :: methlis, nbrpass, iordre
  REAL                               :: weight
  REAL,   DIMENSION(:,:),ALLOCATABLE :: bathy, bathylis
  CHARACTER(len=140)                 :: chfilin,chanswer,chfilout,chfilroot
  CHARACTER(len=80 )                 :: chtreatment,chcomment

  INTEGER                            :: ncid_in ,idvar3_in, idlon, idlat
  INTEGER                            :: istatus, idx

  INTEGER                            :: ji,jj,ipass,jpi,jpj,nperio, newper

  INTEGER :: narg, iargc
  CHARACTER(LEN=80) :: cdum


  narg=iargc()
  IF ( narg /= 2 ) THEN
    PRINT *,' USAGE: batsmooth input_file nbrpass'
    PRINT *,'   Specify the number of passes for shapiro filter'
    PRINT *,'   ORCA025 uses 2 '
    PRINT *,'   Output file name will be build acording to settings'
    STOP
  ENDIF
  !  nbrpass=16 ! make test with 0,1,2,3,4
  CALL getarg(1,chfilin) 
  CALL getarg(2,cdum) ; READ(cdum,*) nbrpass
  chfilroot=chfilin
  idx = INDEX(TRIM(chfilin),'.nc') 
  IF ( idx /= 0 ) chfilroot=chfilin(1:idx-1)
  
  !-------------------------------------------------------------------
  ! 0 - DEFINITION of the treatment ot be applied to the bathymetry 
  ! ------------------------------------------------------------------
  !           Median filter
  !  methlis=3
  !  nbrpass=20
  !           Shapiro filter
  methlis=1
  !
  !      iordre must be even for shapiro filter 
  !
  iordre=2
  weight=0.6
  !
  chtreatment = 'Bathy after smoothing'
  IF (methlis == 1) THEN
     WRITE(chcomment,'("Shapiro filter order ",i2,", weight:",f4.1,", nb passes:",i2)' ) iordre, weight,nbrpass
  ELSEIF (methlis == 2) THEN
     WRITE(chcomment,'("Hanning filter, nb passes:",i2)' ),nbrpass
  ELSEIF (methlis == 3) THEN
     WRITE(chcomment,'("Median filter, nb passes:",i2)' ),nbrpass
  ELSE
     PRINT *, 'smoothing method unknown, we stop'
     STOP
  ENDIF
  !
  ! ------------------------------------------------------------------
  ! 1 -  open bathymetry file
  ! ------------------------------------------------------------------
  !                            default value for file name: 
  IF (methlis == 1) THEN
     !  case with bug in indian ocean
     WRITE(chfilout,'(a,"_shapiro_",i2.2,"p_",f3.1,"w.nc")' ) TRIM(chfilroot),nbrpass,weight
  ELSE IF (methlis == 3) THEN
     WRITE(chfilout,'(a,"_median_",i2.2,"p.nc")' ) chfilroot, nbrpass
  ENDIF
   PRINT *,'SMOOTHED FILE IS ', TRIM(chfilout)
  !                             read actual bathymetry filename
  PRINT *,' file name for bathy+coastline, ''y'' for default:',TRIM(chfilin)

  !  READ dimensions in netcdf file
  !
  istatus = NF90_OPEN(TRIM(chfilin), NF90_NOWRITE, ncid_in)
  IF (istatus .NE. NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(istatus)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilin)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  istatus = NF90_INQ_DIMID(ncid_in, 'x', idlon)      
  istatus = NF90_INQUIRE_DIMENSION(ncid_in, idlon, len=jpi)
  istatus = NF90_INQ_DIMID(ncid_in, 'y', idlat)      
  istatus = NF90_INQUIRE_DIMENSION(ncid_in, idlat, len=jpj)
  PRINT *,' reading netcdf file, dimension jpi=',jpi,' jpj=',jpj
  IF (jpi <= 0 .OR. jpj <= 0) THEN
     PRINT *, ' dimensions read in .nc file is zero'
     STOP
  ENDIF
  istatus=NF90_INQ_VARID(ncid_in, 'Bathymetry'   , idvar3_in)
  !
  ! Allocate arrays
  ! ------------------------------------
  ALLOCATE (bathy(jpi,jpj)  ,  bathylis(jpi,jpj) )
  !
  ! ------------------------------------------------------------------
  ! 2 -  read bathymetry file
  ! ------------------------------------------------------------------
  !   Reading the bathymetry in  meters
  !     note that although ncdump says it is batg(y,x), 
  !     we recover batg(x,y) as expected.
  !
  istatus = NF90_GET_VAR(ncid_in, idvar3_in, bathy  )
  newper=2  ! default orca value

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
  IF (nperio==1) THEN
     PRINT *, 'how many overlap (2 for ORCA type configuration) '
     READ (5,*) newper
  END IF
  !
  !  re-apply initial mask 
  !          the bathymetry must be positive for OPA model.
  !          if it is not, convert it to positive values
  bathylis = ABS(bathy)

  PRINT *, ' Bathy read ok, max =',MAXVAL(bathylis)
  PRINT *, ' Bathy read ok, min =',MINVAL(bathylis)
  !

  ! ------------------------------------------------------------------
  !  3 -   Bathymetry  smoothing
  ! -------------------------------------------------------------------
  PRINT *, 'Bathymetry smoothing'

  IF (methlis == 1) THEN
     DO ipass = 1, nbrpass
        PRINT *,'Smoothing with shapiro filter, pass number :',ipass
        CALL shapiro(bathylis, jpi, jpj, iordre, weight, newper)
     END DO
  ELSEIF   (methlis == 2) THEN
     DO ipass = 1, nbrpass
        PRINT*,'Lissage de la bathymetrie par hanning, passage :',ipass
        CALL hanning(bathylis, 1, jpi, 1, jpj)
     END DO
  ELSEIF   (methlis == 3) THEN
     CALL median_filt(bathylis,jpi,jpj,nbrpass)
  ENDIF
  ! -------------------------------------------------------------------
  !
  !   Apply lateral boundary conditions
  !   calling lbc ( ptab,kjpi,kjpj, kjpk, ktype,ksgn,nperio)
  !   The bathymetry file is at "T" points, no sign change: ktype = 1, ksgn = 1.
  !   it is a 2D array so kjpk = 1
  ! -------------------------------------------------------------------
  !
  CALL lbc ( bathylis,jpi,jpj, 1, 1 ,1, nperio)   
  !
  !  re-apply initial mask 
  !
  WHERE(bathy == 0) bathylis = 0 

  ! -------------------------------------------------------------------
  !  Write new bathy file
  ! -------------------------------
  CALL rewrite_bathy_nc(bathylis, jpi,jpj,chfilin,chfilout, chtreatment,chcomment)

  STOP

  !*********************************************************

CONTAINS
  !*********************************************************
  SUBROUTINE shapiro(bat, jpi, jpj, N, ww, ewper)
    !*********************************************************

    !   uniform weight=ww
    !  with zeroth order 
    ! This amounts to an average w*(1/4), 1-(w*1/2), w*1/4    
    !  for ww= 1  , average (1/4    1/2    1/4 )
    !  for ww= 0.8, average (0.2    0.6    0.2 )
    !  for ww= 0.6, average (0.15   0.7    0.15)
    !  for ww= 0.4, average (0.1    0.8    0.1 )
    IMPLICIT NONE

    INTEGER, INTENT(in) :: jpi, jpj, N, ewper
    REAL,    INTENT(in) :: ww

    REAL, DIMENSION(1:jpi,1:jpj), INTENT(inout) :: bat
    REAL, DIMENSION(:,:),ALLOCATABLE                :: tmp, tmp2
    REAL, DIMENSION(1:jpi,1:jpj) :: bat1
    ! Variables locales
    ! -----------------
    INTEGER :: d
    INTEGER :: i,j,k, ip1, im1

    ALLOCATE ( tmp(jpi,jpj), tmp2(jpi,jpj))

    IF (MOD(N,2).NE.0) THEN
       PRINT *,'N must be even in the shapiro filter'
       CALL EXIT(1)
    END IF

    !  Do the first y pass to initialize the temporary array
    bat1=bat
    DO j=2,jpj-1
       DO i=  1,jpi
          !IF (bat(i,j+1) == 0) bat1(i,j+1)=bat(i,j) !PIERRE
          !IF (bat(i,j-1) == 0) bat1(i,j+1)=bat(i,j) !PIERRE
          tmp(i,j) = 0.25 * (bat1(i,j-1) + bat1(i,j+1) - 2*bat1(i,j))
       END DO
    END DO

    !  Other passes in the y direction.
    !   for order 2, two passes, with d=2 and 1
    DO k=4,N,2
       d = k/2
       DO j=1+d,jpj-d
          DO i=1,jpi
             tmp2(i,j) = - 0.25 * (tmp(i,j-1) + tmp(i,j+1) - 2*tmp(i,j))
          END DO
       END DO
       DO  j=1+d,jpj-d
          DO i=1,jpi
             tmp(i,j) = tmp2(i,j)
          END DO
       END DO
    END DO

    !  Add the changes to u
    DO  j=2,jpj-1
       DO  i=1,jpi
          bat(i,j) = bat(i,j) + ww*tmp(i,j)
       END DO
    END DO

    !  Initialize tmp to filter in the x direction.
    PRINT *, 'Initialize tmp to filter in the x direction.'
    DO j=1,jpj
       !DO i=2,jpi-1
       DO i=1,jpi
          im1=i-1
          ip1=i+1
          IF (i==1) im1=jpi-ewper
          IF (i==jpi) ip1=1+ewper
          IF ( im1 < 0 ) print *, i,j
          IF ( ip1 < 0 ) print *, i,j
          tmp(i,j) = 0.25 * (bat(im1,j) + bat(ip1,j) - 2*bat(i,j))
       END DO
    END DO

    !  Other x passes
    DO k=4,N,2
       d = k/2
       DO j=1,jpj
          DO i=1+d,jpi-d
             im1=i-1
             ip1=i+1
             IF (i==1) im1=jpi-ewper
             IF (i==1) ip1=1+ewper
             tmp2(i,j) = - 0.25 * (tmp(im1,j) + tmp(ip1,j) - 2*tmp(i,j))
          END DO
          DO i=1+d,jpi-d
             tmp(i,j) = tmp2(i,j)
          END DO
       END DO
    END DO

    !  Add changes to u
    DO j=1,jpj
       !DO i=2,jpi-1
       DO i=1,jpi
          bat(i,j) = bat(i,j) + ww*tmp(i,j)
       END DO
    END DO
    RETURN
  END SUBROUTINE shapiro

  !*********************************************************

  SUBROUTINE hanning(bathy, jpimin, jpimax, jpjmin, jpjmax)

    IMPLICIT NONE
    INTEGER, INTENT(in) :: jpimin, jpjmin,  jpimax, jpjmax

    REAL, DIMENSION(jpimin:jpimax,jpjmin:jpjmax), INTENT(inout) :: bathy
    REAL, DIMENSION( jpimin:jpimax,jpjmin:jpjmax)               :: tmp

    INTEGER :: ii, jj

    tmp = bathy
    DO ii= jpimin+1, jpimax-1
       DO jj= jpjmin+1, jpjmax-1
          IF(bathy(ii+1,jj)<=0 .OR. bathy(ii-1,jj)<=0 .OR. & 
               bathy(ii,jj+1)<=0 .OR. bathy(ii,jj-1)<=0) THEN
             CYCLE
          ELSE
             tmp(ii,jj)=0.5*bathy(ii,jj)+0.125*bathy(ii+1,jj) &
                  +0.125*bathy(ii-1,jj)+0.125*bathy(ii,jj+1)+0.125*bathy(ii,jj-1)
          ENDIF
       END DO
    END DO
    bathy = tmp

    RETURN
  END SUBROUTINE hanning

  !*********************************************************

  SUBROUTINE median_filt(batin,jpi,jpj,niter)
    !*********************************************************
    ! Median filter based on Julio Candela's matlab routine.
    !
    !*********************************************************



    IMPLICIT NONE
    INTEGER, PARAMETER 		:: jp2 = 2, jp5= 5, jp17 = jp5+4*jp2+4
    INTEGER, INTENT(in) 	:: jpi, jpj,niter
    REAL, DIMENSION(jpi,jpj), INTENT(inout)  :: batin

    REAL, DIMENSION(jpi,jpj)              ::    batout
    REAL, DIMENSION(jp17)                 :: wk
    INTEGER, DIMENSION(jp17)              :: ind
    INTEGER				  :: ji,jj,jk,iter,kflag,jp9,ipoint

    !   integer division
    jp9 = (jp17+1)/2
    PRINT *, '        jp9=',jp9
    batout = batin
    DO iter=1,niter
       PRINT *, ' Median filter pass number:',iter
       PRINT *, '-----------------------------------'
       ipoint = 0
       DO jj=2,jpj-1
          DO ji=2,jpi-1
             DO jk = 1,jp5
                wk(jk) = batout(ji,jj)
             ENDDO
             DO jk=1,jp2
                wk(jp5+jk)       = batin(ji-1,jj  )
                wk(jp5+jp2+jk)   = batin(ji  ,jj-1)
                wk(jp5+2*jp2+jk) = batin(ji+1,jj  )
                wk(jp5+3*jp2+jk) = batin(ji  ,jj+1)
             ENDDO
             wk(jp5+4*jp2+1) = batin(ji-1,jj-1)
             wk(jp5+4*jp2+2) = batin(ji+1,jj-1)
             wk(jp5+4*jp2+3) = batin(ji-1,jj+1)
             wk(jp5+4*jp2+4) = batin(ji+1,jj+1)
             !   need a sorting routine here
             !          SUBROUTINE SSORT (X, IY, N, KFLAG)
             CALL ssort(wk,ind,jp17,kflag)
             !
             IF (wk(jp9).NE.batout(ji,jj)) ipoint = ipoint+1
             IF (MOD(ji,400) == 0 .AND. MOD(jj,200)==0) THEN
                PRINT '("ji,jj=",i5,1x,i5," bathy:",f6.1," modified :",f6.1)', ji,jj,batout(ji,jj),wk(jp9)
                PRINT '(17(1x,f6.1))', wk
             ENDIF
             batout(ji,jj) = wk(jp9 )
          ENDDO
       ENDDO
       batin = batout
       !                       end of loop on iterations
       PRINT *, ' Number of points changed for this loop:',ipoint
       PRINT *, '------------------------------------------------'
    ENDDO
    RETURN

  END SUBROUTINE median_filt

END PROGRAM batsmooth

