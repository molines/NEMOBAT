!-------------------------------------------------------------------
PROGRAM batinterp
  !-------------------------------------------------------------------
  !  interpolates a bathy onto an ORCA grid.
  !  The method is by averaging all the original bathy points that fall into 
  !  an orca grid cell. This is very good for relatively coarse grids 
  !  (like ORCA025).
  !---------------------------------------------------------
  !   Two averaging methods are provided: 
  !     l_interp = 0 : arithmetic average;
  !     l_interp = 1 : median average.
  !--------------------------------------------------------
  !
  !  if defined etopo2:
  !  original bathy is ETOPO2+bedmap.
  !  input file (from elisabeth remy, MERCATOR): netcdf, 
  !  lon varies from  [ -180:1./30:180-1./30 : 10800 points 
  !  lat varies from  [ 90:-1./30:-90+1./30 ]: 5400 points 
  !  the  variable is depth and it is an integer
  !  (note that for obscure reasons one row and column of grid points have 
  !  been left out the original etopo2 file).
  !  
  !   if defined gebco1 (default):
  !  interpolates gebco 1mn bathy onto ORCA grid, from 88N to 88S only.
  !  the gebco file is netdcf (grd). 
  !  It has latitude from -180 to 180 (jpibat = 21601)
  !  and latitude from  -88, 88. (jpjbat = 10561)
  !  the array (z) is a short integer and one dimensional.
  !
  !
  !  Characteristics of the ORCA025 grid: 
  !  note that near the poles of the grid the longitude can vary considerably
  !  between two adjacent grid points. 
  !
  !   At grid point jj = 1019, there is a difference of 70 degrees lon. 
  !   At grid point jj = 1015, there is a difference of 14 degrees lon.
  !   With Gebco there is a max of 1219 gebco points in each ORCA cell (more with 
  !   ETOPO2 because we interpolate all the way to the pole)
  !   The cells having more points are likely to fall on land anyway 
  !   (the grid poles are on land)
  !-------------------------------------------------------------------
  !  History: derived by Elisabeth Remy from a version of opabat provided by 
  !  Maurice Imbard (LODYC) in 1995.
  !  Modified by A.M. Treguier and  Olivier Le Galloudeg (Mercator)
  !  June 2005:
  !  - support of netcdf files; 
  !  - support of both GEBCO and ETOPO2,
  !  - follow the IDL code developed by G. Madec and students (bathymean.pro)
  !    (much simpler)
  !---------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-07-24 08:23:41 +0200 (Fri, 24 Jul 2009) $
   !!  $Rev: 235 $
   !!  $Id: batinterpnew.F90 235 2009-07-24 06:23:41Z forge $
   !! -------------------------------------------------------------------------
  USE netcdf
  IMPLICIT NONE
  INTEGER                            :: jpi,jpj,jpibat,jpjbat,jpijbat
  !
  REAL                               ::  bat_min,bat_max
  REAL,     DIMENSION(:,:) ,ALLOCATABLE  :: bato
  INTEGER  ,DIMENSION(:,:) ,ALLOCATABLE  :: itmp

  REAL*8,   DIMENSION(:,:,:,:) ,ALLOCATABLE  :: coord
  REAL  ,   DIMENSION(:,:) ,ALLOCATABLE  :: glamf, gphif,glamt,gphit
  INTEGER,  DIMENSION(:,:) ,ALLOCATABLE  :: mpoi,moce
  CHARACTER(len=140)                 :: chfilout
  CHARACTER(len=140)                 :: chfilbatin, chfilgrid
  INTEGER                            :: NCID_IN, NCID_OUT, XID, YID
  INTEGER                            :: IDVAR1, IDVAR2, IDVAR3, IDVAR1_IN, IDVAR2_IN
  INTEGER                            :: IDVAR4, IDVAR5, IDVAR6, IDVAR7
  INTEGER                            :: LONID, LATID
  INTEGER                            :: STATUS
  INTEGER, DIMENSION (2)             :: DIM_MSK
  INTEGER                            :: ji,jj, jk
  INTEGER                            :: l_interp, kperio
!!!!!!!!!!!!!!!!!!!! ajout OPABAT !!!!!!!!!!!!!!!!!!!!
  REAL,DIMENSION(:,:)   ,ALLOCATABLE  :: vardep,zdep
  REAL,DIMENSION(:,:)   ,ALLOCATABLE  :: glamhr,gphihr,pb
  REAL,DIMENSION(:,:)   ,ALLOCATABLE  :: tabtestsw,tabtestse,tabtestne,tabtestnw
  REAL,DIMENSION(:)     ,ALLOCATABLE  :: flo, fla, vardep1d
  INTEGER,DIMENSION(:,:) ,ALLOCATABLE :: mask_box,mask_oce
  INTEGER,DIMENSION(:)   ,ALLOCATABLE :: istmp
  REAL lonmin, lonmax, latmin, latmax, latsw, latse, latne, latnw, gphifmax
  REAL slatmin, slatmax              ! min and max latitudes in input file 
  REAL resolution                    ! resolution in min for the input grid
  REAL nbpd                          ! number of input point per deg
  INTEGER jimin,jimax,nxhr,nyhr, i_calc
  INTEGER :: i,j,jjmin,jjmax
  INTEGER knumax, nperio

  !-------------------------------------------------------------------
  !  Choose averaging method: l_interp=1: median average, 
  !  l_interp=0 : arithmetic average.
  !-------------------------------------------------------------------
  l_interp=1 
  !-------------------------------------------------------------------
  !   sets input/ output file names
  !-------------------------------------------------------------------
  chfilbatin = 'bathy_in.nc'
  chfilout = 'bathy_out.nc'
  !-------------------------------------------------------------------
  ! 0 - Read grid file 
  ! ------------------------------------------------------------------
  chfilgrid = 'coordinates.nc'
  PRINT *, 'Opening  the input file:' ,TRIM(chfilgrid)
  STATUS = NF90_OPEN(TRIM(chfilgrid), NF90_NOWRITE, NCID_IN)
  IF (STATUS /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilgrid)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF

  STATUS = NF90_INQ_DIMID(NCID_IN, 'x', LONID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpi)
  STATUS = NF90_INQ_DIMID(NCID_IN, 'y', LATID)      
  STATUS = NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpj)
  PRINT *,'--    Reading input file, dimension jpi=',jpi,' jpj=',jpj
  !
  !   for a global grid we need to construct a criterion 
  !
  !   ALLOCATE arrays for longitude and latitude 
  ! 
  ALLOCATE(glamf(jpi,jpj),gphif(jpi,jpj),glamt(jpi,jpj),gphit(jpi,jpj))
  ALLOCATE(zdep(jpi,jpj))
  !
  !            read longitude and latitude 
  !
  STATUS = NF90_INQ_VARID(NCID_IN, 'glamf', IDVAR1_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR1_IN, glamf(:,:), start=(/1,1/), count=(/jpi,jpj/) )
  STATUS = NF90_INQ_VARID(NCID_IN, 'gphif', IDVAR1_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR1_IN, gphif(:,:), start=(/1,1/), count=(/jpi,jpj/) )
  STATUS = NF90_INQ_VARID(NCID_IN, 'glamt', IDVAR2_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR2_IN, glamt(:,:), start=(/1,1/), count=(/jpi,jpj/) )
  STATUS = NF90_INQ_VARID(NCID_IN, 'gphit', IDVAR2_IN)
  STATUS = NF90_GET_VAR(NCID_IN, IDVAR2_IN, gphit(:,:), start=(/1,1/), count=(/jpi,jpj/) )
  bat_max= MAXVAL(glamf)
  bat_min= MINVAL(glamf)
  PRINT *, ' read  glamf  OK max=',bat_max,' min =', bat_min 
  bat_max= MAXVAL(glamt)
  bat_min= MINVAL(glamt)
  PRINT *, ' read  glamt  OK max=',bat_max,' min =', bat_min 
  gphifmax= MAXVAL(gphif)
  bat_min= MINVAL(gphif)
  PRINT *, ' read  gphif  OK max=',gphifmax,' min =', bat_min 
  bat_max= MAXVAL(gphit)
  bat_min= MINVAL(gphit)
  PRINT *, ' read  gphit  OK max=',bat_max,' min =', bat_min 

  STATUS=NF90_CLOSE(NCID_IN)

  !--------------------------------------------------------------
  ! Boundary conditions  
  !--------------------------------------------------------------
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
  !--------------------------------------------------------------
  !
  !    Allocate arrays for bathymetry 
  ALLOCATE ( pb(jpi,jpj), moce(jpi,jpj), mpoi(jpi,jpj))

  !-------------------------------------------------------------------
  ! 1 - Read original bathy file 
  ! ------------------------------------------------------------------
  PRINT *, 'Opening  the input file:' ,TRIM(chfilbatin)
  STATUS = NF90_OPEN(TRIM(chfilbatin), NF90_NOWRITE, NCID_IN)
  IF (STATUS /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(STATUS)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(chfilbatin)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  ! We suppose that all input files are netcdf map. If grd, preprocess with grd2nc.
  ! etopo1 and gebco
  CALL hdlerr(NF90_INQ_DIMID(NCID_IN, 'lon', LONID))      
  CALL hdlerr( NF90_INQUIRE_DIMENSION(NCID_IN, LONID, len=jpibat))
  CALL hdlerr(NF90_INQ_DIMID(NCID_IN, 'lat', LATID))      
  CALL hdlerr( NF90_INQUIRE_DIMENSION(NCID_IN, LATID, len=jpjbat))
  ALLOCATE (itmp(jpibat,jpjbat))
  PRINT *,'--   Reading etopox file, dimension ',jpibat,jpjbat

  

  CALL hdlerr(NF90_INQ_VARID(NCID_IN,'z',IDVAR1_IN))
  CALL hdlerr(NF90_GET_VAR(NCID_IN,IDVAR1_IN,itmp))
  !                                 ALLOCATE  bato array 
  ! 
  ALLOCATE(bato(jpibat,jpjbat),flo(jpibat),fla(jpjbat))
  CALL hdlerr(NF90_INQ_VARID(NCID_IN, 'lon', LONID))
  CALL hdlerr(NF90_GET_VAR(NCID_IN, LONID, flo))

  CALL hdlerr(NF90_INQ_VARID(NCID_IN, 'lat', LATID))
  CALL hdlerr(NF90_GET_VAR(NCID_IN, LATID, fla))

  ! determine the resolution in minutes
  resolution=nint( (flo(20)-flo(19) )*60 )
  nbpd=nint(60/resolution)
  PRINT *,'-- Input file resolution is ', INT(resolution),  ' minutes'

  bato = itmp

  slatmin = -99               !  99 values for a global file
  slatmax = 99
  !   calculate max and min value of input array bathy 
  ! 
  bat_max= MAXVAL(itmp)
  bat_min= MINVAL(itmp)
  PRINT *, ' read original bathy OK max=',bat_max,' min =', bat_min
  PRINT *, ' read original lon   OK max=',MAXVAL(flo),' min =',MINVAL(flo) 
  PRINT *, ' read original lat   OK max=',MAXVAL(fla),' min =',MINVAL(fla) 
  ! test 
  slatmin=MINVAL(fla)
  slatmax=MAXVAL(fla)

  DEALLOCATE (itmp)
  CALL hdlerr(NF90_CLOSE(NCID_IN))

  !-------------------------------------------------------------------------
  ! 2 -  INTERPOLATION    : LOOP OVER ORCA GRID  
  !-------------------------------------------------------------------------
  knumax = 0
  DO jj=2,jpj-1
     IF (MOD(jj,100) == 0) THEN
        PRINT *, jj,'/',jpj-1
     END IF
     DO ji=2,jpi-1

        !     Defining grid parameters for T-point (ji,jj)
        lonmin=MIN(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))  
        lonmax=MAX(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))
        latmin=MIN(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))
        latmax=MAX(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))

        jimin=FLOOR(nbpd*(lonmin+180))+1
        jimax=FLOOR(nbpd*(lonmax+180))+2
        jjmin=FLOOR(nbpd*(latmin+90-1./nbpd))+1
        jjmax=FLOOR(nbpd*(latmax+90-1./nbpd))+2
        !
        !   first verifies that the ORCA cell is within the domain of the original bathymetry 
        !   if not, moce is set to a value of -1 as an indicator. 

        IF (latmin < slatmin .OR. latmax >  slatmax ) THEN
           moce(ji,jj)=-1
           mpoi(ji,jj)=0
           zdep(ji,jj) = 0.
        ELSE

           !    THERE ARE THREE CASES:
           !    - grid point near north pole
           !    - grid point near dateline
           !    - normal case 
           !  NOTE that the following test has been modified from abs(lonmax-lonmin) .GE. 180
           !       to 120 only in order to catch the North pole point in ORCA025.
           !
           IF ( ABS(lonmax-lonmin) >=  120) THEN 
              !    -------------- first north pole case. 
              !                   just take the northernmost value of etopo2. 
              IF (latmax >=  gphifmax ) THEN
                 zdep(ji,jj) = bato(1,jpjbat)
                 i_calc = 0
              ELSE
                 !     ------------------ dateline case
                 !     caution: cyclic east west at lon=-180 or 180 in low resolution grid   
                 !     need to fetch data from first and last i-index of bathy file. 
                 !
                 !     If the box is not aligned along a lon,lat line, then we need to 
                 !     redefine lonmin and latmin. Assume there are two corners on each side 
                 !     of the dateline.
                 !     lonmax = the min of the two positive values  (close to 180   = JPIBAT)
                 !     lonmin is the max of the two negative values (close to -180  = 1)
                 IF ( glamf(ji-1,jj-1) > 0.AND.glamf(ji-1,jj) > 0) & 
                      lonmax=MIN(glamf(ji-1,jj-1), glamf(ji-1,jj) )  
                 IF (glamf(ji,jj-1) < 0 .AND. glamf(ji,jj) < 0)  &
                      lonmin=MAX( glamf(ji,jj-1),glamf(ji,jj) )
                 i_calc = 1

#if defined etopo2
                 jimin=FLOOR(nbpd*(lonmin+180))+1
                 jimax=FLOOR(nbpd*(lonmax+180))+2
#else
                 jimin=FLOOR(nbpd*(lonmin+180))+2
                 jimax=FLOOR(nbpd*(lonmax+180))+3
#endif

                 jimax = MIN(jimax, jpibat)
                 jjmax = MIN(jjmax, jpjbat)
                 nxhr = jpibat-jimax+1 + jimin+1
                 nyhr=jjmax-jjmin+1

                 ALLOCATE(vardep(nxhr,nyhr))
                 ! pastes the two extremities in the correct order
                 vardep(1:jpibat-jimax+1,:)= bato(jimax:jpibat,jjmin:jjmax)
                 vardep(jpibat-jimax+2:nxhr,:)=bato(1:jimin+1 ,jjmin:jjmax)

                 ALLOCATE(glamhr(nxhr,nyhr),gphihr(nxhr,nyhr))
                 DO j=1,nyhr
                    glamhr(1:jpibat-jimax+1,j)=flo(jimax:jpibat)-360
                    glamhr(jpibat-jimax+2:nxhr,j)=flo(1:jimin+1)
                 ENDDO
                 DO i=1,nxhr
                    gphihr(i,:)=fla(jjmin:jjmax)
                 ENDDO

                 latsw=glamf(ji-1,jj-1)
                 latse=glamf(ji  ,jj-1)
                 latne=glamf(ji  ,jj  )
                 latnw=glamf(ji-1,jj  )
                 IF( latsw >=  0 ) latsw =latsw-360
                 IF( latse >=  0 ) latse =latse-360
                 IF( latne >=  0 ) latne =latne-360
                 IF( latnw >=  0 ) latnw =latnw-360
                 ALLOCATE(tabtestSW(nxhr,nyhr))
                 ALLOCATE(tabtestSE(nxhr,nyhr))
                 ALLOCATE(tabtestNE(nxhr,nyhr))
                 ALLOCATE(tabtestNW(nxhr,nyhr))

                 tabtestSW=(latse-latsw)*(gphihr-gphif(ji-1,jj-1)) - (gphif(ji  ,jj-1)-gphif(ji-1,jj-1))*(glamhr-latsw)
                 tabtestSE=(latne-latse)*(gphihr-gphif(ji  ,jj-1)) - (gphif(ji  ,jj  )-gphif(ji  ,jj-1))*(glamhr-latse)
                 tabtestNE=(latnw-latne)*(gphihr-gphif(ji  ,jj  )) - (gphif(ji-1,jj  )-gphif(ji  ,jj  ))*(glamhr-latne)
                 tabtestNW=(latsw-latnw)*(gphihr-gphif(ji-1,jj  )) - (gphif(ji-1,jj-1)-gphif(ji-1,jj  ))*(glamhr-latnw)
              ENDIF
              !
              !--------------------------------- Normal case, no problem with periodicity.
           ELSE
              ! reading the bathymetry square surrounding  the orca T-cell ji,jj
              i_calc = 1
              jimax = MIN(jimax, jpibat)
              jjmax = MIN(jjmax, jpjbat)
              nxhr=jimax-jimin+1
              nyhr=jjmax-jjmin+1
              ALLOCATE(vardep(nxhr,nyhr))
              vardep=bato(jimin:jimax,jjmin:jjmax)

              ALLOCATE(glamhr(nxhr,nyhr),gphihr(nxhr,nyhr))
              DO j=1,nyhr
                 glamhr(:,j)=flo(jimin:jimax)
              ENDDO
              DO i=1,nxhr
                 gphihr(i,:)=fla(jjmin:jjmax)
              ENDDO

              ! Define the 4 'equations de droite'
              ALLOCATE(tabtestSW(nxhr,nyhr))
              ALLOCATE(tabtestSE(nxhr,nyhr))
              ALLOCATE(tabtestNE(nxhr,nyhr))
              ALLOCATE(tabtestNW(nxhr,nyhr))
              tabtestSW=(glamf(ji  ,jj-1)-glamf(ji-1,jj-1))*(gphihr-gphif(ji-1,jj-1)) - (gphif(ji  ,jj-1)-gphif(ji-1,jj-1))*(glamhr-glamf(ji-1,jj-1))
              tabtestSE=(glamf(ji  ,jj  )-glamf(ji  ,jj-1))*(gphihr-gphif(ji  ,jj-1)) - (gphif(ji  ,jj  )-gphif(ji  ,jj-1))*(glamhr-glamf(ji  ,jj-1))
              tabtestNE=(glamf(ji-1,jj  )-glamf(ji  ,jj  ))*(gphihr-gphif(ji  ,jj  )) - (gphif(ji-1,jj  )-gphif(ji  ,jj  ))*(glamhr-glamf(ji  ,jj  ))
              tabtestNW=(glamf(ji-1,jj-1)-glamf(ji-1,jj  ))*(gphihr-gphif(ji-1,jj  )) - (gphif(ji-1,jj-1)-gphif(ji-1,jj  ))*(glamhr-glamf(ji-1,jj  ))
           ENDIF
           !--------------------------------------------------------------------------------------- 
           !  NOW we have the original bathymetry in the domain containing the ORCA cell.
           !  We have the 4 line equations to find which points are inside or outside the cell
           !  Find those points and do interpolation
           !--------------------------------------------------------------------------------------- 
           !                      does not do it if we are at the pole 
           IF  ( i_calc == 1 ) THEN
              ALLOCATE(mask_oce(nxhr,nyhr),mask_box(nxhr,nyhr))
              mask_oce=0
              mask_box = 0
              WHERE ( (tabtestSW >=  0.) .AND. &
                   (tabtestSE >=  0.) .AND.  &
                   (tabtestNE >=  0.) .AND.  &
                   (tabtestNW >=  0.)) mask_box = 1
              mpoi(ji,jj) = SUM(mask_box)
              ! 
              !    problem: does not find any ETOPO2 point in the ORCA grid box
              !
              IF (mpoi(ji,jj) == 0. ) THEN
                 IF ( MINVAL(vardep) >= 0 ) THEN
                    !          we are probably close to the pole of the grid, and on land. no problem
                    pb(ji,jj) = -1
                 ELSE
                    pb(ji,jj) = 1
                 ENDIF
              ENDIF   ! mpoi == 0
              !                     find total number of ocean points 
              WHERE ( (vardep <  0.) .AND. (tabtestSW >=  0.) .AND. &
                   (tabtestSE >=  0.) .AND.  &
                   (tabtestNE >=  0.) .AND.  &
                   (tabtestNW >=  0.) ) mask_oce = 1
              moce(ji,jj) = SUM(mask_oce)
              knumax = MAX(knumax,moce(ji,jj))
              !                     find bathymetry value by interpolation
              IF (moce(ji,jj) == 0 )  THEN
                 zdep(ji,jj) = 0
              ELSE
                 IF (l_interp == 0 ) THEN
                    !                Arithmetic average
                    zdep(ji,jj) = SUM (vardep*mask_oce)/float(moce(ji,jj))
                 ELSE
                    !                Median average 
                    !       put very large values on land to eliminate those points when taking the median 
                    vardep = vardep*mask_oce  - 20000*(1-mask_oce)
                    ALLOCATE(vardep1d(nxhr*nyhr),istmp(nxhr*nyhr))
                    vardep1d = RESHAPE(vardep,(/ nxhr*nyhr /) )
                    CALL ssort(vardep1d,istmp,nxhr*nyhr)
                    !                  Calculate median
                    IF (MOD(moce(ji,jj),2) /=  0) THEN
                       zdep(ji,jj) = vardep1d( moce(ji,jj)/2 + 1)
                    ELSE
                       zdep(ji,jj)= ( vardep1d(moce(ji,jj)/2) + vardep1d(moce(ji,jj)/2+1) )/2.0
                    END IF
                    DEALLOCATE (vardep1d, istmp)
                 ENDIF
              ENDIF  ! moce == 0
              DEALLOCATE(glamhr,gphihr)
              DEALLOCATE(tabtestSW,tabtestSE,tabtestNE,tabtestNW,vardep,mask_oce,mask_box)
           ENDIF   !  end i_calc = 1
        ENDIF   !  ORCA point in the domain of the original bathymetry 
        !
        ! --------------------------------  END OF LOOP ON JI, JJ
     ENDDO
  ENDDO
  !-----------------------------------------------------------------
  PRINT *,MAXVAL(zdep),MINVAL(zdep)  

  PRINT *, ' MAX num of original bathy points in ORCA025 MESH:',knumax

  !--------------------------------------------------------------
  ! Boundary conditions  
  !--------------------------------------------------------------
  !          the bathymetry must be positive for OPA model.
  !          if it is not, convert it to positive values
  zdep  = ABS(zdep)

  !
  !   apply boundary condition on  bathy  
  CALL lbc ( zdep ,jpi,jpj, 1, 1 ,1, nperio)   

  !--------------------------------------------------------------
  !  NOW deal with output file 
  !--------------------------------------------------------------
  PRINT *, 'Opening  the output file:' ,TRIM(chfilout)
  STATUS = NF90_CREATE(TRIM(chfilout), NF90_CLOBBER, NCID_OUT)
  IF (STATUS /=  NF90_NOERR) THEN
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
  !
  STATUS = NF90_DEF_DIM(NCID_OUT, 'x', jpi, XID)
  STATUS = NF90_DEF_DIM(NCID_OUT, 'y', jpj, YID)
  !
  !  Definition of variables
  !
  STATUS=NF90_DEF_VAR(NCID_OUT,'nav_lon',NF90_FLOAT,(/XID,YID/),IDVAR1) 
  STATUS=NF90_DEF_VAR(NCID_OUT,'nav_lat',NF90_FLOAT,(/XID,YID/),IDVAR2) 
  !                                     ---------------  
  !                                     Dim_msk has 2 dimensions  XID and YID
  STATUS=NF90_DEF_VAR(NCID_OUT,'Bathymetry' ,NF90_FLOAT,(/XID,YID/),IDVAR3) 
  STATUS=NF90_DEF_VAR(NCID_OUT,'Pb'       ,NF90_FLOAT,(/XID,YID/),IDVAR4) 
  STATUS=NF90_DEF_VAR(NCID_OUT,'moce'       ,NF90_INT,(/XID,YID/),IDVAR5)
  STATUS=NF90_DEF_VAR(NCID_OUT,'mpoi'       ,NF90_INT,(/XID,YID/),IDVAR6)
  STATUS=NF90_ENDDEF(NCID_OUT)
  !
  !  Save variables
  !
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR1, glamt)
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR2, gphit)
  !
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR3, zdep  )
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR4, pb  )
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR5, moce  )
  STATUS=NF90_PUT_VAR(NCID_OUT, IDVAR6, mpoi  )

  !
  !  Close files
  !
  STATUS=NF90_CLOSE(NCID_OUT)
CONTAINS
!-------------------------------------------------------------------
SUBROUTINE hdlerr(istatus)
  !--------------------------------------------------------------------
  INTEGER                 :: istatus
  IF (istatus /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(istatus)
     STOP 'stopped here'
  ENDIF
  RETURN
END SUBROUTINE hdlerr

!-------------------------------------------------------------------
SUBROUTINE prihre( chtab, tab,nx,ny)
  !-------------------------------------------------------------------
  INTEGER          :: nx, ny
  REAL, DIMENSION(nx,ny) :: tab
  CHARACTER (len=10) :: chtab
  PRINT *, '-------------------', chtab, '--------------------'
  DO jj= ny,1,-1
     PRINT *, 'jj=',jj
     PRINT 100, tab(:,jj)
  ENDDO
100 FORMAT (10(1x,1pg10.3))
END SUBROUTINE prihre

END PROGRAM batinterp

