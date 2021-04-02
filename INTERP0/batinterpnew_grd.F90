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
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: batinterpnew_grd.F90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER                            :: jpi,jpj,jpibat,jpjbat,jpijbat
!
  REAL                               ::  bat_min,bat_max
  REAL,     DIMENSION(:,:) ,ALLOCATABLE  :: bato
#if defined etopo2
  INTEGER  ,DIMENSION(:,:)   ,ALLOCATABLE  :: itmp
#else
  INTEGER*2,DIMENSION(:)   ,ALLOCATABLE  :: itmp
#endif
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
      real,DIMENSION(:,:)   ,ALLOCATABLE  :: vardep,zdep
      real,DIMENSION(:,:)   ,ALLOCATABLE  :: glamhr,gphihr,pb
      real,DIMENSION(:,:)   ,ALLOCATABLE  :: tabtestsw,tabtestse,tabtestne,tabtestnw
      real,DIMENSION(:)     ,ALLOCATABLE  :: flo, fla, vardep1d
      integer,DIMENSION(:,:) ,ALLOCATABLE :: mask_box,mask_oce
      integer,DIMENSION(:)   ,ALLOCATABLE :: istmp
      real lonmin, lonmax, latmin, latmax, latsw, latse, latne, latnw, gphifmax
      real slatmin, slatmax              ! min and max latitudes in input file 
      integer jimin,jimax,nxhr,nyhr, i_calc
      integer :: i,j,jjmin,jjmax
      integer knumax, nperio
 
!
 include "netcdf.inc"

!-------------------------------------------------------------------
!  Choose averaging method: l_interp=1: median average, 
!  l_interp=0 : arithmetic average.
!-------------------------------------------------------------------
   l_interp=1 
!-------------------------------------------------------------------
!   sets input/ output file names
!-------------------------------------------------------------------

   chfilbatin = '/home/mathiot/Bathymetrie/GRD/1236703069735.grd';
   chfilout = '../PERIANT8_GMRT_mars2009.nc'
!---------------------------------------------------------------------------------------------
!lesommer@meolx06:~/test_grd$ grdinfo 1236703069735.grd
!1236703069735.grd: Title: Bathymetry Grid
!1236703069735.grd: Command:
!        Projection: Cylindrical Equidistant
!        this grid created by the program MapApp 1.0b
!1236703069735.grd: Remark:
!1236703069735.grd: Gridline node registration used
!1236703069735.grd: Grid file format: cf (# 10)
!1236703069735.grd: x_min: 180 x_max: 539.859 x_inc: 0.140625 name: Longitude nx: 2560
!1236703069735.grd: y_min: -78.062 y_max: 81.0057 y_inc: 0.0813229 name: Latitude ny: 1957
!1236703069735.grd: z_min: -10078.7 z_max: 6586.1 name: Eleveation (m)
!1236703069735.grd: scale_factor: 1 add_offset: 0
!-------------------------------------------------------------------
! 0 - Read grid file 
! ------------------------------------------------------------------
 chfilgrid = '../PERIANT8-coordinates.nc';
       PRINT *, 'Opening  the input file:' ,trim(chfilgrid)
       STATUS = NF_OPEN(trim(chfilgrid), NF_NOWRITE, NCID_IN)
       IF (STATUS .NE. NF_NOERR) THEN
           PRINT *, NF_STRERROR(STATUS)
           PRINT *, ' '
           PRINT *, 'Could not open the file:' ,trim(chfilgrid)
           PRINT *, ' '
           STOP 'We stop here'
       ENDIF

      STATUS = NF_INQ_DIMID(NCID_IN, 'x', LONID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LONID, jpi)
      STATUS = NF_INQ_DIMID(NCID_IN, 'y', LATID)      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LATID, jpj)
      print *,'--    Reading input file, dimension jpi=',jpi,' jpj=',jpj
!
!   for a global grid we need to construct a criterion 
!
!   ALLOCATE arrays for longitude and latitude 
! 
      ALLOCATE(glamf(jpi,jpj),gphif(jpi,jpj),glamt(jpi,jpj),gphit(jpi,jpj))
      ALLOCATE(coord(jpi,jpj,1,1),zdep(jpi,jpj))
!
!            read longitude and latitude 
!
      STATUS=NF_INQ_VARID(NCID_IN, 'glamf', IDVAR1_IN)
      STATUS = NF_GET_VAR_DOUBLE(NCID_IN, IDVAR1_IN, coord)
      glamf(:,:)=coord(:,:,1,1)
      STATUS=NF_INQ_VARID(NCID_IN, 'gphif', IDVAR1_IN)
      STATUS = NF_GET_VAR_DOUBLE(NCID_IN, IDVAR1_IN, coord)
      gphif(:,:)=coord(:,:,1,1)
      STATUS=NF_INQ_VARID(NCID_IN, 'glamt', IDVAR2_IN)
      STATUS = NF_GET_VAR_DOUBLE(NCID_IN, IDVAR2_IN, coord)
      glamt(:,:)=coord(:,:,1,1)
      STATUS=NF_INQ_VARID(NCID_IN, 'gphit', IDVAR2_IN)
      STATUS = NF_GET_VAR_DOUBLE(NCID_IN, IDVAR2_IN, coord)
      gphit(:,:)=coord(:,:,1,1)
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

      STATUS=NF_CLOSE(NCID_IN)

      deallocate(coord)
!--------------------------------------------------------------
! Boundary conditions  
!--------------------------------------------------------------
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
!--------------------------------------------------------------
!
!    Allocate arrays for bathymetry 
     ALLOCATE ( pb(jpi,jpj), moce(jpi,jpj), mpoi(jpi,jpj))

!-------------------------------------------------------------------
! 1 - Read original bathy file 
! ------------------------------------------------------------------
       PRINT *, 'Opening  the input file:' ,trim(chfilbatin)
       STATUS = NF_OPEN(trim(chfilbatin), NF_NOWRITE, NCID_IN)
       IF (STATUS .NE. NF_NOERR) THEN
           PRINT *, NF_STRERROR(STATUS)
           PRINT *, ' '
           PRINT *, 'Could not open the file:' ,trim(chfilbatin)
           PRINT *, ' '
           STOP 'We stop here'
       ENDIF
! .....................; Case ETOPO1
      STATUS = NF_INQ_DIMID(NCID_IN, 'xysize', LONID)
      jpibat=2560
      jpjbat=1957      
      STATUS = NF_INQ_DIMLEN(NCID_IN, LONID, jpijbat)
      PRINT *,'--    Reading gebco file, dimension jpijbat=',jpijbat
      IF (jpijbat .NE. jpibat*jpjbat) THEN
           PRINT *, ' error in dimensions for gebco file'
           STOP
      ENDIF
!
      ALLOCATE(itmp(jpijbat),bato(jpibat,jpjbat),flo(jpibat),fla(jpjbat))
      call hdlerr(NF_INQ_VARID(NCID_IN, 'z', IDVAR1_IN))
      call hdlerr(NF_GET_VAR_INT2(NCID_IN, IDVAR1_IN, itmp))
      jk=1
      DO jj=jpjbat,1,-1
         DO ji=1,jpibat
           bato(ji,jj) = -itmp(jk)
           jk=jk+1
         ENDDO
      ENDDO
      WHERE (bato<=0) bato=0
      
      DO jj=1,jpjbat
         fla(jj)=-78.062+0.0813229*(jj-1)
      ENDDO
      DO ji=1,jpibat
         flo(ji)=180.+0.140625*float(ji-1)
      ENDDO
      flo=flo-360

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
       chfilout='GMRT.nc'
       PRINT *, 'Opening  the output file:' ,trim(chfilout)
       STATUS = NF_CREATE(trim(chfilout), NF_CLOBBER, NCID_OUT)
!
       STATUS = NF_DEF_DIM(NCID_OUT, 'lon', jpibat, XID)
       STATUS = NF_DEF_DIM(NCID_OUT, 'lat', jpjbat, YID)
!  Definition of variables
!
      STATUS=NF_DEF_VAR(NCID_OUT,'lon',NF_FLOAT,1,XID,IDVAR1) 
      STATUS=NF_DEF_VAR(NCID_OUT,'lat',NF_FLOAT,1,YID,IDVAR2) 
!                   
!                  --------------- 
      DIM_msk(1)=XID
      DIM_msk(2)=YID
!                                     Dim_msk has 2 dimensions  XID and YID
      STATUS=NF_DEF_VAR(NCID_OUT,'Bathymetry' ,NF_FLOAT,2,DIM_msk(1:2),IDVAR3)
      STATUS=NF_PUT_ATT_REAL (NCID_OUT, IDVAR3, 'missing_value', NF_REAL, 1, 0)
!     STATUS=NF_DEF_VAR(NCID_OUT,'mcoast'     ,NF_INT,2,DIM_msk(1:2),IDVAR7)
      STATUS=NF_ENDDEF(NCID_OUT)
!
!  Save variables
!
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR1, flo)
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR2, fla)
!
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR3, bato  )



!
!  Close files
!
      STATUS=NF_CLOSE(NCID_OUT)
      STOP
!      slatmax =  88-0.125    !  this file is restricted in latitude 
!      slatmin = -88+0.125    !     (not truly global) 
      slatmax =  81.0057    !  this file is restricted in latitude 
      slatmin = -78.062    !     (not truly global) 

!
!   calculate max and min value of input array bathy 
! 
      bat_max= MAXVAL(itmp)
      bat_min= MINVAL(itmp)
      PRINT *, ' read original bathy OK max=',bat_max,' min =', bat_min
      PRINT *, ' read original lon   OK max=',MAXVAL(flo),' min =',MINVAL(flo) 
      PRINT *, ' read original lat   OK max=',MAXVAL(fla),' min =',MINVAL(fla) 
      deallocate (itmp)
      call hdlerr(NF_CLOSE(NCID_IN))

!-------------------------------------------------------------------------
! 2 -  INTERPOLATION    : LOOP OVER ORCA GRID  
!-------------------------------------------------------------------------
  knumax = 0
  do jj=2,jpj-1
       IF (mod(jj,100).EQ.0) THEN
          PRINT *, jj,'/',jpj-1
       END IF
    do ji=2,jpi-1
!

   
!     Defining grid parameters for T-point (ji,jj)
      lonmin=MIN(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))  
      lonmax=MAX(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))
      latmin=MIN(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))
      latmax=MAX(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))


      jimin=floor(1/0.140625*(lonmin+180))+1 !floor(60*(lonmin+180))+1
      jimax=floor(1/0.140625*(lonmax+180))+2 !floor(60*(lonmax+180))+2
      jjmin=floor(1/0.0813229*(latmin+78.062))+1 !floor(60*(latmin+90))+1
      jjmax=floor(1/0.0813229*(latmax+78.062))+2 !floor(60*(latmax+90))+2

!    ----------------- Control print for debugging
!     IF (mod(jj,200).EQ.1 .AND. mod(ji,200) .EQ. 1) THEN
!        PRINT 102, '---- ji,jj:',ji,jj,' jimin,max:',jimin,jimax,' jjmin,max:',jjmin,jjmax
!        PRINT 103,'      glamf:',glamf(ji-1,jj-1),glamf(ji-1,jj  ), &
!                                 glamf(ji  ,jj ),glamf(ji  ,jj-1  )
!       PRINT 103,'      gphif:',gphif(ji-1,jj-1),gphif( ji-1,jj ), &
!                                 gphif(ji  ,jj ),gphif(ji  ,jj-1 )
!      ENDIF

 103  format(a25,4(1x,f11.3))
 102  format(3(a12,i5,1x,i5))
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
      IF ( abs(lonmax-lonmin) .GE. 120) THEN 
!    -------------- first north pole case. 
!                   just take the northernmost value of etopo2. 
        IF (latmax .GE. gphifmax ) THEN
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
         IF ( glamf(ji-1,jj-1).GT.0.AND.glamf(ji-1,jj).GT.0) & 
              lonmax=MIN(glamf(ji-1,jj-1), glamf(ji-1,jj) )  
         IF (glamf(ji,jj-1).LT.0 .AND. glamf(ji,jj).LT.0)  &
                   lonmin=MAX( glamf(ji,jj-1),glamf(ji,jj) )
         i_calc = 1


         jimin=floor(60*(lonmin+180))+2
         jimax=floor(60*(lonmax+180))+3

         jimax = min(jimax, jpibat)
         jjmax = min(jjmax, jpjbat)
         nxhr = jpibat-jimax+1 + jimin+1
         nyhr=jjmax-jjmin+1

         allocate(vardep(nxhr,nyhr))
! pastes the two extremities in the correct order
         vardep(1:jpibat-jimax+1,:)= bato(jimax:jpibat,jjmin:jjmax)
         vardep(jpibat-jimax+2:nxhr,:)=bato(1:jimin+1 ,jjmin:jjmax)

         allocate(glamhr(nxhr,nyhr),gphihr(nxhr,nyhr))
         do j=1,nyhr
            glamhr(1:jpibat-jimax+1,j)=flo(jimax:jpibat)-360
            glamhr(jpibat-jimax+2:nxhr,j)=flo(1:jimin+1)
         enddo
         do i=1,nxhr
            gphihr(i,:)=fla(jjmin:jjmax)
         enddo

         latsw=glamf(ji-1,jj-1)
         latse=glamf(ji  ,jj-1)
         latne=glamf(ji  ,jj  )
         latnw=glamf(ji-1,jj  )
         if( latsw .GE. 0 ) latsw =latsw-360
         if( latse .GE. 0 ) latse =latse-360
         if( latne .GE. 0 ) latne =latne-360
         if( latnw .GE. 0 ) latnw =latnw-360
         allocate(tabtestSW(nxhr,nyhr))
         allocate(tabtestSE(nxhr,nyhr))
         allocate(tabtestNE(nxhr,nyhr))
         allocate(tabtestNW(nxhr,nyhr))

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
         jimax = min(jimax, jpibat)
         jjmax = min(jjmax, jpjbat)
         nxhr=jimax-jimin+1
         nyhr=jjmax-jjmin+1
         allocate(vardep(nxhr,nyhr))
         vardep=bato(jimin:jimax,jjmin:jjmax)

         allocate(glamhr(nxhr,nyhr),gphihr(nxhr,nyhr))
         do j=1,nyhr
            glamhr(:,j)=flo(jimin:jimax)
         enddo
         do i=1,nxhr
            gphihr(i,:)=fla(jjmin:jjmax)
         enddo

! Define the 4 'equations de droite'
      allocate(tabtestSW(nxhr,nyhr))
      allocate(tabtestSE(nxhr,nyhr))
      allocate(tabtestNE(nxhr,nyhr))
      allocate(tabtestNW(nxhr,nyhr))
      tabtestSW=(glamf(ji  ,jj-1)-glamf(ji-1,jj-1))*(gphihr-gphif(ji-1,jj-1)) - (gphif(ji  ,jj-1)-gphif(ji-1,jj-1))*(glamhr-glamf(ji-1,jj-1))
      tabtestSE=(glamf(ji  ,jj  )-glamf(ji  ,jj-1))*(gphihr-gphif(ji  ,jj-1)) - (gphif(ji  ,jj  )-gphif(ji  ,jj-1))*(glamhr-glamf(ji  ,jj-1))
      tabtestNE=(glamf(ji-1,jj  )-glamf(ji  ,jj  ))*(gphihr-gphif(ji  ,jj  )) - (gphif(ji-1,jj  )-gphif(ji  ,jj  ))*(glamhr-glamf(ji  ,jj  ))
      tabtestNW=(glamf(ji-1,jj-1)-glamf(ji-1,jj  ))*(gphihr-gphif(ji-1,jj  )) - (gphif(ji-1,jj-1)-gphif(ji-1,jj  ))*(glamhr-glamf(ji-1,jj  ))
  endif
!--------------------------------------------------------------------------------------- 
!  NOW we have the original bathymetry in the domain containing the ORCA cell.
!  We have the 4 line equations to find which points are incide or outside the cell
!  Find those points and do interpolation
!--------------------------------------------------------------------------------------- 
!                      does not do it if we are at the pole 
  IF  ( i_calc == 1 ) THEN
     allocate(mask_oce(nxhr,nyhr),mask_box(nxhr,nyhr))
     mask_oce=0
     mask_box = 0
     where ( (tabtestSW .GE. 0.) .and. &
                (tabtestSE .GE. 0.) .and.  &
                (tabtestNE .GE. 0.) .and.  &
                (tabtestNW .GE. 0.)) mask_box = 1
     mpoi(ji,jj) = sum(mask_box)
! 
!    problem: does not find any ETOPO2 point in the ORCA grid box
!
     IF (mpoi(ji,jj) == 0. ) THEN
        IF ( MINVAL(vardep) >= 0 ) THEN
!          we are probably close to the pole of the grid, and on land. no problem
           pb(ji,jj) = -1
        ELSE
           pb(ji,jj) = 1
           !PRINT 102, 'pblm ji,jj:',ji,jj,' jimin,max:',jimin,jimax,' jjmin,max:',jjmin,jjmax
           !PRINT 103,'      glamf:',glamf(ji-1,jj-1),glamf(ji-1,jj  ), &
           !                   glamf(ji  ,jj ),glamf(ji  ,jj-1  )
           !PRINT 103,'      gphif:',gphif(ji-1,jj-1),gphif( ji-1,jj ), &
           !                   gphif(ji  ,jj ),gphif(ji  ,jj-1 )
           !PRINT *, ' nxhr:',nxhr,'nyhr:',nyhr
           !print *, ' glamhr:',glamhr(:,1)
           !print *, ' gphihr:',gphihr(1,:)
           !call prihre('vardep   :',vardep   ,nxhr,nyhr)
           !call prihre('tabtestSW:',tabtestSW,nxhr,nyhr)
           !call prihre('tabtestSE:',tabtestSE,nxhr,nyhr)
           !call prihre('tabtestNE:',tabtestNE,nxhr,nyhr)
           !   call prihre('tabtestNW:',tabtestNW,nxhr,nyhr)
        ENDIF
     ENDIF   ! mpoi == 0
!                     find total number of ocean points 
     where ( (vardep .LT. 0.) .and. (tabtestSW.GE. 0.) .and. &
                 (tabtestSE .GE. 0.) .and.  &
                 (tabtestNE .GE. 0.) .and.  &
                 (tabtestNW .GE. 0.) ) mask_oce = 1
     moce(ji,jj) = sum(mask_oce)
     knumax = max(knumax,moce(ji,jj))
!                     find bathymetry value by interpolation
     if (moce(ji,jj) == 0 )  THEN
        zdep(ji,jj) = 0
     else
        IF (l_interp == 0 ) THEN
!                Arithmetic average
           zdep(ji,jj) = sum (vardep*mask_oce)/float(moce(ji,jj))
        ELSE
!                Median average 
!       put very large values on land to eliminate those points when taking the median 
           vardep = vardep*mask_oce  - 20000*(1-mask_oce)
           allocate(vardep1d(nxhr*nyhr),istmp(nxhr*nyhr))
           vardep1d = reshape(vardep,(/ nxhr*nyhr /) )
           call ssort(vardep1d,istmp,nxhr*nyhr)
!                  Calculate median
           IF (mod(moce(ji,jj),2) .NE. 0) THEN
               zdep(ji,jj) = vardep1d( moce(ji,jj)/2 + 1)
           ELSE
               zdep(ji,jj)= ( vardep1d(moce(ji,jj)/2) + vardep1d(moce(ji,jj)/2+1) )/2.0
           END IF
           deallocate (vardep1d, istmp)
        ENDIF
     endif  ! moce == 0
     deallocate(glamhr,gphihr)
     deallocate(tabtestSW,tabtestSE,tabtestNE,tabtestNW,vardep,mask_oce,mask_box)
   ENDIF   !  end i_calc = 1
   ENDIF   !  ORCA point in the domain of the original bathymetry 
!
! --------------------------------  END OF LOOP ON JI, JJ
enddo
enddo
!-----------------------------------------------------------------
    print *,MAXVAL(zdep),MINVAL(zdep)  

PRINT *, ' MAX num of original bathy points in ORCA025 MESH:',knumax

!--------------------------------------------------------------
! Boundary conditions  
!--------------------------------------------------------------
!          the bathymetry must be positive for OPA model.
!          if it is not, convert it to positive values
  zdep  = abs(zdep)

!
!   apply boundary condition on  bathy  
  CALL lbc ( zdep ,jpi,jpj, 1, 1 ,1, nperio)   
     
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
!
       STATUS = NF_DEF_DIM(NCID_OUT, 'x', JPI, XID)
       STATUS = NF_DEF_DIM(NCID_OUT, 'y', JPJ, YID)
       DIM_msk(1)=XID
       DIM_msk(2)=YID
!
!  Definition of variables
!
      STATUS=NF_DEF_VAR(NCID_OUT,'nav_lon',NF_FLOAT,2,DIM_msk(1:2),IDVAR1) 
      STATUS=NF_DEF_VAR(NCID_OUT,'nav_lat',NF_FLOAT,2,DIM_msk(1:2),IDVAR2) 
!                                     ---------------  
!                                     Dim_msk has 2 dimensions  XID and YID
      STATUS=NF_DEF_VAR(NCID_OUT,'Bathymetry' ,NF_FLOAT,2,DIM_msk(1:2),IDVAR3) 
      STATUS=NF_DEF_VAR(NCID_OUT,'Pb'     ,NF_FLOAT,2,DIM_msk(1:2),IDVAR4) 
      STATUS=NF_DEF_VAR(NCID_OUT,'moce'       ,NF_INT,2,DIM_msk(1:2),IDVAR5)
      STATUS=NF_DEF_VAR(NCID_OUT,'mpoi'       ,NF_INT,2,DIM_msk(1:2),IDVAR6)
!      STATUS=NF_DEF_VAR(NCID_OUT,'mcoast'     ,NF_INT,2,DIM_msk(1:2),IDVAR7)
      STATUS=NF_ENDDEF(NCID_OUT)
!
!  Save variables
!
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR1, glamt)
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR2, gphit)
!
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR3, zdep  )
      STATUS=NF_PUT_VAR_REAL(NCID_OUT, IDVAR4, pb  )
      STATUS=NF_PUT_VAR_INT (NCID_OUT, IDVAR5, moce  )
      STATUS=NF_PUT_VAR_INT (NCID_OUT, IDVAR6, mpoi  )



!
!  Close files
!
      STATUS=NF_CLOSE(NCID_OUT)

      stop
      end

!-------------------------------------------------------------------
      subroutine hdlerr(istatus)
!--------------------------------------------------------------------
      integer                 :: istatus
      include 'netcdf.inc'
      if (istatus .ne. NF_NOERR) then
         print *, NF_STRERROR(istatus)
         stop 'stopped here'
      endif
      return
      end

!-------------------------------------------------------------------
      subroutine prihre( chtab, tab,nx,ny)
!-------------------------------------------------------------------
      integer          :: nx, ny
      real, dimension(nx,ny) :: tab
      character (len=10) :: chtab
      print *, '-------------------', chtab, '--------------------'
      DO jj= ny,1,-1
        print *, 'jj=',jj
        print 100, tab(:,jj)
      ENDDO
 100  format (10(1x,1pg10.3))
      return
      end

