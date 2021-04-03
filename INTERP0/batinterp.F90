!-------------------------------------------------------------------
PROGRAM batinterp
  !-------------------------------------------------------------------
  !  interpolates a bathy onto an ORCA grid.
  !  The method is by averaging all the original bathy points that fall into 
  !  an orca grid cell. This is very good for relatively coarse grids 
  !  (like ORCA025).
  !---------------------------------------------------------
  !   Two averaging methods are provided: 
  !     nn_interp = 0 : arithmetic average;
  !     nn_interp = 1 : median average.
  !--------------------------------------------------------
  !-------------------------------------------------------------------
  !  History: derived by Elisabeth Remy from a version of opabat provided by 
  !  Maurice Imbard (LODYC) in 1995.
  !  Modified by A.M. Treguier and  Olivier Le Galloudeg (Mercator)
  !  June 2005:
  !  Generalized by Pierre Mathiot 2009-2020
  !  Cleaning and good practice Jean-Marc Molines 2021
  !---------------------------------------------------------------------
  USE netcdf
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-07-24 08:23:41 +0200 (ven. 24 juil. 2009) $
  !!  $Rev: 235 $
  !!  $Id: batinterpnew.F90 235 2009-07-24 06:23:41Z forge $
  !! -------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER(KIND=4)                 :: npiglo,npjglo    ! size of NEMO grid
  INTEGER(KIND=4)                 :: npibat,npjbat    ! size of external Bathy file
  INTEGER(KIND=4)                 :: numnam=99        ! namelist logical unit
  INTEGER(KIND=4)                 :: iargc, narg, ijarg
  INTEGER(KIND=4)                 :: iimin, iimax, ijmin, ijmax  
  !
  INTEGER(KIND=4)                 :: ierr
  INTEGER(KIND=4)                 :: ncid, ncio, idx, idy, id
  INTEGER(KIND=4)                 :: idv_lon, idv_lat, idv_bat
  INTEGER(KIND=4)                 :: idv_moce, idv_mpoi

  INTEGER(KIND=4)                 :: npiloc, npjloc     ! size of local subdomain
  INTEGER(KIND=4)                 :: ji,jj,  jii, jij   ! dummy loop index
  INTEGER(KIND=4)                 :: ik                 ! index counter
  INTEGER(KIND=4)                 :: iiloc1, iiloc2, iiloc3, iiloc4
  INTEGER(KIND=4)                 :: ijloc1, ijloc2, ijloc3, ijloc4 
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mpoi, moce
  INTEGER(KIND=4), DIMENSION(:,:), ALLOCATABLE :: mask_box, mask_oce

  REAL(KIND=4)                               :: median          ! external function from scilib
  REAL(KIND=4)                               :: zvmin, zvmax, zcrit
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: glamf, gphif,glamt,gphit
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: vardep, zdep, vardeploc
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE  :: tabtestsw,tabtestse,tabtestne,tabtestnw
  REAL(KIND=4), DIMENSION(:),   ALLOCATABLE  :: vardep1d

  REAL(KIND=8)                               :: dlonmin, dlonmax, dlatmin, dlatmax
  REAL(KIND=8)                               :: dlatmin_in, dlatmax_in
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: dflo, dfla, dfloloc, dflaloc, dzflo, dzfla

  
  CHARACTER(LEN=80)                          :: cldum         ! dummy char variable
  CHARACTER(LEN=80)                          :: cf_tplnamlst  ! template namelist file name created with -p
  CHARACTER(LEN=80)                          :: cf_namlist    ! namelist file name

  ! Namelist variable
  INTEGER(KIND=4)     :: nn_interp   ! interpolation method : O: Arithmetic Average   1: Median average
  INTEGER(KIND=4)     :: nn_perio    ! NEMO Periodic condition. If 99 try to guess the correct condition
                                     !           from the size of the domain (case of ORCA xx grids).

  CHARACTER(len=140)  :: cn_fout     ! NEMO bathymetry
  CHARACTER(len=140)  :: cn_varout   ! NEMO bathymetry variable name
  CHARACTER(len=140)  :: cn_fgrid    ! NEMO coordinates/mesh_hgr/domain_cfg file
  CHARACTER(len=140)  :: cn_fbatin   ! External Bathymetry file
  CHARACTER(len=140)  :: cn_varin    ! External Bathymetry variable name
  CHARACTER(len=140)  :: cn_xdim     ! External Bathymetry file X dimension name
  CHARACTER(len=140)  :: cn_ydim     ! External Bathymetry file Y dimension name
  CHARACTER(len=140)  :: cn_lonv     ! External Bathymetry file longitude name
  CHARACTER(len=140)  :: cn_latv     ! External Bathymetry file Latitude name

  LOGICAL             :: ln_sign     ! If True change the sign of surface elevation (case of Ice Shelves)
  LOGICAL             :: ln_regin    ! True : External Bathymetry file regular, False:  irregular

  NAMELIST /naminterpo/ nn_interp, ln_sign, cn_fbatin, cn_fout, cn_fgrid, nn_perio,      &
         &              cn_varin, cn_varout, cn_xdim, cn_ydim, ln_regin, cn_lonv, cn_latv
  !-------------------------------------------------------------------
  narg = iargc()

  IF ( narg == 0 ) THEN
     PRINT *,' usage : batinterp.exe -f NAMELIST-name'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'        This program is part of the NEMOBAT package used to create a '
     PRINT *,'      bathymetry for NEMO configurations. All parameters are set  '
     PRINT *,'      in a namelist file, whose name is passed to the program with -f option.'
     PRINT *,'      '
     PRINT *,'        This program is used to interpolate an external bathymetry'
     PRINT *,'      data base (such as GEBCO, ETOPO, BedMachine etc... ) to a NEMO'
     PRINT *,'      configuration defined by its horizontal grid. Supported interpolation'
     PRINT *,'      techniques are presently arithmetic averaging or median average.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f NAMELIST-name : Pass the name of the namelist file to use.'
     PRINT *,'              A template namelist can be created using option -p '
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -p TEMPLATE_namelist: create a template namelist whose name is'
     PRINT *,'             TEMPLATE_namelist for further editing before use.' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       All needed files are indicated in the namelist' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : guess bathymetry, whose name is given in the namelist.'
     PRINT *,'         variables : the name of the bathymetry is specified in namelist.'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'       '
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_namlist   ) ; ijarg=ijarg+1
        ! option
     CASE ( '-p'   ) ; CALL getarg(ijarg, cf_tplnamlst ) ; ijarg=ijarg+1 
                     ; CALL PrintTplNamelist(cf_tplnamlst) 
                     ; STOP 'Template namelist created'
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! Set default namelist values
  nn_interp=1 
  ln_sign = .FALSE.   ! to change the sign of data to deal with surface elevation (isf for exemple)
  !
  cn_fbatin = 'RTOPO1/RTopo105b_50S.nc'
  cn_varin  = 'bathy' 
  cn_xdim   = 'lon'
  cn_ydim   = 'lat'
  cn_lonv   = 'lon'
  cn_latv   = 'lat'
  ln_regin  = .TRUE.
  !
  cn_fout   = 'bathy_PIG025_rtopo1.nc'
  cn_varout = 'Bathymetry'
  !
  cn_fgrid  = 'PIG025-I/PIG025_coordinates.nc'
  !
  nn_perio  = 0
  !
  !
  !-------------------------------------------------------------------
  ! READ NAMELIST
  !-------------------------------------------------------------------
  OPEN  ( numnam, FILE=cf_namlist, STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL')
  READ  ( numnam, naminterpo )
  !-------------------------------------------------------------------
  ! 0 - Read grid file 
  ! ------------------------------------------------------------------
  PRINT *, 'Opening  the input file: ' ,TRIM(cn_fgrid)
  ierr = NF90_OPEN(TRIM(cn_fgrid), NF90_NOWRITE, ncid)
  IF (ierr /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(ierr)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(cn_fgrid)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF

  ierr = NF90_INQ_DIMID(ncid, 'x', id)      
  ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npiglo)
  ierr = NF90_INQ_DIMID(ncid, 'y', id)      
  ierr = NF90_INQUIRE_DIMENSION(ncid, id, len=npjglo)
  PRINT *,'--    Reading input file, dimension npiglo=',npiglo,' npjglo=',npjglo
  !
  !   for a global grid we need to construct a criterion 
  !
  !   ALLOCATE arrays for longitude and latitude 
  ! 
  ALLOCATE(glamf(npiglo,npjglo),gphif(npiglo,npjglo),glamt(npiglo,npjglo),gphit(npiglo,npjglo))
  ALLOCATE(zdep(npiglo,npjglo))
  !
  !            read longitude and latitude  of F and T points
  !
  ierr = NF90_INQ_VARID(ncid, 'glamf', id)
  ierr = NF90_GET_VAR(ncid, id, glamf(:,:), start=(/1,1/), count=(/npiglo,npjglo/) )
  ierr = NF90_INQ_VARID(ncid, 'gphif', id)
  ierr = NF90_GET_VAR(ncid, id, gphif(:,:), start=(/1,1/), count=(/npiglo,npjglo/) )
  ierr = NF90_INQ_VARID(ncid, 'glamt', id)
  ierr = NF90_GET_VAR(ncid, id, glamt(:,:), start=(/1,1/), count=(/npiglo,npjglo/) )
  ierr = NF90_INQ_VARID(ncid, 'gphit', id)
  ierr = NF90_GET_VAR(ncid, id, gphit(:,:), start=(/1,1/), count=(/npiglo,npjglo/) )

  ! check (??)
  zvmax= MAXVAL(glamf)  ; zvmin= MINVAL(glamf)
  PRINT *, ' read  glamf  OK max=',zvmax,' min =', zvmin 
  zvmax= MAXVAL(glamt)  ; zvmin= MINVAL(glamt)
  PRINT *, ' read  glamt  OK max=',zvmax,' min =', zvmin 
  zvmax= MAXVAL(gphif)  ; zvmin= MINVAL(gphif)
  PRINT *, ' read  gphif  OK max=',zvmax,' min =', zvmin 
  zvmax= MAXVAL(gphit)  ; zvmin= MINVAL(gphit)
  PRINT *, ' read  gphit  OK max=',zvmax,' min =', zvmin 

  ierr=NF90_CLOSE(ncid)

  !--------------------------------------------------------------
  ! Boundary conditions  
  !--------------------------------------------------------------
  IF (nn_perio == 99) THEN
     IF (npiglo== 182 .AND.npjglo==149) THEN 
        nn_perio= 4
        PRINT *, 'This is an ORCA R2 configuration, nn_perio=4'
     ELSEIF (npiglo== 722 .AND.npjglo==511) THEN
        nn_perio= 6
        PRINT *, 'This is an ORCA R05 configuration, nn_perio=6'
     ELSEIF (npiglo== 1442 .AND.npjglo==1021) THEN
        nn_perio= 4
        PRINT *, 'This is an ORCA R025 configuration, nn_perio=4'
     ELSE
        PRINT *, 'your configuration is unknown: you need to choose nn_perio '
        PRINT *, '     jperio= 0, closed,      jperio= 1, cyclic east-west'
        PRINT *, '     jperio= 2, equatorial symmetric'
        PRINT *, '     jperio= 3, north fold with T-point pivot; jperio= 4, the same + cyclic east-west'   
        PRINT *, '     jperio= 5, north fold with F-point pivot; jperio= 6, the same + cyclic east-west'   
        READ (5,*) nn_perio
     ENDIF
  END IF
  !--------------------------------------------------------------
  !
  !    Allocate arrays for bathymetry 
  ALLOCATE ( moce(npiglo,npjglo), mpoi(npiglo,npjglo))

  !-------------------------------------------------------------------
  ! 1 - Read original bathy file 
  ! ------------------------------------------------------------------
  PRINT *, 'Opening  the input file: ' ,TRIM(cn_fbatin)
  ierr = NF90_OPEN(TRIM(cn_fbatin), NF90_NOWRITE, ncid)
  IF (ierr /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(ierr)
     PRINT *, ' '
     PRINT *, 'Could not open the file:' ,TRIM(cn_fbatin)
     PRINT *, ' '
     STOP 'We stop here'
  ENDIF
  ! We suppose that all input files are netcdf map. If grd, preprocess with grd2nc.
  ! etopo1 and gebco
  ! We also assume that original bathy file have longitude and latitude variables (cn_lonv and cn_latv)
  ! If this is not the case, preprocess the file with NSIDC_map_trf (ask JMM .. )
  CALL hdlerr(NF90_INQ_DIMID(ncid, cn_xdim, id))      
  CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid, id, len=npibat))
  CALL hdlerr(NF90_INQ_DIMID(ncid, cn_ydim, id))      
  CALL hdlerr( NF90_INQUIRE_DIMENSION(ncid, id, len=npjbat))
  ALLOCATE (vardep(npibat,npjbat))
  PRINT *,'--   Reading ',TRIM(cn_fbatin),' file, dimension ',npibat,npjbat

  CALL hdlerr(NF90_INQ_VARID(ncid,cn_varin,id))
  CALL hdlerr(NF90_GET_VAR(ncid,id,vardep))

  ! 
  IF ( ln_regin ) THEN
     ALLOCATE(dflo(npibat,npjbat),dfla(npibat,npjbat), dzflo(npibat,1), dzfla(npjbat,1))
     PRINT *, 'input data are on a case regular grid (1d)'
     CALL hdlerr(NF90_INQ_VARID(ncid, cn_lonv, id))
     CALL hdlerr(NF90_GET_VAR(ncid, id, dzflo))

     CALL hdlerr(NF90_INQ_VARID(ncid, cn_latv, id))
     CALL hdlerr(NF90_GET_VAR(ncid, id, dzfla))

     DO ji=1,npjbat
        dflo(:,ji)=dzflo(:,1)
     END DO
     DO ji=1,npibat
        dfla(ji,:)=dzfla(:,1)
     END DO
  ELSE
     ALLOCATE(dflo(npibat,npjbat),dfla(npibat,npjbat))
     CALL hdlerr(NF90_INQ_VARID(ncid, cn_lonv, id))
     CALL hdlerr(NF90_GET_VAR(ncid, id, dflo))

     CALL hdlerr(NF90_INQ_VARID(ncid, cn_latv, id))
     CALL hdlerr(NF90_GET_VAR(ncid, id, dfla))
  END IF

  dlatmin_in = -99               !  99 values for a global file
  dlatmax_in =  99
  !   calculate max and min value of input array bathy 
  ! 
  zvmax= MAXVAL(vardep)
  zvmin= MINVAL(vardep)
  PRINT *, ' read original bathy OK max=',zvmax,' min =', zvmin
  PRINT *, ' read original lon   OK max=',MAXVAL(dflo),' min =',MINVAL(dflo) 
  PRINT *, ' read original lat   OK max=',MAXVAL(dfla),' min =',MINVAL(dfla) 
  ! test 
  dlatmin_in=MINVAL(dfla)
  dlatmax_in=MAXVAL(dfla)

  CALL hdlerr(NF90_CLOSE(ncid))

  !-------------------------------------------------------------------------
  ! 2 -  INTERPOLATION    : LOOP OVER ORCA GRID  
  !-------------------------------------------------------------------------
  iiloc1=INT(npibat/2); ijloc1=INT(npjbat/2)
  zdep(:,:)=-9999.99
  DO jj=2,npjglo-1
     IF (MOD(jj,1) == 0) PRINT *, jj,'/',npjglo-1
     DO ji=2,npiglo-1
        !
        !     Defining grid parameters for T-point (ji,jj)
        dlonmin=MIN(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))
        dlonmax=MAX(glamf(ji-1,jj-1),glamf(ji,jj-1),glamf(ji,jj),glamf(ji-1,jj))
        dlatmin=MIN(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))
        dlatmax=MAX(gphif(ji-1,jj-1),gphif(ji,jj-1),gphif(ji,jj),gphif(ji-1,jj))

        ! remove point using missing data
        IF (ABS(glamf(ji-1,jj)) + ABS(glamf(ji,jj)) + ABS(glamf(ji,jj-1)) + ABS(glamf(ji-1, jj-1)) >=  1.0e10) THEN
           moce(ji,jj) = 0 ; zdep(ji,jj) = -9999.99 ; mpoi (ji,jj) = 0 ;
           ! skip part out of input data
        ELSE IF ( dlatmax_in < dlatmin .OR. dlatmin_in > dlatmax ) THEN
           moce(ji,jj) = 0 ; zdep(ji,jj) = -9999.99 ; mpoi (ji,jj) = 0 ;
        ELSE
           IF ( ABS(glamf(ji-1,jj) - glamf(ji+1,jj)) >= 180 ) THEN
              IF (MAX(glamf(ji-1,jj),glamf(ji+1,jj)) >  180.) THEN
                 WHERE (glamf > 180.) 
                    glamf = glamf - 360. 
                 END WHERE
                 WHERE (dflo  > 180.d0 )
                    dflo = dflo - 360.d0
                 END WHERE
              ELSE
                 WHERE (glamf < 0.)
                    glamf = glamf + 360.d0
                 END WHERE
                 WHERE (dflo < 0.)
                    dflo = dflo + 360.d0
                 END WHERE
              END IF
           END IF

           ! position in input grid of dlatmin/dlonmin, dlatmin/dlonmax, dlatmax/dlonmin, dlatmax/dlonmax
           CALL NearestPoint(dlonmin, dlatmax, npibat, npjbat, dflo, dfla, iiloc1, ijloc1)
           iiloc2=iiloc1; ijloc2=ijloc1; CALL NearestPoint(dlonmin, dlatmin, npibat, npjbat, dflo, dfla, iiloc2, ijloc2)
           iiloc3=iiloc2; ijloc3=ijloc2; CALL NearestPoint(dlonmax, dlatmax, npibat, npjbat, dflo, dfla, iiloc3, ijloc3)
           iiloc4=iiloc3; ijloc4=ijloc3; CALL NearestPoint(dlonmax, dlatmin, npibat, npjbat, dflo, dfla, iiloc4, ijloc4)

           iimin=MIN(iiloc1, iiloc2, iiloc3, iiloc4)
           ijmin=MIN(ijloc1, ijloc2, ijloc3, ijloc4)
           iimax=MAX(iiloc1, iiloc2, iiloc3, iiloc4)
           ijmax=MAX(ijloc1, ijloc2, ijloc3, ijloc4)
           npiloc=iimax-iimin+1
           npjloc=ijmax-ijmin+1

           ! allocate loc data
           ALLOCATE(tabtestSW(npiloc,npjloc), tabtestSE(npiloc,npjloc), tabtestNE(npiloc,npjloc), tabtestNW(npiloc,npjloc))
           ALLOCATE(dfloloc(npiloc,npjloc), dflaloc(npiloc,npjloc))
           ALLOCATE(mask_box(npiloc,npjloc), mask_oce(npiloc,npjloc))
           ALLOCATE(vardeploc(npiloc,npjloc))
           dfloloc(:,:)=dflo(iimin:iimax,ijmin:ijmax)
           dflaloc(:,:)=dfla(iimin:iimax,ijmin:ijmax)
           vardeploc(:,:)=vardep(iimin:iimax,ijmin:ijmax)

           ! Define the 4 'equations de droite' to remove all data not in the cell
           tabtestSW=  (glamf(ji  ,jj-1)-glamf(ji-1,jj-1))*(dflaloc-gphif(ji-1,jj-1)) & 
                &     - (gphif(ji  ,jj-1)-gphif(ji-1,jj-1))*(dfloloc-glamf(ji-1,jj-1))
           tabtestSE=  (glamf(ji  ,jj  )-glamf(ji  ,jj-1))*(dflaloc-gphif(ji  ,jj-1)) &
                &     - (gphif(ji  ,jj  )-gphif(ji  ,jj-1))*(dfloloc-glamf(ji  ,jj-1))
           tabtestNE=  (glamf(ji-1,jj  )-glamf(ji  ,jj  ))*(dflaloc-gphif(ji  ,jj  )) &
                - (gphif(ji-1,jj  )-gphif(ji  ,jj  ))*(dfloloc-glamf(ji  ,jj  ))
           tabtestNW=  (glamf(ji-1,jj-1)-glamf(ji-1,jj  ))*(dflaloc-gphif(ji-1,jj  )) &
                - (gphif(ji-1,jj-1)-gphif(ji-1,jj  ))*(dfloloc-glamf(ji-1,jj  ))

           !--------------------------------------------------------------------------------------- 
           !  NOW we have the original bathymetry in the domain containing the ORCA cell.
           !  We have the 4 line equations to find which points are inside or outside the cell
           !  Find those points and do interpolation
           !--------------------------------------------------------------------------------------- 
           mask_box(:,:)=0
           WHERE ( (tabtestSW >=  0.) .AND.  &
                (tabtestSE   >=  0.) .AND.  &
                (tabtestNE   >=  0.) .AND.  &
                (tabtestNW   >=  0.)) 
              mask_box = 1
           END WHERE
           mpoi(ji,jj) = MAX(1,SUM(mask_box))  ! to avoid division by 0 (0/0)

           ! remove all the data not in the cell and on land
           IF (ln_sign ) vardeploc = - vardeploc
           mask_oce(:,:)=1
           WHERE (vardeploc <= -9000)
              mask_oce = 0
           END WHERE
           mask_oce = mask_oce * mask_box
           moce(ji,jj) = SUM(mask_oce)

           ALLOCATE(vardep1d(moce(ji,jj)))
           IF (FLOAT(moce(ji,jj))/FLOAT(mpoi(ji,jj)) > 0.5) THEN
              IF (nn_interp == 0 ) THEN
                 !                Arithmetic average
                 zdep(ji,jj) = SUM (vardeploc*mask_oce)/float(moce(ji,jj))
              ELSE
                 !                Median average 
                 ! If the compiler does not like the RESHAPE statement
                 ! use the commented code instead
                 vardep1d=RESHAPE(vardeploc, (/ npiloc*npjloc /) )
!                ik = 0
!                DO jii=1,npiloc
!                   DO jij=1,npjloc
!                      IF (mask_oce(jii,jij) == 1 ) THEN 
!                         ik = ik + 1
!                         vardep1d(ik) = vardeploc(jii,jij) 
!                      END IF
!                   END DO
!                END DO
                 zdep(ji,jj) = median(vardep1d,moce(ji,jj))
              ENDIF
           END IF
           DEALLOCATE(vardep1d)
           DEALLOCATE(tabtestSW, tabtestSE, tabtestNE, tabtestNW)
           DEALLOCATE(vardeploc, dfloloc, dflaloc, mask_oce, mask_box)
           !
        END IF
        ! --------------------------------  END OF LOOP ON JI, JJ
     ENDDO
  ENDDO
  !-----------------------------------------------------------------

  !--------------------------------------------------------------
  ! Boundary conditions  
  !--------------------------------------------------------------
  !          the bathymetry must be positive for OPA model.
  !          if it is not, convert it to positive values
  !  zdep  = ABS(zdep)

  !
  !   apply boundary condition on  bathy  
  CALL lbc ( zdep ,npiglo,npjglo, 1, 1 ,1, nn_perio)   

  !--------------------------------------------------------------
  !  NOW deal with output file 
  !--------------------------------------------------------------
  PRINT *, 'Opening  the output file:' ,TRIM(cn_fout)
  ierr = NF90_CREATE(TRIM(cn_fout), NF90_CLOBBER, ncio)
  IF (ierr /=  NF90_NOERR) THEN
     PRINT *, NF90_STRERROR(ierr)
     PRINT *, ' '
     PRINT *, 'Could not open the output file ', TRIM(cn_fout)
     PRINT *, ' '
     STOP 'We stop here'
  ELSE
     PRINT *,'--   Output file open OK'
  ENDIF
  !
  !  Definition of dimensions.
  !
  ierr = NF90_DEF_DIM(ncio, 'x', npiglo, idx)
  ierr = NF90_DEF_DIM(ncio, 'y', npjglo, idy)
  !
  !  Definition of variables
  !
  ierr=NF90_DEF_VAR(ncio,'nav_lon',NF90_FLOAT,(/idx,idy/),idv_lon) 
  ierr=NF90_DEF_VAR(ncio,'nav_lat',NF90_FLOAT,(/idx,idy/),idv_lat) 
  !                                     ---------------  
  !                                     Dim_msk has 2 dimensions  idx and idy
  ierr=NF90_DEF_VAR(ncio,cn_varout ,NF90_FLOAT,(/idx,idy/),idv_bat) 
  ierr=NF90_DEF_VAR(ncio,'moce'       ,NF90_INT,(/idx,idy/),idv_moce)
  ierr=NF90_DEF_VAR(ncio,'mpoi'       ,NF90_INT,(/idx,idy/),idv_mpoi)
  ierr=NF90_ENDDEF(ncio)
  !
  !  Save variables
  !
  ierr=NF90_PUT_VAR(ncio, idv_lon, glamt)
  ierr=NF90_PUT_VAR(ncio, idv_lat, gphit)
  !
  ierr=NF90_PUT_VAR(ncio, idv_bat, zdep  )
  ierr=NF90_PUT_VAR(ncio, idv_moce, moce  )
  ierr=NF90_PUT_VAR(ncio, idv_mpoi, mpoi  )

  !
  !  Close files
  !
  ierr=NF90_CLOSE(ncio)
  CLOSE(numnam)

CONTAINS

  SUBROUTINE PrintTplNamelist( cd_template)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PrintTplNamelist ***
    !!
    !! ** Purpose :  Print a template of the namelist used in this program
    !!               to be edited before use. 
    !!
    !! ** Method  :  Just write statement 
    !!
    !!----------------------------------------------------------------------
    CHARACTER(LEN=*), INTENT(in) :: cd_template
    
    INTEGER(KIND=4) :: inum =36
    !!----------------------------------------------------------------------
    OPEN(inum, file= cd_template )
    WRITE(inum,'(a)') '! Namelist file to be used by batinterp.exe'
    WRITE(inum,'(a)') '! This file was produced by the command batinterp.exe -p '//TRIM(cd_template)
    WRITE(inum,'(a)') '!    EDIT this file to fit your needs/wish before use.'
    WRITE(inum,'(a)') '  '
    WRITE(inum,'(a)') '&naminterpo'
    WRITE(inum,'(a)') '   nn_interp = 1       ! interpolation method : O: Arithmetic Average'
    WRITE(inum,'(a)') '                       !                        1: Median average'
    WRITE(inum,'(a)') '   nn_perio  = 0       ! NEMO periodic conditions (99 : automatic inference)'
    WRITE(inum,'(a)') ''
    WRITE(inum,'(a)') '   ! NEMO bathymetry output file'
    WRITE(inum,'(a)') "   cn_fout   = 'bathy_meter.nc'   ! name of the output filer"
    WRITE(inum,'(a)') "   cn_varout = 'Bathymetry'       ! name of output variable"
    WRITE(inum,'(a)') ''
    WRITE(inum,'(a)') '   ! NEMO coordinate file'
    WRITE(inum,'(a)') "   cn_fgrid  = 'mesh_hgr.nc'      ! name of horizontal grid file"
    WRITE(inum,'(a)') ''
    WRITE(inum,'(a)') '   ! External Bathymetric file'
    WRITE(inum,'(a)') "   cn_fbatin = 'BedmachineGreenland...'  ! name of external baty file"
    WRITE(inum,'(a)') '   ln_regin  = .FALSE.  ! True  : Regular external bathy file'
    WRITE(inum,'(a)') '                        ! False :Irregular external bathy file'
    WRITE(inum,'(a)') "   cn_varin  = 'bed'    ! name of bathymetry in file"
    WRITE(inum,'(a)') "   cn_xdim   = 'x'      ! name of X dimension"
    WRITE(inum,'(a)') "   cn_ydim   = 'y'      ! name of Y dimension"
    WRITE(inum,'(a)') "   cn_lonv   = 'lon'    ! name of longitude variable"
    WRITE(inum,'(a)') "   cn_latv   = 'lat'    ! name of latitude variable"
    WRITE(inum,'(a)') ''
    WRITE(inum,'(a)') '   ! change sign  of bathymetry'
    WRITE(inum,'(a)') '   ln_sign   = .FALSE.  ! change sign for bathymetry (ISF related)'
    WRITE(inum,'(a)') '/'
    CLOSE(inum)
  
  END SUBROUTINE PrintTplNamelist
  !-------------------------------------------------------------------
  SUBROUTINE hdlerr(istatus)
    !--------------------------------------------------------------------
    INTEGER(KIND=4)                 :: istatus
    IF (istatus /=  NF90_NOERR) THEN
       PRINT *, NF90_STRERROR(istatus)
       STOP 'stopped here'
    ENDIF
    RETURN
  END SUBROUTINE hdlerr

  !-------------------------------------------------------------------
  SUBROUTINE prihre( chtab, tab,nx,ny)
    !-------------------------------------------------------------------
    INTEGER(KIND=4)          :: nx, ny
    REAL, DIMENSION(nx,ny) :: tab
    CHARACTER (len=10) :: chtab
    PRINT *, '-------------------', chtab, '--------------------'
    DO jj= ny,1,-1
       PRINT *, 'jj=',jj
       PRINT 100, tab(:,jj)
    ENDDO
100 FORMAT (10(1x,1pg10.3))
  END SUBROUTINE prihre

  SUBROUTINE NearestPoint(ddlon, ddlat, kpi, kpj, ddlam, ddphi, kpiloc, kpjloc)
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE NearestPoint  ***
    !!
    !! ** Purpose : Computes the positions of the nearest i,j in the grid
    !!              from the given longitudes and latitudes
    !!
    !! ** Method  : Starts on the middle of the grid, search in a 20x20 box,
    !!              and move the box in the direction where the distance 
    !!              between the box and the point is minimum.
    !!              Iterates ...
    !!              Stops when the point is outside the grid.
    !!
    !! References : P.A. Darbon and A. de Miranda acknowledged for this 
    !!              clever algorithm developped in CLIPPER.
    !!----------------------------------------------------------------------
    REAL(KIND=8),                     INTENT(in) :: ddlon, ddlat        !: lon and lat of target point
    INTEGER(KIND=4),                 INTENT (in) :: kpi, kpj          !: grid size
    REAL(KIND=8), DIMENSION(kpi,kpj), INTENT(in) :: ddlam, ddphi        !: model grid layout
    INTEGER(KIND=4),              INTENT (inout) :: kpiloc, kpjloc    !: nearest point location
    LOGICAL                                      :: ld_bnd            !: reach boundary flag

    INTEGER(KIND=4)                              :: ji, jj
    INTEGER(KIND=4), PARAMETER                   :: jp_blk=10
    INTEGER(KIND=4)                              :: ii0, ij0
    INTEGER(KIND=4)                              :: ii1, ij1, kpiloctry1, kpjloctry1, nreachbnd
    REAL(KIND=4)                                 :: zdist
    REAL(KIND=4)                                 :: zdistmin, zdistmin0, disttry1
    LOGICAL, SAVE                                :: ll_bndcell, ll_first=.TRUE.
    !!----------------------------------------------------------------------
    !   IF ( ll_first ) THEN
    !      kpiloc = kpi/2 ; kpjloc = kpj/2    ! seek from the middle of domain
    !      ll_first=.FALSE.
    !   ENDIF

    zdistmin=1000000. ; zdistmin0=1000000.
    ii0 = kpiloc      ; ij0 = kpjloc
    ll_bndcell=.TRUE. ; ld_bnd=.FALSE. ; nreachbnd=1

    ! loop until found or boundary reach
    DO  WHILE ( ll_bndcell .AND. .NOT. ld_bnd )
       ii0 = kpiloc - jp_blk ;  ii1 = kpiloc + jp_blk
       ij0 = kpjloc - jp_blk ;  ij1 = kpjloc + jp_blk

       ! search only the inner domain
       IF (ii0 <= 0 ) ii0 = 2
       IF (ii1 > kpi) ii1 = kpi - 1
       IF (ij0 <= 0 ) ij0 = 2
       IF( ij1 > kpj) ij1 = kpj - 1

       ! within a block jp_blk+1 x jp_blk+1:
       DO jj=ij0,ij1
          DO ji=ii0,ii1
             ! compute true distance (orthodromy) between target point and grid point
             zdist    = dist(ddlon, ddlam(ji,jj), ddlat, ddphi(ji,jj) )
             zdistmin = MIN(zdistmin, zdist)
             ! update kpiloc, kpjloc if distance decreases
             IF (zdistmin /=  zdistmin0 ) THEN
                kpiloc=ji
                kpjloc=jj
             ENDIF
             zdistmin0=zdistmin
          END DO
       END DO

       ll_bndcell=.FALSE.
       ! if kpiloc, kpjloc belong to block boundary proceed to next block, centered on kpiloc, kpjloc
       IF (kpiloc == ii0 .OR. kpiloc == ii1) ll_bndcell=.TRUE.
       IF (kpjloc == ij0 .OR. kpjloc == ij1) ll_bndcell=.TRUE.

       ! boundary reach ---> not found
       IF (kpiloc == 2  .OR. kpiloc ==kpi-1) ld_bnd=.TRUE.
       IF (kpjloc == 2  .OR. kpjloc ==kpj-1) ld_bnd=.TRUE.
       IF (ld_bnd == .TRUE. .AND. nreachbnd == 1) THEN
          disttry1 = zdistmin0 ; kpiloctry1=kpiloc ; kpjloctry1=kpjloc
          ! E/W periodicity :
          IF (kpiloc == 2    ) kpiloc = kpi-4
          IF (kpiloc == kpi-1) kpiloc = 4
          nreachbnd=2
          !PRINT *, 'reach bnd 1st time'
          ld_bnd=.FALSE.
       END IF
    END DO
    IF (nreachbnd == 2) THEN
       !PRINT *, ' cmp dist ', zdistmin0, disttry1, kpiloc, kpiloctry1, kpjloc, kpjloctry1
       IF (zdistmin0 >= disttry1) THEN
          kpiloc = kpiloctry1
          kpjloc = kpjloctry1
       END IF
    END IF

  END SUBROUTINE  NEARESTPOINT

  REAL(KIND=8) FUNCTION dist(ddlona, ddlonb, ddlata, ddlatb)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION dist  ***
    !!
    !! ** Purpose : Compute the distance (km) between
    !!              point A (lona, lata) and B (lonb, latb)  
    !!
    !! ** Method  : Use of double precision is important. Compute the 
    !!              distance along the orthodromy
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=8), INTENT(in) :: ddlata, ddlona, ddlatb, ddlonb

    REAL(KIND=8), SAVE :: dl_latar, dl_latbr, dl_lonar, dl_lonbr
    REAL(KIND=8)       :: dl_pds
    REAL(KIND=8), SAVE :: dl_ux, dl_uy, dl_uz
    REAL(KIND=8)       :: dl_vx, dl_vy, dl_vz
    REAL(KIND=8), SAVE :: dl_prevlat=-1000.d0
    REAL(KIND=8), SAVE :: dl_prevlon=-1000.d0
    REAL(KIND=8), SAVE :: dl_r, dl_pi, dl_conv

    LOGICAL :: ll_first=.TRUE.
    !!----------------------------------------------------------------------
    ! initialise some values at first call
    IF ( ll_first ) THEN
       ll_first = .FALSE.
       ! constants
       dl_pi   = ACOS(-1.d0)
       dl_conv = dl_pi/180.d0  ! for degree to radian conversion
       ! Earth radius
       dl_r    = (6378.137d0+6356.7523d0)/2.0d0 ! km
    ENDIF

    ! compute these term only if they differ from previous call
    IF ( ddlata /= dl_prevlat .OR. ddlona /= dl_prevlon) THEN
       dl_latar   = ddlata*dl_conv
       dl_lonar   = ddlona*dl_conv
       dl_ux      = COS(dl_lonar)*COS(dl_latar)
       dl_uy      = SIN(dl_lonar)*COS(dl_latar)
       dl_uz      = SIN(dl_latar)
       dl_prevlat = ddlata
       dl_prevlon = ddlona
    ENDIF

    dl_latbr = ddlatb*dl_conv
    dl_lonbr = ddlonb*dl_conv
    dl_vx    = COS(dl_lonbr)*COS(dl_latbr)
    dl_vy    = SIN(dl_lonbr)*COS(dl_latbr)
    dl_vz    = SIN(dl_latbr)

    dl_pds   = dl_ux*dl_vx + dl_uy*dl_vy + dl_uz*dl_vz

    IF (dl_pds >= 1.) THEN
       dist = 0.
    ELSE
       dist = dl_r*ACOS(dl_pds)
    ENDIF

  END FUNCTION dist


END PROGRAM batinterp

