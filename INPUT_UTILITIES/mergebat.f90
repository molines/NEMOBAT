PROGRAM mergebat
  !!======================================================================
  !!                     ***  PROGRAM  mergebat  ***
  !!=====================================================================
  !!  ** Purpose : Produce a merged bathymetry between 2 files, using
  !!               a distance as weight
  !!
  !!  ** Method  : use namelist to pass the name of the files
  !!
  !! History :  4.0  : 03/2017  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------
   USE netcdf
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------
  IMPLICIT NONE 
  INTEGER(KIND=4) :: ji,jj
  INTEGER(KIND=4) :: npiglo, npjglo
  INTEGER(KIND=4) :: narg, ijarg, iargc
  

  INTEGER(KIND=4) :: ierr, ncid, id

  REAL(KIND=4)    :: alfa
  REAL(KIND=4), DIMENSION(:,:), ALLOCATABLE ::  Vorig, Vpatch, Vfinal, tdist

  CHARACTER(LEN=80)  :: cmd
  CHARACTER(LEN=80)  :: cldum

  LOGICAL            :: l_print = .FALSE.

  ! Namelist
  INTEGER(KIND=4)    :: numnam
  CHARACTER(LEN=80)  :: cf_nam

  REAL(KIND=4)       :: rn_limit=80 ! km
  CHARACTER(LEN=80)  :: cn_orig="eGREENLAND025.L75_bathymetry.nc"
  CHARACTER(LEN=80)  :: cn_patch="eGREENLAND025.L75_bathy_meter_002.nc2"
  CHARACTER(LEN=80)  :: cn_final="eGREENLAND025.L75_bathymetry_patched.nc"
  CHARACTER(LEN=80)  :: cn_dist="dist.coast"
  
  CHARACTER(LEN=80)  :: cn_xdim="x"
  CHARACTER(LEN=80)  :: cn_ydim="y"
  CHARACTER(LEN=80)  :: cn_vbat="Bathymetry"
  CHARACTER(LEN=80)  :: cn_vbatp="Bathymetry"
  CHARACTER(LEN=80)  :: cn_vdist="Tcoast"

  NAMELIST /nammerge/ rn_limit,cn_orig, cn_xdim, cn_ydim,cn_vbat, &
       &                       cn_patch,cn_vbatp,                 &
       &                       cn_final,                          &
       &                       cn_dist, cn_vdist
  !!----------------------------------------------------------------------
  narg=iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  mergebat -f NAMELIST [-p ]'
     PRINT *,'      '
     PRINT *,'     PURPOSE :' 
     PRINT *,'          Merge two bathymetric files using a distance file giving the'
     PRINT *,'          distance to the boundaries for the patch file (meters).' 
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -f NAMELIST : namelist file name used to describe the file names '   
     PRINT *,'            and  other information'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       -p : print a template file name for use with mergebat' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none except the input files specified in the namelist.' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : Name is specified in the namelist'
     PRINT *,'         variables : Name is specified in the Bathymetry'
     PRINT *,'      '
     PRINT *,'     SEE ALSO :'
     PRINT *,'      ' 
     PRINT *,'      '
     STOP
  ENDIF

  ijarg = 1 
  DO WHILE ( ijarg <= narg )
     CALL getarg(ijarg, cldum ) ; ijarg=ijarg+1
     SELECT CASE ( cldum )
     CASE ( '-f'   ) ; CALL getarg(ijarg, cf_nam ) ; ijarg=ijarg+1
        ! option
     CASE ( '-p'   ) ; l_print = .true.            ; ijarg = ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO


  IF (l_print ) THEN
    CALL PrintNamelist
    STOP
  ENDIF

  OPEN(numnam,file=cf_nam)
  READ(numnam,nammerge)
  CLOSE(numnam)


  cmd="cp "//TRIM(cn_orig)//" "//TRIM(cn_final)
  CALL SYSTEM(cmd)

  ierr=NF90_OPEN(cn_orig,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_DIMID(ncid,cn_xdim,id ) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr= NF90_INQ_DIMID(ncid,cn_ydim,id ) ; ierr=NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)
  PRINT *, ' NPIGLO : ', npiglo
  PRINT *, ' NPJGLO : ', npjglo
  ALLOCATE( Vorig(npiglo,npjglo), Vpatch(npiglo,npjglo), Vfinal(npiglo,npjglo), tdist(npiglo,npjglo))

  ierr= NF90_INQ_VARID(ncid,cn_vbat,id)  ; ierr = NF90_GET_VAR(ncid,id,Vorig)
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cn_patch,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,cn_vbatp,id) ; ierr = NF90_GET_VAR(ncid,id,Vpatch)
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cn_dist,NF90_NOWRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,cn_vdist,id) ; ierr = NF90_GET_VAR(ncid,id,tdist)
  tdist=tdist/1000.  ! from meters to km
  ierr = NF90_CLOSE(ncid)

  ierr=NF90_OPEN(cn_final,NF90_WRITE,ncid)
  ierr= NF90_INQ_VARID(ncid,cn_vbat,id) 
  Vfinal=Vorig
  DO jj=2,npjglo-1
     DO ji = 2, npiglo-1
       IF ( Vpatch(ji,jj) > 0 ) THEN
          IF ( tdist(ji,jj) <= rn_limit)  THEN
            alfa =  tdist(ji,jj)/rn_limit
          ELSE
            alfa = 1.
          ENDIF
          Vfinal(ji,jj) = Vpatch(ji,jj) * alfa + Vorig(ji,jj) *(1 - alfa)
       ELSE IF (Vorig(ji,jj) > 0 ) THEN
         Vfinal(ji,jj) = -50.  ! this indicates an original point  in water, now drown
       ENDIF
     ENDDO
  ENDDO
  ierr = NF90_PUT_VAR(ncid,id,Vfinal)

  ierr = NF90_CLOSE(ncid)
  
  CONTAINS
  
  SUBROUTINE PrintNamelist
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE PrintNamelist  ***
    !!
    !! ** Purpose :  Display a template of the namelist used in this program 
    !!
    !! ** Method  :  Just print text. 
    !!
    !!----------------------------------------------------------------------
    PRINT '(a)','!  Template namelist for mergebat program'
    PRINT '(a)','!  MUST BE EDITED TO FIT YOUR NEEDS'
    PRINT '(a)','!'
    PRINT '(a)','&nammerge'
    PRINT '(a)','   rn_limit  = 80    ! (km) Width of the ramp for merging'
    PRINT '(a)','! Original file : The original bathymetry file to be patched '
    PRINT '(a)','   cn_orig  = "N/A" '
    PRINT '(a)','       cn_xdim  = "x"           ! x-dimension name'
    PRINT '(a)','       cn_ydim  = "y"           ! y-dimension name'
    PRINT '(a)','       cn_vbat  = "Bathymetry"  ! name of the bathymetry variable in original file'
    PRINT '(a)',''
    PRINT '(a)','! Patch file :'
    PRINT '(a)','   cn_patch = "N/A"'
    PRINT '(a)','       cn_vbatp = "Bathymetry"  ! name of the bathymetry variable in the patch file'
    PRINT '(a)',''
    PRINT '(a)','! Output file : (from a copy of the original file)'
    PRINT '(a)','   cn_final = "N/A"'
    PRINT '(a)',''
    PRINT '(a)','! Distance file : give the distance to the boudaries'
    PRINT '(a)','   cn_dist = "N/A"'
    PRINT '(a)','       cn_vdist = "Tcoast"     ! name of the distance in the distance file'
    PRINT '(a)','/'

  
  END SUBROUTINE PrintNamelist


END PROGRAM mergebat
