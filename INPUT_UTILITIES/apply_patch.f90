PROGRAM apply_patch
  !!======================================================================
  !!                     ***  PROGRAM  apply_patch  ***
  !!=====================================================================
  !!  ** Purpose : Merge a subdomain file into a wider file
  !!
  !!  ** Method  : Read the full domain, read the subdomain and use the 
  !!               the window coordinates passed in argument for patching
  !!
  !! History :  1.0  : 04/2021  : J.M. Molines : 
  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  !! $Id$
  !! Copyright (c) 2021, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER(KIND=4) :: iimin, iimax, ijmin, ijmax
  INTEGER(KIND=4) :: npiglo, npjglo
  INTEGER(KIND=4) :: npiloc, npjloc
  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: ncid, id, ierr, idv, ncip


  REAL(KIND=4),    DIMENSION(:,:), ALLOCATABLE :: batp
  REAL(KIND=4),    DIMENSION(:,:), ALLOCATABLE :: bathy
  CHARACTER(LEN=256) :: cf_batio 
  CHARACTER(LEN=256) :: cf_patch
  CHARACTER(LEN=256) :: cf_batic 
  CHARACTER(LEN=256) :: cldum
  CHARACTER(LEN=256) :: cline
  !---------------------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  apply_patch -b  ORIGINAL-BATHY-file -p PATCH-file '
     PRINT *,'          -o FINAL-BATHY-file -w imin imax jmin jmax'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Merge the subdomain represented in the patch file into'
     PRINT *,'       the ORIGINAL-BATHY-file using the windows coordinates.'
     PRINT *,'       '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -b ORIGINAL-BATHY-file : give the name of the original bathy file.' 
     PRINT *,'       -p PATCH-file : give the name of the patch file.' 
     PRINT *,'       -o FINAL-BATHY-file : give the name of the final (patched) bathy file.' 
     PRINT *,'       -w imin imax jmin jmax : give the I-J window of the patch into the '
     PRINT *,'                 ORIGINAL-BATHY-file. Note that for this program, array index'
     PRINT *,'                 start at 1 (double check when translating from ncks).'
     PRINT *,'      '
     PRINT *,'     OPTIONS :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     REQUIRED FILES :'
     PRINT *,'       none' 
     PRINT *,'      '
     PRINT *,'     OUTPUT : '
     PRINT *,'       netcdf file : specified by -o option'
     PRINT *,'         variables : Bathymetry'
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
     CASE ( '-b'   ) ; CALL getarg(ijarg, cf_batio ) ; ijarg=ijarg+1
     CASE ( '-p'   ) ; CALL getarg(ijarg, cf_patch ) ; ijarg=ijarg+1
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_batic ) ; ijarg=ijarg+1
     CASE ( '-w'   ) ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1 ; READ(cldum,*) iimin
          ;          ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1 ; READ(cldum,*) iimax
          ;          ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmin
          ;          ; CALL getarg(ijarg, cldum    ) ; ijarg=ijarg+1 ; READ(cldum,*) ijmax
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  ! copy original file to patched file 
  cldum="cp "//TRIM(cf_batio)//" "//TRIM(cf_batic)
  CALL SYSTEM( cldum )

  ! work on the copy
  ierr = NF90_OPEN(cf_batic,NF90_WRITE,ncid)
  ierr = NF90_INQ_DIMID( ncid,"x",id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr = NF90_INQ_DIMID( ncid,"y",id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)

  ALLOCATE( bathy(npiglo,npjglo) )
  ! read original bathy
  ierr = NF90_INQ_VARID( ncid,"Bathymetry",id) ; ierr= NF90_GET_VAR(ncid,id,bathy)
  ! keep id it will be used for write
  idv = id

  ! Read patch
  ierr = NF90_OPEN(cf_patch,NF90_NOWRITE,ncip)
  ierr = NF90_INQ_DIMID( ncip,"x",id) ; ierr = NF90_INQUIRE_DIMENSION( ncip,id,len=npiloc)
  ierr = NF90_INQ_DIMID( ncip,"y",id) ; ierr = NF90_INQUIRE_DIMENSION( ncip,id,len=npjloc)

  ALLOCATE( batp(npiloc,npjloc) )
  ierr = NF90_INQ_VARID( ncip,"Bathymetry",id) ; ierr= NF90_GET_VAR(ncip,id,batp)
  ierr = NF90_CLOSE(ncip)

  ! some quick check on the size of the patch
  ierr =0
  IF ( iimax - iimin +1 /= npiloc ) THEN
     PRINT *,' ERROR : inconsistent size of patch in I direction vs imin imax'
     ierr = 1
  ENDIF
  IF ( ijmax - ijmin +1 /= npjloc ) THEN
     PRINT *,' ERROR : inconsistent size of patch in J direction vs jmin jmax'
     ierr = ierr +1
  ENDIF
  IF ( ierr /= 0 ) STOP

  ! patch the original bathy
  bathy(iimin:iimax, ijmin:ijmax) = batp(:,:)

  ! write it back to the file
  ierr = NF90_PUT_VAR(ncid,idv,bathy)
  ierr = NF90_CLOSE(ncid)

END PROGRAM apply_patch


