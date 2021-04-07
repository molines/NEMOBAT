PROGRAM apply_history
  !!======================================================================
  !!                     ***  PROGRAM  apply_history  ***
  !!=====================================================================
  !!  ** Purpose : Automatically apply the hand corrections  done with
  !!              BMGtools, saved in an history file. 
  !!
  !!  ** Method  : takes both original bathy file and history file to
  !!              produce a corrected file (similar to the one created
  !!              by BMGtools).
  !!
  !! History :  1.0  : 04/2021  : J.M. Molines : 
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!----------------------------------------------------------------------


  !!----------------------------------------------------------------------
  USE netcdf
  !!----------------------------------------------------------------------
  !! CDFTOOLS_4.0 , MEOM 2017
  !! $Id$
  !! Copyright (c) 2012, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)

  IMPLICIT NONE

  INTEGER(KIND=4) :: ji
  INTEGER(KIND=4) :: il, ii, ij
  INTEGER(KIND=4) :: inum , nl, nlskip
  INTEGER(KIND=4) :: npiglo, npjglo
  INTEGER(KIND=4) :: narg, ijarg, iargc
  INTEGER(KIND=4) :: ncid, id, ierr

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ni,nj

  REAL(KIND=4),    DIMENSION(:)  , ALLOCATABLE :: v0,v1
  REAL(KIND=4),    DIMENSION(:,:), ALLOCATABLE :: bathy
  CHARACTER(LEN=80) :: cf_batio 
  CHARACTER(LEN=80) :: cf_histo = "eGREENLAND025.L75_bathy_meter_002.5_history"
  CHARACTER(LEN=80) :: cf_batic 
  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=256) :: cline
  !---------------------------------------------------------------------------
  narg = iargc()
  IF ( narg == 0 ) THEN
     PRINT *,' usage :  apply_history -b  ORIGINAL-BATHY-file -h HISTORY-file '
     PRINT *,'          -o CORRECTED-BATHY-file'
     PRINT *,'      '
     PRINT *,'     PURPOSE :'
     PRINT *,'       Apply the modifications done on the original bathy file using'
     PRINT *,'       BMGtools  from the history file.'
     PRINT *,'      '
     PRINT *,'     ARGUMENTS :'
     PRINT *,'       -b ORIGINAL-BATHY-file : give the name of the original bathy file.' 
     PRINT *,'       -h HISTORY-file : give the name of the history file.' 
     PRINT *,'       -o CORRECTED-BATHY-file : give the name of the corrected bathy file.' 
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
     CASE ( '-o'   ) ; CALL getarg(ijarg, cf_batic ) ; ijarg=ijarg+1
     CASE ( '-h'   ) ; CALL getarg(ijarg, cf_histo ) ; ijarg=ijarg+1
     CASE DEFAULT    ; PRINT *, ' ERROR : ', TRIM(cldum),' : unknown option.'; STOP 1
     END SELECT
  ENDDO

  cldum="cp "//TRIM(cf_batio)//" "//TRIM(cf_batic)
  CALL SYSTEM( cldum )
  ! work on the copy
  ierr = NF90_OPEN(cf_batic,NF90_WRITE,ncid)
  ierr = NF90_INQ_DIMID( ncid,"x",id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npiglo)
  ierr = NF90_INQ_DIMID( ncid,"y",id) ; ierr = NF90_INQUIRE_DIMENSION(ncid,id,len=npjglo)

  ALLOCATE( bathy(npiglo,npjglo) )
  ierr = NF90_INQ_VARID( ncid,"Bathymetry",id) ; ierr= NF90_GET_VAR(ncid,id,bathy)
  ! keep id it will be used for write

  CALL ReadHisto
  DO ji =1 , nl
   ii=ni(ji) ; ij=nj(ji)
   IF ( ABS( bathy(ii,ij) - v0(ji)) > 1 ) THEN
     PRINT *, 'Old value different from original file ', ii, ij
   ELSE
     bathy(ii,ij) = v1(ji)
   ENDIF
  ENDDO
  ierr = NF90_PUT_VAR(ncid,id,bathy)
  ierr = NF90_CLOSE(ncid)


 CONTAINS
 SUBROUTINE ReadHisto
    !!---------------------------------------------------------------------
    !!                  ***  ROUTINE ReadHisto ***
    !!
    !! ** Purpose :   Read history file from BMGtools and
    !!                fill in ni, nj, v0 and v1
    !!
    !! ** Method  :  read ascii file 
    !!
    !!----------------------------------------------------------------------

  nl=0 ; nlskip=0

  OPEN(inum, file=cf_histo)

  DO 
     READ(inum,'(a)',END=99) cline

     IF ( cline(1:1) == '#') THEN
        nlskip = nlskip + 1
     ELSE
        nl=nl+1
     ENDIF

  ENDDO
99 CONTINUE
  PRINT *, 'NL : ', nl
  ALLOCATE( ni(nl), nj(nl), v0(nl), v1(nl))
  REWIND( inum)
  il = 0
  DO
     READ(inum,'(a)',END=98) cline

     IF ( cline(1:1) /= '#') THEN
        il=il+1
        READ(cline,*) cldum,cldum,cldum,cldum,cldum,ni(il), cldum,cldum,nj(il),cldum,cldum, v0(il), cldum, v1(il)
     ENDIF
  ENDDO
98 CONTINUE
!  DO ji = 1, nl
!    print '( 2i5, 2f8.1)', ni(ji),nj(ji),v0(ji), v1(ji)
!  ENDDO

 END SUBROUTINE ReadHisto



END PROGRAM apply_history


