PROGRAM tread_histo
  USE netcdf

  IMPLICIT NONE

  INTEGER(KIND=4) :: ji
  INTEGER(KIND=4) :: il
  INTEGER(KIND=4) :: inum , nl, nlskip

  INTEGER(KIND=4), DIMENSION(:), ALLOCATABLE :: ni,nj

  REAL(KIND=4),    DIMENSION(:), ALLOCATABLE :: v0,v1
  CHARACTER(LEN=80) :: cf_histo = "eGREENLAND025.L75_bathy_meter_002.5_history"
  CHARACTER(LEN=80) :: cldum
  CHARACTER(LEN=256) :: cline
  !---------------------------------------------------------------------------

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
  DO ji = 1, nl
    print '( 2i5, 2f8.1)', ni(ji),nj(ji),v0(ji), v1(ji)
  ENDDO



END PROGRAM tread_histo


