PROGRAM revertetopo2
 !-----------------------------------------------------------------------------
 !         ***   PROGRAM revertetopo2  ***
 !  Purpose : rewrite etopo2 nc file so that (1,1) corresponds to SW corner
 !
 !  Method: read an rewrite !
 !          put dimension lon lat and add corresponding variable
 !          layout is -180W TO 180 E, ni=o overlap
 !
 ! history : original : J.M. Molines 25/05/2009
 !-----------------------------------------------------------------------------
  !! -------------------------------------------------------------------------
  !!  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
  !!  $Rev: 228 $
  !!  $Id: revertetopo2.f90 228 2009-04-29 14:19:25Z forge $
  !! -------------------------------------------------------------------------
 USE netcdf
 IMPLICIT NONE
 INTEGER :: ncid, idvar, idx, idy, idlon, idlat, istatus
 INTEGER :: ni,nj
 INTEGER :: ji,jj

 INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: idept
 REAL(KIND=4),    DIMENSION(:),   ALLOCATABLE :: rlon, rlat
 REAL(KIND=8) :: dx

 CHARACTER(LEN=80) :: cfilin='etopo2bedmap.nc' ! hard coded because unique !
 CHARACTER(LEN=80) :: cfilout='reformatted.nc' ! user will change this name later

 dx=1./30.d0

 istatus = NF90_OPEN(cfilin,NF90_NOWRITE,ncid)
 istatus = NF90_INQ_DIMID(ncid,'jpio',idx)
 istatus = NF90_INQUIRE_DIMENSION(ncid,idx,len=ni)
 istatus = NF90_INQ_DIMID(ncid,'jpjo',idy)
 istatus = NF90_INQUIRE_DIMENSION(ncid,idy,len=nj)
 ALLOCATE ( rlon(ni), rlat(nj), idept(ni,nj) )

 istatus = NF90_INQ_VARID(ncid,'depth',idvar)
 istatus = NF90_GET_VAR(ncid,idvar,idept)
 istatus = NF90_CLOSE(ncid)

 ! create rlon, rlat
 DO ji=1,ni
  rlon(ji)=-180.+dx/2. + (ji-1) * dx
 END DO
 DO jj=1,nj
  rlat(jj)=-90.+dx/2. + (jj-1) * dx
 ENDDO

 PRINT *, 'RLON RANGE :', rlon(1),' ---> ',rlon(ni)
 PRINT *, 'RLAT RANGE :', rlat(1),' ---> ',rlat(nj)
 
 ! Create output file
 istatus = NF90_CREATE(cfilout,NF90_CLOBBER,ncid)
 istatus = NF90_DEF_DIM(ncid,'lon',ni,idx)
 istatus = NF90_DEF_DIM(ncid,'lat',nj,idy)
 istatus = NF90_DEF_VAR(ncid,'lon',NF90_FLOAT,(/idx/),idlon)
 istatus = NF90_DEF_VAR(ncid,'lat',NF90_FLOAT,(/idy/),idlat)
 istatus = NF90_DEF_VAR(ncid,'z',NF90_SHORT,(/idx,idy/),idvar)
 istatus = NF90_ENDDEF(ncid)

 ! put vars
 istatus = NF90_PUT_VAR(ncid,idlon,rlon)
 istatus = NF90_PUT_VAR(ncid,idlat,rlat)
 istatus = NF90_PUT_VAR(ncid,idvar,idept(:,nj:1:-1))
 istatus = NF90_CLOSE(ncid)
 
END PROGRAM revertetopo2
