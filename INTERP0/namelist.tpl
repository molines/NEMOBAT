! Namelist file to be used by batinterp.exe
! This file was produced by the command batinterp.exe -p namelist.tpl
!    EDIT this file to fit your needs/wish before use.
  
&naminterpo
   nn_interp = 1       ! interpolation method : O: Arithmetic Average
                       !                        1: Median average
                       !                        2: Median average modified
                       ! In case of even sorted list do not use a mean if one
                       ! of the 2 central values is 0.
   nn_perio  = 0       ! NEMO periodic conditions (99 : automatic inference)

   ! NEMO bathymetry output file
   cn_fout   = 'bathy_meter.nc'   ! name of the output filer
   cn_varout = 'Bathymetry'       ! name of output variable
   ln_time   = .FALSE.            ! Creation of a time dim/var

   ! NEMO coordinate file
   cn_fgrid  = 'mesh_hgr.nc'      ! name of horizontal grid file

   ! External Bathymetric file
   cn_fbatin = 'BedmachineGreenland...'  ! name of external baty file
   ln_regin  = .FALSE.  ! True  : Regular external bathy file
                        ! False :Irregular external bathy file
   cn_varin  = 'bed'    ! name of bathymetry in file
   cn_xdim   = 'x'      ! name of X dimension
   cn_ydim   = 'y'      ! name of Y dimension
   cn_lonv   = 'lon'    ! name of longitude variable
   cn_latv   = 'lat'    ! name of latitude variable

   ! change sign  of bathymetry
   ln_sign   = .FALSE.  ! change sign for bathymetry (ISF related)
/
