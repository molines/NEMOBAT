&naminterpo
   nn_interp = 1       ! interpolation method : O: Arithmetic Average
                       !                        1: Median average
   nn_perio  = 0       ! NEMO periodic conditions (99 : automatic inference)

   cn_fout   = 'bathy_meter.nc'   ! name of the output file
   cn_varout = 'Bathymetry'       ! name of output variable

   cn_fgrid  = 'mesh_hgr.nc'      ! name of horizontal grid file

   cn_batin  = 'BedmachineGreenland...'  ! name of external baty file
   cn_varin  = 'bed'    ! name of bathymetry in file
   cn_xdim   = 'x'      ! name of X dimension
   cn_ydim   = 'y'      ! name of Y dimension

   ln_sign   = .FALSE.  ! change sign for bathymetry (ISF related)
   ln_regin  = .FALSE.  ! True  : Regular external bathy file
                        ! False :Irregular external bathy file
/
