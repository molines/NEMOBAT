#!/bin/ksh

# mkbat script to build a bathymetry for a coordinates.nc file, from gebco1.nc and etopo1.nc
# this can be used as template but requires carefull checking to fit your own requirements

#  $Date: 2009-04-29 22:55:31 +0200 (Wed, 29 Apr 2009) $
#  $Rev: 234 $
#  $Id: mkbat.ksh 234 2009-04-29 20:55:31Z forge $

CONFIG=ORCA025
IDIR=/home/users/molines/DATA/ORCA025-I
MSK_ID=${CONFIG}-G70
GEBCO1=/home/users/molines/DATA/BATHY/GEBCO_1min.nc
ETOPO1=/home/users/molines/DATA/BATHY/ETOPO1_Bed_c_gmt4.grd.nc
COORDINATES=/home/users/molines/DATA/ORCA025-I/coordinates_ORCA_R025_lombok+ombai_v2.nc

# make covenient links
ln -sf $COORDINATES coordinates.nc
ln -sf $GEBCO1 gebco1.nc
ln -sf $ETOPO1 etopo1.nc

#(1) INTERP0
#-----------
interp_gebco1() { ln -sf gebco1.nc bathy_in.nc ;
             ../INTERP0/batinterp1 
             mv bathy_out.nc ${CONFIG}_bathy_gebco1.nc ;}

interp_etopo1() { ln -sf etopo1.nc bathy_in.nc ;
             ../INTERP0/batinterp2
             mv bathy_out.nc ${CONFIG}_bathy_etopo1.nc ;}

if [ ! -f ${CONFIG}_bathy_gebco1.nc ] ; then interp_gebco1 ; fi
if [ ! -f ${CONFIG}_bathy_etopo1.nc ] ; then interp_etopo1 ; fi

#(2) Smoothing etopo1
#--------------------
smooth_bathy() { ../SMOOTHING/batsmooth $1 ; }

  n_pass_shapiro=2   #  use 2 shapiro w0.6 in this case
  smooth_bathy ${CONFIG}_bathy_etopo1.nc $n_pass_shapiro  > smooth.log 
  etopo_smooth=$( cat smooth.log | grep SMOOTHED | awk '{print $NF}' )

#(3) Merge gebco1 and etop1 on model grid
#----------------------------------------

merge_bathy() { ../COMBINE/combine $1 $2 ; }

  merge_bathy ${CONFIG}_bathy_gebco1.nc  $etopo_smooth
  mv  combine_etopo1_gebco1.nc ${CONFIG}_bathy_etopo1_gebco1_smoothed.nc

#(4) Apply coast line
#---------------------
apply_coast () { ../APPLYCOAST/batcoast -if $1 -of $2 -msk $3 ; }

  if [ ! -f mask.nc ] ; then ln -s $IDIR/${MSK_ID}_byte_mask.nc mask.nc ; fi
  
  if [ ! -f ${CONFIG}_bathy_etopo1_gebco1_smoothed_coast.nc ] ; then
   apply_coast ${CONFIG}_bathy_etopo1_gebco1_smoothed.nc ${CONFIG}_bathy_etopo1_gebco1_smoothed_coast.nc mask.nc
  fi

#(5) Apply hand correction saves in h90 file
#---------------------------------------
  here=$(pwd)
  cd ../CORRECT
  # check that your correction are in this file (and name of output file)
  ln -sf H90/orca025_hand_modif_april09.h90 hand_modif.h90
  make correct
  cd $here
  ../CORRECT/correct ${CONFIG}_bathy_etopo1_gebco1_smoothed_coast.nc ${CONFIG}_bathy_etopo1_gebco1_smoothed_coast_corrected_apr09.nc
