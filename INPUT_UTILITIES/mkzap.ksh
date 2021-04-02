#!/bin/ksh
## -------------------------------------------------------------------------
##  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 228 $
##  $Id: mkzap.ksh 228 2009-04-29 14:19:25Z forge $
## -------------------------------------------------------------------------
DATADIR=..
GEBCO=GEBCO_1min
ETOPO1=ETOPO1_Bed_c_gmt4.grd
ETOPO2=ETOPO2_Bedmap
SM_SW=topo_11.1

ZAPIOLA_ZOOM="-52 -34 -48 -40"
INDIAN_RIDGE_ZOOM="63 110 -55 -30"
ROMANCHE_ZOOM="-20 -16 -2 2"
AMUNDSEN_ZOOM="-140 -100 -76 -68"
BEAUFORT_ZOOM="-175 -100 65 80"


for prefix in  beaufort ; do
 case $prefix in
  'zapiola') ZOOM=$ZAPIOLA_ZOOM ;;
  'indian_ridge') ZOOM=$INDIAN_RIDGE_ZOOM ;;
  'romanche') ZOOM=$ROMANCHE_ZOOM ;;
  'amundsen') ZOOM=$AMUNDSEN_ZOOM ;;
  'beaufort') ZOOM=$BEAUFORT_ZOOM ;;
  * ) echo $prefix not pre-defined ... skip   ;
      prefix='none' ;;
 esac
 
 if [ $prefix != 'none' ] ; then
  for bathy in $GEBCO $ETOPO1 $ETOPO2 $SM_SW ; do
    batfil=${DATADIR}/${bathy}.nc
    zoom=$(../SRC/find_zoom $batfil $ZOOM)
    LON=$(echo $zoom | awk '{ print "-d lon,"$1","$2}')
    LAT=$(echo $zoom | awk '{ print "-d lat,"$3","$4}')
    ncks -F -O $LON $LAT $batfil ${prefix}_$bathy.nc
    ncview ${prefix}_$bathy.nc &
    # build pseudo coordinates
    ../SRC/mk_coor ${prefix}_$bathy.nc
    mv pseudo_coordinates.nc ${prefix}_${bathy}_pseudo_coordinates.nc
  done
 fi
done

