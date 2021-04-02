#!/bin/ksh
## -------------------------------------------------------------------------
##  $Date: 2009-04-29 22:55:31 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 234 $
##  $Id: look_detail.ksh 234 2009-04-29 20:55:31Z forge $
## -------------------------------------------------------------------------
#set -x
. ../bin/chk_zoom.dat
set +x
bathy=ORCA025_bathy_etopo1_gebco1_smoothed_coast_corrected_apr09.nc

for zone in $zoom ; do

  name=$( echo $zone | tr -d '$' | sed -e 's/_ZOOM//' )
  window=$(eval "echo $zone")
  echo $name
  ijzoom=$(cdffindij $window  coordinates.nc T | tail -2 | head -1)
  # cdffindij had some problem in the arctic ... hard coded value (ORCA025)
  if [ $name = 'ARCTIC_SEA' ] ; then
   ijzoom='150 698 704 1020'
  fi
  imin=$(echo $ijzoom | awk '{print $1}')
  imax=$(echo $ijzoom | awk '{print $2}')
  jmin=$(echo $ijzoom | awk '{print $3}')
  jmax=$(echo $ijzoom | awk '{print $4}')

  lonlim='-d x,'$imin','$imax
  latlim='-d y,'$jmin','$jmax
  echo $ijzoom
  ncks -F -O $lonlim $latlim $bathy ${name}_bathy.nc
  cdfbathy -f $bathy -z $ijzoom -d $name.txt
  echo
done
