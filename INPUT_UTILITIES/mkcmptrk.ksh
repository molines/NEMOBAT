#!/bin/ksh
## -------------------------------------------------------------------------
##  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 228 $
##  $Id: mkcmptrk.ksh 228 2009-04-29 14:19:25Z forge $
## -------------------------------------------------------------------------
#set -x 
TRKDIR=/home/users/molines/GEODAS/trk55200_data/trk55200

#for trk in georgia elt08 hydr03mv rc1214 ; do
for trkf in $TRKDIR/*.xyz  ; do
  trk1=$(basename $trkf) ; trk=${trk1%.xyz}
  echo WORKING for $trk
  # convert file xyz to iyxz 
  ln -s $TRKDIR/$trk.xyz .
  ./xyz2iyxz.ksh $trk.xyz
  \rm $trk.xyz
  #  for each bathy 
   \rm -f tmp
   cat $trk.iyxz | awk '{printf"% 9.4f % 9.4f % 7.1f\n",$3,$2,$4 }' > tmp
   header="longitude latitude hydro"

  for coord in  *pseudo_coordinates.nc ; do
    bathy=${coord%_pseudo_coordinates.nc}
    # compute weight
    cdfweight2D $trk.iyxz $coord T
    mv weight_T.bin ${trk}_${bathy}_weight_T.bin
    # colocolize
    cdfcoloc2D ${trk}_${bathy}_weight $bathy.nc
    mv izb.txt  ${trk}_${bathy}_izb.txt
    name=$(echo $bathy | awk -F_ '{ print $2 }')
    if [ $name = 'topo' ] ; then name='Smith_Sandwell_v11.1' ; fi
    header="$header $name"
    cat ${trk}_${bathy}_izb.txt | awk '{ print $NF}' > tmp2 # take only last column
    paste tmp tmp2 > ztmp
    mv -f ztmp tmp
  done
    echo  ' '$header > ${trk}_izb.csv
    cat tmp >> ${trk}_izb.csv
    \rm *.bin
    \rm *.iyxz
    \rm -f tmp ztmp tmp2 *.txt
done
