#!/bin/ksh
# This script is to convert GEODAS bathy track files (*.xyz) into
# ascii usable files *.iyxz for cdfweight2D and further cdfcoloc2D
## -------------------------------------------------------------------------
##  $Date: 2009-04-29 16:19:25 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 228 $
##  $Id: xyz2iyxz.ksh 228 2009-04-29 14:19:25Z forge $
## -------------------------------------------------------------------------

if [ $# != 1 ] ; then
  echo usage: xyz2iyxz.ksh file.xyz
  exit
fi

fxyz=$1
fiyxz=${fxyz%.xyz}.iyxz

cat $fxyz |  awk '{ print NR " " $2 " " $1 " " $3 }' > $fiyxz
