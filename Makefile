# Makefile to build the whole set of executable used in NEMOBAT
## -------------------------------------------------------------------------
##  $Date: 2009-04-29 18:40:59 +0200 (Wed, 29 Apr 2009) $
##  $Rev: 233 $
##  $Id: Makefile 233 2009-04-29 16:40:59Z forge $
## -------------------------------------------------------------------------
DIRS = APPLYCOAST  COMBINE  INPUT_UTILITIES  SMOOTHING   CORRECT  INTERP0 


all:
	for dir in $(DIRS) ; do \
           cd $$dir ; echo doing make in $$dir ; make ; \
           cd ../ ; \
        done

clean:
	\rm -f *.o *~ svn-commit.tmp*
