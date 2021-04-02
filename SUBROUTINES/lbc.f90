   SUBROUTINE lbc ( ptab,jpi,jpj, jpk, ktype,ksgn,nperio)
      !!---------------------------------------------------------------------
      !!                       ROUTINE lbc
      !!                     ***************
      !! ** Purpose :
      !!      Insure the lateral boundary conditions for a bathymetry file
      !!
      !! ** Method :
      !!
      !! History :
      !!        !  97-06  (G. Madec)  Original code OPA
      !!   8.5  !  02-09  (G. Madec)  F90: Free form and module
      !!   adapted to take it out of OPA code: A.M. Treguier
      !!----------------------------------------------------------------------
   !! -------------------------------------------------------------------------
   !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
   !!  $Rev: 225 $
   !!  $Id: lbc.f90 225 2009-04-29 13:52:46Z molines $
   !! -------------------------------------------------------------------------
      IMPLICIT NONE

!!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT( in )                         ::          jpi,jpj,jpk
         ! define dimension of input array
      REAL, dimension(jpi,jpj,jpk)                  ::          ptab
      INTEGER, INTENT( in )                         ::          ktype,ksgn,nperio 
         ! ktype     define the nature of the grid-pointat which ptab is defined
         !                 = 1 ,  T- and W-points
         !                 = 2 ,  U-point
         !                 = 3 ,  V-point
         !                 = 4 ,  F-point  
         ! ksgn,      control of the sign change
         !                 = 0 , the sign is modified following the type of b.c. used
         !                 = 1 , the sign of the field is unchanged at the boundaries
         ! nperio,    control of  boundary conditions
         !     jperio= 0, closed'
         !     jperio= 1, cyclic east-west'
         !     jperio= 2, equatorial symmetric'
         !     jperio= 3, north fold with T-point pivot'
         !     jperio= 4, same as above, cyclic east-west',   
         !     jperio= 5, north fold with F-point pivot'
         !     jperio= 6, same as above,cyclic east-west',   
 

      !! * Local declarations
      INTEGER ji, jj,jk,jpim1,jpjm1
      INTEGER ijt, iju
      REAL       ::   zsgn
      !!----------------------------------------------------------------------
      !!  OPA 8.5, LODYC-IPSL (2002)
      !!----------------------------------------------------------------------
      
      jpim1=jpi-1
      jpjm1= jpj-1

      ! 0. Sign setting
      ! ---------------
      
      IF( ksgn == 0 ) THEN
         zsgn = -1.
      ELSE
         zsgn =  1.
      ENDIF


         DO jk = 1,jpk


            ! 1. East-West boundary conditions
            ! --------------------------------

            IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
               ! cyclic
               ptab( 1 ,:,jk) = ptab(jpim1,:,jk)
               ptab(jpi,:,jk) = ptab(  2  ,:,jk)
            ELSE
               ! closed
               IF( ktype == 4 ) THEN
                  ptab(jpi,:,jk) = 0.e0
               ELSE
                  ptab( 1 ,:,jk) = 0.e0
                  ptab(jpi,:,jk) = 0.e0
               ENDIF
            ENDIF

            ! 2. North-South boundary conditions
            ! ----------------------------------

            IF( nperio == 2 ) THEN
               ! south symmetric
               IF( ktype == 1 .OR. ktype == 2 ) THEN
                  ptab(:, 1 ,jk) = ptab(:,3,jk)
                  ptab(:,jpj,jk) = 0.e0
               ELSEIF( ktype == 3 .OR. ktype == 4 ) THEN
                  ptab(:, 1 ,jk) = zsgn * ptab(:,2,jk)
                  ptab(:,jpj,jk) = 0.e0
               ENDIF
            ELSEIF( nperio == 3 .OR. nperio == 4 ) THEN
               ! north fold ( Grid defined with a T-point pivot) ORCA 2 degre
               ptab( 1 ,jpj,jk) = 0.e0
               ptab(jpi,jpj,jk) = 0.e0
               IF( ktype == 1 ) THEN
                  DO ji = 2, jpi
                     ijt = jpi-ji+2
                     ptab(ji, 1 ,jk) = 0.e0
                     ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-2,jk)
                  END DO
                  DO ji = jpi/2+1, jpi
                     ijt = jpi-ji+2
                     ptab(ji,jpjm1,jk) = zsgn * ptab(ijt,jpjm1,jk)
                  END DO
               ELSEIF( ktype == 2 ) THEN
                  DO ji = 1, jpi-1
                     iju = jpi-ji+1
                     ptab(ji, 1 ,jk) = 0.e0
                     ptab(ji,jpj,jk) = zsgn * ptab(iju,jpj-2,jk)
                  END DO
                  DO ji = jpi/2, jpi-1
                     iju = jpi-ji+1
                     ptab(ji,jpjm1,jk) = zsgn * ptab(iju,jpjm1,jk)
                  END DO
               ELSEIF( ktype == 3 ) THEN
                  DO ji = 2, jpi
                     ijt = jpi-ji+2
                     ptab(ji,  1  ,jk) = 0.e0
                     ptab(ji,jpj-1,jk) = zsgn * ptab(ijt,jpj-2,jk)
                     ptab(ji,jpj  ,jk) = zsgn * ptab(ijt,jpj-3,jk)
                  END DO
               ELSEIF( ktype == 4 ) THEN
                  DO ji = 1, jpi-1
                     iju = jpi-ji+1
                     ptab(ji,jpj-1,jk) = ptab(iju,jpj-2,jk)
                     ptab(ji,jpj  ,jk) = ptab(iju,jpj-3,jk)
                  END DO
               ENDIF
            ELSEIF( nperio == 5 .OR. nperio == 6 ) THEN
               ! north fold ( Grid defined with a F-point pivot) ORCA 0.5 degre
               ptab( 1 ,jpj,jk) = 0.e0
               ptab(jpi,jpj,jk) = 0.e0
               IF( ktype == 1 ) THEN
                  DO ji = 1, jpi
                     ijt = jpi-ji+1
                     ptab(ji, 1 ,jk) = 0.e0
                     ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-1,jk)
                  END DO
               ELSEIF( ktype == 2 ) THEN
                  DO ji = 1, jpi-1
                     iju = jpi-ji
                     ptab(ji, 1 ,jk) = 0.e0
                     ptab(ji,jpj,jk) = zsgn * ptab(iju,jpj-1,jk)
                  END DO
               ELSEIF( ktype == 3 ) THEN
                  DO ji = 1, jpi
                     ijt = jpi-ji+1
                     ptab(ji, 1 ,jk) = 0.e0
                     ptab(ji,jpj,jk) = zsgn * ptab(ijt,jpj-2,jk)
                  END DO
                  DO ji = jpi/2+1, jpi
                     ijt = jpi-ji+1
                     ptab(ji,jpjm1,jk) = zsgn * ptab(ijt,jpjm1,jk)
                  END DO
               ELSEIF( ktype == 4 ) THEN
                  DO ji = 1, jpi-1
                     iju = jpi-ji
                     ptab(ji,jpj  ,jk) = ptab(iju,jpj-2,jk)
                  END DO
                  DO ji = jpi/2+1, jpi-1
                     iju = jpi-ji
                     ptab(ji,jpjm1,jk) = ptab(iju,jpjm1,jk)
                  END DO
               ENDIF
            ELSE
               ! closed
               IF( ktype == 4 ) THEN
                  ptab(:,jpj,jk) = 0.e0
               ELSE
                  ptab(:, 1 ,jk) = 0.e0
                  ptab(:,jpj,jk) = 0.e0
               ENDIF
            ENDIF
 
            !     End of slab
            !     ===========

         END DO

   END SUBROUTINE lbc
