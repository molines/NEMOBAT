MODULE dommsk
   !!==============================================================================
   !!                       ***  MODULE dommsk   ***
   !! Ocean initialization : domain land/sea mask 
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_msk        : compute land/ocean mask
   !!   dom_msk_nsa    : update land/ocean mask when no-slip accurate
   !!                    option is used.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE obc_oce         ! ocean open boundary conditions
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp
   USE solisl          ! ???
   USE dynspg_fsc      !
   USE dynspg_fsc_atsk !

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dom_msk        ! routine called by inidom.F90

   !! * Module variables
   REAL(wp) ::   &
      shlat = 2.   ! type of lateral boundary condition on velocity (namelist namlbc)
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!---------------------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/dommsk.F90,v 1.5 2005/03/27 18:34:57 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!---------------------------------------------------------------------------------

CONTAINS
   
   SUBROUTINE dom_msk
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk  ***
      !!
      !! ** Purpose :   Compute land/ocean mask arrays at tracer points, hori-
      !!      zontal velocity points (u & v), vorticity points (f) and baro-
      !!      tropic stream function  points (b).
      !!        Set mbathy to the number of non-zero w-levels of a water column
      !!      (if island in the domain (lk_isl=T), this is done latter in
      !!      routine solver_init)
      !!
      !! ** Method  :   The ocean/land mask is computed from the basin bathy-
      !!      metry in level (mbathy) which is defined or read in dommba.
      !!      mbathy equals 0 over continental T-point, -n over the nth 
      !!      island T-point, and the number of ocean level over the ocean.
      !!
      !!      At a given position (ji,jj,jk) the ocean/land mask is given by:
      !!      t-point : 0. IF mbathy( ji ,jj) =< 0
      !!                1. IF mbathy( ji ,jj) >= jk
      !!      u-point : 0. IF mbathy( ji ,jj)  or mbathy(ji+1, jj ) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy(ji+1, jj ) >= jk.
      !!      v-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1) >= jk.
      !!      f-point : 0. IF mbathy( ji ,jj)  or mbathy( ji ,jj+1)
      !!                   or mbathy(ji+1,jj)  or mbathy(ji+1,jj+1) =< 0
      !!                1. IF mbathy( ji ,jj) and mbathy( ji ,jj+1)
      !!                and mbathy(ji+1,jj) and mbathy(ji+1,jj+1) >= jk.
      !!      b-point : the same definition as for f-point of the first ocean
      !!                level (surface level) but with 0 along coastlines.
      !!
      !!        The lateral friction is set through the value of fmask along
      !!      the coast and topography. This value is defined by shlat, a
      !!      namelist parameter:
      !!         shlat = 0, free slip  (no shear along the coast)
      !!         shlat = 2, no slip  (specified zero velocity at the coast)
      !!         0 < shlat < 2, partial slip   | non-linear velocity profile
      !!         2 < shlat, strong slip        | in the lateral boundary layer
      !!
      !!      N.B. If nperio not equal to 0, the land/ocean mask arrays
      !!      are defined with the proper value at lateral domain boundaries,
      !!      but bmask. indeed, bmask defined the domain over which the
      !!      barotropic stream function is computed. this domain cannot
      !!      contain identical columns because the matrix associated with
      !!      the barotropic stream function equation is then no more inverti-
      !!      ble. therefore bmask is set to 0 along lateral domain boundaries
      !!      even IF nperio is not zero.
      !!
      !!      In case of open boundaries (lk_obc=T):
      !!        - tmask is set to 1 on the points to be computed bay the open
      !!          boundaries routines.
      !!        - bmask is  set to 0 on the open boundaries.
      !!
      !!      Set mbathy to the number of non-zero w-levels of a water column
      !!                  mbathy = min( mbathy, 1 ) + 1
      !!      (note that the minimum value of mbathy is 2).
      !!
      !! ** Action :
      !!                     tmask    : land/ocean mask at t-point (=0. or 1.)
      !!                     umask    : land/ocean mask at u-point (=0. or 1.)
      !!                     vmask    : land/ocean mask at v-point (=0. or 1.)
      !!                     fmask    : land/ocean mask at f-point (=0. or 1.)
      !!                          =shlat along lateral boundaries
      !!                     bmask    : land/ocean mask at barotropic stream
      !!                          function point (=0. or 1.) and set to
      !!                          0 along lateral boundaries
      !!                   mbathy   : number of non-zero w-levels 
      !!
      !! History :
      !!        !  87-07  (G. Madec)  Original code
      !!        !  91-12  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
      !!        !  96-01  (G. Madec)  suppression of common work arrays
      !!        !  96-05  (G. Madec)  mask computed from tmask and sup-
      !!                 pression of the double computation of bmask
      !!        !  97-02  (G. Madec)  mesh information put in domhgr.F
      !!        !  97-07  (G. Madec)  modification of mbathy and fmask
      !!        !  98-05  (G. Roullet)  free surface
      !!        !  00-03  (G. Madec)  no slip accurate
      !!        !  01-09  (J.-M. Molines)  Open boundaries
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! *Local declarations
      INTEGER  ::   ji, jj, jk, ii     ! dummy loop indices
      INTEGER  ::   iif, iil, ijf, ijl
      INTEGER  ::   ii0, ii1, ij0, ij1
      INTEGER, DIMENSION(jpi,jpj) ::  imsk

      REAL(wp), DIMENSION(jpi,jpj) ::   zwf

      NAMELIST/namlbc/ shlat
      !!---------------------------------------------------------------------
      

      ! Namelist namlbc : lateral momentum boundary condition
      REWIND( numnam )
      READ  ( numnam, namlbc )
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dommsk : ocean mask '
         WRITE(numout,*) '~~~~~~'
         WRITE(numout,*) '         Namelist namlbc'
         WRITE(numout,*) '            lateral momentum boundary cond. shlat = ',shlat
      ENDIF

      IF ( shlat == 0. ) THEN
          IF(lwp) WRITE(numout,*) '         ocean lateral free-slip '
        ELSEIF ( shlat  ==  2. ) THEN
          IF(lwp) WRITE(numout,*) '         ocean lateral  no-slip '
        ELSEIF ( 0. < shlat .AND. shlat < 2. ) THEN
          IF(lwp) WRITE(numout,*) '         ocean lateral  partial-slip '
        ELSEIF ( 2. < shlat ) THEN
          IF(lwp) WRITE(numout,*) '         ocean lateral  strong-slip '
        ELSE
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) ' shlat is negative = ', shlat
          nstop = nstop + 1
      ENDIF

      ! 1. Ocean/land mask at t-point (computed from mbathy)
      ! -----------------------------
      ! Tmask has already the right boundary conditions since mbathy is ok

      tmask(:,:,:) = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( FLOAT( mbathy(ji,jj)-jk )+.1 >= 0.e0 ) tmask(ji,jj,jk) = 1.e0
            END DO  
         END DO  
      END DO  


      ! Interior domain mask (used for global sum)
      ! --------------------

      tmask_i(:,:) = tmask(:,:,1)
      iif = jpreci                         ! ???
      iil = nlci - jpreci + 1
      ijf = jprecj                         ! ???
      ijl = nlcj - jprecj + 1

      tmask_i( 1 :iif,   :   ) = 0.e0      ! first columns
      tmask_i(iil:jpi,   :   ) = 0.e0      ! last  columns (including mpp extra columns)
      tmask_i(   :   , 1 :ijf) = 0.e0      ! first rows
      tmask_i(   :   ,ijl:jpj) = 0.e0      ! last  rows (including mpp extra rows)

      ! north fold mask
      tpol(1:jpiglo) = 1.e0 
      fpol(1:jpiglo) = 1.e0
      IF( jperio == 3 .OR. jperio == 4 ) THEN      ! T-point pivot
         tpol(jpiglo/2+1:jpiglo) = 0.e0
         fpol(     1    :jpiglo) = 0.e0
      ENDIF
      IF( jperio == 5 .OR. jperio == 6 ) THEN      ! F-point pivot
         tpol(     1    :jpiglo) = 0.e0
         fpol(jpiglo/2+1:jpiglo) = 0.e0
      ENDIF

      IF( nperio == 3 .OR. nperio == 4 ) THEN      ! T-point pivot: only half of the nlcj-1 row
         DO ji = iif+1, iil-1
            tmask_i(ji,ijl-1) = tmask_i(ji,ijl-1) * tpol(mig(ji))
         END DO
      ENDIF 


      ! 2. Ocean/land mask at u-,  v-, and z-points (computed from tmask)
      ! -------------------------------------------
      
      ! Computation
      DO jk = 1, jpk
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector loop
               umask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)
               vmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji  ,jj+1,jk)
               fmask(ji,jj,jk) = tmask(ji,jj  ,jk) * tmask(ji+1,jj  ,jk)   &
                  &            * tmask(ji,jj+1,jk) * tmask(ji+1,jj+1,jk)
            END DO
         END DO
      END DO

      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
         !                                           ! =======================
         ! modified vmask value in                   !  ORCA_R2 configuration
         ! the vicinity of some straits              ! =======================

         IF( n_cla == 1 ) THEN 
            !                                ! vmask = 0. on Gibraltar zonal section
            ij0 = 101   ;   ij1 = 101
            ii0 = 138   ;   ii1 = 139   ;   vmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 19:jpk ) = 0.e0
            !                                ! vmask = 0. on Bab el Mandeb zonal section
            ij0 =  87   ;   ij1 =  87
            ii0 = 161   ;   ii1 = 163   ;   vmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 18:jpk ) = 0.e0
         ENDIF

      ENDIF
      ! akl
         IF( cp_cfg == "itf" .AND. jp_cfg == 25 ) THEN  ! ORCA R025 configuration
            !                                             ! =======================
            ii0 = 125    ;   ii1 = 125        ! East of Ombai strait
            !! ij0 = 71   P18 -  ;   ij1 = 72   ;   fmask( ii0:ii1 , ij0:ij1, 1:jpk ) =  2.0
            ij0 = 72     ;   ij1 = 72   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the East Ombai Strait'
            !! ii0 = 123  P18 -   ;   ii1 = 124        ! West of Ombai strait
            ii0 = 124    ;   ii1 = 124        ! West of Ombai strait
            ij0 = 73     ;   ij1 = 73   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the West Ombai Strait '
            ! ii0 = 123  P18 -   ;   ii1 = 123        ! exit of Ombai strait
            ! ij0 = 71     ;   ij1 = 72   ;   fmask( ii0:ii1 , ij0:ij1, 1:jpk ) =  2.0
            ! IF(lwp) WRITE(numout,*)
            ! IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the exit of Ombai Strait '
            ii0 = 87   ;   ii1 = 87        ! Lombok strait
            ij0 = 70   ;   ij1 = 70   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the Lombok Strait'
            ii0 = 100   ;   ii1 = 100        ! Sape strait
            ij0 = 71   ;   ij1 = 72   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the Sape Strait'
            ii0 = 118   ;   ii1 = 118        ! Alor strait
            ij0 = 72   ;   ij1 = 72   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r025: fmask = 2 at the Alor Strait'
            !
         ENDIF

         IF( cp_cfg == "orca" .AND. jp_cfg == 05 ) THEN   ! ORCA R05 configuration
            !
            ii0 =  93   ;   ii1 =  94        ! Sumba Strait (e2u = 40 km)
            ij0 = 232   ;   ij1 = 232   ;    fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r05: Reduced e2u at the Sumba Strait'
            !
            ii0 = 103   ;   ii1 = 103        ! Ombai Strait (e2u = 15 km)
            ij0 = 232   ;   ij1 = 232   ;    fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r05: Reduced e2u at the Ombai Strait'
            !
            ii0 =  87   ;   ii1 =  87        ! Lombok Strait (e1v = 10 km)
            ij0 = 232   ;   ij1 = 233   ;    fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1), 1:jpk ) =  2.0
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '             orca_r05: Reduced e1v at the Lombok Strait'
            !
            !
         ENDIF

      ! Lateral boundary conditions
      CALL lbc_lnk( umask, 'U', 1. )
      CALL lbc_lnk( vmask, 'V', 1. )
      CALL lbc_lnk( fmask, 'F', 1. )


      ! 4. ocean/land mask for the elliptic equation
      ! --------------------------------------------
      
      ! Computation
      IF( lk_dynspg_fsc .OR. lk_dynspg_fsc_tsk ) THEN
         bmask(:,:) = tmask(:,:,1)       ! elliptic equation is written at t-point
      ELSE
         bmask(:,:) = fmask(:,:,1)       ! elliptic equation is written at f-point
      ENDIF
      
      ! Boundary conditions
      !   cyclic east-west : bmask must be set to 0. on rows 1 and jpi
      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) THEN
         bmask( 1 ,:) = 0.e0
         bmask(jpi,:) = 0.e0
      ENDIF
      
      !   south symmetric :  bmask must be set to 0. on row 1
      IF( nperio == 2 ) THEN
         bmask(:, 1 ) = 0.e0
      ENDIF
      
      !   north fold : 
      IF( nperio == 3 .OR. nperio == 4 ) THEN
         IF( lk_dynspg_fsc .OR. lk_dynspg_fsc_tsk ) THEN
            ! T-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
            DO ji = 1, jpi
               ii = ji + nimpp - 1
               bmask(ji,jpj-1) = bmask(ji,jpj-1) * tpol(ii)
               bmask(ji,jpj  ) = 0.e0
            END DO
         ELSE
            ! T-pt pivot and F-pt elliptic eq. : bmask set to 0. on rows jpj-1 and jpj
            bmask(:,jpj-1) = 0.e0
            bmask(:,jpj  ) = 0.e0
         ENDIF
      ENDIF
      IF( nperio == 5 .OR. nperio == 6 ) THEN
         IF( lk_dynspg_fsc .OR. lk_dynspg_fsc_tsk ) THEN
            ! F-pt pivot and T-pt elliptic eq. : bmask set to 0. on row jpj
            bmask(:,jpj) = 0.e0
         ELSE
            ! F-pt pivot and F-pt elliptic eq. : bmask set to 0. on row jpj and on half jpjglo-1 row
            DO ji = 1, jpi
               ii = ji + nimpp - 1
               bmask(ji,jpj-1) = bmask(ji,jpj-1) * fpol(ii)
               bmask(ji,jpj  ) = 0.e0
            END DO
         ENDIF
      ENDIF

      ! Mpp boundary conditions: bmask is set to zero on the overlap
      ! region for all elliptic solvers

      IF( lk_mpp ) THEN
         IF( nbondi /= -1 .AND. nbondi /= 2 )   bmask(  1 :jpreci,:) = 0.e0
         IF( nbondi /=  1 .AND. nbondi /= 2 )   bmask(nlci:jpi   ,:) = 0.e0
         IF( nbondj /= -1 .AND. nbondj /= 2 )   bmask(:,  1 :jprecj) = 0.e0
         IF( nbondj /=  1 .AND. nbondj /= 2 )   bmask(:,nlcj:jpj   ) = 0.e0
      
         ! north fold : bmask must be set to 0. on rows jpj-1 and jpj 
         IF( npolj == 3 .OR. npolj == 4 ) THEN
            IF( lk_dynspg_fsc .OR. lk_dynspg_fsc_tsk ) THEN
               DO ji = 1, nlci
                  ii = ji + nimpp - 1
                  bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * tpol(ii)
                  bmask(ji,nlcj  ) = 0.e0
               END DO
            ELSE
               DO ji = 1, nlci
                  bmask(ji,nlcj-1) = 0.e0
                  bmask(ji,nlcj  ) = 0.e0
               END DO
            ENDIF
         ENDIF
         IF( npolj == 5 .OR. npolj == 6 ) THEN
            IF( lk_dynspg_fsc .OR. lk_dynspg_fsc_tsk ) THEN
               DO ji = 1, nlci
                  bmask(ji,nlcj  ) = 0.e0
               END DO
            ELSE
               DO ji = 1, nlci
                  ii = ji + nimpp - 1
                  bmask(ji,nlcj-1) = bmask(ji,nlcj-1) * fpol(ii)
                  bmask(ji,nlcj  ) = 0.e0
               END DO
            ENDIF
         ENDIF
      ENDIF


      ! mask for second order calculation of vorticity
      ! ----------------------------------------------
      
      CALL dom_msk_nsa

      
      ! Lateral boundary conditions on velocity (modify fmask)
      ! ---------------------------------------
      
      DO jk = 1, jpk

         zwf(:,:) = fmask(:,:,jk)
         
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
               IF( fmask(ji,jj,jk) == 0. ) THEN
                  fmask(ji,jj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jj), zwf(ji,jj+1),   &
                     &                                    zwf(ji-1,jj), zwf(ji,jj-1)  )  )
               ENDIF
            END DO
         END DO
         
         DO jj = 2, jpjm1
            IF( fmask(1,jj,jk) == 0. ) THEN
               fmask(1  ,jj,jk) = shlat * MIN( 1., MAX( zwf(2,jj), zwf(1,jj+1), zwf(1,jj-1) ) )
            ENDIF
            IF( fmask(jpi,jj,jk) == 0. ) THEN
               fmask(jpi,jj,jk) = shlat * MIN( 1., MAX( zwf(jpi,jj+1), zwf(jpim1,jj), zwf(jpi,jj-1) ) )
            ENDIF
         END DO
         
         DO ji = 2, jpim1
            IF( fmask(ji,1,jk) == 0. ) THEN
               fmask(ji, 1 ,jk) = shlat * MIN( 1., MAX( zwf(ji+1,1), zwf(ji,2), zwf(ji-1,1) ) )
            ENDIF
            IF( fmask(ji,jpj,jk) == 0. ) THEN
               fmask(ji,jpj,jk) = shlat * MIN( 1., MAX( zwf(ji+1,jpj), zwf(ji-1,jpj), zwf(ji,jpjm1) ) )
            ENDIF
         END DO
      END DO
      

      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN
         !                                           ! =======================
         ! Increased lateral friction in             !  ORCA_R2 configuration
         ! the vicinity of some straits              ! =======================
         !
         IF( n_cla == 0 ) THEN
            !                                ! Sound  strait
            ij0 = 116   ;   ij1 = 117
            ii0 = 147   ;   ii1 = 148   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 10.0e0
         ELSE
            !                                ! Gibraltar strait and Gulf of Cadiz
            ij0 = 102   ;   ij1 = 102
            ii0 = 137   ;   ii1 = 139   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.0e0
            ij0 = 101   ;   ij1 = 101
            ii0 = 139   ;   ii1 = 139   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.0e0
            ij0 = 100   ;   ij1 = 100
            ii0 = 137   ;   ii1 = 139   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) =  0.0e0
            !                                ! Sound  strait
            ij0 = 116   ;   ij1 = 117
            ii0 = 147   ;   ii1 = 148   ;   fmask( mi0(ii0):mi1(ii1) , mj0(ij0):mj1(ij1) , 1:jpk ) = 10.0e0
         ENDIF
         !
      ENDIF
      
      ! Lateral boundary conditions on fmask
      CALL lbc_lnk( fmask, 'F', 1. )
      
      ! Mbathy set to the number of w-level (minimum value 2)
      ! -----------------------------------
      IF( lk_isl ) THEN
         ! this is done at the end of solver_init routine
      ELSE
         DO jj = 1, jpj
            DO ji = 1, jpi
               mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
            END DO
         END DO
      ENDIF
      
      ! Control print
      ! -------------
      IF( nprint == 1 .AND. lwp ) THEN
         imsk(:,:) = INT( tmask_i(:,:) )
         WRITE(numout,*) ' tmask_i : '
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                           1, jpj, 1, 1, numout)
         WRITE (numout,*)
         WRITE (numout,*) ' dommsk: tmask for each level'
         WRITE (numout,*) ' ----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( tmask(:,:,jk) )

            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: vmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( vmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: fmask for each level'
         WRITE(numout,*) ' -----------------------------'
         DO jk = 1, jpk
            imsk(:,:) = INT( fmask(:,:,jk) )
            WRITE(numout,*)
            WRITE(numout,*) ' level = ',jk
            CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                              1, jpj, 1, 1, numout )
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' dom_msk: bmask '
         WRITE(numout,*) ' ---------------'
         WRITE(numout,*)
         imsk(:,:) = INT( bmask(:,:) )
         CALL prihin( imsk(:,:), jpi, jpj, 1, jpi, 1,   &
               &                           1, jpj, 1, 1, numout )
      ENDIF

   END SUBROUTINE dom_msk

#if defined key_noslip_accurate
   !!----------------------------------------------------------------------
   !!   'key_noslip_accurate' :         accurate no-slip boundary condition
   !!----------------------------------------------------------------------
   
   SUBROUTINE dom_msk_nsa
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dom_msk_nsa  ***
      !! 
      !! ** Purpose :
      !!
      !! ** Method  :
      !!
      !! ** Action :
      !!
      !! History :
      !!        !  00-03  (G. Madec)  no slip accurate
      !!----------------------------------------------------------------------
      !! *Local declarations
      INTEGER  :: ji, jj, jk, ii     ! dummy loop indices
      INTEGER ::   ine, inw, ins, inn, itest, ierror, iind, ijnd
      INTEGER, DIMENSION(jpi*jpj*jpk,3) ::  icoord
      REAL(wp) ::   zaa
      !!---------------------------------------------------------------------
      !!  OPA 9.0 , LOCEAN-IPSL (2005)
      !!---------------------------------------------------------------------
      

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) 'dom_msk_nsa : noslip accurate boundary condition'
      IF(lwp)WRITE(numout,*) '~~~~~~~~~~~   using Schchepetkin and O Brian scheme'
      IF( lk_mpp ) THEN
         IF(lwp)WRITE(numout,cform_err)
         IF(lwp)WRITE(numout,*) ' mpp version is not yet implemented'
         nstop = nstop + 1
      ENDIF

      ! mask for second order calculation of vorticity
      ! ----------------------------------------------
      ! noslip boundary condition: fmask=1  at convex corner, store
      ! index of straight coast meshes ( 'west', refering to a coast,
      ! means west of the ocean, aso)
      
      DO jk = 1, jpk
         DO jl = 1, 4
            npcoa(jl,jk) = 0
            DO ji = 1, 2*(jpi+jpj)
               nicoa(ji,jl,jk) = 0
               njcoa(ji,jl,jk) = 0
            END DO
         END DO
      END DO
      
      IF( jperio == 2 ) THEN
         WRITE(numout,*) ' '
         WRITE(numout,*) ' symetric boundary conditions need special'
         WRITE(numout,*) ' treatment not implemented. we stop.'
         STOP
      ENDIF
      
      ! convex corners
      
      DO jk = 1, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               zaa = tmask(ji  ,jj,jk) + tmask(ji  ,jj+1,jk)   &
                  &+ tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk)
               IF( ABS(zaa-3.) <= 0.1 )   fmask(ji,jj,jk) = 1.
            END DO
         END DO
      END DO

      ! north-south straight coast

      DO jk = 1, jpkm1
         inw = 0
         ine = 0
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zaa = tmask(ji+1,jj,jk) + tmask(ji+1,jj+1,jk)
               IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
                  inw = inw + 1
                  nicoa(inw,1,jk) = ji
                  njcoa(inw,1,jk) = jj
                  IF( nprint == 1 ) WRITE(numout,*) ' west  : ', jk, inw, ji, jj
               ENDIF
               zaa = tmask(ji,jj,jk) + tmask(ji,jj+1,jk)
               IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
                  ine = ine + 1
                  nicoa(ine,2,jk) = ji
                  njcoa(ine,2,jk) = jj
                  IF( nprint == 1 ) WRITE(numout,*) ' east  : ', jk, ine, ji, jj
               ENDIF
            END DO
         END DO
         npcoa(1,jk) = inw
         npcoa(2,jk) = ine
      END DO

      ! west-east straight coast

      DO jk = 1, jpkm1
         ins = 0
         inn = 0
         DO jj = 2, jpjm1
            DO ji =2, jpim1
               zaa = tmask(ji,jj+1,jk) + tmask(ji+1,jj+1,jk)
               IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
                  ins = ins + 1
                  nicoa(ins,3,jk) = ji
                  njcoa(ins,3,jk) = jj
                  IF( nprint == 1 ) WRITE(numout,*) ' south : ', jk, ins, ji, jj
               ENDIF
               zaa = tmask(ji+1,jj,jk) + tmask(ji,jj,jk)
               IF( ABS(zaa-2.) <= 0.1 .AND. fmask(ji,jj,jk) == 0 ) THEN
                  inn = inn + 1
                  nicoa(inn,4,jk) = ji
                  njcoa(inn,4,jk) = jj
                  IF( nprint == 1 ) WRITE(numout,*) ' north : ', jk, inn, ji, jj
               ENDIF
            END DO
         END DO
         npcoa(3,jk) = ins
         npcoa(4,jk) = inn
      END DO

      itest = 2 * ( jpi + jpj )
      DO jk = 1, jpk
         IF( npcoa(1,jk) > itest .OR. npcoa(2,jk) > itest .OR.   &
             npcoa(3,jk) > itest .OR. npcoa(4,jk) > itest ) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' level jk = ',jk
            WRITE(numout,*) ' straight coast index arraies are too small.:'
            WRITE(numout,*) ' npe, npw, nps, npn = ', npcoa(1,jk), npcoa(2,jk),   &
                &                                     npcoa(3,jk), npcoa(4,jk)
            WRITE(numout,*) ' 2*(jpi+jpj) = ',itest,'. we stop.'
            STOP   !!bug nstop to be used
        ENDIF
      END DO

      ierror = 0
      iind = 0
      ijnd = 0
      IF( nperio == 1 .OR. nperio == 4 .OR. nperio == 6 ) iind = 2
      IF( nperio == 3 .OR. nperio == 4 .OR. nperio == 5 .OR. nperio == 6 ) ijnd = 2
      DO jk = 1, jpk
         DO jl = 1, npcoa(1,jk)
            IF( nicoa(jl,1,jk)+3 > jpi+iind ) THEN
               ierror = ierror+1
               icoord(ierror,1) = nicoa(jl,1,jk)
               icoord(ierror,2) = njcoa(jl,1,jk)
               icoord(ierror,3) = jk
            ENDIF
         END DO
         DO jl = 1, npcoa(2,jk)
            IF(nicoa(jl,2,jk)-2 < 1-iind ) THEN
               ierror = ierror + 1
               icoord(ierror,1) = nicoa(jl,2,jk)
               icoord(ierror,2) = njcoa(jl,2,jk)
               icoord(ierror,3) = jk
            ENDIF
         END DO
         DO jl = 1, npcoa(3,jk)
            IF( njcoa(jl,3,jk)+3 > jpj+ijnd ) THEN
               ierror = ierror + 1
               icoord(ierror,1) = nicoa(jl,3,jk)
               icoord(ierror,2) = njcoa(jl,3,jk)
               icoord(ierror,3) = jk
            ENDIF
         END DO
         DO jl=1,npcoa(4,jk)
            IF( njcoa(jl,4,jk)-2 < 1) THEN
               ierror=ierror+1
               icoord(ierror,1)=nicoa(jl,4,jk)
               icoord(ierror,2)=njcoa(jl,4,jk)
               icoord(ierror,3)=jk
            ENDIF
         END DO
      END DO
      
      IF( ierror > 0 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '              Problem on lateral conditions'
         IF(lwp) WRITE(numout,*) '                 Bad marking off at points:'
         DO jl = 1, ierror
            IF(lwp) WRITE(numout,*) 'Level:',icoord(jl,3),   &
               &                  '  Point(',icoord(jl,1),',',icoord(jl,2),')'
         END DO
         IF(lwp) WRITE(numout,*) 'We stop...'   !!cr print format to be used
         nstop = nstop + 1
      ENDIF

   END SUBROUTINE dom_msk_nsa

#else
   !!----------------------------------------------------------------------
   !!   Default option :                                      Empty routine
   !!----------------------------------------------------------------------
   SUBROUTINE dom_msk_nsa       
   END SUBROUTINE dom_msk_nsa
#endif
   
   !!======================================================================
END MODULE dommsk
