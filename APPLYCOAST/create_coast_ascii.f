      program create_coast_ascii
c
c    reads an ascii (old-fashioned)  bathymetry opa file 
c    creates an ascii file of zeroes and ones (mask)
c  
c
c   

c  !! -------------------------------------------------------------------------
c  !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
c  !!  $Rev: 225 $
c  !!  $Id: create_coast_ascii.f 225 2009-04-29 13:52:46Z molines $
c  !! -------------------------------------------------------------------------
      implicit none
      
      integer jpi,jpj    
C        ORCA05
      PARAMETER(JPI = 722,JPJ = 511 )
      integer mbathy(jpi,jpj)
      character filnam*80,filout*80
      character clexp*15
      
      integer ifreq,il1,il2
      integer ji,jj,ij,jn,ii
      integer iim,ijm,numin,numout
C
C
c        
      numin=10
      numout = 20
c
c  1.  OPEN FILES
c
      filnam= '/home/planier/ORCA_R05/BATHY_LEVEL/bathymetry_ORCA_R05.ascii'
      filout = 'coast_original_ORCA_R05.a'
      OPEN(UNIT=numin,FILE=filnam,STATUS='OLD',
     $     FORM='FORMATTED')
      open(unit=numout,file=filout)
      print *, ' output file will be ',filout
c
c  2. Reads data in the files
c
c
c           read bathymetry
c      
      READ(numin,9101) clexp,iim,ijm
      if (iim.ne.jpi.or.ijm.ne.jpj) then
        print *, ' modify jpi and jpj in coasline.f'
        print *, ' bathy read: jpi=',iim,' jpj=',ijm
        print *, ' differs from jpi=',jpi,' jpj=',jpj
        stop
      endif
      READ(numin,'(/)')
          ifreq=40
          il1=1
          DO jn=1,jpi/ifreq+1
            READ(numin,'(/)')
            il2=min0(jpi,il1+ifreq-1)
                READ(numin,9201) (ii,ji=il1,il2,5)
            READ(numin,'(/)')
            DO jj=jpj,1,-1
              READ(numin,9202) ij,(mbathy(ji,jj),ji=il1,il2)
            END DO
            il1=il1+ifreq
          END DO
C
 9101     FORMAT(1x,a15,2i8)
 9201     FORMAT(3x,13(i3,12x))
 9202     FORMAT(i3,41i3)
        print *, ' read bathymetry OK'
c
c  3.  transforms the bathymetry into a mask if it is not.
c       values of zero (land), lt 0 (islands) or 1 (ocean) only.
c
      do jj=1,jpj
        do ji=1,jpi
          mbathy(ji,jj)=min(mbathy(ji,jj),1)
        enddo
      enddo
c
C  now write coast mask
C
        do jj=jpj,1,-1
           write(numout,900) (mbathy(ji,jj),ji=1,jpi)
        enddo
 900    format(2000i1)

      end
