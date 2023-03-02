
c  !! -------------------------------------------------------------------------
c  !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
c  !!  $Rev: 225 $
c  !!  $Id: scilib.F 225 2009-04-29 13:52:46Z molines $
c  !! -------------------------------------------------------------------------

c************************************************************************
      real function sdot(n,sx,incx,sy,incy)
c************************************************************************
      real sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     +   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end

c************************************************************************
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c************************************************************************
      real sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end

c************************************************************************
      function ssum(n,fsum,nstep)
c************************************************************************
      dimension fsum(n)

      tmp=0.
      do 100 ij=1,n,nstep
        tmp=tmp+fsum(ij)
 100  continue
      ssum=tmp
      return
      end

c************************************************************************
      subroutine scopy(n,sx,incx,sy,incy)                                       
c************************************************************************
clll. optimize                                                                  
c                                                                               
c     copies a vector, x, to a vector, y.                                       
c     uses unrolled loops for increments equal to 1.                            
c     jack dongarra, linpack, 3/11/78.                                          
c                                                                               
      real sx(1),sy(1)                                                          
      integer i,incx,incy,ix,iy,m,mp1,n                                         
c                                                                               
      if(n.le.0)return                                                          
      if(incx.eq.1.and.incy.eq.1)go to 20                                       
c                                                                               
c        code for unequal increments or equal increments                        
c          not equal to 1                                                       
c                                                                               
      ix = 1                                                                    
      iy = 1                                                                    
      if(incx.lt.0)ix = (-n+1)*incx + 1                                         
      if(incy.lt.0)iy = (-n+1)*incy + 1                                         
      do 10 i = 1,n                                                             
        sy(iy) = sx(ix)                                                         
        ix = ix + incx                                                          
        iy = iy + incy                                                          
   10 continue                                                                  
      return                                                                    
c                                                                               
c        code for both increments equal to 1                                    
c                                                                               
c                                                                               
c        clean-up loop                                                          
c                                                                               
   20 m = mod(n,7)                                                              
      if( m .eq. 0 ) go to 40                                                   
      do 30 i = 1,m                                                             
        sy(i) = sx(i)                                                           
   30 continue                                                                  
      if( n .lt. 7 ) return                                                     
   40 mp1 = m + 1                                                               
      do 50 i = mp1,n,7                                                         
        sy(i) = sx(i)                                                           
        sy(i + 1) = sx(i + 1)                                                   
        sy(i + 2) = sx(i + 2)                                                   
        sy(i + 3) = sx(i + 3)                                                   
        sy(i + 4) = sx(i + 4)                                                   
        sy(i + 5) = sx(i + 5)                                                   
        sy(i + 6) = sx(i + 6)                                                   
   50 continue                                                                  
      return                                                                    
      end 
                                                                      
c************************************************************************
      subroutine sscal(n,sa,sx,incx)                                            
c************************************************************************
clll. optimize                                                                  
c                                                                               
c     scales a vector by a constant.                                            
c     uses unrolled loops for increment equal to 1.                             
c     jack dongarra, linpack, 6/17/77.                                          
c                                                                               
      real sa,sx(1)                                                             
      integer i,incx,m,mp1,n,nincx                                              
c                                                                               
      if(n.le.0)return                                                          
      if(incx.eq.1)go to 20                                                     
c                                                                               
c        code for increment not equal to 1                                      
c                                                                               
      nincx = n*incx                                                            
      do 10 i = 1,nincx,incx                                                    
        sx(i) = sa*sx(i)                                                        
   10 continue                                                                  
      return                                                                    
c                                                                               
c        code for increment equal to 1                                          
c                                                                               
c                                                                               
c        clean-up loop                                                          
c                                                                               
   20 m = mod(n,5)                                                              
      if( m .eq. 0 ) go to 40                                                   
      do 30 i = 1,m                                                             
        sx(i) = sa*sx(i)                                                        
   30 continue                                                                  
      if( n .lt. 5 ) return                                                     
   40 mp1 = m + 1                                                               
      do 50 i = mp1,n,5                                                         
        sx(i) = sa*sx(i)                                                        
        sx(i + 1) = sa*sx(i + 1)                                                
        sx(i + 2) = sa*sx(i + 2)                                                
        sx(i + 3) = sa*sx(i + 3)                                                
        sx(i + 4) = sa*sx(i + 4)                                                
   50 continue                                                                  
      return                                                                    
      end  

c***********************************************************
      REAL FUNCTION median(list,num)
c***********************************************************

      INTEGER num
      REAL list(num),listord(num)

! re-order the vector in a monoton way
cc      call qcksrt(list,listord,num)     
      call ssort(list,listord,num)  

! Calculate median
      IF (mod(num,2) .NE. 0) THEN
         median = list(num/2 + 1)
      ELSE
         median = ( list(num/2) + list(num/2+1) )/2.0
      END IF
      RETURN
      END FUNCTION median

c***********************************************************
      REAL FUNCTION median_no0(list,num)
c  modif by JMM in order not to make average with 0 value
c***********************************************************

      INTEGER num
      REAL list(num),listord(num)

! re-order the vector in a monoton way
cc      call qcksrt(list,listord,num)     
      call ssort(list,listord,num)

! Calculate median
      IF (mod(num,2) .NE. 0) THEN
         median_no0 = list(num/2 + 1)
      ELSE
         median_no0 = ( list(num/2) + list(num/2+1) )/2.0
         IF ( list(num/2) .EQ. 0 .OR. list(num/2+1) .EQ. 0 ) THEN 
           median_no0 = 0.
         ENDIF
      END IF
      RETURN
      END FUNCTION median_no0


c***********************************************************
      subroutine qcksrt(arr1,arr2,n)
c***********************************************************

      implicit none
      integer m,nstack
      real fm,fa,fc,fmi
c      parameter (m=7,nstack=500,fm=7875.,fa=211.,fc=1663.,fmi=1./fm)
      parameter (m=60,nstack=100,fm=7875.,fa=211.,fc=1663.,fmi=1./fm)
      integer l,ir,j,i,iq,n,jstack
      integer istack(nstack)
      real fx,a,b
      real arr1(n),arr2(n)

      jstack=0
      l=1
      ir=n
      fx=0.
 10   if(ir-l .lt.m)then
        do j=l+1,ir
           a=arr1(j)
           b=arr2(j)
           do i=j-1,1,-1
              if (arr1(i).le.a) goto 12
              arr1(i+1)=arr1(i)
              arr2(i+1)=arr2(i)
           enddo
           i=0
 12        arr1(i+1)=a
           arr2(i+1)=b
        enddo
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
        fx=MOD(fx*fa+fc,fm)
        iq=l+(ir-l+1)*(fx*fmi)
        a=arr1(iq)
        arr1(iq)=arr1(l)
        b=arr2(iq)
        arr2(iq)=arr2(l)
 20     continue
 21     if (j.gt.0) then
           if(a.lt.arr1(j))then
             j=j-1
             goto 21
           endif
        endif
        if (j.le.i)then
           arr1(i)=a
           arr2(i)=b
           goto 30
        endif
        arr1(i)=arr1(j)
        arr2(i)=arr2(j)
        i=i+1
 22     if (i.le.n) then
           if (a.gt.arr1(i)) then
              i=i+1
              goto 22
           endif
        endif
        if (j.le.i)then
           arr1(j)=a
           arr2(j)=b
           i=j
           goto 30
        endif
        arr1(j)=arr1(i)
        arr2(j)=arr2(i)
        j=j-1
        goto 20
 30     jstack=jstack+2
        if(jstack.gt.nstack)pause 'nstack trop petit'
        if (ir-i.ge.i-1) then
           istack(jstack)=ir
           istack(jstack-1)=i+1
           ir=i-1
        else
           istack(jstack)=i-1
           istack(jstack-1)=l
           l=i+1
        endif
      endif
      goto 10
      end

c***********************************************************    
      SUBROUTINE SSORT (X, IY, N)
c***********************************************************
c
c    Example of a Bubble Sort
c
C***BEGIN PROLOGUE  SSORT
C***PURPOSE  Sort an array and make the same interchanges in
C            an auxiliary array.  The array is sorted in
C            decreasing order.
C***TYPE      SINGLE PRECISION
C***KEYWORDS  SORT, SORTING
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      IY - array to be carried with X (all swaps of X elements are
C          matched in IY .  After the sort IY(J) contains the original
C          postition of the value X(J) in the unsorted X array.
C      N - number of values in array X to be sorted
C      KFLAG - Not used in this implementation
C
C***REVISION HISTORY  (YYMMDD)
C   950310  DATE WRITTEN
C   John Mahaffy
C***END PROLOGUE  SSORT
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL X(*)
      INTEGER IY(*)
C     .. Local Scalars ..
      REAL TEMP
      INTEGER I, J, JMAX, ITEMP
C     .. External Subroutines ..
C     None
C     .. Intrinsic Functions ..
C     None
C
C***FIRST EXECUTABLE STATEMENT  SSORT
C
      JMAX=N-1
      DO 200 I=1,N-1
         TEMP=1.E38
         DO 100 J=1,JMAX
            IF(X(J).GT.X(J+1)) GO TO 100
              TEMP=X(J)
              X(J)=X(J+1)
              X(J+1)=TEMP
              ITEMP=IY(J)
              IY(J)=IY(J+1)
              IY(J+1)=ITEMP
  100    CONTINUE
         IF(TEMP.EQ.1.E38) GO TO 300
         JMAX=JMAX-1
  200 CONTINUE
  300 RETURN
      END
c***********************************************************    
      SUBROUTINE MYSSORT (X, IY, N)
c***********************************************************
c
c    Example of a Bubble Sort
c
C***BEGIN PROLOGUE  MYSSORT
C
C   SSORT is changed by suppressing the go to command 
C   to optimize the code
C
C***REVISION HISTORY  (YYMMDD)
C   950310  DATE WRITTEN
C   John Mahaffy
C
C   030425 Elisabeth REMY
C
C***END PROLOGUE  SSORT
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL X(*)
      INTEGER IY(*)
C     .. Local Scalars ..
      REAL TEMP
      INTEGER I, J, JMAX, ITEMP
C     .. External Subroutines ..
C     None
C     .. Intrinsic Functions ..
C     None
C
C***FIRST EXECUTABLE STATEMENT  SSORT
C
      JMAX=N-1
      DO 200 I=1,N-1
         DO 100 J=1,JMAX
C            IF(X(J).GT.X(J+1)) GO TO 100
            IF(X(J).LE.X(J+1)) THEN
              TEMP=X(J)
              X(J)=X(J+1)
              X(J+1)=TEMP

              ITEMP=IY(J)
              IY(J)=IY(J+1)
              IY(J+1)=ITEMP
            ENDIF
  100    CONTINUE
         JMAX=JMAX-1
  200 CONTINUE

      RETURN
      END
