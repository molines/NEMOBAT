      SUBROUTINE SSORT (X, IY, N, KFLAG)
C  C! -------------------------------------------------------------------------
C  !!  $Date: 2009-04-29 15:52:46 +0200 (Wed, 29 Apr 2009) $
C  !!  $Rev: 225 $
C  !!  $Id: sort2.f 225 2009-04-29 13:52:46Z molines $
C  !! -------------------------------------------------------------------------
      IMPLICIT NONE
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
