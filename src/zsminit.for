C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SMINIT(LA)
C
C **  LAST MODIFIED BY JOHN HAMRICK AND MIKE MORTON ON 10 april 1999
C **  Last modified by AEE 3/18/07
C**********************************************************************C
C
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'
C
cxh      INSMICI=40
CXH      INSMRST=40
CXH      ISMORST=45
CXH      ISMOZB=46
C
      SMTSNAME(1) = 'som'
      SMTSNAME(2) = 'sim'
      SMTSNAME(3) = 'sbf'
C
      DO L=2,LA
        SMHYST(L)=.FALSE.
      END DO
C
      RETURN
      END
