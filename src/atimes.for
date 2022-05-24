C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE ATIMES(PCGM,PTMPM)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION PCGM(LCM),PTMPM(LCM)
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PTMPM(L)=CC(L)*PCGM(L)+CS(L)*PCGM(LS)+CW(L)*PCGM(L-1)
     $       +CE(L)*PCGM(L+1)+CN(L)*PCGM(LN)
      END DO
      PTMPM(1)=0.
      PTMPM(LC)=0.
C
      RETURN
      END
