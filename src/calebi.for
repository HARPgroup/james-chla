C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEBI
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  CALEBI CALCULATES THE EXTERNAL BUOYANCY INTEGRALS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C     DIMENSION CH(LCM,KCM)
C
C**********************************************************************C
C
      DO L=2,LA
      BI1(L)=0.
      BI2(L)=0.
      CH(L,KC)=DZC(KC)*B(L,KC)
      BE(L)=GP*DZC(KC)*B(L,KC)
      END DO
C
      DO K=KS,1,-1
      DO L=2,LA
      CH(L,K)=CH(L,K+1)+DZC(K)*B(L,K)
      BE(L)=BE(L)+GP*DZC(K)*B(L,K)
      END DO
      END DO
C
      DO K=1,KC
      DO L=2,LA
      BI1(L)=BI1(L)+GP*DZC(K)*(CH(L,K)-0.5*DZC(K)*B(L,K))
      BI2(L)=BI2(L)+GP*DZC(K)*(CH(L,K)+Z(K-1)*B(L,K))
      END DO
      END DO
C
C**********************************************************************C
C
      RETURN
      END
