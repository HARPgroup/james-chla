C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALBAL4
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C **  CALCULATE MOMENTUM AND ENERGY DISSIPATION
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
c     UUEOUT=UUEOUT+SPB(L)*SPB(L+1)*DXYU(L)
c    $      *(U(L,1)*TBX(L)-U(L,KC)*TSX(L))
c     VVEOUT=VVEOUT+SPB(L)*SPB(LN)*DXYV(L)
c    $      *(V(L,1)*TBY(L)-V(L,KC)*TSY(L))
      UUEOUT=UUEOUT+0.5*SPB(L)*DXYP(L)*(U(L,1)*TBX(L)+U(L+1,1)*TBX(L+1)
     $       -U(L,KC)*TSX(L)-U(L+1,KC)*TSX(L+1))
      VVEOUT=VVEOUT+0.5*SPB(L)*DXYP(L)*(V(L,1)*TBY(L)+V(LN,1)*TBX(LN)
     $       -V(L,KC)*TSY(L)-V(LN,KC)*TSX(LN))
      END DO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      DUTMP=0.5*( U(L,K+1)+U(L+1,K+1)-U(L,K)-U(L+1,K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUT=UUEOUT+SPB(L)*2.0*DXYP(L)*AV(L,K)
     $      *( DUTMP*DUTMP )/(DZC(K+1)+DZC(K))
      VVEOUT=VVEOUT+SPB(L)*2.0*DXYP(L)*AV(L,K)
     $      *( DVTMP*DVTMP )/(DZC(K+1)+DZC(K))
      BBEOUT=BBEOUT+SCB(L)*DXYP(L)*HP(L)
     $      *GP*AB(L,K)*(B(L,K+1)-B(L,K))
      END DO
      END DO
C 
C**********************************************************************C
C
      RETURN
      END
