C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALHDMF
C
C **  SUBROUTINE CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND
C **  DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING
C **  SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      AHMAX=AHO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE HORIZONTAL VELOCITY SHEARS
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      DXU1(L,K)=(U1(L+1,K)-U1(L,K))/DXP(L)
      DYV1(L,K)=(V1(LN,K)-V1(L,K))/DYP(L)
      DYU1(L,K)=2.*SUB(LS)*(U1(L,K)-U1(LS,K))/(DYU(L)+DYU(LS))
      DXV1(L,K)=2.*SVB(L-1)*(V1(L,K)-V1(L-1,K))/(DXV(L)+DXV(L-1))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE HORIZONTAL VISCOSITY
C
      IF(AHD.GT.0.0) THEN
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LNE=LNEC(L)
      LSW=LSWC(L)
      LSE=LSEC(L)
      AH(L,K)=AHO+AHD*DXP(L)*DYP(L)*SQRT( 2.*DXU1(L,K)*DXU1(L,K)
     $                                   +2.*DYV1(L,K)*DYV1(L,K)
     $  +0.0625*(DYU1(L,K)+DYU1(LN,K)+DYU1(L+1,K)+DYU1(LNE,K)
     $          +DXV1(L,K)+DXV1(LN,K)+DXV1(L+1,K)+DXV1(LNE,K))**2 )
      AHC(L,K)=AHO+AHD*0.0625*( (DXP(L)+DXP(L-1)+DXP(LS)+DXP(LSW))**2 )
     $   *SQRT( 0.125*(DXU1(L,K)+DXU1(L-1,K)+DXU1(LS,K)+DXU1(LSW,K))**2
     $         +0.125*(DYV1(L,K)+DYV1(L-1,K)+DYV1(LS,K)+DYV1(LSW,K))**2
     $         +DYU1(L,K)*DYU1(L,K)+2.*DYU1(L,K)*DXV1(L,K)
     $         +DXV1(L,K)*DXV1(L,K) )
      END DO
      END DO
C
      ELSE
C
      DO K=1,KC
      DO L=2,LA
       AH(L,K)=AHO
       AHC(L,K)=AHO
      END DO
      END DO
C
      END IF
C
C----------------------------------------------------------------------C
C
C **  CALCULATE HORIZONTAL DIFFUSION DUE TO WAVE BREAKING
C
      IF (ISWAVE.GE.1) THEN
      IF (WVLSH.GT.0.0.OR.WVLSX.GT.0.0) THEN
C
       IF (N.LT.NTSWV) THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
        ELSE
         WVFACT=1.0
       END IF
C
c     DO K=1,KC
c     DO L=2,LA
c     DTMP=(WVFACT*WVDISP(L,K)/HMP(L))**0.3333
c     RLTMP=WVLSH*( HMP(L)**1.3333)+WVLSX*( DXYP(L)**0.6667 )
c     AH(L,K)=AH(L,K)+DTMP*RLTMP
c     END DO
c     END DO
C
c     DO K=2,KC
c     DO L=2,LA
c     LS=LSC(L)
c     LSW=LSWC(L)
c     DTMP=0.25*WVFACT*(WVDISP(L,K)+WVDISP(L-1,K)
c    $          +WVDISP(LS,K)+WVDISP(LSW,K))
c     HCTMP=0.25*(HMP(L)+HMP(L-1)
c    $          +HMP(LS)+HMP(LSW))
c     DTMP=(DTMP/HCTMP)**0.3333
c     RLTMP=WVLSH*(HCTMP**1.3333)+WVLSX*( DXYP(L)**0.6667 )
c     AHC(L,K)=AHC(L,K)+DTMP*RLTMP
c     END DO
c     END DO
C
C
      AHWVX=WVFACT*WVLSX*WVPRD*WVPRD
C
      DO K=1,KC
      DO L=2,LA
      DTMPH=(WVFACT*WVDISP(L,K))**0.3333
      DTMPX=WVDISP(L,K)/HMP(L)
      AH(L,K)=AH(L,K)+WVLSH*DTMPH*HMP(L)+AHWVX*DTMPX
      END DO
      END DO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LSW=LSWC(L)
      DTMP=0.25*(WVDISP(L,K)+WVDISP(L-1,K)
     $          +WVDISP(LS,K)+WVDISP(LSW,K))
      DTMPH=(WVFACT*DTMP)**0.3333
      DTMPX=DTMP/HMC(L)
      AHC(L,K)=AHC(L,K)+WVLSH*DTMPH*HMC(L)+AHWVX*DTMPX
      END DO
      END DO
C
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C **  CALCULATE DIFFUSIVE MOMENTUM FLUXES
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      FMDUX(L,K)=2.*DYP(L)*H1P(L)*AH(L,K)*DXU1(L,K)
      FMDUY(L,K)=0.5*(DXU(L)+DXU(LS))*H1C(L)*AHC(L,K)
     $              *(DYU1(L,K)+DXV1(L,K))
      FMDVY(L,K)=2.*DXP(L)*H1P(L)*AH(L,K)*DYV1(L,K)
      FMDVX(L,K)=0.5*(DYV(L)+DYV(L-1))*H1C(L)*AHC(L,K)
     $              *(DYU1(L,K)+DXV1(L,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      RETURN
      END
