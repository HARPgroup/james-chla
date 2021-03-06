C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP1D (ISTL)
C
C **  SUBROUTINE CALEXP1D CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
C **  FOR 1D CHANNEL MODE
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
       DELT=DT2
       DELTD2=DT
      IF (ISTL.EQ.2) THEN
       DELT=DT
       DELTD2=0.5*DT
      END IF
C
      DELTI=1./DELT
C
C**********************************************************************C
C
C **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
C
C----------------------------------------------------------------------C
C
      FXE(1)=0.
      FYE(1)=0.
C
      FXE(LC)=0.
      FYE(LC)=0.
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        FXE(L)=0.
        FYE(L)=0.
       END DO
      END DO
C
C**********************************************************************C
C
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.2) GO TO 200
      IF (ISTL.EQ.3) THEN
        GO TO 300
c       IF (ISCDMA.EQ.0) GO TO 300
c       IF (ISCDMA.EQ.1) GO TO 400
c       IF (ISCDMA.EQ.2) GO TO 350
      END IF
C
C**********************************************************************C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N 
C
  200 CONTINUE
C
C----------------------------------------------------------------------C
C
C      DO K=1,KC
C       DO L=2,LA
C        UHDY2(L,K)=UHDY(L,K)
C        VHDX2(L,K)=VHDX(L,K)
C        W2(L,K)=W(L,K)
C       END DO
C      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.5*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.5*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.5*(VHDX(L,K)+VHDX(LN,K))
C
       FUHU(L,K)=MAX(UHB,0.)*U(L,K)
     $         +MIN(UHB,0.)*U(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U(LS,K)
     $         +MIN(VHC,0.)*U(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V(L-1,K)
     $         +MIN(UHC,0.)*V(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V(L,K)
     $         +MIN(VHB,0.)*V(LN,K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C      DO K=1,KS
C       DO L=2,LA
C       LS=LSC(L)
C       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
C       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
C       FWU(L,K)=MAX(WU,0.)*U1(L,K)
C     $        +MIN(WU,0.)*U1(L,K+1)
C       FWV(L,K)=MAX(WV,0.)*V1(L,K)
C     $        +MIN(WV,0.)*V1(L,K+1)
C       END DO
C      END DO
C
C----------------------------------------------------------------------C
C
      GO TO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.5*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.5*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.5*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.5*(VHDX(L,K)+VHDX(LN,K))
       FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
     $         +MIN(UHB,0.)*U1(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
     $         +MIN(VHC,0.)*U1(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
     $         +MIN(UHC,0.)*V1(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
     $         +MIN(VHB,0.)*V1(LN,K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C      DO K=1,KS
C       DO L=2,LA
C       LS=LSC(L)
C       WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
C       WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
C       FWU(L,K)=MAX(WU,0.)*U1(L,K)
C     $        +MIN(WU,0.)*U1(L,K+1)
C       FWV(L,K)=MAX(WV,0.)*V1(L,K)
C     $        +MIN(WV,0.)*V1(L,K+1)
C       END DO
C      END DO
C
C----------------------------------------------------------------------C
C
      GO TO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  FIRST HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
C **  WITH TRANSPORT AT (N-1/2) AND TRANSPORTED FIELD AT (N-1)
C **  SECOND HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE  
C **  WITH TRANSPORT AT (N+1/2) AND TRANSPORTED FIELD AT (N)
C
  350 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       U2(L,K)=U1(L,K)+U(L,K)
       V2(L,K)=V1(L,K)+V(L,K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)               
       UHC=0.25*(UHDY(L,K)+UHDY(LS,K))
       UHB=0.25*(UHDY(L,K)+UHDY(L+1,K))
       VHC=0.25*(VHDX(L,K)+VHDX(L-1,K))
       VHB=0.25*(VHDX(L,K)+VHDX(LN,K))
       FUHU(L,K)=MAX(UHB,0.)*U2(L,K)
     $         +MIN(UHB,0.)*U2(L+1,K)
       FVHU(L,K)=MAX(VHC,0.)*U2(LS,K)
     $         +MIN(VHC,0.)*U2(L,K)
       FUHV(L,K)=MAX(UHC,0.)*V2(L-1,K)
     $         +MIN(UHC,0.)*V2(L,K)
       FVHV(L,K)=MAX(VHB,0.)*V2(L,K)
     $         +MIN(VHB,0.)*V2(LN,K)
       END DO
      END DO
C
c     DO K=1,KC
c     DO L=2,LA
c     LN=LNC(L)
c     LS=LSC(L)               
c     UHC=0.125*(UHDY(L,K)+UHDY(LS,K)+UHDY1(L,K)+UHDY1(LS,K))
c     UHB=0.125*(UHDY(L,K)+UHDY(L+1,K)+UHDY1(L,K)+UHDY1(L+1,K))
c     VHC=0.125*(VHDX(L,K)+VHDX(L-1,K)+VHDX1(L,K)+VHDX1(L-1,K))
c     VHB=0.125*(VHDX(L,K)+VHDX(LN,K)+VHDX1(L,K)+VHDX1(LN,K))
c     FUHU(L,K)=MAX(UHB,0.)*U1(L,K)
c    $         +MIN(UHB,0.)*U1(L+1,K)
c     FVHU(L,K)=MAX(VHC,0.)*U1(LS,K)
c    $         +MIN(VHC,0.)*U1(L,K)
c     FUHV(L,K)=MAX(UHC,0.)*V1(L-1,K)
c    $         +MIN(UHC,0.)*V1(L,K)
c     FVHV(L,K)=MAX(VHB,0.)*V1(L,K)
c    $         +MIN(VHB,0.)*V1(LN,K)
c     END DO
c     END DO
C
c     DO K=1,KC
c     DO L=2,LA
c     LN=LNC(L)
c     LS=LSC(L)               
c     UHC=0.125*(3.*UHDY(L,K)+3.*UHDY(LS,K)-UHDY1(L,K)-UHDY1(LS,K))
c     UHB=0.125*(3.*UHDY(L,K)+3.*UHDY(L+1,K)-UHDY1(L,K)-UHDY1(L+1,K))
c     VHC=0.125*(3.*VHDX(L,K)+3.*VHDX(L-1,K)-VHDX1(L,K)-VHDX1(L-1,K))
c     VHB=0.125*(3.*VHDX(L,K)+3.*VHDX(LN,K)-VHDX1(L,K)-VHDX1(LN,K))
c     FUHU(L,K)=FUHU(L,K)+MAX(UHB,0.)*U(L,K)
c    $                   +MIN(UHB,0.)*U(L+1,K)
c     FVHU(L,K)=FVHU(L,K)+MAX(VHC,0.)*U(LS,K)
c    $                   +MIN(VHC,0.)*U(L,K)
c     FUHV(L,K)=FUHV(L,K)+MAX(UHC,0.)*V(L-1,K)
c    $                   +MIN(UHC,0.)*V(L,K)
c     FVHV(L,K)=FVHV(L,K)+MAX(VHB,0.)*V(L,K)
c    $                   +MIN(VHB,0.)*V(LN,K)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
C      DO K=1,KS
C       DO L=2,LA
C       LS=LSC(L)
C       WU=0.25*DXYU(L)*(W(L,K)+W(L-1,K))
C       WV=0.25*DXYV(L)*(W(L,K)+W(LS,K))
C       FWU(L,K)=MAX(WU,0.)*U2(L,K)
C     $        +MIN(WU,0.)*U2(L,K+1)
C       FWV(L,K)=MAX(WV,0.)*V2(L,K)
C     $        +MIN(WV,0.)*V2(L,K+1)
C       END DO
C      END DO
C
c     DO K=1,KS
c     DO L=2,LA
c     LS=LSC(L)
c     WU=0.125*DXYU(L)*(W(L,K)+W(L-1,K)+W1(L,K)+W1(L-1,K))
c     WV=0.125*DXYV(L)*(W(L,K)+W(LS,K)+W1(L,K)+W1(LS,K))
c     FWU(L,K)=MAX(WU,0.)*U1(L,K)
c    $        +MIN(WU,0.)*U1(L,K+1)
c     FWV(L,K)=MAX(WV,0.)*V1(L,K)
c    $        +MIN(WV,0.)*V1(L,K+1)
c     END DO
c     END DO
c
c     DO K=1,KS
c     DO L=2,LA
c     LS=LSC(L)
c     WU=0.125*DXYU(L)*(3.*W(L,K)+3.*W(L-1,K)-W1(L,K)-W1(L-1,K))
c     WV=0.125*DXYV(L)*(3.*W(L,K)+3.*W(LS,K)-W1(L,K)-W1(LS,K))
c     FWU(L,K)=FWU(L,K)+MAX(WU,0.)*U(L,K)
c    $                 +MIN(WU,0.)*U(L,K+1)
c     FWV(L,K)=FWV(L,K)+MAX(WV,0.)*V(L,K)
c    $                 +MIN(WV,0.)*V(L,K+1)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
      GO TO 500
C
C**********************************************************************C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT 
C **  AT (N) AND TRANSPORTED FIELD AT (N)
C
  400 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)  
      FUHU(L,K)=0.25*(UHDY(L+1,K)+UHDY(L,K))*(U(L+1,K)+U(L,K))
      FVHU(L,K)=0.25*(VHDX(L,K)+VHDX(L-1,K))*(U(L,K)+U(LS,K))
      FUHV(L,K)=0.25*(UHDY(L,K)+UHDY(LS,K))*(V(L,K)+V(L-1,K))
      FVHV(L,K)=0.25*(VHDX(L,K)+VHDX(LN,K))*(V(LN,K)+V(L,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C      DO K=1,KS
C      DO L=2,LA
C      LS=LSC(L)    
C      FWU(L,K)=0.25*DXYU(L)*(W(L,K)+W(L-1,K))*(U(L,K+1)+U(L,K))
C      FWV(L,K)=0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
C      END DO
C      END DO
C
C**********************************************************************C
C
  500 CONTINUE
C
      DO K=1,KC
       DO L=1,LA
       FUHU(L,K)=STCUV(L)*FUHU(L,K)
       FVHV(L,K)=STCUV(L)*FVHV(L,K)
       END DO
      END DO
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE ACCELERATIONS 
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)      
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K))
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K))
       END DO
      END DO
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      IF (ISHDMF.GE.1) THEN
C
      DO K=1,KC
       DO L=2,LA
       LS=LSC(L)
       LN=LNC(L)
       FX(L,K)=FX(L,K)
     $        -SDX(L)*(FMDUX(L,K)-FMDUX(L-1,K)+FMDUY(LN,K)-FMDUY(L,K))
       FY(L,K)=FY(L,K)
     $        -SDY(L)*(FMDVX(L+1,K)-FMDVX(L,K)+FMDVY(L,K)-FMDVY(LS,K))
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
C
      IF (ISWAVE.GE.2) THEN
C
      IF (N.LT.NTSWV) THEN
         TMPVAL=FLOAT(N)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVAL)
       ELSE
        WVFACT=1.0
      END IF        
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FXE(L)=FXE(L)+WVFACT*SAAX(L)*FXWAVE(L,K)*DZC(K)
         FYE(L)=FYE(L)+WVFACT*SAAY(L)*FYWAVE(L,K)*DZC(K)
        END DO
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION OF INTERNAL ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C 
C      DO K=1,KC
C       DO L=2,LA
C         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
C         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
C        END DO
C      END DO
C
C----------------------------------------------------------------------C
C 
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF (ISWAVE.GE.2) THEN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         FX(L,K)=FX(L,K)+WVFACT*SAAX(L)*FXWAVE(L,K)
         FY(L,K)=FY(L,K)+WVFACT*SAAY(L)*FYWAVE(L,K)
        END DO
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
C **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FBBX(L,K)=ROLD*FBBX(L,K)+RNEW*SBX(L)*GP*HU(L)*
     $            ( HU(L)*( (B(L,K+1)-B(L-1,K+1))*DZC(K+1)
     $                     +(B(L,K)-B(L-1,K))*DZC(K) )
     $           +(B(L,K+1)-B(L,K)+B(L-1,K+1)-B(L-1,K))*
     $            (HMP(L)-HMP(L-1)-Z(K)*(HP(L)-HP(L-1))) )
      FBBY(L,K)=ROLD*FBBY(L,K)+RNEW*SBY(L)*GP*HV(L)*
     $            ( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(K+1)
     $                     +(B(L,K)-B(LS,K))*DZC(K) )
     $           +(B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*
     $            (HMP(L)-HMP(LS)-Z(K)*(HP(L)-HP(LS))) )
      END DO
      END DO
C
C     IF(N.EQ.1) THEN
C       OPEN(1,FILE='buoy.diag',STATUS='UNKNOWN')
C       DO L=2,LA
C        DO K=1,KS
C        TMP3D(K)=SUBO(L)*FBBX(L,K)
C        END DO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C        DO K=1,KS
C        TMP3D(K)=SVBO(L)*FBBY(L,K)
C        END DO
C       WRITE(1,1111)IL(L),JL(L),(TMP3D(K),K=1,KS)
C       END DO
C       CLOSE(1)
C     END IF
C
 1111 FORMAT(2I5,2X,8E12.4)       
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR	EQUATION TERMS
C
C----------------------------------------------------------------------C
C
C      DO K=1,KS
C       RCDZF=CDZF(K)
C       DO L=2,LA
C        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
C     $           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
C     $           +SNLT*(FX(L,K)-FX(L,K+1))) )
C        DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
C     $           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
C     $           +SNLT*(FY(L,K)-FY(L,K+1))) )
C       END DO
C      END DO
C
C      IF(ISTL.EQ.2) THEN
C 
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO L=LF,LL
C        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
C        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
C       END DO
C      END DO
C
C      END IF
C
C**********************************************************************C
C
      RETURN
      END
