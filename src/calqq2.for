C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN SPLIT ADVECTION SCHEME
C
C1000 CONTINUE
C
C**********************************************************************C
C
C **  FIRST SPLIT HORIZONTAL ADVECTIVE FLUX CALCULATION
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     FUHU(L,K)=0.5*(UHDY2(L,K)*(CON1(L,K)+CON1(L-1,K))
C    $          +ABS(UHDY2(L,K))*(CON1(L-1,K)-CON1(L,K)))
C    $          +0.5*SUB(L)*DYU(L)*H1U(L)*(AH(L,K)+AH(L-1,K))*
C    $          (CON1(L-1,K)-CON1(L,K))*DXIU(L)
C     FVHU(L,K)=0.5*(VHDX2(L,K)*(CON1(L,K)+CON1(LS,K))
C    $          +ABS(VHDX2(L,K))*(CON1(LS,K)-CON1(L,K)))
C    $          +0.5*SVB(L)*DXV(L)*H1V(L)*(AH(L,K)+AH(LS,K))*
C    $          (CON1(LS,K)-CON1(L,K))*DYIV(L)
C     END DO
C     END DO
C
C     IF (ISTL.EQ.2) THEN
C      DO L=2,LA
C      HTMP(L)=0.5*(HP(L)+H1P(L))
C      END DO
C     ELSE
C      DO L=2,LA
C      HTMP(L)=0.5*(HP(L)+H2P(L))
C      END DO
C     END IF
C
C**********************************************************************C
C
C **  FIRST SPLIT HORIZONTAL ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     CH(L,K)=CON1(L,K)*(S2TL*H1P(L)+S3TL*H2P(L))+DELTD2*(FUHU(L,K)
C    $       -FUHU(L+1,K)+FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
C     END DO
C     END DO
C
C     IF (ISUD.EQ.1) THEN
C     DO K=1,KC
C     DO L=2,LA
C     CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
C     END DO
C     END DO
C     END IF
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)/HTMP(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING SALINITY OR SPECIFY INFLOW SALINITY
C **  AT OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBS
C     L=LCBS(LL)
C     LN=LNC(L)
C
C     IF (VHDX2(LN,K).LT.0.) THEN
C      CTMP=CON1(L,K)+DELTD2*(VHDX2(LN,K)*CON1(L,K)
C    $                     -FVHU(LN,K))*DXYIP(L)/HTMP(L)
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
C      CLOS(LL,K,M)=CON(L,K)
C      NLOS(LL,K,M)=N
C     ELSE
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
C      NMNLO=N-NLOS(LL,K,M)
C      IF (NMNLO.GE.NTSCRS(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOS(LL,K,M)
C    $         +(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBW
C     L=LCBW(LL)
C
C     IF (UHDY2(L+1,K).LT.0.) THEN
C      CTMP=CON1(L,K)+DELTD2*(UHDY2(L+1,K)*CON1(L,K)
C    $                     -FUHU(L+1,K))*DXYIP(L)/HTMP(L)
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
C      CLOW(LL,K,M)=CON(L,K)
C      NLOW(LL,K,M)=N
C     ELSE
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
C      NMNLO=N-NLOW(LL,K,M)
C      IF (NMNLO.GE.NTSCRW(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOW(LL,K,M)
C    $        +(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBE
C     L=LCBE(LL)C
C
C     IF (UHDY2(L,K).GT.0.) THEN
C      CTMP=CON1(L,K)+DELTD2*(FUHU(L,K)
C    $C    C    C    C     -UHDY2(L,K)*CON1(L,K))*DXYIP(L)/HTMP(L)
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
C      CLOE(LL,K,M)=CON(L,K)
C      NLOE(LL,K,M)=N
C     ELSE
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
C      NMNLO=N-NLOE(LL,K,M)
C      IF (NMNLO.GE.NTSCRE(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOE(LL,K,M)
C    $C        +(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBN
C     L=LCBN(LL)
C     LS=LSC(L)
C
C     IF (VHDX2(L,K).GT.0.) THEN
C      CTMP=CON1(L,K)+DELTD2*(FVHU(L,K)
C    $C    C    C    C     -VHDX2(L,K)*CON1(L,K))*DXYIP(L)/HTMP(L)
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
C      CLON(LL,K,M)=CON(L,K)
C      NLON(LL,K,M)=N
C     ELSE
C      IF (ISUD.EQ.1) CON1(L,K)=CON(L,K)
C      CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
C      NMNLO=N-NLON(LL,K,M)
C      IF (NMNLO.GE.NTSCRN(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLON(LL,K,M)
C    $C        +(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  FIRST SPLIT HORIZONTAL ANTI-DIFFUSION FLUX CALCULATION
C
C     IF (ISADAH(M).EQ.0) GO TO 1100
C
C----------------------------------------------------------------------C
C
C      DO K=1,KC
C      DO L=2,LA
C      LS=LSC(L)C    C
C      UUU(L,K)=U2(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
C      VVV(L,K)=V2(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
C      END DO
C      END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)C    C
C     LNW=LNWC(L)
C     LSE=LSEC(L)
C     AUHU=ABS(UHDY2(L,K))
C     AVHV=ABS(VHDX2(L,K))
C     UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
C     VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
CC    AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
CC    AHU(L,K)=MAX(AHTMP,0.)
CC    AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
CC    AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
C     VTERM=VTERM*(CON(L,K)-CON(LS,K))
C     UTERM=UTERM-0.25*DELTD2*UHDY2(L,K)*
C    $C     (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K))
C     VTERM=VTERM-0.25*DELTD2*VHDX2(L,K)*
C    $C     (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K))
C     UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
C     VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
CC    AUHU=ABS(UHU)
CC    AVHV=ABS(VHV)
CC    UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
CC    VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
CC    AHTMP=UTERM*DXU(L)/(DYU(L)*HU(L))
CC    AHU(L,K)=MAX(AHTMP,0.)
CC    AHTMP=VTERM*DYV(L)/(DXV(L)*HV(L))
CC    AHV(L,K)=MAX(AHTMP,0.)
C     FUHU(L,K)=0.5*(UHU*(CON(L,K)+CON(L-1,K))
C    $C    C    +ABS(UHU)*(CON(L-1,K)-CON(L,K)))
C     FVHU(L,K)=0.5*(VHV*(CON(L,K)+CON(LS,K))
C    $C    C    +ABS(VHV)*(CON(LS,K)-CON(L,K)))
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  FIRST SPLIT HORIZONTAL ANTI-DIFFUSIVE ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C
C     DO LL=1,NCBS
C     L=LCBS(LL)
C     LN=LNC(L)
C     FVHU(LN,K)=0.0
C     END DO
C
C     DO LL=1,NCBW
C     L=LCBW(LL)
C     FUHU(L+1,K)=0.0
C     END DO
C
C     DO LL=1,NCBE
C     L=LCBE(LL)
C     FUHU(L,K)=0.0
C     END DO
C
C     DO LL=1,NCBN
C     L=LCBN(LL)
C     FVHU(L,K)=0.0
C     END DO
C
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     CH(L,K)=CON(L,K)*HTMP(L)+DELTD2*(FUHU(L,K)-FUHU(L+1,K)
C    $C      +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
C     END DO
C     END DO
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)/HTMP(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  VERTICAL ADVECTIVE FLUX CALCULATION
C
C1100 CONTINUE
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY. IF ISCDCA
C **  EQUALS 1, CALCULATE 3 TIME LEVEL FLUXES BY CENTRAL DIFFERENCE
C
C----------------------------------------------------------------------C
C
C     IF (ISTL.EQ.2) THEN
C      DO K=1,KS
C      DO L=2,LA
C      FWU(L,K)=0.5*(W2(L,K)*(CON(L,K+1)+CON(L,K))
C    $C    C    +ABS(W2(L,K))*(CON(L,K)-CON(L,K+1)))
C      END DO
C      END DO
C     ELSE
C     IF (ISCDCA(MVAR).EQ.0) THEN
C      DO K=1,KS
C      DO L=2,LA
C      FWU(L,K)=0.5*(W2(L,K)*(CON(L,K+1)+CON(L,K))
C    $C    C    +ABS(W2(L,K))*(CON(L,K)-CON(L,K+1)))
C      END DO
C      END DO
C     ELSE
C      DO K=1,KS
C      DO L=2,LA
C      FWU(L,K)=0.5*W2(L,K)*(CON(L,K+1)+CON(L,K))
C      END DO
C      END DO
C     END IF
C     END IF
C
C**********************************************************************C
C
C **  VERTICAL ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)C
C     CH(L,K)=CON(L,K)*HTMP(L)
C    $    +DELT*(FWU(L,K-1)-FWU(L,K))*DZIC(K)
C     END DO
C     END DO
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)/HTMP(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  VERTICAL ANTI-DIFFUSION FLUX CALCULATION
C
C     IF (ISADAV(M).EQ.0) GO TO 1200
C
C----------------------------------------------------------------------C
C
C     DO K=1,KS
C     DO L=2,LA
C     AWW=ABS(W2(L,K))
C     WTERM=AWW*(1.-DELT*AWW*DZIG(K)*HPI(L))*(CON(L,K+1)-CON(L,K))
C     WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
C     FWU(L,K)=0.5*(WW*(CON(L,K+1)+CON(L,K))
C    $C        +ABS(WW)*(CON(L,K)-CON(L,K+1)))
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  VERTICAL ANTI-DIFFUSIVE ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     CH(L,K)=CON(L,K)*HTMP(L)
C    $C     +DELT*(FWU(L,K-1)-FWU(L,K))*DZIC(K)
C     END DO
C     END DO
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)/HTMP(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  SECOND SPLIT HORIZONTAL ADVECTIVE FLUX CALCULATION
C
C1200 CONTINUE
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LS=LSC(L)
C     FUHU(L,K)=0.5*(UHDY2(L,K)*(CON(L,K)+CON(L-1,K))
C    $C    C    +ABS(UHDY2(L,K))*(CON(L-1,K)-CON(L,K)))
C     FVHU(L,K)=0.5*(VHDX2(L,K)*(CON(L,K)+CON(LS,K))
C    $C    C    +ABS(VHDX2(L,K))*(CON(LS,K)-CON(L,K)))
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  SECOND SPLIT HORIZONTAL ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)C
C     CH(L,K)=CON(L,K)*HTMP(L)+DELTD2*(FUHU(L,K)-FUHU(L+1,K)
C    $C      +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
C     END DO
C     END DO
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING SALINITY OR SPECIFY INFLOW SALINITY
C **  AT OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBS
C     L=LCBS(LL)
C     LN=LNC(L)
C
C     IF (VHDX2(LN,K).LT.0.) THEN
C      CON(L,K)=CON(L,K)+DELTD2*(VHDX2(LN,K)*CON(L,K)
C    $C    C    C    C      -FVHU(LN,K))*DXYIP(L)*HPI(L)
C      IF (CON(L,K).GT.CBS(LL,1,M)) CON(L,K)=CBS(LL,1,M)
C      CLOS(LL,K,M)=CON(L,K)
C      NLOS(LL,K,M)=N
C     ELSE
C      CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)
C      NMNLO=N-NLOS(LL,K,M)
C      IF (NMNLO.GE.NTSCRS(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOS(LL,K,M)
C    $C        +(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBW
C     L=LCBW(LL)C
C
C     IF (UHDY2(L+1,K).LT.0.) THEN
C      CON(L,K)=CON(L,K)+DELTD2*(UHDY2(L+1,K)*CON(L,K)
C    $C    C    C    C      -FUHU(L+1,K))*DXYIP(L)*HPI(L)
C      IF (CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
C      CLOW(LL,K,M)=CON(L,K)
C      NLOW(LL,K,M)=N
C     ELSE
C      CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)
C      NMNLO=N-NLOW(LL,K,M)
C      IF (NMNLO.GE.NTSCRW(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOW(LL,K,M)
C    $C        +(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBE
C     L=LCBE(LL)C
C
C     IF (UHDY2(L,K).GT.0.) THEN
C      CON(L,K)=CON(L,K)+DELTD2*(FUHU(L,K)
C    $C    C    C    C      -UHDY2(L,K)*CON(L,K))*DXYIP(L)*HPI(L)
C      IF (CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
C      CLOE(LL,K,M)=CON(L,K)
C      NLOE(LL,K,M)=N
C     ELSE
C      CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)
C      NMNLO=N-NLOE(LL,K,M)
C      IF (NMNLO.GE.NTSCRE(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLOE(LL,K,M)
C    $C        +(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO LL=1,NCBN
C     L=LCBN(LL)
C     LS=LSC(L)
C
C     IF (VHDX2(L,K).GT.0.) THEN
C      CON(L,K)=CON(L,K)+DELTD2*(FVHU(L,K)
C    $C    C    C    C    C     -VHDX2(L,K)*CON(L,K))*DXYIP(L)*HPI(L)
C      IF (CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
C      CLON(LL,K,M)=CON(L,K)
C      NLON(LL,K,M)=N
C     ELSE
C      CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)
C      NMNLO=N-NLON(LL,K,M)
C      IF (NMNLO.GE.NTSCRN(LL)) THEN
C       CON(L,K)=CBT
C      ELSE
C       CON(L,K)=CLON(LL,K,M)
C    $C        +(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
C      END IF
C     END IF
C
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  SECOND SPLIT HORIZONTAL ANTI-DIFFUSIVE FLUX CALCULATION
C
C     IF (ISADAH(M).EQ.0) RETURN
C
C----------------------------------------------------------------------C
C
C      DO K=1,KC
C      DO L=2,LA
C      LS=LSC(L)C    C
C      UUU(L,K)=U2(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
C      VVV(L,K)=V2(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
C      END DO
C      END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)C    C
C     LNW=LNWC(L)
C     LSE=LSEC(L)
C     AUHU=ABS(UHDY2(L,K))
C     AVHV=ABS(VHDX2(L,K))
C     UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
C     VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
CC    AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
CC    AHU(L,K)=MAX(AHTMP,0.)
CC    AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
CC    AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(CON(L,K)-CON(L-1,K))
C     VTERM=VTERM*(CON(L,K)-CON(LS,K))
C     UTERM=UTERM-0.25*DELTD2*UHU*
C    $C     (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K))
C     VTERM=VTERM-0.25*DELTD2*VHV*
C    $C     (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K))
C     UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
C     VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
CC    AUHU=ABS(UHU)
CC    AVHV=ABS(VHV)
CC    UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
CC    VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
CC    AHTMP=UTERM*DXU(L)/(DYU(L)*HU(L))
CC    AHU(L,K)=MAX(AHTMP,0.)
CC    AHTMP=VTERM*DYV(L)/(DXV(L)*HV(L))
CC    AHV(L,K)=MAX(AHTMP,0.)
C     FUHU(L,K)=0.5*(UHU*(CON(L,K)+CON(L-1,K))
C    $C    C    +ABS(UHU)*(CON(L-1,K)-CON(L,K)))
C     FVHU(L,K)=0.5*(VHV*(CON(L,K)+CON(LS,K))
C    $C    C    +ABS(VHV)*(CON(LS,K)-CON(L,K)))
C     END DO
C     END DO
C
C**********************************************************************C
C
C **  SECOND SPLIT HORIZONTAL ANTI-DIFFUSIVE ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C
C     DO LL=1,NCBS
C     L=LCBS(LL)
C     LN=LNC(L)
C     FVHU(LN,K)=0.0
C     END DO
C
C     DO LL=1,NCBW
C     L=LCBW(LL)
C     FUHU(L+1,K)=0.0
C     END DO
C
C     DO LL=1,NCBE
C     L=LCBE(LL)
C     FUHU(L,K)=0.0
C     END DO
C
C     DO LL=1,NCBN
C     L=LCBN(LL)
C     FVHU(L,K)=0.0
C     END DO
C
C     END DO
C
C----------------------------------------------------------------------C
C
C     DO K=1,KC
C     DO L=2,LA
C     LN=LNC(L)
C     CH(L,K)=CON(L,K)*HP(L)+DELTD2*(FUHU(L,K)-FUHU(L+1,K)
C    $C      +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
C     END DO
C     END DO
C
C     DO K=1,KC
C     DO L=2,LA
C     CON(L,K)=SCB(L)*CH(L,K)*HPI(L)+(1.-SCB(L))*CON(L,K)
C     END DO
C     END DO
C
C**********************************************************************C
C
c     RETURN
c     END
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALQQ2 (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C **  SUBROUTINE CALQQ2 CALCULATES THE TURBULENT INTENSITY SQUARED AT
C **  TIME LEVEL (N+1).  THE VALUE OF ISTL INDICATES THE NUMBER OF
C **  TIME LEVELS INVOLVED.  THIS VERSION USES A SEPARATE ADVECTIVE
C **  TRANSPORT SUBROUTINE CALTRANQ
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF (ISTL.EQ.2) THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      END IF
      BSMALL=1.E-12
C
C**********************************************************************C
C
C **  MOVE UHDY2, VHDX2 AND W2 TO QQ TRANSPORT LOCATIONS
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      U2(L,K)=0.5*(U2(L,K)+U2(L,K+1))
      V2(L,K)=0.5*(V2(L,K)+V2(L,K+1))
      UHDY2(L,K)=0.5*(UHDY2(L,K)+UHDY2(L,K+1))
      VHDX2(L,K)=0.5*(VHDX2(L,K)+VHDX2(L,K+1))
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      W2(L,K)=0.5*(W2(L,K)+W2(L,K+1))
      END DO
      END DO
C
C**********************************************************************C
C
C **  ELIMINATE INFLOWS ACROSS OPEN BOUNDARIES
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF (VHDX2(LN,K).GT.0.) VHDX2(LN,K)=0.
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBW
      L=LCBW(LL)
      IF (UHDY2(L+1,K).GT.0.) UHDY2(L+1,K)=0.
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBE
      L=LCBE(LL)
      IF (UHDY2(L,K).LT.0.) UHDY2(L,K)=0.
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      IF (VHDX2(L,K).LT.0.) VHDX2(L,K)=0.
      END DO
      END DO
C
C**********************************************************************C
C
C **  CALL ADVECTIVE TRANSPORT SUBROUTINE FOR QQ AND QQL
C
C----------------------------------------------------------------------C
C
      CALL CALTRANQ (ISTL,0,QQ,QQ1)
C
      CALL CALTRANQ (ISTL,0,QQL,QQL1)
C
      DO L=2,LA
      W2(L,0)=0.
      W2(L,KC)=0.
      END DO
C
C**********************************************************************C
C
C **  CALCULATE PRODUCTION, LOAD BOUNDARY CONDITIONS AND SOLVE
C **  TRANSPORT EQUATIONS
C
C **  CU1=CUQ, CU2=CUQL, UUU=QQH, VVV=QQLH
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      PQQ=DELT*(AB(L,K)*GP*HP(L)*DZIG(K)*(B(L,K+1)-B(L,K))
     $    +AV(L,K)*DZIGSD4(K)*(U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $    +AV(L,K)*DZIGSD4(K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2)
      UUU(L,K)=QQ(L,K)*HP(L)+2.*PQQ
      VVV(L,K)=QQL1(L,K)*HP(L)+CTE1*DML(L,K)*PQQ
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      CLQTMP=-DELT*CDZKK(1)*AQ(L,1)*HPI(L)
      CUQTMP=-DELT*CDZKKP(1)*AQ(L,2)*HPI(L)
      CMQTMP=1.-CLQTMP-CUQTMP
     $       +2.*DELT*SQRT(QQ(L,1))/(CTURB*DML(L,1)*HP(L))
      CMQLTMP=1.-CLQTMP-CUQTMP
     $       +DELT*(SQRT(QQ(L,1))/(CTURB*DML(L,1)*HP(L)))*(1.
     $       +CTE2*DML(L,1)*DML(L,1)*FPROX(1))
      EQ=1./CMQTMP
      EQL=1./CMQLTMP
      CU1(L,1)=CUQTMP*EQ
      CU2(L,1)=CUQTMP*EQL
      UUU(L,1)=(UUU(L,1)-CLQTMP*HP(L)*QQ(L,0))*EQ
      VVV(L,1)=VVV(L,1)*EQL
      CUQTMP=-DELT*CDZKKP(KS)*AQ(L,KC)*HPI(L)
      UUU(L,KS)=UUU(L,KS)-CUQTMP*HP(L)*QQ(L,KC)
      END DO
C
      DO K=2,KS
      DO L=2,LA
      CLQTMP=-DELT*CDZKK(K)*AQ(L,K)*HPI(L)
      CUQTMP=-DELT*CDZKKP(K)*AQ(L,K+1)*HPI(L)
      CMQTMP=1.-CLQTMP-CUQTMP
     $       +2.*DELT*SQRT(QQ(L,K))/(CTURB*DML(L,K)*HP(L))
      CMQLTMP=1.-CLQTMP-CUQTMP
     $       +DELT*(SQRT(QQ(L,K))/(CTURB*DML(L,K)*HP(L)))*(1.
     $       +CTE2*DML(L,K)*DML(L,K)*FPROX(K))
      EQ=1./(CMQTMP-CLQTMP*CU1(L,K-1))
      EQL=1./(CMQLTMP-CLQTMP*CU2(L,K-1))
      CU1(L,K)=CUQTMP*EQ
      CU2(L,K)=CUQTMP*EQL
C     IF (EQ.0) PAUSE
      UUU(L,K)=(UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
      VVV(L,K)=(VVV(L,K)-CLQTMP*VVV(L,K-1))*EQL
      END DO
      END DO
C
      DO K=KS-1,1,-1
      DO L=2,LA
      UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
      VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      QQ1(L,K)=S2TL*QQ1(L,K)+S3TL*QQ(L,K)
      QQHDH=UUU(L,K)*HPI(L)
      QQ(L,K)=MAX(QQHDH,QQMIN)
      QQ(L,K)=SPB(L)*QQ(L,K)+(1.-SPB(L))*QQMIN
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      QQL1(L,K)=S2TL*QQL1(L,K)+S3TL*QQL(L,K)
      QQHDH=VVV(L,K)*HPI(L)
      QQL(L,K)=MAX(QQHDH,QQLMIN)
      QQL(L,K)=SPB(L)*QQL(L,K)+(1.-SPB(L))*QQLMIN
      DMLTMP=QQL(L,K)/QQ(L,K)
      DMLTMP=MAX(DMLTMP,DMLMIN)
      DELB=B(L,K)-B(L,K+1)
      IF(DELB.GT.BSMALL) THEN
       DMLMAX=CTE3*SQRT(QQ(L,K)/(G*HP(L)*DZIG(K)*DELB))
       DML(L,K)=MIN(DMLMAX,DMLTMP)
      ELSE
       DML(L,K)=DMLTMP
      END IF
      DML(L,K)=SPB(L)*DML(L,K)+(1.-SPB(L))*DMLMIN
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
c     DO K=1,KS
c     DO LL=1,NCBS
c     LN=LNC(L)
c     QQ(L,K)=QQ(LN,K)
c     QQL(L,K)=QQL(LN,K)
c     DML(L,K)=DML(LN,K)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
c     DO K=1,KS
c     DO LL=1,NCBW
c     L=LCBW(LL)
c     QQ(L,K)=QQ(L+1,K)
c     QQL(L,K)=QQL(L+1,K)
c     DML(L,K)=DML(L+1,K)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
c     DO K=1,KS
c     DO LL=1,NCBE
c     L=LCBE(LL)
c     QQ(L,K)=QQ(L-1,K)
c     QQL(L,K)=QQL(L-1,K)
c     DML(L,K)=DML(L-1,K)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
c     DO K=1,KS
c     DO LL=1,NCBN
c     L=LCBN(LL)
c     LS=LSC(L)
c     QQ(L,K)=QQ(LS,K)
c     QQL(L,K)=QQL(LS,K)
c     DML(L,K)=DML(LS,K)
c     END DO
c     END DO
C
C----------------------------------------------------------------------C
C
      RETURN
      END
