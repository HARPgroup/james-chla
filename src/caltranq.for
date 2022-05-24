C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTRANQ (ISTL,M,QCON,QCON1)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE
C **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
C **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES
C **  THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION QCON(LCM,0:KCM), QCON1(LCM,0:KCM)
C     DIMENSION  QCH(LCM,0:KCM), QCMAX(LCM,0:KCM), QCMIN(LCM,0:KCM),
C    $         QCONT(LCM,0:KCM), QCON2(LCM,0:KCM),
C
C**********************************************************************C
C
      BSMALL=1.0E-6
C
      DELT=DT2
      DELTA=DT2
      IF(ISCDCA(M).EQ.2) DELTA=DT
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF (ISTL.NE.3) THEN
       DELT=DT
       DELTA=DT
       DELTD2=0.5*DT
       S3TL=0.0
       S2TL=1.0
       ISUD=0
      END IF
C
      IF (ISLTMT.GE.1) ISUD=1
C
C**********************************************************************C
C
      DO K=0,KC
      DO L=1,LC
      QCONT(L,K)=0.
      QCMAX(L,K)=0.
      QCMIN(L,K)=0.
      END DO
      END DO
C
C**********************************************************************C
C
C **  ADVECTIVE FLUX CALCULATION
C
      IF (ISTL.EQ.2) GO TO 300
      IF (ISCDCA(M).EQ.0) GO TO 300
      IF (ISCDCA(M).EQ.1) GO TO 400
      IF (ISCDCA(M).EQ.2) GO TO 350
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*QCON1(L-1,K)
     $         +MIN(UHDY2(L,K),0.)*QCON1(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*QCON1(LS,K)
     $         +MIN(VHDX2(L,K),0.)*QCON1(L,K)
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*QCON1(L,K)
     $        +MIN(W2(L,K),0.)*QCON1(L,K+1)
      END DO
      END DO
C
      GO TO 500
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN  (N-1) AND (N+1) AND ADVECTED FIELD AVERAGED
C **  BETWEEN AT (N-1) AND (N) IF ISTL 3 ONLY
C
  350 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=0,KC
      DO L=2,LA
      QCON2(L,K)=0.5*(QCON(L,K)+QCON1(L,K))
      END DO
      END DO
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MAX(UHDY2(L,K),0.)*QCON2(L-1,K)
     $         +MIN(UHDY2(L,K),0.)*QCON2(L,K)
      FVHU(L,K)=MAX(VHDX2(L,K),0.)*QCON2(LS,K)
     $         +MIN(VHDX2(L,K),0.)*QCON2(L,K)
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      FWU(L,K)=MAX(W2(L,K),0.)*QCON2(L,K)
     $        +MIN(W2(L,K),0.)*QCON2(L,K+1)
      END DO
      END DO
C
      GO TO 500
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT
C **  AVERAGED BETWEEN (N+1) AND (N-1) AND TRANSPORTED FIELD AT (N)
C
  400 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=0.5*UHDY2(L,K)*(QCON(L,K)+QCON(L-1,K))
      FVHU(L,K)=0.5*VHDX2(L,K)*(QCON(L,K)+QCON(LS,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF (VHDX2(LN,K).LT.0.) FVHU(LN,K)=VHDX2(LN,K)*QCON1(LN,K)
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(UHDY2(L+1,K).LT.0.) FUHU(L+1,K)=UHDY2(L+1,K)*QCON1(L+1,K)
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(UHDY2(L,K).GT.0.) FUHU(L,K)=UHDY2(L,K)*QCON1(L-1,K)
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS =LSC(L)
      IF (VHDX2(L,K).GT.0.) FVHU(L,K)=VHDX2(L,K)*QCON1(LS,K)
      END DO
C
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=0,KS
      DO L=2,LA
      FWU(L,K)=0.5*W2(L,K)*(QCON(L,K+1)+QCON(L,K))
      END DO
      END DO
C
C**********************************************************************C
C
C **  ADVECTION CALCULATION
C
  500 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.2)THEN
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      QCH(L,K)=QCON1(L,K)*H1P(L)
     $       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     $              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $             +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO
C
      IF (ISFCT(M).GE.1) THEN
       DO K=0,KC
       DO L=2,LA
       QCON2(L,K)=QCON1(L,K)
       END DO
       END DO
      END IF
C
      ELSE
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      QCH(L,K)=QCON1(L,K)*H2P(L)
     $       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     $              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $             +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO
C
      IF (ISFCT(M).GE.1) THEN
       DO K=0,KC
       DO L=2,LA
       QCON2(L,K)=QCON(L,K)
       END DO
       END DO
      END IF
C
      END IF
C
      IF (ISUD.EQ.1) THEN
      DO K=1,KS
      DO L=2,LA
      QCON1(L,K)=SCB(L)*QCON(L,K)+(1.-SCB(L))*QCON1(L,K)
      END DO
      END DO
      END IF
C
      DO K=1,KS
      DO L=2,LA
      QCON(L,K)=SCB(L)*QCH(L,K)*HPI(L)+(1.-SCB(L))*QCON(L,K)
      END DO
      END DO
C
C**********************************************************************C
C
C **  ASSIGN OPEN BOUNDARY CELL VALUES NOT TO INTERFER WITH FCT SCHEME
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO LL=1,NCBS
      LN=LNC(L)
      QCON(L,K)=QCON(LN,K)
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBW
      L=LCBW(LL)
      QCON(L,K)=QCON(L+1,K)
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBE
      L=LCBE(LL)
      QCON(L,K)=QCON(L-1,K)
      END DO
      END DO
C
      DO K=1,KS
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      QCON(L,K)=QCON(LS,K)
      END DO
      END DO
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION
C
      IF (ISADAC(M).EQ.0) RETURN
      IF (ISCDCA(M).EQ.1) RETURN
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      UUU(L,K)=U2(L,K)*(QCON(L,K)-QCON(L-1,K))*DXIU(L)
      VVV(L,K)=V2(L,K)*(QCON(L,K)-QCON(LS,K))*DYIV(L)
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      LS=LSC(L)
      WWW(L,K)=W2(L,K)*(QCON(L,K+1)-QCON(L,K))*DZIG(K)*HPI(L)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      LNW=LNWC(L)
      LSE=LSEC(L)
      AUHU=ABS(UHDY2(L,K))
      AVHV=ABS(VHDX2(L,K))
      UTERM=AUHU*(1.-DELTA*AUHU/(DXYU(L)*HU(L)))
     $          *(QCON(L,K)-QCON(L-1,K))
      VTERM=AVHV*(1.-DELTA*AVHV/(DXYV(L)*HV(L)))
     $          *(QCON(L,K)-QCON(LS,K))
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
C     UTERM=UTERM*(QCON(L,K)-QCON(L-1,K))
C     VTERM=VTERM*(QCON(L,K)-QCON(LS,K))
      UTERM=UTERM-0.25*DELTA*UHDY2(L,K)*
     $      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
     $      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
      VTERM=VTERM-0.25*DELTA*VHDX2(L,K)*
     $      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
     $      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
      UHU=UTERM/(QCON(L,K)+QCON(L-1,K)+BSMALL)
      VHV=VTERM/(QCON(L,K)+QCON(LS,K)+BSMALL)
C     AUHU=ABS(UHU)
C     AVHV=ABS(VHV)
C     UTERM=AUHU*(1.-DELTD2*AUHU/(DXYU(L)*HU(L)))
C     VTERM=AVHV*(1.-DELTD2*AVHV/(DXYV(L)*HV(L)))
C     AHTMP=0.5*UTERM*DXU(L)/(DYU(L)*HU(L))
C     AHU(L,K)=MAX(AHTMP,0.)
C     AHTMP=0.5*VTERM*DYV(L)/(DXV(L)*HV(L))
C     AHV(L,K)=MAX(AHTMP,0.)
      FUHU(L,K)=MAX(UHU,0.)*QCON(L-1,K)
     $         +MIN(UHU,0.)*QCON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*QCON(LS,K)
     $         +MIN(VHV,0.)*QCON(L,K)
      END DO
      END DO
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(W2(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)*HPI(L))*(QCON(L,K+1)-QCON(L,K))
      WTERM=WTERM-0.25*DELTA*W2(L,K)*
     $      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
     $      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
      WW=WTERM/(QCON(L,K+1)+QCON(L,K)+BSMALL)
      FWU(L,K)=MAX(WW,0.)*QCON(L,K)
     $        +MIN(WW,0.)*QCON(L,K+1)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DO K=1,KS
      FVHU(LN,K)=0.0
      END DO
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DO K=1,KS
      FUHU(L+1,K)=0.0
      END DO
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DO K=1,KS
      FUHU(L,K)=0.0
      END DO
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      DO K=1,KS
      FVHU(L,K)=0.0
      END DO
      END DO
C
C**********************************************************************C
C
C **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS
C
      IF (ISFCT(M).EQ.0) GO TO 600
C
C----------------------------------------------------------------------C
C
C **  DETERMINE MAX AND MIN CONCENTRATIONS
C
      DO K=0,KC
      DO L=2,LA
      QCONT(L,K)=MAX(QCON(L,K),QCON2(L,K))
      END DO
      END DO
C
      DO K=1,KS
      DO L=2,LA
      QCMAX(L,K)=MAX(QCONT(L,K-1),QCONT(L,K),QCONT(L,K+1))
      END DO
      END DO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB(L)*QCONT(L-1,K)
      CEMAX=SUB(L+1)*QCONT(L+1,K)
      CSMAX=SVB(L)*QCONT(LS,K)
      CNMAX=SVB(LN)*QCONT(LN,K)
      QCMAX(L,K)=MAX(QCMAX(L,K),CNMAX,CEMAX,CSMAX,CWMAX)
      END DO
      END DO
C
      DO L=2,LA
      DO K=1,KC
      QCONT(L,K)=MIN(QCON(L,K),QCON2(L,K))
      END DO
      END DO
C
      DO L=2,LA
      QCMIN(L,1)=MIN(QCONT(L,1),QCONT(L,2))
      QCMIN(L,KC)=MIN(QCONT(L,KS),QCONT(L,KC))
      END DO
C
      DO K=2,KS
      DO L=2,LA
      QCMIN(L,K)=MIN(QCONT(L,K-1),QCONT(L,K),QCONT(L,K+1))
      END DO
      END DO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMIN=SUB(L)*QCONT(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB(L+1)*QCONT(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB(L)*QCONT(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB(LN)*QCONT(LN,K)+1.E+6*(1.-SVB(LN))
      QCMIN(L,K)=MIN(QCMIN(L,K),CNMIN,CEMIN,CSMIN,CWMIN)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SEPARATE POSITIVE AND NEGATIVE FLUXES PUTTING NEGATIVE FLUXES
C **  INTO FUHV, FVHV, AND FWV
C
      DO K=1,KS
      DO L=2,LA
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE INFLUX AND OUTFLUX IN  CONCENTRATION UNITS AND LOAD
C **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      DU(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L,K)-FUHV(L+1,K)
     $                               +FVHU(L,K)-FVHV(LN,K))
     $                      +DZIC(K)*(FWU(L,K-1)-FWV(L,K)) )*HPI(L)
      DV(L,K)=DELT*SCB(L)*( DXYIP(L)*(FUHU(L+1,K)-FUHV(L,K)
     $                               +FVHU(LN,K)-FVHV(L,K))
     $                      +DZIC(K)*(FWU(L,K)-FWV(L,K-1)) )*HPI(L)
      END DO
      END DO
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DO K=1,KS
      DU(LN,K)=0.
      DV(LN,K)=0.
      END DO
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DO K=1,KS
      DU(L+1,K)=0.
      DV(L+1,K)=0.
      END DO
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DO K=1,KS
      DU(L-1,K)=0.
      DV(L-1,K)=0.
      END DO
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      DO K=1,KS
      DU(LS,K)=0.
      DV(LS,K)=0.
      END DO
      END DO
C
      DO NS=1,NQSIJ
       L=LQS(NS)
       NQSTMP=NQSERQ(NS)
       DO K=1,KC
        QQQTMP=ABS(QSS(K,NS)+QSERT(K,NQSTMP))
        IF(QQQTMP.GE.1.E-12) THEN
          DU(L,K)=0.
          DV(L,K)=0.
C         DU(L+1,K)=0.
C         DV(LNC(L),K)=0.
        END IF
       END DO
      END DO
C
      DO NCTL=1,NQCTL
       IU=IQCTLU(NCTL)
       JU=JQCTLU(NCTL)
       LU=LIJ(IU,JU)
       ID=IQCTLD(NCTL)
       JD=JQCTLD(NCTL)
       IF (ID.EQ.0.AND.JD.EQ.0) THEN
         LD=LC
        ELSE
         LD=LIJ(ID,JD)
       END IF
       DO K=1,KC
        QQQTMP=ABS(QCTLT(K,NCTL))
        IF(QQQTMP.GE.1.E-12) THEN
          DU(LU,K)=0.
          DV(LU,K)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
          DU(LD,K)=0.
          DV(LD,K)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
        END IF
       END DO
      END DO
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
        NQSTMP=NQWRSERQ(NWR)
        QQQTMP=ABS(QWR(NWR)+QWRSERT(NQSTMP))
        IF(QQQTMP.GE.1.E-12) THEN
          DU(LU,KU)=0.
          DV(LU,KU)=0.
C         DU(LU+1,K)=0.
C         DV(LNC(LU),K)=0.
          DU(LD,KD)=0.
          DV(LD,KD)=0.
C         DU(LD+1,K)=0.
C         DV(LNC(LD),K)=0.
        END IF
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DU AND DV
C
      DO K=1,KS
      DO L=2,LA
      IF (DU(L,K).GT.0.) DU(L,K)=(QCMAX(L,K)-QCON(L,K))/(DU(L,K)+BSMALL)
      IF (DV(L,K).GT.0.) DV(L,K)=(QCON(L,K)-QCMIN(L,K))/(DV(L,K)+BSMALL)
      END DO
      END DO
C
      DO K=1,KS
      DO L=2,LA
      DU(L,K)=MIN(DU(L,K),1.)
      DV(L,K)=MIN(DV(L,K),1.)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  LIMIT FLUXES
C
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)
     $         +MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)
     $         +MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      END DO
      END DO
C
      DO K=0,KS
      DO L=2,LA
      FWU(L,K)=MIN(DV(L,K),DU(L,K+1))*FWU(L,K)
     $        +MIN(DU(L,K),DV(L,K+1))*FWV(L,K)
      END DO
      END DO
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTION CALCULATION
C
  600 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      QCH(L,K)=QCON(L,K)*HP(L)
     $       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     $              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $             +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO
C
      DO K=0,KC
      DO L=2,LA
      QCON(L,K)=SCB(L)*QCH(L,K)*HPI(L)+(1.-SCB(L))*QCON(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  DIAGNOSE FCT SCHEME
C
      IF (ISFCT(M).EQ.2) THEN
      WRITE(6,6010)N
C
      DO K=1,KS
      DO L=2,LA
      CQCMAX=SCB(L)*(QCON(L,K)-QCMAX(L,K))
      IF (CQCMAX.GT.0.) THEN
       WRITE(6,6011)QCON(L,K),QCMAX(L,K),IL(L),JL(L),K
      END IF
      CQCMIN=SCB(L)*(QCMIN(L,K)-QCON(L,K))
      IF (CQCMIN.GT.0.) THEN
       WRITE(6,6012)QCMIN(L,K),QCON(L,K),IL(L),JL(L),K
      END IF
      END DO
      END DO
C
      END IF
C
 6010 FORMAT(1X,'FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT(1X,'QCON = ',E12.4,3X,'QCMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT(1X,'QCMIN = ',E12.4,3X,'QCON = ',E12.4,3X,'I,J,K=',(3I10))
C
C----------------------------------------------------------------------C
C
      RETURN
      END
