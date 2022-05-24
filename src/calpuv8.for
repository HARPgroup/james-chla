C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV8 (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
C ** SUBROUTINE CALPUV8 CALCULATES THE EXTERNAL SOLUTION FOR HP, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOWS WITH PROVISIONS FOR WETTING
C ** AND DRYING OF CELLS, USING THE SMOLARKIEWICZ AND MARGOLIN SCHEME
C ** (MONTHLY WEATHER REVIEW, VOL 121, 1847-1859, JUNE 1993)
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION QSUMTMP(LCM)
C
C**********************************************************************C
C
       DELT=DT2
       DELTD2=DT
      IF (ISTL.EQ.2) THEN
       DELT=DT
       DELTD2=0.5*DT
      END IF
      DELTI=1./DELT
C
C**********************************************************************C
C
C **  ADVANCE THE CONTINUITY EQUATION 1/2 OF A TIME STEP FOR VOLUME
C **  SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
       DO L=2,LA
       IF (QSUME(L).GE.0.) THEN
         HPPTMP=H1P(L)+DELTD2*SPB(L)*DXYIP(L)*QSUME(L)
         QSUMTMP(L)=QSUME(L)
        ELSE
         IF (H1P(L).LE.HDRY) THEN
           HPPTMP=H1P(L)
           QSUMTMP(L)=0.
          ELSE
           QSUMTMP(L)=-(H1P(L)-HDRY)*DXYP(L)*DELTI
           QSUMTMP(L)=MAX(QSUMTMP(L),QSUME(L))
           HPPTMP=H1P(L)+DELTD2*SPB(L)*DXYIP(L)*QSUMTMP(L)
         END IF
       END IF
       HPTMP(L)=HPPTMP
       END DO
C
C**********************************************************************C
C
C **  ADVANCE CONTIUNITY EQUATION FOR ADVECTIVE FLOWS
C
C----------------------------------------------------------------------C
C
C **  ADVANCE THE CONTINUITY EQUATION A FULL TIME STEP USING UPWIND
C **  DIFFERENCING
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.EQ.9.OR.ISDRY.EQ.10) THEN
C
      DO L=1,LC
      FUHU(L,1)=0.
      FVHU(L,1)=0.
      FVHU(L,1)=0.
      FVHV(L,1)=0.
      UHDY2E(L)=0.
      VHDX2E(L)=0.
      END DO
C
      IF(ISTL.EQ.2) THEN
        DO L=2,LA
        UHDY2E(L)=0.5*(UHDY1E(L)*H1UI(L)+UHDYE(L)*HUI(L))
        VHDX2E(L)=0.5*(VHDX1E(L)*H1VI(L)+VHDXE(L)*HVI(L))
        END DO
       ELSE
        DO L=2,LA
        UHDY2E(L)=UHDYE(L)*HUI(L)
        VHDX2E(L)=VHDXE(L)*HVI(L)
        END DO
      END IF
C
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,1)=MAX(UHDY2E(L),0.)*HPTMP(L-1)
     $         +MIN(UHDY2E(L),0.)*HPTMP(L  )
      FVHU(L,1)=MAX(VHDX2E(L),0.)*HPTMP(LS )
     $         +MIN(VHDX2E(L),0.)*HPTMP(L  )
      END DO
C
      DO L=2,LA
      LN=LNC(L)
      HPTMP(L)=HPTMP(L)
     $        +DELT*SPB(L)*DXYIP(L)*( FUHU(L,1)-FUHU(L+1,1)
     $                               +FVHU(L,1)-FVHU(LN ,1) )
      END DO
C
      END IF
C
C----------------------------------------------------------------------C
C
C **  ADVANCE THE CONTINUITY EQUATION A FULL TIME STEP USING ACTUAL
C **  VOLUME FLUXES
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.EQ.11.OR.ISDRY.EQ.12) THEN
C
      DO L=1,LC
      UHDY2E(L)=0.
      VHDX2E(L)=0.
      END DO
C
      IF(ISTL.EQ.2) THEN
        DO L=2,LA
        UHDY2E(L)=0.5*(UHDY1E(L)+UHDYE(L))
        VHDX2E(L)=0.5*(VHDX1E(L)+VHDXE(L))
        END DO
       ELSE
        DO L=2,LA
        UHDY2E(L)=UHDYE(L)
        VHDX2E(L)=VHDXE(L)
        END DO
      END IF
C
      DO L=2,LA
      LN=LNC(L)
      HPTMP(L)=HPTMP(L)
     $        +DELT*SPB(L)*DXYIP(L)*( UHDY2E(L)-UHDY2E(L+1)
     $                               +VHDX2E(L)-VHDX2E(LN ) )
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  ADVANCE THE CONTINUITY EQUATION 1/2 OF A TIME STEP FOR VOLUME
C **  SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      HPPTMP=HPTMP(L)+DELTD2*SPB(L)*DXYIP(L)*QSUMTMP(L)
      IF(HPPTMP.LE.0.0) THEN
        HPPTMP=HPTMP(L)
        QSUMTMP(L)=0.5*QSUMTMP(L)
      END IF
      HPTMP(L)=HPPTMP
      END DO
C
C**********************************************************************C
C
C **  SET ELEVATION BOUNDARY CONDITIONS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      FP(L)=0.
      END DO
C
      TN=DT*FLOAT(N)+TCON*TBEGIN
C
      DO M=1,MTIDE
      TM=MOD(TN,TCP(M))
      TM=PI2*TM/TCP(M)
      CCCOS(M)=COS(TM)
      SSSIN(M)=SIN(TM)
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBW
      L=LPBW(LL)
      FP(L)=PSERT(NPSERW(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
      END DO
      IF(ISPBW(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
       CET=(1.-TMP)/(1.+TMP)
       CE(L)=CET
       CCE(L)=CET
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CER(LR)=CET
        CCER(LR)=CET
       ELSE
        LB=LLBC(L)
        CEB(LB)=CET
        CCEB(LB)=CET
       END IF
       FP(L)=(4.*FP(L)-2.*SQRT(G/HMU(L+1))*FUHDYE(L+1)*DYIU(L+1))/
     $       (1.+TMP)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)
      FP(L)=PSERT(NPSERE(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
      END DO
      IF(ISPBE(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
       CWT=(1.-TMP)/(1.+TMP)
       CW(L)=CWT
       CCW(L)=CWT
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CWR(LR)=CWT
        CCWR(LR)=CWT
       ELSE
        LB=LLBC(L)
        CWB(LB)=CWT
        CCWB(LB)=CWT
       END IF
       FP(L)=(4.*FP(L)+2.*SQRT(G/HMU(L))*FUHDYE(L)*DYIU(L))/(1.+TMP)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      FP(L)=PSERT(NPSERS(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
      END DO
      IF(ISPBS(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
       CNT=(1.-TMP)/(1.+TMP)
       CN(L)=CNT
       CCN(L)=CNT
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CNR(LR)=CNT
        CCNR(LR)=CNT
       ELSE
        LB=LLBC(L)
        CNB(LB)=CNT
        CCNB(LB)=CNT
       END IF
       FP(L)=(4.*FP(L)-2.*SQRT(G/HMV(LN))*FVHDXE(LN)*DXIV(LN))/(1.+TMP)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
      LS=LSC(L)
      FP(L)=PSERT(NPSERN(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP(L)=FP(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
      END DO
      IF(ISPBN(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
       CST=(1.-TMP)/(1.+TMP)
       CS(L)=CST
       CCS(L)=CST
       IF (ISRED(L).EQ.1) THEN
        LR=LLRC(L)
        CSR(LR)=CST
        CCSR(LR)=CST
       ELSE
        LB=LLBC(L)
        CSB(LB)=CST
        CCSB(LB)=CST
       END IF
       FP(L)=(4.*FP(L)+2.*SQRT(G/HMV(L))*FVHDXE(L)*DXIV(L))/(1.+TMP)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      HPTMP(L)=SPB(L)*HPTMP(L)+(1.-SPB(L))*(GI*FP(L)-BELV(L))
      END DO
C
C**********************************************************************C
C
C **  ADJUST VERTICAL SOURCE-SINK DISTRIBUTION
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
       DIFQVOL=QSUME(L)-QSUMTMP(L)
       DO K=1,KC
       QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
       END DO
       QSUME(L)=QSUMTMP(L)
      END DO
C
C**********************************************************************C
C
C **  SET TEMPORARY NEW TIME LEVEL DEPTHS AT VELOCITY POINTS
C
C----------------------------------------------------------------------C

      DO L=2,LA
      LS=LSC(L)
      HUTMP(L)=0.5*(HPTMP(L)+HPTMP(L-1))
      HVTMP(L)=0.5*(HPTMP(L)+HPTMP(LS ))
      END DO
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N)
C
      CALL CALEBI
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL PRESSURE GRADIENTS
C **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
C **  SNLPX=SNLPX*GID2*DYU & SNLPY=SNLPY*GID2*DXV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)
      FPGXE(L)=ROLD*FPGXE(L)+RNEW*(
     $           -SBX(L)*HU(L)*((BI2(L)+BI2(L-1))*(HP(L)-HP(L-1))
     $           +2.*HU(L)*(BI1(L)-BI1(L-1))
     $           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1))) )
      FPGYE(L)=ROLD*FPGYE(L)+RNEW*(
     $           -SBY(L)*HV(L)*((BI2(L)+BI2(LS))*(HP(L)-HP(LS))
     $           +2.*HV(L)*(BI1(L)-BI1(LS))
     $           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS))) )
      END DO
C
C**********************************************************************C
C
C **  CALCULATE 1/2 OF THE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION
C **  TERMS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)
      FUHDYE(L)=UHDY1E(L)-DELTD2*SUB(L)*HRUO(L)*H1U(L)*(P1(L)-P1(L-1))
     $         +SUB(L)*DELTD2*DXIU(L)*(DXYU(L)*(TSX1(L)-TBX1(L))
     $         +FCAXE(L)+FPGXE(L))
      FUHDYE(L)=DYIU(L)*H1UI(L)*FUHDYE(L)
      FVHDXE(L)=VHDX1E(L)-DELTD2*SVB(L)*HRVO(L)*H1V(L)*(P1(L)-P1(LS))
     $         +SVB(L)*DELTD2*DYIV(L)*(DXYV(L)*(TSY1(L)-TBY1(L))
     $         -FCAYE(L)+FPGYE(L))
      FVHDXE(L)=DXIV(L)*H1VI(L)*FVHDXE(L)
      END DO
C
C**********************************************************************C
C
C **  ADVANCE THE ADVECTIVE TERMS BY A FULL TIME STEP
C **  TERMS
C
C----------------------------------------------------------------------C
C
      DO L=1,LC
      FUHU(L,1)=0.
      FVHU(L,1)=0.
      FVHU(L,1)=0.
      FVHV(L,1)=0.
      UHDY2E(L)=0.
      VHDX2E(L)=0.
      FXE(L)=0.
      FYE(L)=0.
      END DO
C
      IF(ISTL.EQ.3) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  TWO TIME LEVEL STEP
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      UHDY2E(L)=UHDY1E(L)+UHDYE(L)
      VHDX2E(L)=VHDX1E(L)+VHDXE(L)
      END DO
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      UHC=0.25*(UHDY2E(L)+UHDY2E(LS ))
      UHB=0.25*(UHDY2E(L)+UHDY2E(L+1))
      VHC=0.25*(VHDX2E(L)+VHDX2E(L-1))
      VHB=0.25*(VHDX2E(L)+VHDX2E(LN ))
      FUHU(L,1)=MAX(UHB,0.)*FUHDYE(L  )
     $         +MIN(UHB,0.)*FUHDYE(L+1)
      FVHU(L,1)=MAX(VHC,0.)*FUHDYE(LS )
     $         +MIN(VHC,0.)*FUHDYE(L  )
      FUHV(L,1)=MAX(UHC,0.)*FVHDXE(L-1)
     $         +MIN(UHC,0.)*FVHDXE(L  )
      FVHV(L,1)=MAX(VHB,0.)*FVHDXE(L  )
     $         +MIN(VHB,0.)*FVHDXE(LN )
      END DO
C
      GO TO 400
C
C----------------------------------------------------------------------C
C
C **  THREE TIME LEVEL (LEAP-FROG) STEP
C **  WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      UHC=0.5*(UHDYE(L)+UHDYE(LS ))
      UHB=0.5*(UHDYE(L)+UHDYE(L+1))
      VHC=0.5*(VHDXE(L)+VHDXE(L-1))
      VHB=0.5*(VHDXE(L)+VHDXE(LN ))
      FUHU(L,1)=MAX(UHB,0.)*FUHDYE(L  )
     $         +MIN(UHB,0.)*FUHDYE(L+1)
      FVHU(L,1)=MAX(VHC,0.)*FUHDYE(LS )
     $         +MIN(VHC,0.)*FUHDYE(L  )
      FUHV(L,1)=MAX(UHC,0.)*FVHDXE(L-1)
     $         +MIN(UHC,0.)*FVHDXE(L  )
      FVHV(L,1)=MAX(VHB,0.)*FVHDXE(L  )
     $         +MIN(VHB,0.)*FVHDXE(LN )
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE  NET ADVECTIVE FLUXES
C
C----------------------------------------------------------------------C
C
  400 CONTINUE
C
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      FXE(L)=FUHU(L  ,1)-FUHU(L-1,1)+FVHU(LN,1)-FVHU(L ,1)
      FYE(L)=FUHV(L+1,1)-FUHV(L  ,1)+FVHV(L ,1)-FVHV(LS,1)
      END DO
C
C----------------------------------------------------------------------C
C
C **  COMPLETE ADVECTIVE FLUX ADVANCE
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      FUHDYE(L)=DYU(L)*H1U(L)*FUHDYE(L)
      FVHDXE(L)=DXV(L)*H1V(L)*FVHDXE(L)
      END DO
C
      DO L=2,LA
      FUHDYE(L)=FUHDYE(L)-SNLT*SUB(L)*DELT*DXIU(L)*FXE(L)
      FVHDXE(L)=FVHDXE(L)-SNLT*SVB(L)*DELT*DYIV(L)*FYE(L)
      END DO
C
C**********************************************************************C
C
C **  ADVANCE EXTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.3) THEN
        DO L=2,LA
        UHDY2E(L)=UHDY1E(L)
        VHDX2E(L)=VHDX1E(L)
        UHDY1E(L)=UHDYE(L)
        VHDX1E(L)=VHDXE(L)
        U1V(L)=UV(L)
        V1U(L)=VU(L)
        P1(L)=P(L)
        H1U(L)=HU(L)
        H1V(L)=HV(L)
        H1UI(L)=HUI(L)
        H1VI(L)=HVI(L)
        H2P(L)=H1P(L)
        H1P(L)=HP(L)
        END DO
       ELSE
        DO L=2,LA
        UHDY2E(L)=0.5*UHDY2E(L)
        VHDX2E(L)=0.5*VHDX2E(L)
        END DO
      END IF
C
C**********************************************************************C
C
C **  UPDATE BOUNDARY CONDITION SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      SUB1(L)=SUB(L)
      SVB1(L)=SVB(L)
      SUB(L)=SUBO(L)
      SVB(L)=SVBO(L)
      SBX(L)=SBXO(L)
      SBY(L)=SBYO(L)
      END DO
C
C**********************************************************************C
C
C **  ADVANCE MOMENTUM EQUATIONS REMAINING 1/2 TIME STEP
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HPTMP(L)+BELV(L))
      END DO
C
  500 CONTINUE
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      UHDYE(L)=SUB(L)*( FUHDYE(L)-DELTD2*HRUO(L)*HUTMP(L)*(P(L)-P(L-1))
     $         +DELTD2*SUB1(L)*DXIU(L)*(DXYU(L)*(TSX1(L)-TBX1(L))
     $         +FCAXE(L)+FPGXE(L)) )
      VHDXE(L)=SVB(L)*( FVHDXE(L)-DELTD2*HRVO(L)*HVTMP(L)*(P(L)-P(LS ))
     $         +DELTD2*SVB1(L)*DYIV(L)*(DXYV(L)*(TSY1(L)-TBY1(L))
     $         -FCAYE(L)+FPGYE(L)) )
      END DO
C
C**********************************************************************C
C
C **  COMPLETE CALCULATION IF CONTINUITY ADJUSTMENT IS NOT ENFORCED
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.EQ.9.OR.ISDRY.EQ.11) THEN
        DO L=2,LA
        HP(L)=HPTMP(L)
        END DO
        GO TO 600
       END IF
C
C**********************************************************************C
C
C **  CALCULATE HP AT NEW TIME LEVEL ENFORCING CONTINUITY
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.3) THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HPTMP(L)
        END DO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HPTMP(L)
        END DO
      END IF
C
C**********************************************************************C
C
C **  CHECK FOR UNDETECTED DRYING
C
C----------------------------------------------------------------------C
C
      OPEN(1,FILE='drywet.log',ACCESS='APPEND',STATUS='UNKNOWN')
C
      ICORDRY=0
      DO L=2,LA
      LN=LNC(L)
      IF (HP(L).LE.0.0) THEN
        ICORDRY=1
        WRITE(1,6945) N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
        WRITE(6,6945) N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
      END IF
      END DO
C
      CLOSE(1)
C
      IF(ICORDRY.EQ.1) THEN
        IF (NCORDRY.GT.IDRYCK) THEN
C         WRITE(6,6961)NCORDRY
CTMP          WRITE(8,6961)NCORDRY
          STOP
        END IF
        NCORDRY=NCORDRY+1
        GO TO 500
      END IF
C
C     WRITE(6,6960)NCORDRY
CTMP      WRITE(8,6960)NCORDRY
C
 6960 FORMAT(' NCORDRY =', I5)
 6961 FORMAT(' UNSTABLE, NCORDRY =', I5)
C
C**********************************************************************C
C
C **  FINAL UPDATE OF SURFACE ELEVATION AND DEPTHS
C
C----------------------------------------------------------------------C
C
  600 CONTINUE
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      END DO
C
      DO L=2,LA
      LS=LSC(L)
      HU(L)=0.5*(HP(L)+HP(L-1))
      HV(L)=0.5*(HP(L)+HP(LS))
      END DO
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      END DO
C
C**********************************************************************C
C
C **  CHECK FOR NEGATIVE DEPTHS
C
C----------------------------------------------------------------------C
C
      IF (ISNEGH.EQ.1) THEN
C
      DO L=2,LA
      IF (HP(L).LT.0.) THEN
      LN=LNC(L)
      WRITE (6,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (6,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (6,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (6,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (6,6064)IL(L),JL(L),HV(LN),H1V(LN)
      WRITE (8,6060)IL(L),JL(L),HP(L),H1P(L),H2P(L)
      WRITE (8,6061)IL(L),JL(L),HU(L),H1U(L)
      WRITE (8,6062)IL(L),JL(L),HU(L+1),H1U(L+1)
      WRITE (8,6063)IL(L),JL(L),HV(L),H1V(L)
      WRITE (8,6064)IL(L),JL(L),HV(LN),H1V(LN)
      END IF
      END DO
C
      END IF
C
 6060 FORMAT(1X,'NEG DEPTH AT I,J =',2I4,'  HP,H1P,H2P =',3(2X,E12.4))
 6061 FORMAT(1X,'NEG DEPTH AT I,J =',2I4,'  HUW,H1UW =',2(2X,E12.4))
 6062 FORMAT(1X,'NEG DEPTH AT I,J =',2I4,'  HUE,H1UE =',2(2X,E12.4))
 6063 FORMAT(1X,'NEG DEPTH AT I,J =',2I4,'  HVS,H1VS =',2(2X,E12.4))
 6064 FORMAT(1X,'NEG DEPTH AT I,J =',2I4,'  HVN,H1VN =',2(2X,E12.4))
C
C**********************************************************************C
C
 6910 FORMAT(1X,'DRYING AT N,I,J =',I10,2I6,'  HP,H1P,H2P ='
     $         ,3(2X,E12.4))
 6911 FORMAT(1X,'DRY W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6912 FORMAT(1X,'DRY E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6913 FORMAT(1X,'DRY S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6914 FORMAT(1X,'DRY N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6920 FORMAT(1X,'WETTING AT N,I,J =',I10,2I6,' HP,H1P,H2P ='
     $         ,3(2X,E12.4))
 6921 FORMAT(1X,'WET S FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
 6922 FORMAT(1X,'WET W FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6923 FORMAT(1X,'WET E FACE N,I,J =',I10,2I6,' HU,H,H1 =',3(2X,E12.4))
 6924 FORMAT(1X,'WET N FACE N,I,J =',I10,2I6,' HV,H,H1 =',3(2X,E12.4))
C
 6930 FORMAT(1X,'WET BY VOL  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     $         ,3(2X,E12.4))
 6940 FORMAT(1X,'RESOLVE,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     $         ,3(2X,E12.4))
 6941 FORMAT(1X,'RESOLVE,  N,I,J =',I10,2I6,' HUE,HP,H1P ='
     $         ,3(2X,E12.4))
 6942 FORMAT(1X,'RESOLVE,  N,I,J =',I10,2I6,' HUW,HP,H1P ='
     $         ,3(2X,E12.4))
 6943 FORMAT(1X,'RESOLVE,  N,I,J =',I10,2I6,' HVS,HP,H1P ='
     $         ,3(2X,E12.4))
 6944 FORMAT(1X,'RESOLVE,  N,I,J =',I10,2I6,' HVN,HP,H1P ='
     $         ,3(2X,E12.4))
 6945 FORMAT(1X,'RESOLVE NEG,  N,I,J =',I10,2I6,' HP,H1P,H2P ='
     $         ,3(2X,E12.4))
 6950 FORMAT(1X,'RESOLVE, NEG DEP N,I,J =',I10,2I6,' HP,H1P,H2P ='
     $         ,3(2X,E12.4))
C
C**********************************************************************C
C
C **  CALCULATE THE EXTERNAL DIVERGENCE
C
C----------------------------------------------------------------------C
C
      IF(ISDIVEX.EQ.1) THEN
C
      DIVEXMX=0.
      DIVEXMN=1000000.
C
C----------------------------------------------------------------------C
C
      IF(ISTL.EQ.3) THEN
C
      DO L=2,LA
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI
     $     +0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $     +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))-QSUME(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      END IF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      END IF
      END DO
C
      ELSE
C
      DO L=2,LA
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI
     $     +0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      END IF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      END IF
      END DO
C
      END IF
C
      IMAX=IL(LMAX)
      JMAX=JL(LMAX)
      IMIN=IL(LMIN)
      JMIN=JL(LMIN)
C
      WRITE(6,6628)DIVEXMX,IMAX,JMAX
      WRITE(6,6629)DIVEXMN,IMIN,JMIN
C
C----------------------------------------------------------------------C
C
      END IF
C
C----------------------------------------------------------------------C
C
  566 FORMAT(1X,'I=',I5,3X,'J=',I5,3X,'HP=',F12.4)
 6628 FORMAT(1X,'DIVEXMX=',E12.4,5X,2I10)
 6629 FORMAT(1X,'DIVEXMN=',E12.4,5X,2I10)
C
C**********************************************************************C
C
      RETURN
      END
