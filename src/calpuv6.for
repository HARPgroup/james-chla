C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV6 (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
C ** SUBROUTINE CALPUV6 CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR FREE SURFACE FLOWS WITH PROVISIONS FOR WETTING
C ** AND DRYING OF CELLS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION QSUMTMP(LCM)
C
C**********************************************************************C
C
      IF (ISDSOLV.EQ.1) THEN
      OPEN(1,FILE='eqcoef.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='eqterm.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='fp.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      IF (N.EQ.1) THEN
        OPEN(1,FILE='eqcoef1.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      END IF
      END IF
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
C **  SET SWITCHES FOR DRYING AND WETTING
C
C----------------------------------------------------------------------C
C
      ITERHP=0
      NCORDRY=0
      ICORDRY=0
      DO L=1,LC
      ISCDRY(L)=0
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
       HUTMP(L)=MAX(HU(L),0.)
       HVTMP(L)=MAX(HV(L),0.)
       TVAR3N(L)=MAX(HP(L),0.)
      END DO
C
      DO L=2,LA
      LS=LSC(L)
      FPGXE(L)=ROLD*FPGXE(L)+RNEW*(
     $  -SBX(L)*HUTMP(L)*((BI2(L)+BI2(L-1))*(TVAR3N(L)-TVAR3N(L-1))
     $           +2.*HUTMP(L)*(BI1(L)-BI1(L-1))
     $           +(BE(L)+BE(L-1))*(BELV(L)-BELV(L-1))) )
      FPGYE(L)=ROLD*FPGYE(L)+RNEW*(
     $  -SBY(L)*HVTMP(L)*((BI2(L)+BI2(LS))*(TVAR3N(L)-TVAR3N(LS))
     $           +2.*HVTMP(L)*(BI1(L)-BI1(LS))
     $           +(BE(L)+BE(LS))*(BELV(L)-BELV(LS))) )
      END DO
C
C**********************************************************************C
C
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
       HUTMP(L)=MAX(H1U(L),0.)
       HVTMP(L)=MAX(H1V(L),0.)
       TVAR3N(L)=G*BELV(L)
      END DO
C
      DO L=2,LA
       TVAR3N(L)=MAX(TVAR3N(L),P1(L))
      END DO
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDY1E(L)
     $      -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(TVAR3N(L)-TVAR3N(L-1))
     $         +SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX1(L)-RITB1*TBX1(L))
     $         +FCAXE(L)+FPGXE(L)-SNLT*FXE(L))
      FVHDXE(L)=VHDX1E(L)
     $      -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(TVAR3N(L)-TVAR3N(LS ))
     $         +SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY1(L)-RITB1*TBY1(L))
     $         -FCAYE(L)+FPGYE(L)-SNLT*FYE(L))
      END DO
C
      DO L=2,LA
      RCX(L)=1.
      RCY(L)=1.
      END DO
C
      RCX(1)=0.
      RCY(1)=0.
      RCX(LC)=0.
      RCY(LC)=0.
C
C**********************************************************************C
C
C **  PROJECT SOLUTION
C
      IF(ISTL.EQ.2) THEN
        DO L=2,LA
         HUTMP(L)=MAX(HU(L),0.)
         HVTMP(L)=MAX(HV(L),0.)
         TVAR3E(L)=P(L)
         TVAR3E(L)=MAX(TVAR3E(L),0.)
        END DO
       ELSE
        DO L=2,LA
         HUTMP(L)=2.*HU(L)-H1U(L)
         HVTMP(L)=2.*HV(L)-H1V(L)
         HUTMP(L)=MAX(HUTMP(L),0.)
         HVTMP(L)=MAX(HVTMP(L),0.)
         TVAR3E(L)=2*P(L)-P1(L)
         TVAR3E(L)=MAX(TVAR3E(L),0.)
        END DO
      END IF
C
      DO L=2,LA
      LS=LSC(L)
       TVAR3W(L)=SUB(L)*( FUHDYE(L)
     $    -DELTD2*HRUO(L)*RCX(L)*HUTMP(L)*(TVAR3E(L)-TVAR3E(L-1)) )
       TVAR3S(L)=SVB(L)*( FVHDXE(L)
     $    -DELTD2*HRVO(L)*RCY(L)*HVTMP(L)*(TVAR3E(L)-TVAR3E(LS )) )
      END DO
C
C **  LOOP U FACES
C
      DO L=2,LA
       SUB(L)=SUBO(L)
       GBWTMP=G*BELV(L-1)
       PWTMP=TVAR3E(L-1)-G*HWET
       GBETMP=G*BELV(L)
       PETMP=TVAR3E(L)-G*HWET
       IF(PWTMP.GT.GBWTMP) THEN
C**      WEST IS WET
         IF(PETMP.GT.GBETMP) THEN
C**        EAST IS WET
           SUB(L)=SUBO(L)
          ELSE
C**        EAST IS DRY
           IF(TVAR3W(L).LT.0.) SUB(L)=0.
         END IF
        ELSE
C**      WEST IS DRY
         IF(PETMP.GT.GBETMP) THEN
C**        EAST IS WET
           IF(TVAR3W(L).GT.0.) SUB(L)=0.
          ELSE
C**        EAST IS DRY
           SUB(L)=0.
         END IF
       END IF
      END DO
C
C **  LOOP V FACES
C
      DO L=2,LA
       LS=LSC(L)
       SVB(L)=SVBO(L)
       GBSTMP=G*BELV(LS)
       PSTMP=TVAR3E(LS)-G*HWET
       GBNTMP=G*BELV(L)
       PNTMP=TVAR3E(L)-G*HWET
       IF(PSTMP.GT.GBSTMP) THEN
C**      SOUTH IS WET
         IF(PNTMP.GT.GBNTMP) THEN
C**        NORTH IS WET
           SVB(L)=SVBO(L)
          ELSE
C**        NORTH IS DRY
           IF(TVAR3S(L).LT.0.) SVB(L)=0.
         END IF
        ELSE
C**      SOUTH IS DRY
         IF(PNTMP.GT.GBNTMP) THEN
C**        NORTH IS WET
           IF(TVAR3W(L).GT.0.) SVB(L)=0.
          ELSE
C**        NORTH IS DRY
           SVB(L)=0.
         END IF
       END IF
      END DO
C
C**********************************************************************C
C
C **  RESET BOUNDARY CONDITIONS SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=1,LC
      FP(L)=0.
      FP1(L)=0.
      END DO
C
C**********************************************************************C
C
C **  SET OPEN BOUNDARY SURFACE ELEVATIONS
C
C----------------------------------------------------------------------C
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
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      FP1(L)=PSERT(NPSERW(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
      END DO
      CET=-0.5*DELTD2*G*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/(1.-TMP)
       CE(L)=CET
       FP1(L)=CET*(4.*FP1(L)
     $ -2.*SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/(1.-TMP)
      ELSE
       FP1(L+1)=-CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      FP1(L)=PSERT(NPSERE(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
      END DO
      CWT=-0.5*DELTD2*G*HRUO(L  )*RCX(L  )*HUTMP(L  )
      IF(ISPBE(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/(1.-TMP)
       CW(L)=CWT
       FP1(L)=CWT*(4.*FP1(L)
     $ +2.*SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/(1.-TMP)
      ELSE
       FP1(L-1)=-CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LN=LNC(L)
      FP1(L)=PSERT(NPSERS(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
      END DO
      CNT=-0.5*DELTD2*G*HRVO(LN )*RCY(LN )*HVTMP(LN )
      IF(ISPBS(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
       CC(L)=CNT*(1.+TMP)/(1.-TMP)
       CN(L)=CNT
       FP1(L)=CNT*(4.*FP1(L)
     $ -2.*SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/(1.-TMP)
      ELSE
       FP1(LN)=-CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
       CC(L)=DELTI*DXYP(L)
       CS(L)=0.
       CW(L)=0.
       CE(L)=0.
       CN(L)=0.
      LS=LSC(L)
      FP1(L)=PSERT(NPSERN(LL))
      DO M=1,MTIDE
      TC=CCCOS(M)
      TS=SSSIN(M)
      FP1(L)=FP1(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
      END DO
      CST=-0.5*DELTD2*G*HRVO(L  )*RCY(L  )*HVTMP(L  )
      IF(ISPBN(LL).EQ.1)THEN
       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
       CC(L)=CST*(1.+TMP)/(1.-TMP)
       CS(L)=CST
       FP1(L)=CST*(4.*FP1(L)
     $ +2.*SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/(1.-TMP)
      ELSE
       FP1(LS)=-CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C
C**********************************************************************C
C
C **  ADJUST VOLUME SOURCE AND SINKS
C
C----------------------------------------------------------------------C
C
      IF (ISGWIE.EQ.0) THEN
C
      DO L=2,LA
      IF (QSUME(L).LE.0.) THEN
        IF (H1P(L).LE.HDRY) THEN
          QSUMTMP(L)=0.
         ELSE
          QSUMTMP(L)=-(H1P(L)-HDRY)*DXYP(L)*DELTI
          QSUMTMP(L)=MAX(QSUMTMP(L),QSUME(L))
        END IF
       ELSE
        QSUMTMP(L)=QSUME(L)
      END IF
      END DO
C
      DO L=2,LA
       DIFQVOL=QSUME(L)-QSUMTMP(L)
       DO K=1,KC
       QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
       END DO
       QSUME(L)=QSUMTMP(L)
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  ADJUST SOURCES AND SINKS ESTIMATING SURFACE AND GROUNDWATER
C **  AVAILABLE FOR EVAPOTRANSPIRATON AND INFILTRATION
C
C----------------------------------------------------------------------C
C
      IF (ISGWIE.GE.1) THEN
C
      DO L=2,LA
      RIFTR(L)=0.
      EVAPSW(L)=0.
      EVAPGW(L)=0.
      IF (H1P(L).GT.HDRY) THEN
C       APPLY MAXIMUM ET
        IF(EVAPCVT.LT.0.) THEN
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     $              (1.+0.00412*TEM(L,KC))))
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)
        END IF
        EVAPSW(L)=EVAPT(L)*DXYP(L)
        RIFTR(L)=0.
C       CALCULATE DEPTH OF ACTIVE GROUNDWATER ELEV BELOW SURFACE
        DTAGW=BELV(L)-AGWELV(L)
        IF (DTAGW.GT.0.0) THEN
C         INFLITRATION CAN OCCUR, CALCULATE LIMITING RATE TO BRING
C         GW ELEV TO SOIL SURFACE
          RIFTRL=RNPOR*DTAGW*DELTI
C         SET RIFTRL TO MIN OF LIMITING RATE OR ACTUAL RATE
          RIFTRL=MIN(RIFTRM,RIFTRL)
C         ESTIMATE RATE BASED ON AVAILABLE SURFACE WATER
          RAVAIL=(H1P(L)-HDRY)*DELTI-EVAPT(L)
C         SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
          RIFTRL=MIN(RAVAIL,RIFTRL)
C         CONVERT TO VOLUME FLOW UNITS
          RIFTR(L)=RIFTRL*DXYP(L)
        END IF
C       ADJUST VOLUME OUTFLOWS OF WET CELLS
        IF (QSUME(L).LT.0.0) THEN
          QSUMIET=RIFTR(L)+EVAPSW(L)
          QEAVAIL=DXYP(L)*(H1P(L)-HDRY)*DELTI-QSUMIET
          QEAVAIL=MAX(QEAVAIL,0.0)
          QEAVAIL=-QEAVAIL
          QSUMTMP(L)=MAX(QSUME(L),QEAVAIL)
         ELSE
          QSUMTMP(L)=QSUME(L)
        END IF
       ELSE
        RIFTR(L)=0.
        EVAPSW(L)=0.
        QSUMTMP(L)=MAX(QSUME(L),0.0)
      END IF
      END DO
C
      DO L=2,LA
      DIFQVOL=QSUME(L)-QSUMTMP(L)
      DO K=1,KC
      QSUM(L,K)=QSUM(L,K)-DIFQVOL*DZC(K)
      END DO
      QSUME(L)=QSUMTMP(L)
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR
C **  NON BOUNDARY POINTS
C **  HRU=HMU*DYU/DXU & HRV=HMV*DXV/DYV
C **  DXYIP=1/(DXP*DYP)
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*P1(L)
     $      -0.5*G*(UHDY1E(L+1)-UHDY1E(L)
     $             +VHDX1E(LN )-VHDX1E(L)) )
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
C       DELP=P(L)-P1(L)
        P1(L)=P(L)
C       P(L)=P(L)+DELP
        H1U(L)=HU(L)
        H1V(L)=HV(L)
        H1UI(L)=HUI(L)
        H1VI(L)=HVI(L)
        H2P(L)=H1P(L)
        H1P(L)=HP(L)
        AGWELV2(L)=AGWELV1(L)
        AGWELV1(L)=AGWELV(L)
        END DO
      END IF
C
cc     DO L=2,LA
cc      PAM(L)=P(L)
cc      END DO
C
      IF (MDCHHQ.GE.1) THEN
        DO NMD=1,MDCHH
        QCHANU(NMD)=0.
        QCHANV(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANVN(NMD)=0.
        END DO
      END IF
C
C**********************************************************************C
C
C **  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
C **  HOST-GUEST CHANNAL INTERACTION FOR NON BOUNDARY POINTS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LN=LNC(L)
      FP(L)=FP1(L)-0.5*G*SPB(L)*
     $      ( SUB(L+1)*FUHDYE(L+1)-SUB(L)*FUHDYE(L)
     $       +SVB(LN )*FVHDXE(LN )-SVB(L)*FVHDXE(L)
     $       -2.0*QSUME(L) )
cc      P(L)=0.
      END DO
C
      IF (ISGWIE.GE.1) THEN
        DO L=2,LA
        FP(L)=FP(L)-G*SPB(L)*(RIFTR(L)+EVAPSW(L))
        END DO
      END IF
C
      IF (MDCHH.GE.1) THEN
        MDCHITR=0
        DO NMD=1,MDCHH
        IF (MDCHTYP(NMD).EQ.1) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANU(NMD)
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANV(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*(QCHANU(NMD)+QCHANV(NMD))
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        END DO
      END IF
C
      GO TO 1000
C
C**********************************************************************C
C
C **  REENTER AT 3000 FOR DRYING CORRECTION
C **  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
C **  HOST-GUEST CHANNAL INTERACTIONS
C
C----------------------------------------------------------------------C
C
 3000 CONTINUE
C
      DO L=2,LA
      LN=LNC(L)
      IF (ISCDRY(L).GE.1) QSUME(L)=MAX(QSUME(L),0.0)
      FP(L)=FP1(L)-0.5*G*SPB(L)*
     $        ( SUB(L+1)*FUHDYE(L+1)-SUB(L)*FUHDYE(L)
     $        +SVB(LN )*FVHDXE(LN )-SVB(L)*FVHDXE(L)
     $        -2.0*QSUME(L) )
cc      P(L)=0.
      END DO
C
      IF (ISGWIE.GE.1) THEN
        DO L=2,LA
        FP(L)=FP(L)-G*SPB(L)*(RIFTR(L)+EVAPSW(L))
        END DO
      END IF
C
      IF (MDCHHQ.EQ.2) THEN
        MDCHITR=0
        DO NMD=1,MDCHH
        QCHANU(NMD)=0.
        QCHANV(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANVN(NMD)=0.
        END DO
      END IF
C
      IF (MDCHHQ.EQ.3.AND.NCORDRY.EQ.1) THEN
        MDCHITR=0
        DO NMD=1,MDCHH
        QCHANU(NMD)=0.
        QCHANV(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANVN(NMD)=0.
        END DO
      END IF
C
      IF (MDCHH.GE.1) THEN
        MDCHITR=0
        DO NMD=1,MDCHH
        IF (ISCDRY(LMDCHH(NMD)).GE.1) THEN
          QCHANU(NMD)=MAX(QCHANU(NMD),0.)
          QCHANV(NMD)=MAX(QCHANV(NMD),0.)
        END IF
        IF (MDCHTYP(NMD).EQ.1) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANU(NMD)
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANV(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*(QCHANU(NMD)+QCHANV(NMD))
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  RESET EQUATION COEFFICIENTS AND REENTER AT 1000 FOR NONLINEAR
C **  ITERATION ON CELL FACE DEPTHS AND AT 2000 FOR ITERATION ON
C **  HOST CELL SUBGRID CHANNEL CELL EXCHANGE
C
C----------------------------------------------------------------------C
C
 2000 CONTINUE
C
      IF (MDCHH.GE.1.AND.MDCHITR.GE.1) THEN
C
        DO L=2,LA
        LN=LNC(L)
        FP(L)=FP1(L)-0.5*G*SPB(L)*
     $        ( SUB(L+1)*FUHDYE(L+1)-SUB(L)*FUHDYE(L)
     $        +SVB(LN )*FVHDXE(LN )-SVB(L)*FVHDXE(L)
     $        -2.0*QSUME(L) )
cc        P(L)=0.
        END DO
C
        IF (ISGWIE.GE.1) THEN
          DO L=2,LA
          FP(L)=FP(L)-G*SPB(L)*(RIFTR(L)+EVAPSW(L))
          END DO
        END IF
C
        DO NMD=1,MDCHH
        IF (ISCDRY(LMDCHH(NMD)).GE.1) THEN
          QCHANU(NMD)=MAX(QCHANU(NMD),0.)
          QCHANV(NMD)=MAX(QCHANV(NMD),0.)
        END IF
        IF (MDCHTYP(NMD).EQ.1) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANU(NMD)
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*QCHANV(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          FP(LMDCHH(NMD))=FP(LMDCHH(NMD))+G*(QCHANU(NMD)+QCHANV(NMD))
          FP(LMDCHU(NMD))=FP(LMDCHU(NMD))-G*QCHANU(NMD)
          FP(LMDCHV(NMD))=FP(LMDCHV(NMD))-G*QCHANV(NMD)
        END IF
        END DO
C
      END IF
C
 1000 CONTINUE
C
      CCMNM=1.E+18
      DO L=2,LA
      IF (SPB(L).GT.0.) THEN
      LN=LNC(L)
      C1=-0.5*DELTD2*G*SPB(L)
      CS(L)=C1*SVB(L  )*HRVO(L  )*RCY(L  )*HVTMP(L  )
C    $     +(1.-SPB(L))*CS(L)
      CW(L)=C1*SUB(L  )*HRUO(L  )*RCX(L  )*HUTMP(L  )
C    $     +(1.-SPB(L))*CW(L)
      CE(L)=C1*SUB(L+1)*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
C    $     +(1.-SPB(L))*CE(L)
      CN(L)=C1*SVB(LN )*HRVO(LN )*RCY(LN )*HVTMP(LN )
C    $     +(1.-SPB(L))*CN(L)
      CC(L)=SPB(L)*(DELTI*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
C    $     +(1.-SPB(L))*CC(L)
      END IF
      CCMNM=MIN(CCMNM,CC(L))
      FPTMP(L)=FP(L)
      END DO
C
      CCMNMI=1./CCMNM
C     WRITE(6,666) CCMNM,CCMNMI
C 666 FORMAT(' CCMIN, CCMINI = ',E12.4,1X,E12.4)
C
      DO LL=1,NPBW
      IF(ISPBW(LL).EQ.0)THEN
        L=LPBW(LL)
        CW(L+1)=0.
      END IF
      END DO
C
      DO LL=1,NPBE
      IF(ISPBE(LL).EQ.0)THEN
        L=LPBE(LL)
        CE(L-1)=0.
      END IF
      END DO
C
      DO LL=1,NPBS
      IF(ISPBS(LL).EQ.0)THEN
        L=LPBS(LL)
        LN=LNC(L)
        CS(LN)=0.
      END IF
      END DO
C
      DO LL=1,NPBN
      IF(ISPBN(LL).EQ.0)THEN
        L=LPBN(LL)
        LS=LSC(L)
        CN(LS)=0.
      END IF
      END DO
C
      DO L=2,LA
      CCS(L)=CS(L)*CCMNMI
      CCW(L)=CW(L)*CCMNMI
      CCE(L)=CE(L)*CCMNMI
      CCN(L)=CN(L)*CCMNMI
      CCC(L)=CC(L)*CCMNMI
      FPTMP(L)=FPTMP(L)*CCMNMI
      CCCI(L)=1./CCC(L)
      END DO
C
      IF (ISDSOLV.EQ.1) THEN
        OPEN(1,FILE='eqcoef.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     $               FPTMP(L)
        END DO
        CLOSE(1)
        IF(N.EQ.1) THEN
        OPEN(1,FILE='eqcoef1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     $               FPTMP(L)
        END DO
        CLOSE(1)
        END IF
      END IF
C
      IF (ISDSOLV.EQ.1.AND.N.EQ.1) THEN
        OPEN(1,FILE='eqterm.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     $               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
        END DO
        CLOSE(1)
      END IF
      IF (ISDSOLV.EQ.1.AND.N.EQ.2) THEN
        OPEN(1,FILE='eqterm.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,NCORDRY
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),ISCDRY(L),SUB(L),SVB(L),HRUO(L),
     $               HRVO(L),RCX(L),RCY(L),HUTMP(L),HVTMP(L)
         END DO
        CLOSE(1)
      END IF
 1001 FORMAT(2I5,6(1X,E12.4))
 1002 FORMAT(3I4,8(1X,E9.2))
C
C
cc      IF (IRVEC.EQ.2.OR.IRVEC.GE.4) THEN
cc        DO L=2,LA
cc        LN=LNC(L)
cc        LS=LSC(L)
cc        FPTMP(L)=FPTMP(L)-CC(L)*PAM(L)-CS(L)*PAM(LS)-CW(L)*PAM(L-1)
cc     $          -CE(L)*PAM(L+1)-CN(L)*PAM(LN)
cc        END DO
cc        GO TO 5000
cc      END IF
C
C----------------------------------------------------------------------C
C
      IF (IRVEC.LE.1.OR.IRVEC.EQ.3) THEN
C
      IF (ISTL.EQ.3) THEN
C
      DO LR=1,NRC
      L=LRC(LR)
      CSR(LR)=CS(L)*CCI(L)
      CWR(LR)=CW(L)*CCI(L)
      CER(LR)=CE(L)*CCI(L)
      CNR(LR)=CN(L)*CCI(L)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      CSB(LB)=CS(L)*CCI(L)
      CWB(LB)=CW(L)*CCI(L)
      CEB(LB)=CE(L)*CCI(L)
      CNB(LB)=CN(L)*CCI(L)
      END DO
C
      ELSE
C
      DO LR=1,NRC
      L=LRC(LR)
      CCSR(LR)=CS(L)*CCI(L)
      CCWR(LR)=CW(L)*CCI(L)
      CCER(LR)=CE(L)*CCI(L)
      CCNR(LR)=CN(L)*CCI(L)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      CCSB(LB)=CS(L)*CCI(L)
      CCWB(LB)=CW(L)*CCI(L)
      CCEB(LB)=CE(L)*CCI(L)
      CCNB(LB)=CN(L)*CCI(L)
      END DO
C
      END IF
C
      DO LR=1,NRC
      L=LRC(LR)
      FPR(LR)=FPTMP(L)*CCI(L)
cc      PRED(LR)=PAM(L)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      FPB(LB)=FPTMP(L)*CCI(L)
cc      PBLK(LB)=PAM(L)
      END DO
C
cc      IF (ISTL.EQ.3) THEN
C
cc      DO L=1,NRC
cc      LN=LBNRC(L)
cc      LS=LBSRC(L)
cc      LE=LBERC(L)
cc      LW=LBWRC(L)
cc      FPR(L)=FPR(L)-PRED(L)-CSR(L)*PBLK(LS)-CWR(L)*PBLK(LW)
cc     $                     -CER(L)*PBLK(LE)-CNR(L)*PBLK(LN)
cc      END DO
C
cc      DO L=1,NBC
cc      LN=LRNBC(L)
cc      LS=LRSBC(L)
cc      LE=LREBC(L)
cc      LW=LRWBC(L)
cc      FPB(L)=FPB(L)-PBLK(L)-CSB(L)*PRED(LS)-CWB(L)*PRED(LW)
cc     $                     -CEB(L)*PRED(LE)-CNB(L)*PRED(LN)
cc      END DO
C
cc      ELSE
C
cc      DO L=1,NRC
cc      LN=LBNRC(L)
cc      LS=LBSRC(L)
cc      LE=LBERC(L)
cc      LW=LBWRC(L)
cc      FPR(L)=FPR(L)-PRED(L)-CCSR(L)*PBLK(LS)-CCWR(L)*PBLK(LW)
cc     $                     -CCER(L)*PBLK(LE)-CCNR(L)*PBLK(LN)
cc      END DO
C
cc      DO L=1,NBC
cc      LN=LRNBC(L)
cc      LS=LRSBC(L)
cc      LE=LREBC(L)
cc      LW=LRWBC(L)
cc      FPB(L)=FPB(L)-PBLK(L)-CCSB(L)*PRED(LS)-CCWB(L)*PRED(LW)
cc     $                     -CCEB(L)*PRED(LE)-CCNB(L)*PRED(LN)
cc      END DO
C
cc      END IF
C
      END IF
C
C----------------------------------------------------------------------C
C
 5000 CONTINUE
C
      IF (IRVEC.EQ.0) CALL RELAX (ISTL)
      IF (IRVEC.EQ.1) CALL RELAXV (ISTL)
      IF (IRVEC.EQ.2) CALL CONGRAD (ISTL)
      IF (IRVEC.EQ.3) CALL CGRS (ISTL)
      IF (IRVEC.GE.4) THEN
        IF(IRVEC.GE.40.AND.IRVEC.LT.50) THEN
          ITMP=IRVEC-40
         ELSE
          ITMP=IRVEC-50
        END IF
        CALL LINBCG (LC,FPTMP,P,ITMP,RSQM,ITERM,ITER,RSQ,IRVEC)
      END IF
C
      ITERHP=ITERHP+1
C
      IF(ISDRY.LT.2) THEN
        DO L=2,LA
        LS=LSC(L)
        HUTMP(L)=0.5*GI*(P(L)-BELV(L)+P(L-1)-BELV(L-1))
        HVTMP(L)=0.5*GI*(P(L)-BELV(L)+P(LS )-BELV(LS ))
        END DO
      END IF
C
      IF (ITERHP.LE.ITERHPM.AND.ISDRY.LT.2) GO TO 1000
C
C**********************************************************************C
C
C **   UPDATE SURFACE ELEVATION
C
C----------------------------------------------------------------------C
C
cc      DO L=2,LA
cc      P(L)=P(L)+PAM(L)
cc      END DO
C
C**********************************************************************C
C
C **  CALCULATE UHEX AND VHEX AND TOTAL DEPTHS AT TIME LEVEL (N+1)
C **  HRU=SUB*DYU/DXU & HRV=SVB*DXV/DYV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      LS=LSC(L)
      UTMP=SUB(L)*( FUHDYE(L)
     $            -DELTD2*HRUO(L)*RCX(L)*HUTMP(L)*(P(L)-P(L-1)) )
      VTMP=SVB(L)*( FVHDXE(L)
     $            -DELTD2*HRVO(L)*RCY(L)*HVTMP(L)*(P(L)-P(LS )) )
      UHDYE(L)=UTMP
      VHDXE(L)=VTMP
      UHE(L)=UTMP*DYIU(L)
      VHE(L)=VTMP*DXIV(L)
      END DO
C
C**********************************************************************C
C
C **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL
C **  TRANSPORTS AT (N+1)
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.3) THEN
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
        IF (ISGWIE.GE.1) HPPTMP=HPPTMP
     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        END DO
       ELSE
        DO L=2,LA
        LN=LNC(L)
        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
        IF (ISGWIE.GE.1) HPPTMP=HPPTMP
     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
        END DO
      END IF
C
      IF(ISDRY.EQ.99) GO TO 9999
C
C**********************************************************************C
C
C **  ADD CHANNEL INTERACTION EXCHANGES
C
C----------------------------------------------------------------------C
C
      IF (MDCHH.GE.1) THEN
        DO NMD=1,MDCHH
        IF (MDCHTYP(NMD).EQ.1) THEN
          HP(LMDCHH(NMD))=HP(LMDCHH(NMD))
     $                   +DELT*DXYIP(LMDCHH(NMD))*QCHANU(NMD)
          HP(LMDCHU(NMD))=HP(LMDCHU(NMD))
     $                   -DELT*DXYIP(LMDCHU(NMD))*QCHANU(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          HP(LMDCHH(NMD))=HP(LMDCHH(NMD))
     $                   +DELT*DXYIP(LMDCHH(NMD))*QCHANV(NMD)
          HP(LMDCHV(NMD))=HP(LMDCHV(NMD))
     $                   -DELT*DXYIP(LMDCHV(NMD))*QCHANV(NMD)
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          HP(LMDCHH(NMD))=HP(LMDCHH(NMD))
     $                   +DELT*DXYIP(LMDCHH(NMD))*QCHANU(NMD)
     $                   +DELT*DXYIP(LMDCHH(NMD))*QCHANV(NMD)
          HP(LMDCHU(NMD))=HP(LMDCHU(NMD))
     $                   -DELT*DXYIP(LMDCHU(NMD))*QCHANU(NMD)
          HP(LMDCHV(NMD))=HP(LMDCHV(NMD))
     $                   -DELT*DXYIP(LMDCHV(NMD))*QCHANV(NMD)
        END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  REVISE P
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      END DO
C
C**********************************************************************C
C
C **  UPDATE TEMPORARY CELL FACE DEPTHS FOR NONLINEAR ITERATIONS
C
C----------------------------------------------------------------------C
C
      IF(ISDRY.LT.2) THEN
        DO L=2,LA
        LS=LSC(L)
        HUTMP(L)=0.5*(HP(L)+HP(L-1))
        HVTMP(L)=0.5*(HP(L)+HP(LS))
        END DO
      END IF
C
C**********************************************************************C
C
C **  CHECK FOR DRYING IN HOST CELLS BEFORE PROCESSING CHANNELS
C
C----------------------------------------------------------------------C
C
      IF (MDCHH.GE.1.AND.MDCHHD2.EQ.1) THEN
C
      DO NMD=1,MDCHH
      L=LMDCHH(NMD)
      LS=LSC(L)
      LN=LNC(L)
      IDRYTMP=0
      HUTMPP=0.5*(HP(L)+HP(L-1))
      IF (HUTMPP.LE.HUWET(L  ).OR.SUBO(L  ).EQ.0.0) THEN
        SUB(L)=0.
        SBX(L)=0.
        IDRYTMP=IDRYTMP+1
      END IF
      HUTMPP=0.5*(HP(L)+HP(L+1))
      IF (HUTMPP.LE.HUWET(L+1).OR.SUBO(L+1).EQ.0.0) THEN
        SUB(L+1)=0.
        SBX(L+1)=0.
        IDRYTMP=IDRYTMP+1
      END IF
      HVTMPP=0.5*(HP(L)+HP(LS ))
      IF (HVTMPP.LE.HVWET(L  ).OR.SVBO(L  ).EQ.0.0) THEN
        SVB(L)=0.
        SBY(L)=0.
        IDRYTMP=IDRYTMP+1
      END IF
      HVTMPP=0.5*(HP(L)+HP(LN ))
      IF (HVTMPP.LE.HVWET(LN ).OR.SVBO(LN ).EQ.0.0) THEN
        SVB(LN)=0.
        SBY(LN)=0.
        IDRYTMP=IDRYTMP+1
      END IF
      IF (HP(L).LE.HDRY.AND.ISCDRY(L).EQ.0) THEN
        ISCDRY(L)=1
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
       ELSE
        IF (IDRYTMP.EQ.4.AND.ISCDRY(L).EQ.0) ISCDRY(L)=1
      END IF
      IF (ISGWIE.EQ.1) THEN
        IF (ISCDRY(L).GE.1) THEN
          IF (HP(L).LT.HDRY) THEN
            RIADDEV=RIFTR(L)+EVAPSW(L)
            IF (RIADDEV.GT.0.0) THEN
              HREDUCE=DELTI*DXYP(L)*(1.01*HDRY-HP(L))
              RPERCNT=RIFTR(L)/RIADDEV
              EPERCNT=EVAPSW(L)/RIADDEV
              RIFTR(L)=RIFTR(L)-RPERCNT*HREDUCE
              EVAPSW(L)=EVAPSW(L)-EPERCNT*HREDUCE
              RIFTR(L)=MAX(RIFTR(L),0.)
              EVAPSW(L)=MAX(EVAPSW(L),0.)
            END IF
          END IF
        END IF
      END IF
      IF (ISGWIE.EQ.2) THEN
        IF (ISCDRY(L).GE.1) THEN
          RIFTR(L)=0.0
          EVAPSW(L)=0.0
        END IF
      END IF
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  PROCESS SUBGRID SCALE CHANNELS AND UPDATE INTERACTIONS
C **  AND REVISE EXCHANGE FLOWS
C
C----------------------------------------------------------------------C
C
      IF (MDCHH.GE.1) THEN
C
        DO NMD=1,MDCHH
        QCHANUN(NMD)=0.
        QCHANVN(NMD)=0.
        END DO
        QCHUERR=-1.E+10
        QCHVERR=-1.E+10
C
        DO NMD=1,MDCHH
C
C       BEGINNING CHANNEL PROCESS FOR WET HOST CELLS
C
        IF (HP(LMDCHH(NMD)).GT.HDRY.AND.
     $      ISCDRY(LMDCHH(NMD)).EQ.0) THEN
C
C         PROCESS U GUEST CHANNEL ONLY
C
          IF (MDCHTYP(NMD).EQ.1) THEN
            ISUDPC(NMD)=0
            PCNEW=( P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $             +P(LMDCHU(NMD))*DXYP(LMDCHU(NMD)) )
     $            /( DXYP(LMDCHH(NMD))+DXYP(LMDCHU(NMD)) )
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPU=GI*PCNEW-BELV(LMDCHU(NMD))
            IF (HPPTMPH.GT.0.0) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHU(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
                IF (ISGWIE.GE.1) QCHANUN(NMD)=QCHANUN(NMD)
     $                          -RIFTR(L)-EVAPSW(L)
               ELSE
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
                IF (ISGWIE.GE.1) QCHANUN(NMD)=QCHANUN(NMD)
     $                          -RIFTR(L)-EVAPSW(L)
              END IF
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
             ELSE
              ISUDPC(NMD)=1
              QCHANUN(NMD)=0.0
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
            END IF
          END IF
C
C         PROCESS V GUEST CHANNEL ONLY
C
          IF (MDCHTYP(NMD).EQ.2) THEN
            ISUDPC(NMD)=0
            PCNEW=( P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $             +P(LMDCHV(NMD))*DXYP(LMDCHV(NMD)) )
     $            /( DXYP(LMDCHH(NMD))+DXYP(LMDCHV(NMD)) )
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPV=GI*PCNEW-BELV(LMDCHV(NMD))
            IF (HPPTMPH.GT.HDRY) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
             ELSE
              PCNEW=G*(HDRY+BELV(LMDCHH(NMD)))
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
            END IF
          END IF
C
C         PROCESS U AND V GUEST CHANNELS
C
          IF (MDCHTYP(NMD).EQ.3) THEN
            ISUDPC(NMD)=0
            PCNEW=(P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $           +P(LMDCHU(NMD))*DXYP(LMDCHU(NMD))
     $           +P(LMDCHV(NMD))*DXYP(LMDCHV(NMD)))
     $           /(DXYP(LMDCHH(NMD))+DXYP(LMDCHU(NMD))
     $                              +DXYP(LMDCHV(NMD)))
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPU=GI*PCNEW-BELV(LMDCHU(NMD))
            HPPTMPV=GI*PCNEW-BELV(LMDCHV(NMD))
            IF (HPPTMPH.GT.HDRY) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHU(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
             ELSE
              PCNEW=G*(HDRY+BELV(LMDCHH(NMD)))
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHU(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
            END IF
          END IF
C
C       END CHANNEL PROCESSING FOR WET HOST CELLS
C
        END IF
C
C       BEGINNING CHANNEL PROCESSING FOR DRY HOST CELL IF WETTING IS
C       ALLOWED
C
        IF (MDCHHD.EQ.1) THEN
        IF (HP(LMDCHH(NMD)).LE.HDRY) THEN
        IF (ISCDRY(LMDCHH(NMD)).LT.NDRYSTP) THEN
           QCHANU(NMD)=0.
           QCHANV(NMD)=0.
          ELSE
C
C         PROCESS U GUEST CHANNEL ONLY
C
          IF (MDCHTYP(NMD).EQ.1) THEN
            ISUDPC(NMD)=0
            PCNEW=( P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $             +P(LMDCHU(NMD))*DXYP(LMDCHU(NMD)) )
     $            /( DXYP(LMDCHH(NMD))+DXYP(LMDCHU(NMD)) )
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPU=GI*PCNEW-BELV(LMDCHU(NMD))
            IF (HPPTMPH.GT.HP(LMDCHH(NMD))) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHU(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
                IF (ISGWIE.GE.1) QCHANUN(NMD)=QCHANUN(NMD)
     $                          -RIFTR(L)-EVAPSW(L)
               ELSE
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
                IF (ISGWIE.GE.1) QCHANUN(NMD)=QCHANUN(NMD)
     $                          -RIFTR(L)-EVAPSW(L)
              END IF
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
             ELSE
              ISUDPC(NMD)=1
              QCHANUN(NMD)=0.0
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
            END IF
          END IF
C
C         PROCESS V GUEST CHANNEL ONLY
C
          IF (MDCHTYP(NMD).EQ.2) THEN
            ISUDPC(NMD)=0
            PCNEW=( P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $             +P(LMDCHV(NMD))*DXYP(LMDCHV(NMD)) )
     $            /( DXYP(LMDCHH(NMD))+DXYP(LMDCHV(NMD)) )
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPV=GI*PCNEW-BELV(LMDCHV(NMD))
            IF (HPPTMPH.GT.HP(LMDCHH(NMD))) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
             ELSE
              ISUDPC(NMD)=0
            END IF
          END IF
C
C         PROCESS U AND V GUEST CHANNELS
C
          IF (MDCHTYP(NMD).EQ.3) THEN
            ISUDPC(NMD)=0
            PCNEW=(P(LMDCHH(NMD))*DXYP(LMDCHH(NMD))
     $           +P(LMDCHU(NMD))*DXYP(LMDCHU(NMD))
     $           +P(LMDCHV(NMD))*DXYP(LMDCHV(NMD)))
     $           /(DXYP(LMDCHH(NMD))+DXYP(LMDCHU(NMD))
     $                              +DXYP(LMDCHV(NMD)))
            HPPTMPH=GI*PCNEW-BELV(LMDCHH(NMD))
            HPPTMPU=GI*PCNEW-BELV(LMDCHU(NMD))
            HPPTMPV=GI*PCNEW-BELV(LMDCHV(NMD))
            IF (HPPTMPH.GT.HP(LMDCHH(NMD))) THEN
              ISUDPC(NMD)=1
              PMDCH(NMD)=PCNEW
              L=LMDCHU(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANUN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANUN(NMD)-QCHANU(NMD))
              QCHUERR=MAX(QCHUERR,AQTMP)
              L=LMDCHV(NMD)
              LN=LNC(L)
              IF (ISTL.EQ.3) THEN
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H2P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $          +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))
               ELSE
                QCHANVN(NMD)=
     $          -(GI*PCNEW-BELV(L)-H1P(L))*DELTI*DXYP(L)+QSUME(L)
     $          -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $          +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))
              END IF
              AQTMP=ABS(QCHANVN(NMD)-QCHANV(NMD))
              QCHVERR=MAX(QCHVERR,AQTMP)
             ELSE
              ISUDPC(NMD)=0
            END IF
          END IF
C
C       END CHANNEL PROCESSING FOR WETTING HOST CELLS
C
        END IF
        END IF
        END IF
C
        END DO
C
        IF(QCHUERR.LT.QCHERR.AND.QCHVERR.LT.QCHERR) GO TO 4000
C
        IF (MDCHITR.LT.MDCHITM) THEN
          DO NMD=1,MDCHH
          IF (ISUDPC(NMD).EQ.1) THEN
            QCHANU(NMD)=QCHANUN(NMD)
            QCHANV(NMD)=QCHANVN(NMD)
C217           P(LMDCHH(NMD))=PMDCH(NMD)
C217           P(LMDCHU(NMD))=PMDCH(NMD)
          END IF
          END DO
          DO L=2,LA
C217         P(L)=P(L)-PAM(L)
          P(L)=0.
          END DO
          MDCHITR=MDCHITR+1
          GO TO 2000
        END IF
C
      END IF
C
      IF (MDCHH.GE.1) THEN
      OPEN(1,FILE='modchan.dia',STATUS='UNKNOWN')
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      WRITE(1,8001) TIME
      DO NMD=1,MDCHH
      WRITE(1,8002)NMD,QCHANU(NMD),QCHANUN(NMD)
      END DO
      CLOSE(1)
      END IF
C
 4000 CONTINUE
      QCHERM=MAX(QCHUERR,QCHVERR)
      IF (MDCHH.GE.1) WRITE(8,8000)MDCHITR,QCHERM
C
 8000 FORMAT(' MDCHITR = ',I5,'  QCHERM = ',E13.5)
 8001 FORMAT(' TIME = ',F12.5,/)
 8002 FORMAT(2X,I5,2(2X,E12.4))
C
C**********************************************************************C
C
C **  CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY
C
C----------------------------------------------------------------------C
C
      OPEN(1,FILE='drywet.log',ACCESS='APPEND',STATUS='UNKNOWN')
C
      HWMDI=1./(HWET-HDRY)
      ICORDRY=0
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      IDRYTMP=0
      RDRYTMP=0.
      HUTMPP=0.5*(HP(L)+HP(L-1))
      IF (HUTMPP.LE.HUWET(L  ).OR.SUBO(L  ).EQ.0.0) THEN
        IF(SUB(L).EQ.1.) THEN
          ICORDRY=1
C         WRITE(1,6941)N,IL(L),JL(L),HU(L),HP(L),H1P(L)
C         WRITE(6,6941)N,IL(L),JL(L),HU(L),HP(L),H1P(L)
C         WRITE(8,6941)N,IL(L),JL(L),HU(L),HP(L),H1P(L)
        END IF
        SUTMP=HWMDI*(HUTMPP-HUDRY(L))
        SUTMP=MAX(SUTMP,0.)
        SUB(L)=SUBO(L)*SUTMP
        SBX(L)=SBXO(L)*SUTMP
        IDRYTMP=IDRYTMP+1
        RDRYTMP=RDRYTMP+SUBO(L  )
      END IF
      HUTMPP=0.5*(HP(L)+HP(L+1))
      IF (HUTMPP.LE.HUWET(L+1).OR.SUBO(L+1).EQ.0.0) THEN
        IF(SUB(L+1).EQ.1) THEN
          ICORDRY=1
C         WRITE(1,6942)N,IL(L),JL(L),HU(L+1),HP(L),H1P(L)
C         WRITE(6,6942)N,IL(L),JL(L),HU(L+1),HP(L),H1P(L)
C         WRITE(8,6942)N,IL(L),JL(L),HU(L+1),HP(L),H1P(L)
        END IF
        SUTMP=HWMDI*(HUTMPP-HUDRY(L+1))
        SUTMP=MAX(SUTMP,0.)
        SUB(L+1)=SUBO(L+1)*SUTMP
        SBX(L+1)=SBXO(L+1)*SUTMP
        IDRYTMP=IDRYTMP+1
        RDRYTMP=RDRYTMP+SUBO(L+1)
      END IF
      HVTMPP=0.5*(HP(L)+HP(LS))
      IF (HVTMPP.LE.HVWET(L).OR.SVBO(L  ).EQ.0.0) THEN
        IF(SVB(L).EQ.1.) THEN
          ICORDRY=1
C         WRITE(1,6943)N,IL(L),JL(L),HV(L),HP(L),H1P(L)
C         WRITE(6,6943)N,IL(L),JL(L),HV(L),HP(L),H1P(L)
C         WRITE(8,6943)N,IL(L),JL(L),HV(L),HP(L),H1P(L)
        END IF
        SVTMP=HWMDI*(HVTMPP-HVDRY(L))
        SVTMP=MAX(SVTMP,0.)
        SVB(L)=SVBO(L)*SVTMP
        SBY(L)=SBYO(L)*SVTMP
        IDRYTMP=IDRYTMP+1
        RDRYTMP=RDRYTMP+SVBO(L  )
      END IF
      HVTMPP=0.5*(HP(L)+HP(LN))
      IF (HVTMPP.LE.HVWET(LN).OR.SVBO(LN ).EQ.0.0) THEN
        IF(SVB(LN).EQ.1.) THEN
          ICORDRY=1
C         WRITE(1,6944)N,IL(L),JL(L),HV(LN),HP(L),H1P(L)
C         WRITE(6,6944)N,IL(L),JL(L),HV(LN),HP(L),H1P(L)
C         WRITE(8,6944)N,IL(L),JL(L),HV(LN),HP(L),H1P(L)
        END IF
        SVTMP=HWMDI*(HVTMPP-HVDRY(LN))
        SVTMP=MAX(SVTMP,0.)
        SVB(LN)=SVBO(LN)*SVTMP
        SBY(LN)=SBYO(LN)*SVTMP
        IDRYTMP=IDRYTMP+1
        RDRYTMP=RDRYTMP+SVBO(LN )
      END IF
      IF (HP(L).LE.HDRY) THEN
        IF (ISCDRY(L).EQ.0) THEN
          ISCDRY(L)=1
          ICORDRY=1
C         WRITE(1,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(6,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(8,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
        END IF
        SUB(L)=0.
        SVB(L)=0.
        SUB(L+1)=0.
        SVB(LN)=0.
        SBX(L)=0.
        SBY(L)=0.
        SBX(L+1)=0.
        SBY(LN)=0.
       ELSE
        IF(IDRYTMP.EQ.4.AND.RDRYTMP.LT.0.5) THEN
          IF (ISCDRY(L).EQ.0) THEN
          ISCDRY(L)=1
          ICORDRY=1
C         WRITE(1,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(6,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
C         WRITE(8,6945)N,IL(L),JL(L),HP(L),H1P(L),H2P(L)
          END IF
        END IF
      END IF
      IF (ISGWIE.EQ.1) THEN
        IF (ISCDRY(L).GE.1) THEN
          IF (HP(L).LT.HDRY) THEN
            RIADDEV=RIFTR(L)+EVAPSW(L)
            IF (RIADDEV.GT.0.0) THEN
              HREDUCE=DELTI*DXYP(L)*(1.01*HDRY-HP(L))
              RPERCNT=RIFTR(L)/RIADDEV
              EPERCNT=EVAPSW(L)/RIADDEV
              RIFTR(L)=RIFTR(L)-RPERCNT*HREDUCE
              EVAPSW(L)=EVAPSW(L)-EPERCNT*HREDUCE
              RIFTR(L)=MAX(RIFTR(L),0.)
              EVAPSW(L)=MAX(EVAPSW(L),0.)
            END IF
          END IF
        END IF
      END IF
      IF (ISGWIE.EQ.2) THEN
        IF (ISCDRY(L).GE.1) THEN
          RIFTR(L)=0.0
          EVAPSW(L)=0.0
        END IF
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
        ITERHP=0
        GO TO 3000
      END IF
C
C     WRITE(6,6960)NCORDRY
CTMP      WRITE(8,6960)NCORDRY
C
 6960 FORMAT(' NCORDRY =', I5)
 6961 FORMAT(' UNSTABLE, NCORDRY =', I5)
C
 9999 CONTINUE
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      RSMALL=1.E-12
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
      LS=LSC(L)
      HP(L)=MAX(HP(L),0.)
      HU(L)=MAX(HU(L),0.)
      HV(L)=MAX(HV(L),0.)
      END DO
C
      DO L=2,LA
      HPI(L)=1./(HP(L)+RSMALL)
      HUI(L)=1./(HU(L)+RSMALL)
      HVI(L)=1./(HV(L)+RSMALL)
      END DO
C
C**********************************************************************C
C
C **  PERFORM UPDATE ON GROUNDWATER ELEVATION
C
C----------------------------------------------------------------------C
C
      IF (ISGWIE.GE.1) THEN
C
        DO L=2,LA
        QSUM(L,KC)=QSUM(L,KC)-EVAPSW(L)
        QSUM(L,1 )=QSUM(L,1 )-RIFTR(L)
        END DO
C
C       INFILTRATION STEP
C
        RNPORI=1./RNPOR
        IF (ISTL.EQ.3) THEN
          DO L=2,LA
          AGWELV(L)=AGWELV2(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          END DO
         ELSE
          DO L=2,LA
          AGWELV(L)=AGWELV1(L)+RNPORI*DELT*DXYIP(L)*RIFTR(L)
          END DO
        END IF
        DO L=2,LA
        AGWELV(L)=MIN(AGWELV(L),BELV(L))
        END DO
C
C       ET STEP
C
        DO L=2,LA
        IF(EVAPCVT.LT.0.) THEN
          SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     $              (1.+0.00412*TEM(L,KC))))
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L)
        END IF
        ETGWTMP=EVAPT(L)-EVAPSW(L)*DXYIP(L)
        ETGWTMP=MAX(ETGWTMP,0.0)
        ETGWAVL=RNPOR*DELTI*(AGWELV(L)-BELAGW(L))
        ETGWAVL=MAX(ETGWAVL,0.0)
        ETGWTMP=MIN(ETGWTMP,ETGWAVL)
        EVAPGW(L)=ETGWTMP*DXYP(L)
        END DO
        DO L=2,LA
        AGWELV(L)=AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)
        END DO
        DO L=2,LA
        AGWELV(L)=MAX(AGWELV(L),BELAGW(L))
        END DO
C
      END IF
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
      IF(SPB(L).NE.0) THEN
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H2P(L))*DELTI
     $     +0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $     +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L))-QSUME(L)
     $     +RIFTR(L)+EVAPSW(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      END IF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      END IF
      END IF
      END DO
C
      ELSE
C
      DO L=2,LA
      IF(SPB(L).NE.0) THEN
      LN=LNC(L)
      DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI
     $     +0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $     +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L)
     $     +RIFTR(L)+EVAPSW(L))
      IF(DIVEX.GT.DIVEXMX)THEN
       DIVEXMX=DIVEX
       LMAX=L
      END IF
      IF(DIVEX.LT.DIVEXMN)THEN
       DIVEXMN=DIVEX
       LMIN=L
      END IF
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
C **  UPDATE ZERO DIMENSION VOLUME BALANCE
C
C----------------------------------------------------------------------C
C
      IF (ISDRY.GE.1.AND.ISTL.EQ.3) THEN
        VOLADD=0.
        DO L=2,LA
        IF (SPB(L).NE.0) THEN
          VOLADD=VOLADD+QSUME(L)-RIFTR(L)-EVAPSW(L)
        END IF
        END DO
        VOLADD=VOLADD*DT
        VOLZERD=VOLZERD+VOLADD
        VETZERD=VETZERD+VOLADD+DT*EVAPSW(L)
      END IF

 5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C
C**********************************************************************C
C
      RETURN
      END
