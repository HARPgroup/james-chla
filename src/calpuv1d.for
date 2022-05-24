C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALPUV1D (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999  
C
C**********************************************************************C
C
C ** SUBROUTINE CALPUV1D CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
C ** AND VHDXE, FOR 1D CHANNEL NETWORKD FREE SURFACE FLOWS WITH 
C ** PROVISIONS FOR WETTING AND DRYING OF CELLS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION QSUMTMP(LCM)
C
C**********************************************************************C
C
c      WRITE(6,6000)N
c 6000 FORMAT(' CALLED CALPUV9, N = ',I10)
C
      IF (N.EQ.1.AND.ISDSOLV.GT.1) THEN
        OPEN(1,FILE='eqcoef.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqterm.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqterm1.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='fp.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqcoef1.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      END IF
C
      IF (ISDSOLV.EQ.1) THEN
        OPEN(1,FILE='eqcoef.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqterm.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='fp.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
      END IF
C
CDIAG      IF (N.EQ.1.AND.ISDSOLV.GE.1) THEN
CDIAG        OPEN(1,FILE='chanpuv.dia',STATUS='UNKNOWN')
CDIAG        CLOSE(1,STATUS='DELETE')
CDIAG      END IF
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
C **  CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS
C **  HRU=SUB*HMU*DYU/DXU & HRV=SVB*HMV*DXV/DYV 
C
C----------------------------------------------------------------------C
C
CDIAG      OPEN(1,FILE='chanpuv.dia',ACCESS='APPEND')
C
CDIAG      WRITE(1,2111)N,ISTL
CDIAG      WRITE(1,2112)
C
C
      IF(ISITB.GE.1) THEN
        DO L=2,LA
         TMPVALU=WPDYU(L)/FADYU(L)
         TMPVALV=WPDXV(L)/FADXV(L)
C MANNING
         STBX(L)=G*ZBR(L)*ZBR(L)*(TMPVALU**0.333)
         STBY(L)=G*ZBR(L)*ZBR(L)*(TMPVALV**0.333)
C MANNING
         TBX1(L)=TMPVALU*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1))
     $           *UHDYE(L)
         TBY1(L)=TMPVALV*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1))
     $           *VHDXE(L)
        END DO
      END IF
C
      IF(ISTL.EQ.2) THEN
        DO L=2,LA
C        HUTMP(L)=0.5*(HU(L)+H1U(L))
C        HVTMP(L)=0.5*(HV(L)+H1V(L))
        HUTMP(L)=FADYU(L)
        HVTMP(L)=FADXV(L)
        END DO
       ELSE
        DO L=2,LA
C        HUTMP(L)=HU(L)
C        HVTMP(L)=HV(L)
        HUTMP(L)=FADYU(L)
        HVTMP(L)=FADXV(L)
        END DO
      END IF
C
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      FUHDYE(L)=UHDYE(L)
     $         -DELTD2*SUB(L)*HRUO(L)*HUTMP(L)*(P(L)-P(L-1))
     $         +SUB(L)*DELT*SRFYU(L)*TSX(L)
     $         -SUB(L)*DELT*WPDYU(L)*RITB1*TBX1(L)
     $         +SUB(L)*DELT*DXIU(L)*( FPGXE(L)-SNLT*FXE(L) )
      FVHDXE(L)=VHDXE(L)
     $         -DELTD2*SVB(L)*HRVO(L)*HVTMP(L)*(P(L)-P(LS ))
     $         +SVB(L)*DELT*SRFXV(L)*TSY(L)
     $         -SVB(L)*DELT*WPDXV(L)*RITB1*TBY1(L)
     $         +SVB(L)*DELT*DYIV(L)*( FPGYE(L)-SNLT*FYE(L) )
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
      IF(ISITB.GE.1) THEN
        DO L=2,LA
         TMPVALU=WPDYU(L)/FADYU(L)
         TMPVALV=WPDXV(L)/FADXV(L)
C MANNING
         STBX(L)=G*ZBR(L)*ZBR(L)*(TMPVALU**0.333)
         STBY(L)=G*ZBR(L)*ZBR(L)*(TMPVALV**0.333)
C MANNING
c        RCX(L)=1./( 1.
c     $    +RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
         RCX(L)=1./( 1.
     $    +RITB*DELT*TMPVALU*STBX(L)*SQRT(VU(L)*VU(L)+U(L,1)*U(L,1)) )
c        RCY(L)=1./( 1.
c     $    +RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
         RCY(L)=1./( 1.
     $    +RITB*DELT*TMPVALV*STBY(L)*SQRT(UV(L)*UV(L)+V(L,1)*V(L,1)) )
CDIAG         WRITE(1,2113)IL(L),JL(L),FUHDYE(L),RCX(L),FVHDXE(L),RCY(L),
CDIAG     $                TMPVALV
         FUHDYE(L)=FUHDYE(L)*RCX(L)
         FVHDXE(L)=FVHDXE(L)*RCY(L)
        END DO
      END IF
C
CDIAG      CLOSE(1)
C
 2111 FORMAT(2I10)
 2112 FORMAT(' I,J,FUHDYE,RCX,FVHDXE,RCY,TMPVALV')
 2113 FORMAT(2I5,6F10.2)
C
C**********************************************************************C
C
C **  RESET BOUNDARY CONDITIONS SWITCHES
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
        LN=LNC(L)
        SUB(L)=SUBO(L)
        SVB(L)=SVBO(L)
        SBX(L)=SBXO(L)
        SBY(L)=SBYO(L)
        SUB(L+1)=SUBO(L+1)
        SVB(LN)=SVBO(LN)
        SBX(L+1)=SBXO(L+1)
        SBY(LN)=SBYO(LN)
      END DO
c
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
       CC(L)=DELTI*DXYP(L)*DADH(L)
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
      CET=0.5*DELTD2*G*HRUO(L+1)*RCX(L+1)*HUTMP(L+1)
      IF(ISPBW(LL).EQ.1)THEN
       TMP=DELTD2*SQRT(G*HMU(L+1))*DXIU(L+1)
       CC(L)=CET*(1.+TMP)/TMP
       CE(L)=-CET
       FP1(L)=CET*(2.*FP1(L)
     $ -SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/TMP
      ELSE
       FP1(L+1)=CET*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBE
      L=LPBE(LL)
       CC(L)=DELTI*DXYP(L)*DADH(L)
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
      CWT=0.5*DELTD2*G*HRUO(L  )*RCX(L)*HUTMP(L  )
      IF(ISPBE(LL).EQ.1)THEN
       TMP=DELTD2*SQRT(G*HMU(L))*DXIU(L)
       CC(L)=CWT*(1.+TMP)/TMP
       CW(L)=-CWT
       FP1(L)=CWT*(2.*FP1(L)
     $ +SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/TMP
      ELSE
       FP1(L-1)=CWT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBS
      L=LPBS(LL)
       CC(L)=DELTI*DXYP(L)*DADH(L)
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
      CNT=0.5*DELTD2*G*HRVO(LN )*RCY(LN)*HVTMP(LN )
      IF(ISPBS(LL).EQ.1)THEN
       TMP=DELTD2*SQRT(G*HMV(LN))*DYIV(LN)
       CC(L)=CNT*(1.+TMP)/TMP
       CN(L)=-CNT
       FP1(L)=CNT*(2.*FP1(L)
     $ -SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/TMP
      ELSE
       FP1(LN)=CNT*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
      DO LL=1,NPBN
      L=LPBN(LL)
       CC(L)=DELTI*DXYP(L)*DADH(L)
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
      CST=0.5*DELTD2*G*HRVO(L  )*RCY(L)*HVTMP(L  )
      IF(ISPBN(LL).EQ.1)THEN
       TMP=DELTD2*SQRT(G*HMV(L))*DYIV(L)
       CC(L)=CST*(1.+TMP)/TMP
       CS(L)=-CST
       FP1(L)=CST*(2.*FP1(L)
     $ +SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/TMP
      ELSE
       FP1(LS)=CST*FP1(L)
       FP1(L)=CC(L)*FP1(L)
      END IF
      END DO
C
C----------------------------------------------------------------------C
C
c      DO LL=1,NPBW
c      L=LPBW(LL)
c       CC(L)=DELTI*DXYP(L)
c       CS(L)=0.
c       CW(L)=0.
c       CE(L)=0.
c       CN(L)=0.
c      FP1(L)=PSERT(NPSERW(LL))
c      DO M=1,MTIDE
c      TC=CCCOS(M)
c      TS=SSSIN(M)
c      FP1(L)=FP1(L)+PCBW(LL,M)*TC+PSBW(LL,M)*TS
c      END DO
c      CET=-0.5*DELTD2*G*HRUO(L+1)*HUTMP(L+1)
c      IF(ISPBW(LL).EQ.1)THEN
c       TMP=DELT*SQRT(G*HMU(L+1))*DXIU(L+1)
c       CC(L)=CET*(1.+TMP)/(1.-TMP)
c       CE(L)=CET
c       FP1(L)=CET*(4.*FP1(L)
c     $ -2.*SQRT(G*HMU(L+1))*FUHDYE(L+1)*DYIU(L+1)*HUI(L+1))/(1.-TMP)
c      ELSE
c       FP1(L+1)=-CET*FP1(L)
c       FP1(L)=CC(L)*FP1(L)
c      END IF
c      END DO
C
C----------------------------------------------------------------------C
C
c      DO LL=1,NPBE
c      L=LPBE(LL)
c       CC(L)=DELTI*DXYP(L)
c       CS(L)=0.
c       CW(L)=0.
c       CE(L)=0.
c       CN(L)=0.      
c      FP1(L)=PSERT(NPSERE(LL))
c      DO M=1,MTIDE
c      TC=CCCOS(M)
c      TS=SSSIN(M)
c      FP1(L)=FP1(L)+PCBE(LL,M)*TC+PSBE(LL,M)*TS
c      END DO
c      CWT=-0.5*DELTD2*G*HRUO(L  )*HUTMP(L  )
c      IF(ISPBE(LL).EQ.1)THEN
c       TMP=DELT*SQRT(G*HMU(L))*DXIU(L)
c       CC(L)=CWT*(1.+TMP)/(1.-TMP)
c       CW(L)=CWT
c       FP1(L)=CWT*(4.*FP1(L)
c     $ +2.*SQRT(G*HMU(L))*FUHDYE(L)*DYIU(L)*HUI(L))/(1.-TMP)
c      ELSE
c       FP1(L-1)=-CWT*FP1(L)
c       FP1(L)=CC(L)*FP1(L)
c      END IF
c      END DO
C
C----------------------------------------------------------------------C
C
c      DO LL=1,NPBS
c      L=LPBS(LL)
c       CC(L)=DELTI*DXYP(L)
c       CS(L)=0.
c       CW(L)=0.
c       CE(L)=0.
c       CN(L)=0.
c      LN=LNC(L)
c      FP1(L)=PSERT(NPSERS(LL))
c      DO M=1,MTIDE
c      TC=CCCOS(M)
c      TS=SSSIN(M)
c      FP1(L)=FP1(L)+PCBS(LL,M)*TC+PSBS(LL,M)*TS
c      END DO
c      CNT=-0.5*DELTD2*G*HRVO(LN )*HVTMP(LN )
c      IF(ISPBS(LL).EQ.1)THEN
c       TMP=DELT*SQRT(G*HMV(LN))*DYIV(LN)
c       CC(L)=CNT*(1.+TMP)/(1.-TMP)
c       CN(L)=CNT
c       FP1(L)=CNT*(4.*FP1(L)
c     $ -2.*SQRT(G*HMV(LN))*FVHDXE(LN)*DXIV(LN)*HVI(LN))/(1.-TMP)
c      ELSE
c       FP1(LN)=-CNT*FP1(L)
c       FP1(L)=CC(L)*FP1(L)
c      END IF
c      END DO
C
C----------------------------------------------------------------------C
C
c      DO LL=1,NPBN
c      L=LPBN(LL)
c       CC(L)=DELTI*DXYP(L)
c       CS(L)=0.
c       CW(L)=0.
c       CE(L)=0.
c       CN(L)=0.
c      LS=LSC(L)
c      FP1(L)=PSERT(NPSERN(LL))
c      DO M=1,MTIDE
c      TC=CCCOS(M)
c      TS=SSSIN(M)
c      FP1(L)=FP1(L)+PCBN(LL,M)*TC+PSBN(LL,M)*TS
c      END DO
c      CST=-0.5*DELTD2*G*HRVO(L  )*HVTMP(L  )
c      IF(ISPBN(LL).EQ.1)THEN
c       TMP=DELT*SQRT(G*HMV(L))*DYIV(L)
c       CC(L)=CST*(1.+TMP)/(1.-TMP)
c       CS(L)=CST
c       FP1(L)=CST*(4.*FP1(L)
c     $ +2.*SQRT(G*HMV(L))*FVHDXE(L)*DXIV(L)*HVI(L))/(1.-TMP)
c      ELSE
c       FP1(LS)=-CST*FP1(L)
c       FP1(L)=CC(L)*FP1(L)
c      END IF
c      END DO
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
        IF (HP(L).LE.HDRY) THEN
          QSUMTMP(L)=0.
         ELSE
          QSUMTMP(L)=-(HP(L)-HDRY)*DXYP(L)*DELTI
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
          RAVAIL=(HP(L)-HDRY)*DELTI-EVAPT(L)
C         SET RIFTRL TO MIN OF AVAILABLE RATE OR LIMITING RATE
          RIFTRL=MIN(RAVAIL,RIFTRL)
C         CONVERT TO VOLUME FLOW UNITS
          RIFTR(L)=RIFTRL*DXYP(L)         
        END IF
C       ADJUST VOLUME OUTFLOWS OF WET CELLS
        IF (QSUME(L).LT.0.0) THEN
          QSUMIET=RIFTR(L)+EVAPSW(L)
          QEAVAIL=DXYP(L)*(HP(L)-HDRY)*DELTI-QSUMIET
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
c      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*P1(L)
c     $      -0.5*G*(UHDY1E(L+1)-UHDY1E(L)
c     $             +VHDX1E(LN )-VHDX1E(L)) )
      FP1(L)=FP1(L)+SPB(L)*( DELTI*DXYP(L)*DADH(L)*P(L)
     $      -0.5*G*(UHDYE(L+1)-UHDYE(L)
     $             +VHDXE(LN )-VHDXE(L)) )
      END DO
C
C**********************************************************************C
C
C **  ADVANCE EXTERNAL VARIABLES FOR THREE TIME LEVEL STEP
C
C----------------------------------------------------------------------C
C
C      IF (ISTL.EQ.3) THEN
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
        FADXP2(L)=FADXP1(L)
        FADYP2(L)=FADYP1(L)
        FADXP1(L)=FADXP(L)
        FADYP1(L)=FADYP(L)
        FADYU1(L)=FADYU(L)
        WPDYU1(L)=WPDYU(L)
        SRFYU1(L)=SRFYU(L)
        FADXV1(L)=FADXV(L)
        WPDXV1(L)=WPDXV(L)
        SRFXV1(L)=SRFXV(L)
        AGWELV2(L)=AGWELV1(L)
        AGWELV1(L)=AGWELV(L)
        END DO
C      END IF
C
cc     DO L=2,LA
cc      PAM(L)=P(L)
cc      END DO
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
      CC(L)=SPB(L)*(DELTI*DADH(L)*DXYP(L)-CS(L)-CW(L)-CE(L)-CN(L))
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
      CC(1)=1.
      CC(LC)=1.
C
C **  scale by minimum diagonal
C
c      DO L=2,LA
c      CCS(L)=CS(L)*CCMNMI
c      CCW(L)=CW(L)*CCMNMI
c      CCE(L)=CE(L)*CCMNMI
c      CCN(L)=CN(L)*CCMNMI
c      CCC(L)=CC(L)*CCMNMI
c      FPTMP(L)=FPTMP(L)*CCMNMI
c      CCCI(L)=1./CCC(L)
c      END DO
C
C
c **  scale to normal form
c
      DO L=2,LA
      CCS(L)=CS(L)/SQRT( CC(L)*CC(LSC(L)) )
      CCW(L)=CW(L)/SQRT( CC(L)*CC(L-1   ) )
      CCE(L)=CE(L)/SQRT( CC(L)*CC(L+1   ) )
      CCN(L)=CN(L)/SQRT( CC(L)*CC(LNC(L)) )
      CCC(L)=1.
      FPTMP(L)=FPTMP(L)/SQRT( CC(L) )
      CCCI(L)=1.
      END DO
C
      CALL CONGRAD (ISTL)
c
      DO L=2,LA
      P(L)=P(L)/SQRT( CC(L) )
      END DO
C
c **  CALCULATE DEPTHS IN OPEN BOUNDARY CELLS
c
      DO LL=1,NPBW
        L=LPBW(LL)
        HP(L)=GI*P(L)-BELV(L)
      END DO
C
      DO LL=1,NPBE
        L=LPBE(LL)
        HP(L)=GI*P(L)-BELV(L)
      END DO
C
      DO LL=1,NPBS
        L=LPBS(LL)
        HP(L)=GI*P(L)-BELV(L)
      END DO
C
      DO LL=1,NPBN
        L=LPBN(LL)
        HP(L)=GI*P(L)-BELV(L)
      END DO
C
c **  CALCULATE DEPTHS IN OPEN BOUNDARY CELLS
c
C
c **  WRITE DIAGNOSTIC FILES
c
      IF (ISDSOLV.GE.1) THEN
        OPEN(1,FILE='eqcoef.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     $               FPTMP(L),SURFTMP
        END DO
        CLOSE(1)
        IF(N.EQ.1) THEN
          OPEN(1,FILE='eqcoef1.out',ACCESS='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
         SURFTMP=GI*P(L)
         WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     $               FPTMP(L),SURFTMP
         WRITE(1,1001)IL(L),JL(L),CS(L),CW(L),CC(L),CE(L),CN(L),
     $               FP(L),SURFTMP
          END DO
          CLOSE(1)
        END IF
      END IF
C
      IF (ISDSOLV.GE.1) THEN
        OPEN(1,FILE='eqterm.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     $               HRVO(L),HUTMP(L),HVTMP(L)
        END DO
        CLOSE(1)
        IF (N.EQ.1) THEN
          OPEN(1,FILE='eqterm1.out',ACCESS='APPEND',STATUS='UNKNOWN')
          WRITE(1,1001)N,ISTL
          DO L=2,LA
          WRITE(1,1001)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     $               HRVO(L),HUTMP(L),HVTMP(L)
          END DO
          CLOSE(1)
        END IF
      END IF
C
 1001 FORMAT(2I5,7(1X,E12.4))
 1002 FORMAT(3I4,8(1X,E9.2))
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
C      UHE(L)=UTMP*DYIU(L)
C      VHE(L)=VTMP*DXIV(L)
      END DO
C
C**********************************************************************C
C
C **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL 
C **  TRANSPORTS AT (N+1)
C
C----------------------------------------------------------------------C
C
c      IF (ISTL.EQ.3) THEN
c        DO L=2,LA
c        LN=LNC(L)
c        HPPTMP=H2P(L)+DELT*DXYIP(L)*(QSUME(L)
c     $       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
c     $       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
c        IF (ISGWIE.GE.1) HPPTMP=HPPTMP
c     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
c        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
c        END DO
c       ELSE
c        DO L=2,LA
c        LN=LNC(L)
c        HPPTMP=H1P(L)+DELT*DXYIP(L)*(QSUME(L)
c     $       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
c     $       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
c        IF (ISGWIE.GE.1) HPPTMP=HPPTMP
c     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
c        HP(L)=SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))
c        END DO
c      END IF
C
C
      DO L=2
      IF (ISTL.EQ.3) THEN
        DO L=2,LA
         AREAOLD(L)=0.
         IF(LCT(L).EQ.6) AREAOLD(L)=FADYP2(L)
         IF(LCT(L).EQ.7) AREAOLD(L)=FADXP2(L)
        END DO
        DO L=2,LA
        LN=LNC(L)
        AREATMP=AREAOLD(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY2E(L+1)-UHDYE(L)-UHDY2E(L)
     $       +VHDXE(LN)+VHDX2E(LN)-VHDXE(L)-VHDX2E(L)))
        IF (ISGWIE.GE.1) AREATMP=AREATMP
     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
         AREANEW(L)=SPB(L)*AREATMP+(1.-SPB(L))*AREANEW(L)
        END DO
       ELSE
        DO L=2,LA
         AREAOLD(L)=0.
         IF(LCT(L).EQ.6) AREAOLD(L)=FADYP1(L)
         IF(LCT(L).EQ.7) AREAOLD(L)=FADXP1(L)
        END DO
        DO L=2,LA
        LN=LNC(L)
        AREATMP=AREAOLD(L)+DELT*DXYIP(L)*(QSUME(L)
     $       -0.5*(UHDYE(L+1)+UHDY1E(L+1)-UHDYE(L)-UHDY1E(L)
     $       +VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L)))
        IF (ISGWIE.GE.1) AREATMP=AREATMP
     $                         -DELT*DXYIP(L)*(RIFTR(L)+EVAPSW(L))
         AREANEW(L)=SPB(L)*AREATMP+(1.-SPB(L))*AREANEW(L)
        END DO
      END IF
C
      DO L=2,LA
       IF(LCT(L).EQ.6) FADYP(L)=AREANEW(L)
       IF(LCT(L).EQ.7) FADXP(L)=AREANEW(L)
      END DO
C
      IF (ISTL.EQ.3) THEN
        DO L=2,LA
        DADH1(L)=DADH(L)
        END DO
      END IF
C
      CALL CALAREA(4) 
      CALL CALAREA(3) 
CDIAG      CALL SURFPLT
C
      DO L=2,LA
       LS=LSC(L)      
       UHE(L)=UHDYE(L)/FADYU(L)
       VHE(L)=VHDXE(L)/FADXV(L)
      END DO
C
C**********************************************************************C
C
C **  PERFORM FINAL UPDATES OF P,HU, AND HV
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      P(L)=G*(HP(L)+BELV(L))
      END DO
C
      DO L=2,LA
       LS=LSC(L)      
       HU(L)=0.5*(HP(L)+HP(L-1))
       HV(L)=0.5*(HP(L)+HP(LS ))
      END DO
C
      DO L=2,LA
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      END DO
C
      DO L=2,LA
       UHE(L)=UHE(L)*HU(L)
       VHE(L)=VHE(L)*HV(L)
      END DO
C

C**********************************************************************C
C
C **  OUTPUT 1D SOLUTION DIAGNOSTICS
C
C----------------------------------------------------------------------C
C
CDIAG      OPEN(1,FILE='chanpuv.dia',ACCESS='APPEND')
C
CDIAG      WRITE(1,1111)N,ISTL
CDIAG      WRITE(1,1112)
C
CDIAG      DO L=2,LA
CDIAG      SEL=GI*P(L)
CDIAG      UTMP=UHDYE(L)/FADYU(L)
CDIAG      WRITE(1,1113)IL(L),JL(L),HP(L),SEL,FADYP(L),UHDYE(L),
CDIAG     $             UTMP,FADYU(L)
CDIAG      END DO
C
CDIAG      CLOSE(1)
C
 1111 FORMAT(2I10)
 1112 FORMAT(' I,J,H,SEL,AP,QX,U,AU')
 1113 FORMAT(2I5,6F10.2)
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
      IF (ISNEGH.GE.1) THEN
      INEGFLG=0
C
      DO L=2,LA
      IF (HP(L).LT.0.) THEN
      INEGFLG=1
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
      DO L=2,LA
      IF (HU(L).LT.0.) THEN
      INEGFLG=1
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
C
      DO L=2,LA
      IF (HV(L).LT.0.) THEN
      INEGFLG=1
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
      IF (ISNEGH.EQ.2) THEN
      IF (INEGFLG.EQ.1) THEN
C
        CALL RESTOUT(1)
C
        OPEN(1,FILE='eqcoef.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqcoef.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        SURFTMP=GI*P(L)
        WRITE(1,1001)IL(L),JL(L),CCS(L),CCW(L),CCC(L),CCE(L),CCN(L),
     $               FPTMP(L),SURFTMP
        END DO
        CLOSE(1)
C
        OPEN(1,FILE='eqterm.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='eqterm.out',ACCESS='APPEND',STATUS='UNKNOWN')
        WRITE(1,1001)N,ISTL
        DO L=2,LA
        WRITE(1,1002)IL(L),JL(L),SUB(L),SVB(L),HRUO(L),
     $               HRVO(L),HUTMP(L),HVTMP(L)
        END DO
        CLOSE(1)
C
        STOP
C
      END IF
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
 6628 FORMAT(1X,'DIVEXMX=',E13.5,5X,2I10)
 6629 FORMAT(1X,'DIVEXMN=',E13.5,5X,2I10)
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
 
 5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E13.5))
C                         
C**********************************************************************C
C
      RETURN
      END
