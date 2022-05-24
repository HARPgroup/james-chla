C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE CALTRPCB (M,NW,CON,CON1)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE
C **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO
C **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES 
C **  THE NUMBER OF TIME LEVELS IN THE STEP
C
!     Transport every 2 timestep
!     i.e.  N-2, N,  N+2, ....
!     
C     Called at even time-step
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'wq.par'
      INCLUDE 'efdc.cmn'
      DIMENSION CON(LCM,KCM), CON1(LCM,KCM)
      COMMON/WQOBC/
     $ CSERTWQ(KCWM,0:NWQCSRM,NWQVM),
     $ CWQLOS(NBBSWM,KCWM,NWQVM),CWQLOW(NBBWWM,KCWM,NWQVM),
     $ CWQLOE(NBBEWM,KCWM,NWQVM),CWQLON(NBBNWM,KCWM,NWQVM),
     $ NWQLOS(NBBSWM,KCWM,NWQVM),NWQLOW(NBBWWM,KCWM,NWQVM),
     $ NWQLOE(NBBEWM,KCWM,NWQVM),NWQLON(NBBNWM,KCWM,NWQVM),
     $ NWQOBS,NWQOBW,NWQOBE,NWQOBN,IWQOBS(NBBSWM,NWQVM),
     & IWQOBW(NBBWWM,NWQVM),IWQOBE(NBBEWM,NWQVM),
     * IWQOBN(NBBNWM,NWQVM),
     * IWQCBS(NBBSWM),IWQCBW(NBBWWM),IWQCBE(NBBEWM),IWQCBN(NBBNWM),
     * JWQCBS(NBBSWM),JWQCBW(NBBWWM),JWQCBE(NBBEWM),JWQCBN(NBBNWM)
     
  !     FLOAD=1.0-F_Load/100  ! reduce load
       FLOAD=1
C
C**********************************************************************C
C
      BSMALL=1.0E-6
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400./365.0
C
      DELT=DT2 ! Called at even timestep
      DELTA=DT2
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      ISTL=3
C
      IF (ISLTMT.GE.1) ISUD=1
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
       DO L=1,LC
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FVHU(L,K)=0.
        FUHV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        DU(L,K)=0.
        DV(L,K)=0.
       END DO
      END DO

C
C**********************************************************************C
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=1,LC
      CONT(L,K)=0.
      CMAX(L,K)=0.
      CMIN(L,K)=0.
      END DO
      END DO
 
C
C**********************************************************************C
C
C **  CALCULATED EXTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      MVAR=NW
!      CALL CALFQC (ISTL,MVAR,M,CON,CON1)
      IF(N.LE.6) THEN
       OPEN(10,FILE='pcb1.dia',ACCESS='APPEND')
       WRITE(10,*)N,TIME

      DO L=2,LA
        DO k=1,kc
        IF(FQC(L,K).GT.0.001) THEN
          WRITE(10,110)IL(L),JL(L),K,FQC(L,K)
        ENDIF
        ENDDO
      END DO
      CLOSE(10)
      ENDIF
  110 FORMAT(1X,3I4,2X,7E12.4,/,15X,7E12.4,/,15X,7E12.4)
     

C**********************************************************************C
C
C **  BEGIN COMBINED ADVECTION SCHEME
C
C**********************************************************************C
C
C **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
C **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED
C **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY
C
C----------------------------------------------------------------------C
C
      IF (ISCDCA(MVAR).EQ.0) GO TO 300
      IF (ISCDCA(MVAR).EQ.1) GO TO 400
      
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
 300  DO K=1,KC
      DO L=2,LA
      CONT(L,K)=CON1(L,K)
      END DO
      END DO
C
CWQJH     WFQC=0.5
C
C----------------------------------------------------------------------C
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS)
      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
      FUHU(L,K)=MAX(UHDYWQ(L,K),0.)*CONT(L-1,K)
     $         +MIN(UHDYWQ(L,K),0.)*CONT(L,K)
      FVHU(L,K)=MAX(VHDXWQ(L,K),0.)*CONT(LS,K)
     $         +MIN(VHDXWQ(L,K),0.)*CONT(L,K)
      END DO
      END DO

C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(WWQ(L,K),0.)*CONT(L,K)
     $        +MIN(WWQ(L,K),0.)*CONT(L,K+1)
      END DO
      END DO
C
      goto 500
C
  400 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=CON1(L,K)
      END DO
      END DO
C
CX     WFQC=1.0
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)    
      FUHU(L,K)=0.5*UHDYWQ(L,K)*(CON(L,K)+CON(L-1,K))
      FVHU(L,K)=0.5*VHDXWQ(L,K)*(CON(L,K)+CON(LS,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
C
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      IF (VHDXWQ(LN,K).LT.0.) FVHU(LN,K)=VHDXWQ(LN,K)*CON1(LN,K)
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      IF(UHDYWQ(L+1,K).LT.0.) FUHU(L+1,K)=UHDYWQ(L+1,K)*CON1(L+1,K)
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      IF(UHDYWQ(L,K).GT.0.) FUHU(L,K)=UHDYWQ(L,K)*CON1(L-1,K)
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS =LSC(L)
      IF (VHDXWQ(L,K).GT.0.) FVHU(L,K)=VHDXWQ(L,K)*CON1(LS,K)
      END DO
C
      END DO
C
C----------------------------------------------------------------------C
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=0.5*WWQ(L,K)*(CON(L,K+1)+CON(L,K))
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX FOR CENTRAL DIFF
C
      CALL CALHDMF
      CALL CALDIFF (3,M,CON1)
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX
C
C----------------------------------------------------------------------C
C
C     IF (ISTRAN(M).GE.1) CALL CALDIFF (ISTL,M,CON1)
C
C**********************************************************************C
C
C **  ADVECTION CALCULATION
C
C----------------------------------------------------------------------C
C
      IF (ISCDCA(MVAR).EQ.0) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LN)
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)      
      CH(L,K)=CONT(L,K)*H2WQ(L)
     $       +DELT*((FQC(L,K)*FLOAD+FUHU(L,K)-FUHU(L+1,K)
     $                       +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $                       +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO
C

      IF (ISFCT(5).GE.1) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)     
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON1(L,K)
       END DO
       END DO     
      END IF
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON1(L,K)
      CON(L,K)=min(CON(L,K),300.0)
      CONT(L,K)=0.0
      END DO
      END DO
C
      else
c 
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,RDZIC)     
      DO K=1,KC
      RDZIC=DZIC(K)
      DO L=2,LA
cx      CH(L,K)=CONT(L,K)*H2P(L)
cx     $       +DELT*( (FUHU(L,K)-FUHU(L+1,K)
      CH(L,K)=CON1(L,K)*H2WQ(L)
     $       +DELT*( ( FQC(L,K)*FLOAD+FUHU(L,K)-FUHU(L+1,K)
     $               +FVHU(L,K)-FVHU(LNC(L),K))*DXYIP(L)
     $               +(FWU(L,K-1)-FWU(L,K))*RDZIC )
      END DO
      END DO
C
      IF (ISFCT(MVAR).GE.1) THEN
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
       DO K=1,KC
       DO L=2,LA
       CON2(L,K)=CON(L,K)
       END DO
       END DO
      END IF
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)      
      DO K=1,KC
      DO L=2,LA
      CON1(L,K)=SCB(L)*CON(L,K)+(1.-SCB(L))*CON1(L,K)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON(L,K)
      CON(L,K)=max(CON(L,K),0.0)
      CONT(L,K)=0.0
      END DO
      END DO
C  
      endif
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES FOR WQVARS, M=8
C
!      IF(M.EQ.8) THEN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBS
      NSID=NCSERS(LL,M)
      L=LCBS(LL)
      LN=LNC(L)
C
      IF (VHDXWQ(LN,K).LT.0.) THEN
       CTMP=CON1(L,K)+DELT*(VHDXWQ(LN,K)*CON1(L,K)
     $      -FVHU(LN,K))*DXYIP(L)/HWQ(L)
       CON1(L,K)=CON(L,K)
       CON(L,K)=CTMP
       CBSTMP=CBS(LL,1,M)+CSERT(1,NSID,M)
       CLOS(LL,K,M)=CON(L,K)
       NLOS(LL,K,M)=N
      ELSE
       CBT= !WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+
     &   + CSERT(K,NSID,M)
       NMNLO=N-NLOS(LL,K,M)
       IF (NMNLO.GE.NTSCRS(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOS(LL,K,M)
     $         +(CBT-CLOS(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       END IF
      END IF
C
      END DO
      END DO
C
 6001 FORMAT('N,K,CBTS = ',2I10,F12.3)
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBW
      NSID=NCSERW(LL,M)
      L=LCBW(LL)      
C
      IF (UHDYWQ(L+1,K).LT.0.) THEN
       CTMP=CON1(L,K)
     $   +DELT*(UHDYWQ(L+1,K)*CON1(L,K)-FUHU(L+1,K))*DXYIP(L)/HWQ(L)
       CON1(L,K)=CON(L,K)
       CON(L,K)=CTMP
       CBWTMP=CBW(LL,1,M)+CSERT(1,NSID,M)
       CLOW(LL,K,M)=CON(L,K)
       NLOW(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLOW(LL,K,M)
       IF (NMNLO.GE.NTSCRW(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOW(LL,K,M)
     $         +(CBT-CLOW(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       END IF
      END IF
C
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      AATMP=1.0  !A_K**(TIME) implement reduce due to time (day)
      DO K=1,KC
      DO LL=1,NCBE
      NSID=NCSERE(LL,M)
      L=LCBE(LL)
	IF (UHDYWQ(L,K).GT.0.) THEN      
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     $      -UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON1(L,K)=CON(L,K)
       CON(L,K)=CTMP
       CBETMP=CBE(LL,1,M)+CSERT(1,NSID,M)
       CLOE(LL,K,M)=CON(L,K)
       NLOE(LL,K,M)=N
      ELSE
       if(NPCB.EQ.0) then
        CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
       else
        CBT=WTCI(K,1)*CBE(LL,1,M)+WTCI(K,2)*CBE(LL,2,M)+CSERT(K,NSID,M)
        CBT=CBT*AATMP   
       endif
       NMNLO=N-NLOE(LL,K,M)
       IF (NMNLO.GE.NTSCRE(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLOE(LL,K,M)
     $         +(CBT-CLOE(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       END IF
      END IF
C
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NCBN
      NSID=NCSERN(LL,M)
      L=LCBN(LL)
      LS=LSC(L)
C
      IF (VHDXWQ(L,K).GT.0.) THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     $      -VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON1(L,K)=CON(L,K)
       CON(L,K)=CTMP
       CBNTMP=CBN(LL,1,M)+CSERT(1,NSID,M)
       CLON(LL,K,M)=CON(L,K)
       NLON(LL,K,M)=N
      ELSE
       CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)
       NMNLO=N-NLON(LL,K,M)
       IF (NMNLO.GE.NTSCRN(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CLON(LL,K,M)
     $         +(CBT-CLON(LL,K,M))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       END IF
      END IF
C
      END DO
      END DO
C
 6002 FORMAT('N,K,CBTN = ',2I10,F12.3)
C----------------------------------------------------------------------C
C
!      END IF
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
      IF (ISADAC(5).EQ.0) RETURN
      IF (ISCDCA(5).EQ.1) RETURN
C
C----------------------------------------------------------------------C
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS)
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      UUU(L,K)=UWQ(L,K)*(CON(L,K)-CON(L-1,K))*DXIU(L)
      VVV(L,K)=VWQ(L,K)*(CON(L,K)-CON(LS,K))*DYIV(L)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS)
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      WWW(L,K)=WWQ(L,K)*(CON(L,K+1)-CON(L,K))*DZIG(K)/HWQ(L)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
c!$OMP DO PRIVATE(L,LS,LN,LNW,LSE,AUHU,AVHV,UTERM,VTERM,UHU,VHV) 
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)          
      LNW=LNWC(L)
      LSE=LSEC(L)      
      AUHU=ABS(UHDYWQ(L,K))
      AVHV=ABS(VHDXWQ(L,K))
      UTERM=AUHU*(1.-DELTA*AUHU/(DXYU(L)*HU(L)))*(CON(L,K)-CON(L-1,K))
      VTERM=AVHV*(1.-DELTA*AVHV/(DXYV(L)*HV(L)))*(CON(L,K)-CON(LS,K))
      UTERM=UTERM-0.25*DELTA*UHDYWQ(L,K)*
     $      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
     $      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
      VTERM=VTERM-0.25*DELTA*VHDXWQ(L,K)*
     $      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
     $      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
      UHU=UTERM/(CON(L,K)+CON(L-1,K)+BSMALL)
      VHV=VTERM/(CON(L,K)+CON(LS,K)+BSMALL)
      FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)
     $         +MIN(UHU,0.)*CON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)
     $         +MIN(VHV,0.)*CON(L,K)
      END DO
      END DO
C
c!$OMP DO PRIVATE(L,LN,AWW,WTERM,WW) 
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(WWQ(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)/HWQ(L))*(CON(L,K+1)-CON(L,K))
      WTERM=WTERM-0.25*DELTA*WWQ(L,K)*
     $      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
     $      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
      WW=WTERM/(CON(L,K+1)+CON(L,K)+BSMALL)
      FWU(L,K)=MAX(WW,0.)*CON(L,K)
     $        +MIN(WW,0.)*CON(L,K+1)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C

      DO K=1,KC
c
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      FVHU(LN,K)=0.0
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      FUHU(L+1,K)=0.0
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      FUHU(L,K)=0.0
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      FVHU(L,K)=0.0
      END DO
C
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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=MAX(CON(L,K),CON2(L,K))
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO L=2,LA
      CMAX(L,1)=MAX(CONT(L,1),CONT(L,2))
      CMAX(L,KC)=MAX(CONT(L,KS),CONT(L,KC))
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L) 
      DO K=2,KS
      DO L=2,LA
      CMAX(L,K)=MAX(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
      END DO
      END DO
C 
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMAX=SUB(L)*CONT(L-1,K) 
      CEMAX=SUB(L+1)*CONT(L+1,K)
      CSMAX=SVB(L)*CONT(LS,K)
      CNMAX=SVB(LN)*CONT(LN,K)
      CMAX(L,K)=MAX(CMAX(L,K),CNMAX,CEMAX,CSMAX,CWMAX)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CONT(L,K)=MIN(CON(L,K),CON2(L,K))
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO L=2,LA
      CMIN(L,1)=MIN(CONT(L,1),CONT(L,2))
      CMIN(L,KC)=MIN(CONT(L,KS),CONT(L,KC))
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=2,KS
      DO L=2,LA
      CMIN(L,K)=MIN(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
      END DO
      END DO
C
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      LN=LNC(L)
      CWMIN=SUB(L)*CONT(L-1,K)+1.E+6*(1.-SUB(L))
      CEMIN=SUB(L+1)*CONT(L+1,K)+1.E+6*(1.-SUB(L+1))
      CSMIN=SVB(L)*CONT(LS,K)+1.E+6*(1.-SVB(L))
      CNMIN=SVB(LN)*CONT(LN,K)+1.E+6*(1.-SVB(LN))
      CMIN(L,K)=MIN(CMIN(L,K),CNMIN,CEMIN,CSMIN,CWMIN)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SEPARATE POSITIVE AND NEGATIVE FLUXES PUTTING NEGATIVE FLUXES 
C **  INTO FUHV, FVHV, AND FWV
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      FUHV(L,K)=MIN(FUHU(L,K),0.)
      FUHU(L,K)=MAX(FUHU(L,K),0.)
      FVHV(L,K)=MIN(FVHU(L,K),0.)
      FVHU(L,K)=MAX(FVHU(L,K),0.)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KS
      DO L=2,LA
      FWV(L,K)=MIN(FWU(L,K),0.)
      FWU(L,K)=MAX(FWU(L,K),0.)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS AND LOAD
C **  INTO DU AND DV, THEN ADJUCT VALUES AT BOUNDARIES
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
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

      DO K=1,KC
c
      DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      DU(LN,K)=0.
      DV(LN,K)=0.
      END DO
C
      DO LL=1,NCBW
      L=LCBW(LL)
      DU(L+1,K)=0.
      DV(L+1,K)=0.
      END DO
C
      DO LL=1,NCBE
      L=LCBE(LL)
      DU(L-1,K)=0.
      DV(L-1,K)=0.
      END DO
C
      DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      DU(LS,K)=0.
      DV(LS,K)=0.
      END DO
C
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
          DU(LD,K)=0.
          DV(LD,K)=0.
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
          DU(LD,KD)=0.
          DV(LD,KD)=0.
        END IF
      END DO

C
C----------------------------------------------------------------------C
C
C **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DU AND DV
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      IF (DU(L,K).GT.0.) DU(L,K)=(CMAX(L,K)-CON(L,K))/(DU(L,K)+BSMALL)
      IF (DV(L,K).GT.0.) DV(L,K)=(CON(L,K)-CMIN(L,K))/(DV(L,K)+BSMALL)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS)
      DO K=1,KC
      DO L=2,LA
      LS=LSC(L)
      FUHU(L,K)=MIN(DV(L-1,K),DU(L,K))*FUHU(L,K)
     $         +MIN(DU(L-1,K),DV(L,K))*FUHV(L,K)
      FVHU(L,K)=MIN(DV(LS,K),DU(L,K))*FVHU(L,K)
     $         +MIN(DU(LS,K),DV(L,K))*FVHV(L,K)
      END DO
      END DO

C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KS
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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LN)
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      CH(L,K)=CON(L,K)*HWQ(L) 
     $       +DELT*((FUHU(L,K)-FUHU(L+1,K)
     $              +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $             +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON(L,K)
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
      DO K=1,KC
      DO L=2,LA
      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
      IF (CCMAX.GT.0.) THEN
       WRITE(6,6011)CON(L,K),CMAX(L,K),IL(L),JL(L),K
      END IF
      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
      IF (CCMIN.GT.0.) THEN
       WRITE(6,6012)CMIN(L,K),CON(L,K),IL(L),JL(L),K
      END IF
      END DO
      END DO
C
      END IF
C
 6010 FORMAT(1X,'FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT(1X,'CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT(1X,'CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
C
C----------------------------------------------------------------------C
C
      RETURN
      END
