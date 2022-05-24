C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE CALTRWQ(M,NW,CON,CON1)
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
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'wq.par'
      INCLUDE 'efdc.cmn'
!      INTEGER FUNCTION OMP_GET_NUM_THREADS          
      DIMENSION CON(LCM,KCM), CON1(LCM,KCM),ITRC(19)
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

C
C**********************************************************************C
C
  !    CALL OMP_SET_NUM_THREADs(4)
      BSMALL=1.0E-6 
C
      DELT=DT2
      DELTA=DT2
      DELTD2=DT
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      ISTL=3
      
      CMAC_A=-1000
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
      MVAR=8
! For water quality is IWQS<0 this is will bypass IWQS is set to zero on C92 as 0 for al full version  ! JS 2018
! loading is added to kinatic equation in wqske.for

      CALL CALFQC(ISTL,MVAR,NW,CON,CON1)
 !      if(NW.eq.3.or.NW.eq.14) then
 !       write(502,'(2I6,8f12.4)')L,NW,(FQC(LQS(8),5))
 !       endif

!      if(N.LT.5) then
!      do L=2,LA
!      if(FQC(L,KC)>0)
!     $   write(502,'(2I6,8f12.4)')L,NW,(FQC(L,K),K=1,KC)
!      end do
!      endif

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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
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
!      CMAC_A=max(CMAC_A,CONT(L,K)) 
      END DO
      END DO
!      write(*,*)'Conce ',M,IN_WQ, CMAX_A,CMAY_A,CMAC_A
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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)      
      CH(L,K)=CONT(L,K)*H2WQ(L)
     $       +DELT*(( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     $                       +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $                       +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO

C
      if(N<100.and.NW.eq.14)then
      write(991,*)FQC(LQS(8),KC)
      endif

      IF (ISFCT(M).GE.1) THEN
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
      CONT(L,K)=0.0
      END DO
      END DO
C
C
C**********************************************************************C
C
C **  CALCULATE LAST OUTFLOWING CONCENTRATION OR SPECIFY INFLOW 
C **  CONCENTRATION AT OPEN BOUNDARIES FOR WQVARS, M=8
C
      IF(M.EQ.8) THEN
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBS
      NSID=IWQOBS(LL,NW)
      L=LIJ( IWQCBS(LL),JWQCBS(LL) )
      LN=LNC(L)
C
      IF (VHDXWQ(LN,K).LT.0.) THEN
       CTMP=CON1(L,K)+DELT*(VHDXWQ(LN,K)*CON1(L,K)
     $      -FVHU(LN,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
       CWQLOS(LL,K,NW)=CON(L,K)        
       NWQLOS(LL,K,NW)=N
      ELSE
       CBT= ! WTCI(K,1)*WQOBCS(LL,1,NW)+WTCI(K,2)*WQOBCS(LL,2,NW)
     $    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOS(LL,K,NW)
       IF (NMNLO.GE.NTSCRS(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOS(LL,K,NW)
     $         +(CBT-CWQLOS(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRS(LL))
       END IF
      END IF
C
      END DO
      END DO

C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBW
      NSID=IWQOBW(LL,NW)
      L=LIJ( IWQCBW(LL),JWQCBW(LL) )
C
      IF (UHDYWQ(L+1,K).LT.0.) THEN
       CTMP=CON1(L,K)+DELT*(UHDYWQ(L+1,K)*CON1(L,K)
     $      -FUHU(L+1,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBW(LL,1,M)) CON(L,K)=CBW(LL,1,M)
       CWQLOW(LL,K,NW)=CON(L,K)
       NWQLOW(LL,K,NW)=N
      ELSE
       CBT= ! WTCI(K,1)*WQOBCW(LL,1,NW)+WTCI(K,2)*WQOBCW(LL,2,NW)
     $    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOW(LL,K,NW)
       IF (NMNLO.GE.NTSCRW(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOW(LL,K,NW)
     $         +(CBT-CWQLOW(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRW(LL))
       END IF
      END IF
C
      END DO
      END DO

C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NWQOBE
      NSID=IWQOBE(LL,NW)
      L=LIJ( IWQCBE(LL),JWQCBE(LL) )
C
      IF (UHDYWQ(L,K).GT.0.) THEN
       CTMP=CON1(L,K)+DELT*(FUHU(L,K)
     $      -UHDYWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBE(LL,1,M)) CON(L,K)=CBE(LL,1,M)
       CWQLOE(LL,K,NW)=CON(L,K)
       NWQLOE(LL,K,NW)=N
      ELSE
       CBT= ! WTCI(K,1)*WQOBCE(LL,1,NW)+WTCI(K,2)*WQOBCE(LL,2,NW)
     $    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLOE(LL,K,NW)
       IF (NMNLO.GE.NTSCRE(LL)) THEN
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLOE(LL,K,NW)
     $         +(CBT-CWQLOE(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRE(LL))
       END IF
      END IF
C
      END DO
      END DO

C
C----------------------------------------------------------------------C
C

      DO K=1,KC
      DO LL=1,NWQOBN
      NSID=IWQOBN(LL,NW)
      L=LIJ( IWQCBN(LL),JWQCBN(LL) )
      LS=LSC(L)
C
      IF (VHDXWQ(L,K).GT.0.) THEN
       CTMP=CON1(L,K)+DELT*(FVHU(L,K)
     $      -VHDXWQ(L,K)*CON1(L,K))*DXYIP(L)/HWQ(L)
       CON(L,K)=CTMP
C      IF (CON(L,K).GT.CBN(LL,1,M)) CON(L,K)=CBN(LL,1,M)
       CWQLON(LL,K,NW)=CON(L,K)
       NWQLON(LL,K,NW)=N
      ELSE
       CBT= ! WTCI(K,1)*WQOBCN(LL,1,NW)+WTCI(K,2)*WQOBCN(LL,2,NW)
     $    +CSERTWQ(K,NSID,NW)
       NMNLO=N-NWQLON(LL,K,NW)
       IF (NMNLO.GE.NTSCRN(LL)) THEN 
        CON(L,K)=CBT
       ELSE
        CON(L,K)=CWQLON(LL,K,NW)
     $         +(CBT-CWQLON(LL,K,NW))*FLOAT(NMNLO)/FLOAT(NTSCRN(LL))
       END IF
      END IF
C
      END DO
      END DO

C
C----------------------------------------------------------------------C
C
      END IF
C
C**********************************************************************C
C
C **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION 
C
      IF (ISADAC(M).EQ.0) RETURN
      IF (ISCDCA(M).EQ.1) RETURN
      
!       if(ISADAC(M).EQ.3.and.NW.NE.19) RETURN
      if(NW.NE.19) RETURN
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
!$OMP  PARALLEL DO DEFAULT(SHARED)
!$OMP& PRIVATE(L,LS,LNW,LSE,AUHU,AVHV,UTERM,VTERM,UTERM1,VTERM1,UHU,VHV)   
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
      UTERM1=UTERM-0.25*DELTA*UHDYWQ(L,K)*
     $      (VVV(L,K)+VVV(LN,K)+VVV(LNW,K)+VVV(L-1,K)
     $      +WWW(L,K)+WWW(L-1,K)+WWW(L-1,K-1)+WWW(L,K-1))
      VTERM1=VTERM-0.25*DELTA*VHDXWQ(L,K)*
     $      (UUU(L,K)+UUU(LS,K)+UUU(LSE,K)+UUU(L+1,K)
     $      +WWW(L,K)+WWW(LS,K)+WWW(LS,K-1)+WWW(L,K-1))
      UHU=UTERM1/(CON(L,K)+CON(L-1,K)+BSMALL)
      VHV=VTERM1/(CON(L,K)+CON(LS,K)+BSMALL)
      FUHU(L,K)=MAX(UHU,0.)*CON(L-1,K)
     $         +MIN(UHU,0.)*CON(L,K)
      FVHU(L,K)=MAX(VHV,0.)*CON(LS,K)
     $         +MIN(VHV,0.)*CON(L,K)
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LN,AWW,WTERM,WW,WTERM1) 
      DO K=1,KS
      DO L=2,LA
      LN=LNC(L)
      AWW=ABS(WWQ(L,K))
      WTERM=AWW*(1.-DELTA*AWW*DZIG(K)/HWQ(L))*(CON(L,K+1)-CON(L,K))
      WTERM1=WTERM-0.25*DELTA*WWQ(L,K)*
     $      (UUU(L,K)+UUU(L+1,K)+UUU(L+1,K+1)+UUU(L,K+1)
     $      +VVV(L,K)+VVV(LN,K)+VVV(LN,K+1)+VVV(L,K+1))
      WW=WTERM1/(CON(L,K+1)+CON(L,K)+BSMALL)
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
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS,LN,CWMAX,CEMAX,CSMAX,CNMAX)  
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

      DO K=2,KS
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)      
      DO L=2,LA
      CMIN(L,K)=MIN(CONT(L,K-1),CONT(L,K),CONT(L,K+1))
      END DO
      END DO
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L,LS,LN,CWMIN,CEMIN,CSMIN,CNMIN)  
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
!$OMP&         PRIVATE(L,LN)
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
