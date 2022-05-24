C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALEXP (ISTL)
C
C **  SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS
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
      FCAXE(1)=0.
      FCAYE(1)=0.
      FCAX1E(1)=0.
      FCAY1E(1)=0.
      FXE(1)=0.
      FYE(1)=0.
C
      FCAXE(LC)=0.
      FCAYE(LC)=0.
      FCAX1E(LC)=0.
      FCAY1E(LC)=0.
      FXE(LC)=0.
      FYE(LC)=0.
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L) 
       DO L=2,LA
        FCAXE(L)=0.
        FCAYE(L)=0.
        FCAX1E(L)=0.
        FCAY1E(L)=0.
        FXE(L)=0.
        FYE(L)=0.
       END DO
!$OMP END DO 
!$OMP END PARALLEL      
C
C**********************************************************************C
C
C **  SELECT ADVECTIVE FLUX FORM
C
C----------------------------------------------------------------------C
C
      IF (ISTL.EQ.2) GO TO 200
      IF (ISTL.EQ.3) THEN
       IF (ISCDMA.EQ.0) GO TO 300
       IF (ISCDMA.EQ.1) GO TO 400
       IF (ISCDMA.EQ.2) GO TO 350
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
!$OMP PARALLEL
!$OMP DO PRIVATE(L) 
      DO K=1,KC
       DO L=2,LA
        UHDY2(L,K)=UHDY1(L,K)+UHDY(L,K)
        VHDX2(L,K)=VHDX1(L,K)+VHDX(L,K)
        W2(L,K)=W1(L,K)+W(L,K)
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL 
C
C----------------------------------------------------------------------C
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN,LS,UHC, UHB,VHC,VHB) 
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)
       UHC=0.25*(UHDY2(L,K)+UHDY2(LS,K))
       UHB=0.25*(UHDY2(L,K)+UHDY2(L+1,K))
       VHC=0.25*(VHDX2(L,K)+VHDX2(L-1,K))
       VHB=0.25*(VHDX2(L,K)+VHDX2(LN,K))
C
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
!$OMP END DO 
!$OMP END PARALLEL 
C
C----------------------------------------------------------------------C
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LS,WU,WV) 
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W2(L,K)+W2(L-1,K))
       WV=0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)
     $        +MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)
     $        +MIN(WV,0.)*V1(L,K+1)
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL 
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
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN,LS,UHC,UHB,VHC,VHB) 
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
!$OMP END DO 
!$OMP END PARALLEL 
C
C----------------------------------------------------------------------C
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LS,WU,WV)
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.5*DXYU(L)*(W(L,K)+W(L-1,K))
       WV=0.5*DXYV(L)*(W(L,K)+W(LS,K))
       FWU(L,K)=MAX(WU,0.)*U1(L,K)
     $        +MIN(WU,0.)*U1(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V1(L,K)
     $        +MIN(WV,0.)*V1(L,K+1)
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL       
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
      DO K=1,KS
       DO L=2,LA
       LS=LSC(L)
       WU=0.25*DXYU(L)*(W(L,K)+W(L-1,K))
       WV=0.25*DXYV(L)*(W(L,K)+W(LS,K))
       FWU(L,K)=MAX(WU,0.)*U2(L,K)
     $        +MIN(WU,0.)*U2(L,K+1)
       FWV(L,K)=MAX(WV,0.)*V2(L,K)
     $        +MIN(WV,0.)*V2(L,K+1)
       END DO
      END DO
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
      DO K=1,KS
      DO L=2,LA
      LS=LSC(L)
      FWU(L,K)=0.25*DXYU(L)*(W(L,K)+W(L-1,K))*(U(L,K+1)+U(L,K))
      FWV(L,K)=0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
      END DO
      END DO
C
C**********************************************************************C
C
  500 CONTINUE
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L)
      DO K=1,KC
       DO L=1,LA
       FUHU(L,K)=STCUV(L)*FUHU(L,K)
       FVHV(L,K)=STCUV(L)*FVHV(L,K)
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL      
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
C
C----------------------------------------------------------------------C
C
      IF(ISDCCA.EQ.0) THEN
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN)
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       CAC(L,K)=( FCORC(L)*DXYP(L)
     $        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     $        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL    
C
      ELSE
C
      CFMAX=CF
C
cC!$OMP PARALLEL
cC!$OMP DO PRIVATE(L,LN,CFEFF,CFMAX)
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       CAC(L,K)=( FCORC(L)*DXYP(L)
     $        +0.5*SNLT*(V(LN,K)+V(L,K))*(DYU(L+1)-DYU(L))
     $        -0.5*SNLT*(U(L+1,K)+U(L,K))*(DXV(LN)-DXV(L)) )*HP(L)
       CFEFF=ABS(CAC(L,K))*DXYIP(L)*HPI(L)
       CFMAX=MAX(CFMAX,CFEFF)
       END DO
      END DO
cC!$OMP END DO 
cC!$OMP END PARALLEL 
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L)
      DO K=1,KC
       DO L=1,LC
        FCAX1(L,K)=0.
        FCAY1(L,K)=0.
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL 
C
      IF(ISTL.EQ.3) THEN
C      IF(IRVEC.EQ.10.OR.IRVEC.EQ.11) THEN
      IF(IRVEC.EQ.14) THEN
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L)
      DO L=1,LC
       TVAR3E(L)=0.
       TVAR3N(L)=0.
       TVAR3W(L)=0.
       TVAR3S(L)=0.
	   TVAR1S(L,1)=0.
	   TVAR2S(L,1)=0.
      END DO
!$OMP END DO      
C
!$OMP DO PRIVATE(L,DZTMP)
      DO K=1,KC
      DZTMP=DZC(K)
      DO L=2,LA
       TVAR3E(L)=TVAR3E(L)+DZTMP*U(L,K)
       TVAR3N(L)=TVAR3N(L)+DZTMP*V(L,K)
       TVAR3W(L)=TVAR3W(L)+DZTMP*U1(L,K)
       TVAR3S(L)=TVAR3S(L)+DZTMP*V1(L,K)
      END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL

C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN)
      DO L=2,LA
       LN=LNC(L)
       TVAR1S(L,1)=( FCORC(L)*DXYP(L)
     $      +0.5*SNLT*(TVAR3N(LN)+TVAR3N(L))*(DYU(L+1)-DYU(L))
     $      -0.5*SNLT*(TVAR3E(L+1)+TVAR3E(L))*(DXV(LN)-DXV(L)) )*HP(L)
       TVAR2S(L,1)=( FCORC(L)*DXYP(L)
     $      +0.5*SNLT*(TVAR3S(LN)+TVAR3S(L))*(DYU(L+1)-DYU(L))
     $      -0.5*SNLT*(TVAR3W(L+1)+TVAR3W(L))*(DXV(LN)-DXV(L)) )*H1P(L)
      END DO
!$OMP END DO 
!$OMP END PARALLEL      
C
C       FCAX(L,K)=ROLD*FCAX(L,K)
C     $          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
C     $                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
C       FCAY(L,K)=ROLD*FCAY(L,K)
C     $          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
C     $                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
C
C
C       TVAR1S=CAC @ N
C       TVAR2S=CAC @ N-1
C
C       TVAR3E=U @ N
C       TVAR3N=V @ N
C       TVAR3W=U @ N-1
C       TVAR3S=V @ N-1
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN,LS,LNW,LSE)
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=0.0625*SCAX(L)*
     $            ( TVAR2S(L,1)*(TVAR3S(LN)+TVAR3S(L))
     $             +TVAR2S(L-1,1)*(TVAR3S(LNW)+TVAR3S(L-1)) )
       FCAY1E(L)=0.0625*SCAY(L)*
     $            ( TVAR2S(L,1)*(TVAR3W(L+1)+TVAR3W(L))
     $             +TVAR2S(LS,1)*(TVAR3W(LSE)+TVAR3W(LS)) )
      END DO
!$OMP END DO 
C
!$OMP DO PRIVATE(L,LN,LS,LNW,LSE)
      DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX1E(L)=FCAX1E(L)-0.125*SCAX(L)*
     $            ( TVAR1S(L,1)*(TVAR3N(LN)+TVAR3N(L))
     $             +TVAR1S(L-1,1)*(TVAR3N(LNW)+TVAR3N(L-1)) )
       FCAY1E(L)=FCAY1E(L)-0.125*SCAY(L)*
     $            ( TVAR1S(L,1)*(TVAR3E(L+1)+TVAR3E(L))
     $             +TVAR1S(LS,1)*(TVAR3E(LSE)+TVAR3E(LS)) )
      END DO
!$OMP END DO 
!$OMP END PARALLEL      
C
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN,LS,LNW,LSE)
      DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
       LS=LSC(L)
       LNW=LNWC(L)
       LSE=LSEC(L)
       FCAX(L,K)=ROLD*FCAX(L,K)
     $          +0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K))
     $                             +CAC(L-1,K)*(V(LNW,K)+V(L-1,K)))
       FCAY(L,K)=ROLD*FCAY(L,K)
     $          +0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(L+1,K)+U(L,K))
     $                             +CAC(LS,K)*(U(LSE,K)+U(LS,K)))
       FX(L,K)=SAAX(L)*(FUHU(L,K)-FUHU(L-1,K)+FVHU(LN,K)-FVHU(L,K))
       FY(L,K)=SAAY(L)*(FUHV(L+1,K)-FUHV(L,K)+FVHV(L,K)-FVHV(LS,K))
       END DO
      END DO
!$OMP END DO 
!$OMP END PARALLEL        
C
C**********************************************************************C
C
C **  ADD HORIZONTAL MOMENTUN DIFFUSION TO ADVECTIVE ACCELERATIONS
C
C----------------------------------------------------------------------C
C
      IF (ISHDMF.EQ.1) THEN
C
!$OMP PARALLEL
!$OMP DO PRIVATE(L,LN,LS)
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
!$OMP END DO 
!$OMP END PARALLEL 
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE EXTERNAL ACCELERATIONS
C
C----------------------------------------------------------------------C
C
c!$OMP PARALLEL
c!$OMP DO PRIVATE(L)
      DO K=1,KC
       DO L=2,LA
        FCAXE(L)=FCAXE(L)+FCAX(L,K)*DZC(K)
        FCAYE(L)=FCAYE(L)+FCAY(L,K)*DZC(K)
        FXE(L)=FXE(L)+FX(L,K)*DZC(K)
        FYE(L)=FYE(L)+FY(L,K)*DZC(K)
       END DO
      END DO
c!$OMP END DO 
c!$OMP END PARALLEL 
C
C----------------------------------------------------------------------C
C
C **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
C
      IF (ISWAVE.GE.1) THEN
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
c!$OMP PARALLEL
c!$OMP DO PRIVATE(L)
      DO K=1,KC
       DO L=2,LA
         FX(L,K)=FX(L,K)+SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(K)
         FY(L,K)=FY(L,K)+SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(K)
        END DO
      END DO
c!$OMP END DO 
c!$OMP END PARALLEL 
C
C----------------------------------------------------------------------C
C
C **  ADD NET WAVE REYNOLDS STRESSES TO INTERNAL ADVECTIVE ACCEL.
C
      IF (ISWAVE.GE.1) THEN
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
c!$OMP PARALLEL
c!$OMP DO PRIVATE(L,LS)
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
c!$OMP END DO 
c!$OMP END PARALLEL 
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
!$OMP PARALLEL
!$OMP DO PRIVATE(L,RCDZF)
      DO K=1,KS
       RCDZF=CDZF(K)
       DO L=2,LA
        DU(L,K)=RCDZF*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI
     $           +DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)+FBBX(L,K)
     $           +SNLT*(FX(L,K)-FX(L,K+1))) )
        DV(L,K)=RCDZF*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI
     $           +DXYIV(L)*(FCAY(L,K)-FCAY(L,K+1)+FBBY(L,K)
     $           +SNLT*(FY(L,K)-FY(L,K+1))) )
       END DO
      END DO
!$OMP END DO
!$OMP END PARALLEL     
C
      IF(ISTL.EQ.2) THEN
C

      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        DU(L,KS)=DU(L,KS)-CDZU(KS)*TSX(L)
        DV(L,KS)=DV(L,KS)-CDZU(KS)*TSY(L)
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
      RETURN
      END
