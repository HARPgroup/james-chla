C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C 
      SUBROUTINE CALTRWQH(M,NW,CON,CON1)
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
       END DO
      END DO

C
C**********************************************************************C
C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
!      DO K=1,KC
!      DO L=1,LC
!      CONT(L,K)=0.
!      CMAX(L,K)=0.
!      CMIN(L,K)=0.
!      END DO
!      END DO
    
C
C**********************************************************************C
C
C **  CALCULATED EXTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      IN_WQ=NW
      MVAR=8
 
      CALL CALFQC(ISTL,MVAR,NW,CON,CON1)

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
!      DO K=1,KC
!      DO L=2,LA
!      CONT(L,K)=CON1(L,K)
!      END DO
!      END DO
      CMAX_A=-1000
      CMAY_A=-1000
      CMAC_A=-1000
C 
C -X-direction
C
      DO K=1,KC
      DO L=3,LA
        IF (UWQ(L,K).GE.0.) THEN
         DELTM32=CON1(L-1,K)-CON1(L-2,K)
         DELTM12=CON1(L,K)-CON1(L-1,K)
         IF(abs(DELTM12).LT.1.0E-4) THEN
         RL=0
         ELSE
         RL=DELTM32/(DELTM12)
         RKAL=CKAX(L-1,K)/CKAX(L,K)
         ENDIF
         BETAL=BETAX1(L,K)+BETAX2(L,K)*RL
         SL=CON1(L-1,K)+0.5*AMAX1(0.,AMIN1(2.,2*RL*RKAL,BETAL))
     *   *DELTM12*max(0.,CKAX(L,K))  
        
         IF (SUB(L-1).EQ.0) SL=CON1(L-1,K)  ! Reset CON1(L-2) to CON1(L-1) at cell next to boundary
              
         ELSE
         
         DELTM12=CON1(L,K)-CON1(L-1,K)
         DELTP12=CON1(L+1,K)-CON1(L,K) 
         IF( abs(DELTM12).LT.1.0E-4)THEN
          RTR=0
         ELSE
          RTR=DELTP12/DELTM12
          RKAR=CKAX(L+1,K)/CKAX(L,K)         
         ENDIF 
         BETAR=BETAX1(L,K)+BETAX2(L,K)*RTR
         SR=CON1(L,K)-0.5*AMAX1(0.,AMIN1(2.,2*RTR*RKAR,BETAR))
     *   *DELTM12*max(0.,CKAX(L,K))                
         ENDIF
         IF (SUB(L+1).EQ.0) SR=CON1(L,K)  ! Reset CON1
         
         FUHU(L,K)=MAX(UHDYWQ(L,K),0.)*SL
     $         +MIN(UHDYWQ(L,K),0.)*SR    
  
!         CMAX_A=max(CMAX_A,SL)
!         CMAY_A=max(CMAY_A,SR)  
!         CMAC_A=max(CMAC_A,CON1(L,K))   
       ENDDO
       FUHU(2,K)=0
       ENDDO
C
C Y-direction
C  
       DO K=1,KC
       DO L=2,LA
         LS=LSC(L)     ! J-1
         LSS=LSSC(L)   ! J-2
         LN=LNC(L)     ! J+1
         LNN=LNNC(L)   ! J+2   
    
         IF (VWQ(L,K).GT.0.) THEN
         DELTM32=(CON1(LS,K)-CON1(LSS,K)) !*SVB(LS)
         DELTM12=(CON1(L,K)-CON1(LS,K))*SVB(L)
         IF(abs(DELTM12).LT.1.0E-4.OR.LSS.EQ.LC) THEN
         RL=0
         ELSE
         RL=DELTM32/(DELTM12)
         RKAL=CKAY(LS,K)/CKAY(L,K)
         ENDIF
         BETAL=BETAY1(L,K)+BETAY2(L,K)*RL
         SL=CON1(LS,K)+0.5*AMAX1(0.,AMIN1(2.,2*RL*RKAL,BETAL))
     *   *DELTM12*max(0.,CKAY(L,K)) 
        
         IF (SVB(LS).EQ.0) SL=CON1(LS,K)  ! Reset CON1(L-2) to CON1(L-1) at cell next to boundary
              
         ELSE

         DELTM12=(CON1(L,K)-CON1(LS,K)) !*SVB(L)
         DELTP12=(CON1(LN,K)-CON1(L,K))*SVB(LN) 
         IF( abs(DELTM12).LT.1.0E-4)THEN
          RTR=0
         ELSE
          RTR=DELTP12/DELTM12
          RKAR=CKAY(LN,K)/CKAY(L,K)         
         ENDIF 
         BETAR=BETAY1(L,K)+BETAY2(L,K)*RTR
         SR=CON1(L,K)-0.5*AMAX1(0.,AMIN1(2.,2*RTR*RKAR,BETAR))
     *   *DELTM12*max(0.,CKAY(L,K))
                   
          IF (SVB(LN).EQ.0) SR=CON1(L,K)  ! Reset CON1

         ENDIF   
               
         FVHU(L,K)=MAX(VHDXWQ(L,K),0.)*SL
     $         +MIN(VHDXWQ(L,K),0.)*SR 
 !         CMAY_A=max(CMAY_A,SL)
 !         CMAY_A=max(CMAY_A,SR)  
       ENDDO
       ENDDO
 !      write(*,*)'Max C ',M,IN_WQ, CMAX_A,CMAY_A,CMAC_A
C
C--Z- direction
C
      DO L=2,LA
      IF(WWQ(L,1).GE.0) THEN
      SL=CON1(L,1)
      ELSE
       DELTM12=CON1(L,2)-CON1(L,1)
         DELTP12=CON1(L,3)-CON1(L,2) 
         IF( abs(DELTM12).LT.1.0E-4)THEN
          RTR=0
         ELSE
          RTR=DELTP12/DELTM12
          RKAR=CKAW(L,2)/CKAX(L,1)         
         ENDIF 
         BETAR=BETAZ1(L,1)+BETAZ2(L,1)*RTR
         SR=CON1(L,1)-0.5*AMAX1(0.,AMIN1(2.,2*RTR*RKAR,BETAR))
     *   *DELTM12*max(0.,CKAW(L,1))                   
      ENDIF
      ENDDO
      
      DO K=2,KC-2
      IF(WWQ(L,K).GE.0) THEN
         DELTM32=CON1(L,K)-CON1(L,K-1)
         DELTM12=CON1(L,K+1)-CON1(L,K)
         IF(abs(DELTM12).LT.1.0E-4) THEN
         RL=0
         ELSE
         RL=DELTM32/(DELTM12)
         RKAL=CKAW(L,K-1)/CKAW(L,K)
         ENDIF
         BETAL=BETAZ1(L,K)+BETAZ2(L,K)*RL
         SL=CON1(L,K)+0.5*AMAX1(0.,AMIN1(2.,2*RL*RKAL,BETAL))
     *   *DELTM12*max(0.,CKAW(L,K))  
      ELSE
         DELTM12=CON1(L,K+1)-CON1(L,K)
         DELTP12=CON1(L,K+2)-CON1(L,K+1) 
         IF( abs(DELTM12).LT.1.0E-4)THEN
          RTR=0
         ELSE
          RTR=DELTP12/DELTM12
          RKAR=CKAW(L,K+1)/CKAW(L,K)         
         ENDIF 
         BETAR=BETAZ1(L,K)+BETAZ2(L,K)*RTR
         SR=CON1(L,K+1)-0.5*AMAX1(0.,AMIN1(2.,2*RTR*RKAR,BETAR))
     *   *DELTM12*max(0.,CKAW(L,K))                
       ENDIF
      ENDDO

!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KS
      DO L=2,LA
      FWU(L,K)=MAX(WWQ(L,K),0.)*CON1(L,K)
     $        +MIN(WWQ(L,K),0.)*CON1(L,K+1)
      END DO
      END DO

C  
  



      DO K=1,KC    
      DO L=2,LA
      LS=LSC(L)      
  !    FUHU(L,K)=MAX(UHDYWQ(L,K),0.)*CON1(L-1,K)
  !   $         +MIN(UHDYWQ(L,K),0.)*CON1(L,K)
  !    FVHU(L,K)=MAX(VHDXWQ(L,K),0.)*CON1(LS,K)
  !   $         +MIN(VHDXWQ(L,K),0.)*CON1(L,K)
     
 !     CMAC_A=max(CMAC_A,CON1(L,K)) 
      END DO
      END DO
!      write(*,*)'Conce ',M,IN_WQ, CMAX_A,CMAY_A,CMAC_A
c      DO K=1,KS
c      DO L=2,LA
c      FWU(L,K)=MAX(WWQ(L,K),0.)*CON1(L,K)
c     $        +MIN(WWQ(L,K),0.)*CON1(L,K+1)
c      END DO
c      END DO

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
      CH(L,K)=CON1(L,K)*H2WQ(L)
     $       +DELT*(( FQC(L,K)+FUHU(L,K)-FUHU(L+1,K)
     $                       +FVHU(L,K)-FVHU(LN,K))*DXYIP(L)
     $                       +(FWU(L,K-1)-FWU(L,K))*DZIC(K))
      END DO
      END DO

C
!$OMP          PARALLEL DO DEFAULT(SHARED)
!$OMP&         PRIVATE(L)
      DO K=1,KC
      DO L=2,LA
      CON(L,K)=SCB(L)*CH(L,K)/HWQ(L)+(1.-SCB(L))*CON1(L,K)
!      CONT(L,K)=0.0
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
C----------------------------------------------------------------------C
C
C **  DIAGNOSE FCT SCHEME
C
!      IF (ISFCT(M).EQ.2) THEN
!      WRITE(6,6010)N
! C
!      DO K=1,KC
!      DO L=2,LA
!      CCMAX=SCB(L)*(CON(L,K)-CMAX(L,K))
!      IF (CCMAX.GT.0.) THEN
!       WRITE(6,6011)CON(L,K),CMAX(L,K),IL(L),JL(L),K
!      END IF
!      CCMIN=SCB(L)*(CMIN(L,K)-CON(L,K))
!      IF (CCMIN.GT.0.) THEN
!       WRITE(6,6012)CMIN(L,K),CON(L,K),IL(L),JL(L),K
!      END IF
!      END DO
!      END DO
C
!      END IF
C
 6010 FORMAT(1X,'FCT DIAGNOSTICS AT N = ',I5)
 6011 FORMAT(1X,'CON = ',E12.4,3X,'CMAX = ',E12.4,3X,'I,J,K=',(3I10))
 6012 FORMAT(1X,'CMIN = ',E12.4,3X,'CON = ',E12.4,3X,'I,J,K=',(3I10))
C
C----------------------------------------------------------------------C
C
      RETURN 
      END
