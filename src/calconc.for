C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALCONC (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND
C **  SUSPENDED CONSTITUTENTS, INCLUDING SALINITY, TEMPERATURE, DYE AND
C **  AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL
C **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION EEB(LCM),CCLBTMP(LCM)
C
CDHP  DIMENSION CUBTMP(LCM),CMBTMP(LCM),CLBTMP(LCM),EB(LCM),VTMP(LCM),
CDHP $          ABHPI(LCM,KSM)
C
C**********************************************************************C
C
      NTOX=NPCB
      
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      IF (ISTL.EQ.2) THEN
       DELT=DT
       S3TL=0.0
       S2TL=1.0
      END IF
      DELTD2=DELT
      
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% JS 1/31/2014      

C     DELTD2=0.5*DELT
C
C**********************************************************************C
C
C **  OPTIONAL MASS BALANCE CALCULATION
C
C----------------------------------------------------------------------C
C
C     IF (NCTBC.NE.NTSTBC.AND.ISBAL.GE.1) THEN
      IF (ISTL.NE.2.AND.ISBAL.GE.1) THEN
         CALL CALBAL2
         CALL CALBAL3
         NTMP=MOD(N,2)
         IF(NTMP.EQ.0) THEN
           CALL CBALEV2
           CALL CBALEV3
          ELSE
           CALL CBALOD2
           CALL CBALOD3
         END IF
      END IF
C
C**********************************************************************C
C
C **  SEDIMENT BUDGET CALCULATION    (DLK 10/15)
C
C----------------------------------------------------------------------C
C
C      IF (NCTBC.NE.NTSTBC) THEN
      IF (ISTL.NE.2.AND.ISSBAL.GE.1) THEN
         CALL BUDGET2
         CALL BUDGET3
C         NTMP=MOD(N,2)
C         IF(NTMP.EQ.0) THEN
C           CALL BUDGEV2
C           CALL BUDGEV3
C          ELSE
C           CALL BUDGOD2
C           CALL BUDGOD3
C         END IF
      END IF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION EXPLICIT HALF STEP CALCULATION
C
C----------------------------------------------------------------------C
C
C      IF(KC.EQ.1) GO TO 500
C
C      TTMP=SECNDS(0.0)
C
C      K=1
C      RCDZKK=-DELTD2*CDZKK(1)
C      DO L=2,LA
C      CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
C      CCMBTMP=1.-CCUBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCUBTMP*(SAL1(L,K+1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCUBTMP*(TEM1(L,K+1)
C      DU(L,K) =CCMBTMP*DYE1(L,K)+CCUBTMP*(DYE1(L,K+1)
C      DV(L,K) =CCMBTMP*SED1(L,K)+CCUBTMP*(SED1(L,K+1)
C      END DO
C
C      DO K=2,KS
C      RCDZKMK=-DELTD2*CDZKMK(K)
C      RCDZKK=-DELTD2*CDZKK(K)
C      DO L=2,LA
C      CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)
C      CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
C      CCMBTMP=1.-CCUBTMP-CCLBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCUBTMP*SAL1(L,K+1)
C     $                          +CCLBTMP*SAL1(L,K-1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCUBTMP*TEM1(L,K+1)
C     $                          +CCLBTMP*TEM1(L,K-1)
C      DU(L,K) =CCMBTMP*DYE1(L,K)+CCUBTMP*DYE1(L,K+1)
C     $                          +CCLBTMP*DYE1(L,K-1)
C      DV(L,K) =CCMBTMP*SED1(L,K)+CCUBTMP*SED1(L,K+1)
C     $                          +CCLBTMP*SED1(L,K-1)
C      END DO
C      END DO
C
C      K=KC
C      RCDZKMK=-DELTD2*CDZKMK(K)
C      DO L=2,LA
C      CCLBTMP=RCDZKMK*HPI(L)*AB(L,K-1)
C      CCMBTMP=1.-CCLBTMP
C      UUU(L,K)=CCMBTMP*SAL1(L,K)+CCLBTMP*SAL1(L,K-1)
C      VVV(L,K)=CCMBTMP*TEM1(L,K)+CCLBTMP*TEM1(L,K-1)
C      DU(L,K)=CCMBTMP*DYE1(L,K)+CCLBTMP*DYE1(L,K-1)
C      DV(L,K)=CCMBTMP*SED1(L,K)+CCLBTMP*SED1(L,K-1)
C      END DO
C
C      TVDIF=TVDIF+SECNDS(TTMP)
C
C**********************************************************************C
C
  500 CONTINUE
C
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.EQ.0) THEN
C
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
      IF (ISTRAN(1).EQ.1) CALL CALTRAN (ISTL,1,1,SAL,SAL1)
      IF (ISTRAN(2).EQ.1) CALL CALTRAN (ISTL,2,2,TEM,TEM1)
      IF (ISTRAN(3).EQ.1) CALL CALTRAN (ISTL,3,3,DYE,DYE1)
C     IF (ISTRAN(3).EQ.1) CALL CALTRAN (ISTL,4,4,SFL,SFL1)
      IF (ISTRAN(5).EQ.1.and.NPCB.EQ.0) THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         END DO
         END DO
         CALL CALTRAN (ISTL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF   
            
      IF (ISTRAN(6).EQ.1) THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         END DO
         END DO
         CALL CALTRAN (ISTL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF
      IF (ISTRAN(7).EQ.1) THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         END DO
         END DO
         CALL CALTRAN (ISTL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF
C     IF (ISTRAN(1).EQ.1) CALL CALTRAN (ISTL,1,SAL,UUU)
C     IF (ISTRAN(2).EQ.1) CALL CALTRAN (ISTL,2,TEM,VVV)
C     IF (ISTRAN(3).EQ.1) CALL CALTRAN (ISTL,3,DYE,DU)
C     IF (ISTRAN(4).EQ.1) CALL CALTRAN (ISTL,4,SED,DV)
C
C     IF (ISTRAN(1).EQ.2) CALL CALTRANI (ISTL,1,SAL,SAL1)
C     IF (ISTRAN(2).EQ.2) CALL CALTRANI (ISTL,2,TEM,TEM1)
C     IF (ISTRAN(3).EQ.2) CALL CALTRANI (ISTL,3,DYE,DYE1)
C     IF (ISTRAN(4).EQ.2) CALL CALTRANI (ISTL,4,SED,SED1)
C     IF (ISTRAN(1).EQ.2) CALL CALTRANI (ISTL,1,SAL,UUU)
C     IF (ISTRAN(2).EQ.2) CALL CALTRANI (ISTL,2,TEM,VVV)
C     IF (ISTRAN(3).EQ.2) CALL CALTRANI (ISTL,3,DYE,DU)
C     IF (ISTRAN(4).EQ.2) CALL CALTRANI (ISTL,4,SED,DV)
C
      IF(ISCRAY.EQ.0) THEN
        TSADV=TSADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSADV=TSADV+T2TMP-T1TMP
        WTSADV=WTSADV+(WT2TMP-WT1TMP)*0.001
      END IF
C
      END IF
C
C**********************************************************************C
C
C **  3D ADVECTI0N TRANSPORT CALCULATION
C
C----------------------------------------------------------------------C
C
      IF(IS1DCHAN.GE.1) THEN
C
       IF(ISTL.EQ.3) THEN
         DO L=2,LA
          IF(LCT(L).EQ.6) AREAOLD(L)=FADYP2(L)
          IF(LCT(L).EQ.7) AREAOLD(L)=FADXP2(L)
          IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
          IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
         END DO
c       ELSE
         DO L=2,LA
          IF(LCT(L).EQ.6) AREAOLD(L)=FADYP1(L)
          IF(LCT(L).EQ.7) AREAOLD(L)=FADXP1(L)
          IF(LCT(L).EQ.6) AREANEW(L)=FADYP(L)
          IF(LCT(L).EQ.7) AREANEW(L)=FADXP(L)
         END DO
       END IF
C
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
!      IF (ISTRAN(1).EQ.1) CALL CALTRAN1D (ISTL,1,1,SAL,SAL1)
!      IF (ISTRAN(2).EQ.1) CALL CALTRAN1D (ISTL,2,2,TEM,TEM1)
!      IF (ISTRAN(3).EQ.1) CALL CALTRAN1D (ISTL,3,3,DYE,DYE1)
C     IF (ISTRAN(4).EQ.1) CALL CALTRAN1D (ISTL,4,4,SFL,SFL1)
      IF (ISTRAN(5).EQ.1.AND.NOCB.EQ.0) THEN
        DO NT=1,NTOX
         M=MSVTOX(NT)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=TOX1(L,K,NT)
          TVAR2S(L,K)=TOX(L,K,NT)
         END DO
         END DO
  !       CALL CALTRAN1D (ISTL,5,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          TOX1(L,K,NT)=TVAR1S(L,K)
          TOX(L,K,NT)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF
      IF (ISTRAN(6).EQ.1) THEN
        DO NS=1,NSED
         M=MSVSED(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SED1(L,K,NS)
          TVAR2S(L,K)=SED(L,K,NS)
         END DO
         END DO
   !      CALL CALTRAN1D (ISTL,6,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SED1(L,K,NS)=TVAR1S(L,K)
          SED(L,K,NS)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF
      IF (ISTRAN(7).EQ.1) THEN
        DO NS=1,NSND
         M=MSVSND(NS)
         DO K=1,KC
         DO L=1,LC
          TVAR1S(L,K)=SND1(L,K,NS)
          TVAR2S(L,K)=SND(L,K,NS)
         END DO
         END DO
 !        CALL CALTRAN1D (ISTL,7,M,TVAR2S,TVAR1S)
         DO K=1,KC
         DO L=1,LC
          SND1(L,K,NS)=TVAR1S(L,K)
          SND(L,K,NS)=TVAR2S(L,K)
         END DO
         END DO
        END DO
      END IF
C
      IF(ISCRAY.EQ.0) THEN
        TSADV=TSADV+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TSADV=TSADV+T2TMP-T1TMP
        WTSADV=WTSADV+(WT2TMP-WT1TMP)*0.001
      END IF
C
      END IF
C
C**********************************************************************C
C
C **  SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
C
      IF (ISTRAN(2).GE.1) CALL CALHEAT(ISTL)
C
C**********************************************************************C
C
C **  FULL IMPLICIT DYE AND TOXIC CONTAMINANT DECAY CALCULATION
C
      IF (ISTRAN(3).GE.1) THEN
C
      CDYETMP=1./(1.+DELT*RKDYE)
C     CDYETMP=(1.-DELTD2*RKDYE)/(1.+DELTD2*RKDYE)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=1,KC
        DO L=LF,LL
         DYE(L,K)=CDYETMP*DYE(L,K)
        END DO
       END DO
      END DO
C
      END IF
C
C
C      IF (ISTRAN(5).GE.1) THEN
C      DO NT=1,NTOX
C      CDYETMP=1./(1.+DELT*RKTOXW(NT))
C     CDYETMP=(1.-DELTD2*RKTOXW(NT))/(1.+DELTD2*RKTOXW(NT))
C      DO ND=1,NDM
C       LF=2+(ND-1)*LDM
C       LL=LF+LDM-1
C       DO K=1,KC
C        DO L=LF,LL
C         TOX(L,K,NT)=CDYETMP*TOX(L,K,NT)
C        END DO
C       END DO
C      END DO
C      END DO
C      END IF
C
C**********************************************************************C
C
C **  BOTTOM AND INTERNAL SEDIMENT AND TOXIC CONTAMINAT
C **  SOURCE-SINK CALCULATION
C
cjh5/13/77      IF (ISTRAN(6).GE.1) CALL CALSED(ISTL,1.0)
C
cjh5/13/77           IF (ISTRAN(7).GE.1) THEN
cjh5/13/77           NTSWVD=ISWVSD
cjh5/13/77            IF (ISTOPT(7).LE.1) CALL CALSED3(ISTL,1.0)
cjh5/13/77            IF (ISTOPT(7).EQ.2.AND.N.GE.NTSWVD)  CALL CALSED3(ISTL,1.0)
cjh5/13/77           END IF
C
cjh5/13/77           IF (ISTRAN(5).GE.1) CALL CALTOX(ISTL,1.0)
C
cx       NTMP=MOD(N,2)
cx       IF(NTMP.EQ.0) THEN
         IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
          if(NPCB.EQ.0) then
           CALL SSEDTOX(ISTL,1.0)
          else
           CALL SSED(ISTL,1.0)
          endif
         ENDIF           
cx       END IF
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION
C
C----------------------------------------------------------------------C
C
      IF(KC.EQ.1) GO TO 1500
C
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
      ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
      RCDZKK=-DELTD2*CDZKK(1)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCUBTMP=RCDZKK*HPI(L)*AB(L,1)
        CCMBTMP=1.-CCUBTMP
        EEB(L)=1./CCMBTMP
        CU1(L,1)=CCUBTMP*EEB(L)
       END DO
       IF(ISTRAN(1).GE.1) THEN
         DO L=LF,LL
          SAL(L,1)=SAL(L,1)*EEB(L)
         END DO
       END IF
       IF(ISTRAN(2).GE.1) THEN
         DO L=LF,LL
          TEM(L,1)=TEM(L,1)*EEB(L)
         END DO
       END IF
       IF(ISTRAN(3).GE.1) THEN
         DO L=LF,LL
          DYE(L,1)=DYE(L,1)*EEB(L)
         END DO
       END IF
       IF(ISTRAN(5).GE.1.and.NPCB.EQ.0) THEN
         DO NT=1,NTOX
          DO L=LF,LL
           TOX(L,1,NT)=TOX(L,1,NT)*EEB(L)
          END DO
         END DO
       END IF
       IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          DO L=LF,LL
           SED(L,1,NS)=SED(L,1,NS)*EEB(L)
          END DO
         END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
        DO NS=1,NSND
         DO L=LF,LL
          SND(L,1,NS)=SND(L,1,NS)*EEB(L)
         END DO
        END DO
       END IF
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=2,KS
        RCDZKMK=-DELTD2*CDZKMK(K)
        RCDZKK=-DELTD2*CDZKK(K)
        DO L=LF,LL
         CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)
         CCUBTMP=RCDZKK*HPI(L)*AB(L,K)
         CCMBTMP=1.-CCLBTMP(L)-CCUBTMP
         EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
         CU1(L,K)=CCUBTMP*EEB(L)
        END DO
        IF(ISTRAN(1).GE.1)THEN
          DO L=LF,LL
           SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)
          END DO
        END IF
        IF(ISTRAN(2).GE.1)THEN
          DO L=LF,LL
           TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)
          END DO
        END IF
        IF(ISTRAN(3).GE.1)THEN
          DO L=LF,LL
           DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)
          END DO
        END IF
        IF(ISTRAN(5).GE.1.and.NPCB.EQ.0)THEN
          DO NT=1,NTOX
           DO L=LF,LL
            TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
           END DO
          END DO
        END IF
        IF(ISTRAN(6).GE.1)THEN
          DO NS=1,NSED
           DO L=LF,LL
            SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)
           END DO
          END DO
        END IF
        IF(ISTRAN(7).GE.1)THEN
          DO NS=1,NSND
           DO L=LF,LL
            SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)
           END DO
          END DO
        END IF
       END DO
      END DO
C
      K=KC
      RCDZKMK=-DELTD2*CDZKMK(K)
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)
        CCMBTMP=1.-CCLBTMP(L)
        EEB(L)=1./(CCMBTMP-CCLBTMP(L)*CU1(L,K-1))
       END DO
       IF(ISTRAN(1).GE.1) THEN
         DO L=LF,LL
          SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)
         END DO
       END IF
       IF(ISTRAN(2).GE.1) THEN
         DO L=LF,LL
          TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)
         END DO
       END IF
       IF(ISTRAN(3).GE.1) THEN
         DO L=LF,LL
          DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)
         END DO
       END IF
       IF(ISTRAN(5).GE.1.and.NPCB.EQ.0) THEN
         DO NT=1,NTOX
          DO L=LF,LL
           TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)
          END DO
         END DO
       END IF
       IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          DO L=LF,LL
           SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)
          END DO
         END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          DO L=LF,LL
           SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)
          END DO
         END DO
       END IF
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO K=KC-1,1,-1
        IF(ISTRAN(1).GE.1) THEN
          DO L=LF,LL
           SAL(L,K)=SAL(L,K)-CU1(L,K)*SAL(L,K+1)
          END DO
        END IF
        IF(ISTRAN(2).GE.1) THEN
          DO L=LF,LL
           TEM(L,K)=TEM(L,K)-CU1(L,K)*TEM(L,K+1)
          END DO
        END IF
        IF(ISTRAN(3).GE.1) THEN
          DO L=LF,LL
           DYE(L,K)=DYE(L,K)-CU1(L,K)*DYE(L,K+1)
          END DO
        END IF
        IF(ISTRAN(5).GE.1.and.NPCB.EQ.0) THEN
          DO NT=1,NTOX
           DO L=LF,LL
            TOX(L,K,NT)=TOX(L,K,NT)-CU1(L,K)*TOX(L,K+1,NT)
           END DO
          END DO
        END IF
        IF(ISTRAN(6).GE.1) THEN
          DO NS=1,NSED
           DO L=LF,LL
            SED(L,K,NS)=SED(L,K,NS)-CU1(L,K)*SED(L,K+1,NS)
           END DO
          END DO
        END IF
        IF(ISTRAN(7).GE.1) THEN
          DO NS=1,NSND
           DO L=LF,LL
            SND(L,K,NS)=SND(L,K,NS)-CU1(L,K)*SND(L,K+1,NS)
           END DO
          END DO
        END IF
       END DO
      END DO
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
      END DO
      END DO
C
      DO K=1,KC
       DO L=1,LC
        SEDT(L,K)=0.
        SNDT(L,K)=0.
       END DO
      END DO
C
      DO K=1,KB
      DO NS=1,NSED
       DO L=1,LC
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
       END DO
      END DO
      END DO
C
      DO NS=1,NSND
      DO K=1,KB
       DO L=1,LC
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       END DO
      END DO
      END DO
C
      DO NS=1,NSED
       DO K=1,KC
        DO L=1,LC
         SEDT(L,K)=SEDT(L,K)+SED(L,K,NS)
c	write(904,9114) N, L,K,ns,sedt(L,K),sed(L,K,ns)
9114	format(4i6,9e12.4)
        END DO
       END DO
      END DO
C
      DO NS=1,NSND
       DO K=1,KC
        DO L=1,LC
         SNDT(L,K)=SNDT(L,K)+SND(L,K,NS)
        END DO
       END DO
      END DO
C
cjh5/13/97      DO NT=1,NTOX
cjh5/13/97       DO K=1,KC
cjh5/13/97        DO L=1,LC
cjh5/13/97         TOXPFT(L,K,NT)=0.
cjh5/13/97        END DO
cjh5/13/97       END DO
cjh5/13/97      END DO
C
cjh5/13/97      DO NT=1,NTOX
cjh5/13/97       DO NS=1,NSED+NSND
cjh5/13/97        DO K=1,KC
cjh5/13/97         DO L=1,LC
cjh5/13/97          TOXPFT(L,K,NT)=TOXPFT(L,K,NT)+TOXPF(L,K,NS,NT)
cjh5/13/97         END DO
cjh5/13/97        END DO
cjh5/13/97       END DO
cjh5/13/97      END DO
C
      IF(ISCRAY.EQ.0) THEN
        TVDIF=TVDIF+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TVDIF=TVDIF+T2TMP-T1TMP
        WTVDIF=WTVDIF+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
 1500 CONTINUE
C
C**********************************************************************C
C
C **  DATA ASSIMILATION
C
      IF(ISCDA(1).GT.0) THEN
C
 !     open(1,file='assim.out',ACCESS='APPEND',STATUS='UNKNOWN')
      DO K=1,KC
       DO NLC=1,NLCDA
        L=LIJ(ICDA(NLC),JCDA(NLC))
C
        NS=NCSERA(NLC,1)

        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,1).GT.1.0e-8) THEN
            T_lag=abs(CSERTT(NS)-TIME)
            dd_0=SAL(L,K)
            SAL(L,K)=TSCDA*CSERT(K,NS,1)+(1.-TSCDA)*SAL(L,K)
            dd_f=SAL(L,K)-dd_0
            II=IL(L)
            JJ=JL(L)
            do i1=-5,5
             do j1=-10,10
             if(IJCT(II+i1,JJ+j1).EQ.5) then
              L0=LIJ(II+i1,JJ+j1)
              x00=DLON(L0)
              y00=DLAT(L0)
              x11=XCDA(NLC)
              y11=YCDA(NLC) 
              d11=sqrt((x00-x11)*(x00-x11)+(y00-y11)*(y00-y11))
              dd0=SAL(L0,K)
              dec_0=TSCDA*exp(-d11/SLAGCDA)*exp(-T_lag/STLAGCDA/2.0)
              SAL(L0,K)=dec_0*CSERT(K,NS,1)+(1.-dec_0)*SAL(L0,K)
 !             write(1,*)time,NLC,dd0,dd0-SAL(L0,K)
             endif
             enddo
            enddo
          END IF
        END IF
C
        NS=NCSERA(NLC,2)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,2).GT.0) THEN
            TEM(L,K)=TSCDA*CSERT(K,NS,2)+(1.-TSCDA)*TEM(L,K)
          END IF
        END IF
C
        NS=NCSERA(NLC,3)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,3).GT.0) THEN
            DYE(L,K)=TSCDA*CSERT(K,NS,3)+(1.-TSCDA)*DYE(L,K)
          END IF
        END IF
C
        NS=NCSERA(NLC,4)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,4).GT.0) THEN
            SFL(L,K)=TSCDA*CSERT(K,NS,4)+(1.-TSCDA)*SFL(L,K)
          END IF
        END IF
C
        NS=NCSERA(NLC,5)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,5).GT.0) THEN
            TOX(L,K,1)=TSCDA*CSERT(K,NS,5)+(1.-TSCDA)*TOX(L,K,1)
          END IF
        END IF
C
        NS=NCSERA(NLC,6)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,6).GT.0) THEN
            SED(L,K,1)=TSCDA*CSERT(K,NS,6)+(1.-TSCDA)*SED(L,K,1)  ! John hardwired here, Ji, 10/26/00
cji, John used efdc.inp C66A & 66B to force model results = data at certain locations, Ji, 10/26/00
c
          END IF
        END IF
C
        NS=NCSERA(NLC,7)
        IF(NS.GT.0) THEN
          IF(CSERT(K,NS,7).GT.0) THEN
            SND(L,K,1)=TSCDA*CSERT(K,NS,7)+(1.-TSCDA)*SND(L,K,1)
          END IF
        END IF
C
       END DO
      END DO
C
      END IF
      
!      close(1)
C
C**********************************************************************C
C
C **  SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION
C
c     IF (ISTRAN(2).GE.1) CALL CALHEAT(ISTL)
C
C**********************************************************************C
C
C **  DYE DECAY CALCULATION
C
c     CDYETMP=1./(1.+DELTD2*RKDYE)
C     CDYETMP=(1.-DELTD2*RKDYE)/(1.+DELTD2*RKDYE)
c     DO K=1,KC
c     DO L=2,LA
c     DYE(L,K)=CDYETMP*DYE(L,K)
c     END DO
c     END DO
C
C**********************************************************************C
C
C **  BOTTOM AND INTERNAL SEDIMENT SOURCE-SINK CALCULATION
C
C     IF (ISTRAN(4).GE.1) CALL CALSED(ISTL,1.0)
c     IF (ISTRAN(4).EQ.5) CALL CALSED3(ISTL,1.0,SED)
C
C**********************************************************************C
C
      RETURN
C
C**********************************************************************C
C
C **  VERTICAL DIFFUSION CALCULATION OPTIMIZED FOR HP 9000 S700
C
 1000 CONTINUE
C
C----------------------------------------------------------------------C
C
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
CDHP  DO K=1,KS
CDHP  CALL vec_$mult_vector(HPI(1),AB(1,K),LA,ABHPI(1,K))
CDHP  END DO
C
CDHP  K=1
CDHP  RUNIT=1.0
CDHP  RTMP=-DELT*CDZKK(1)
C     DO L=1,LA
C     CUBTMP(L)=RTMP*HPI(L)*AB(L,1)=RTMP*ABHPI(L,1)
C     END DO
CDHP  CALL vec_$mult_constant(ABHPI(1,1),LA,RTMP,CUBTMP(1))
C     DO L=1,LA
C     CMBTMP(L)=1.-CUBTMP(L)
C     END DO
CDHP  CALL vec_$SUB_constant(CUBTMP(1),LA,RUNIT,CMBTMP(1))
CDHP  DO L=1,LA
CDHP  EB(L)=1./CMBTMP(L)
CDHP  END DO
C     DO L=1,LA
C     CU1(L,1)=CUBTMP(L)*EB(L)
C     END DO
CDHP  CALL vec_$mult_vector(CUBTMP(1),EB(1),LA,CU1(1,1))
C     IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SAL(L,1)=SAL(L,1)*EB(L)
C       END DO
C     END IF
CDHP  IF (ISTRAN(1).GE.1)
CDHP $  CALL vec_$mult_vector(SAL(1,1),EB(1),LA,SAL(1,1))
CDHP  IF (ISTRAN(2).GE.1)
CDHP $  CALL vec_$mult_vector(TEM(1,1),EB(1),LA,TEM(1,1))
CDHP  IF (ISTRAN(3).GE.1)
CDHP $  CALL vec_$mult_vector(DYE(1,1),EB(1),LA,DYE(1,1))
CDHP  IF (ISTRAN(4).GE.1)
CDHP $  CALL vec_$mult_vector(SED(1,1),EB(1),LA,SED(1,1))
C
CDHP  DO K=2,KS
CDHP  RTMP1=-DELT*CDZKMK(K)
CDHP  RTMP2=-DELT*CDZKK(K)
C      DO L=1,LA
C      CLBTMP(L)=RTMP1*AB(L,K-1)*HPI(L)=RTMP1*ABHPI(L,K-1)
C      CUBTMP(L)=RTMP2*AB(L,K)*HPI(L)=RTMP2*ABHPI(L,K)
C      END DO
CDHP  CALL vec_$mult_constant(ABHPI(1,K-1),LA,RTMP1,CLBTMP(1))
CDHP  CALL vec_$mult_constant(ABHPI(1,K),LA,RTMP2,CUBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=1.-CLBTMP(L)-CUBTMP(L)
C      END DO
CDHP  CALL vec_$SUB_constant(CUBTMP(1),LA,RUNIT,CMBTMP(1))
CDHP  CALL vec_$SUB_vector(CMBTMP(1),CLBTMP(1),LA,CMBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=CMBTMP(L)-CLBTMP(L)*CU1(L,K-1)
C      END DO
CDHP  CALL vec_$mult_rSUB_vector(CLBTMP(1),CU1(1,K-1),CMBTMP(1),LA,
CDHP $                           CMBTMP(1))
CDHP   DO L=1,LA
CDHP   EB(L)=1./CMBTMP(L)
CDHP   END DO
C      DO L=1,LA
C      CU1(L,K)=CUBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CUBTMP(1),EB(1),LA,CU1(1,K))
C      DO L=1,LA
C      CLBTMP(L)=-1.*CLBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CLBTMP(1),EB(1),LA,CLBTMP(1))
CDHP  CALL vec_$neg(CLBTMP(1),LA,CLBTMP(1))
CDHP  IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SAL(L,K)=EB(L)*SAL(L,K)+CLBTMP(L)*SAL(L,K-1)
C       END DO
CDHP    CALL vec_$mult_vector(CLBTMP(1),SAL(1,K-1),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SAL(1,K),VTMP(1),LA,SAL(1,K))
CDHP  END IF
CDHP  IF (ISTRAN(2).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),TEM(1,K-1),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),TEM(1,K),VTMP(1),LA,TEM(1,K))
CDHP  END IF
CDHP  IF (ISTRAN(3).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),DYE(1,K-1),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),DYE(1,K),VTMP(1),LA,DYE(1,K))
CDHP  END IF
CDHP  IF (ISTRAN(4).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),SED(1,K-1),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SED(1,K),VTMP(1),LA,SED(1,K))
CDHP  END IF
CDHP  END DO
C
CDHP  K=KC
CDHP  RTMP=-DELT*CDZKMK(KC)
C      DO L=1,LA
C      CLBTMP(L)=RTMP*AB(L,KC-1)*HPI(L)=RTMP*ABHPI(L,KS)
C      END DO
CDHP  CALL vec_$mult_constant(ABHPI(1,KS),LA,RTMP,CLBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=1.-CLBTMP(L)
C      END DO
CDHP  CALL vec_$SUB_constant(CLBTMP(1),LA,RUNIT,CMBTMP(1))
C      DO L=1,LA
C      CMBTMP(L)=CMBTMP(L)-CLBTMP(L)*CU1(L,K-1)
C      END DO
CDHP  CALL vec_$mult_rSUB_vector(CLBTMP(1),CU1(1,KS),CMBTMP(1),LA,
CDHP $                           CMBTMP(1))
CDHP   DO L=1,LA
CDHP   EB(L)=1./CMBTMP(L)
CDHP   END DO
C      DO L=1,LA
C      CLBTMP(L)=-1.*CLBTMP(L)*EB(L)
C      END DO
CDHP  CALL vec_$mult_vector(CLBTMP(1),EB(1),LA,CLBTMP(1))
CDHP  CALL vec_$neg(CLBTMP(1),LA,CLBTMP(1))
CDHP  IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SAL(L,K)=EB(L)*SAL(L,K)+CLBTMP(L)*SAL(L,K-1)
C       END DO
CDHP    CALL vec_$mult_vector(CLBTMP(1),SAL(1,KS),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SAL(1,KC),VTMP(1),LA,SAL(1,KC))
CDHP  END IF
CDHP  IF (ISTRAN(2).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),TEM(1,KS),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),TEM(1,KC),VTMP(1),LA,TEM(1,KC))
CDHP  END IF
CDHP  IF (ISTRAN(3).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),DYE(1,KS),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),DYE(1,KC),VTMP(1),LA,DYE(1,KC))
CDHP  END IF
CDHP  IF (ISTRAN(4).GE.1) THEN
CDHP    CALL vec_$mult_vector(CLBTMP(1),SED(1,KS),LA,VTMP(1))
CDHP    CALL vec_$mult_add_vector(EB(1),SED(1,KC),VTMP(1),LA,SED(1,KC))
CDHP  END IF
C
CDHP  DO K=KS,1,-1
C     IF (ISTRAN(1).GE.1) THEN
C       DO L=1,LA
C       SAL(L,K)=SAL(L,K)-CU1(L,K)*SAL(L,K+1)
C       END DO
C     END IF
CDHP  IF (ISTRAN(1).GE.1)
CDHP $  CALL vec_$mult_rSUB_vector(CU1(1,K),SAL(1,K+1),SAL(1,K),LA,
CDHP $                             SAL(1,K))
CDHP  IF (ISTRAN(2).GE.1)
CDHP $  CALL vec_$mult_rSUB_vector(CU1(1,K),TEM(1,K+1),TEM(1,K),LA,
CDHP $                             TEM(1,K))
CDHP  IF (ISTRAN(3).GE.1)
CDHP $  CALL vec_$mult_rSUB_vector(CU1(1,K),DYE(1,K+1),DYE(1,K),LA,
CDHP $                             DYE(1,K))
CDHP  IF (ISTRAN(4).GE.1)
CDHP $  CALL vec_$mult_rSUB_vector(CU1(1,K),SED(1,K+1),SED(1,K),LA,
CDHP $                             SED(1,K))
CDHP  END DO
C
      IF(ISCRAY.EQ.0) THEN
        TVDIF=TVDIF+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TVDIF=TVDIF+T2TMP-T1TMP
        WTVDIF=WTVDIF+(WT2TMP-WT1TMP)*0.001
      END IF
C
C----------------------------------------------------------------------C
C
      RETURN
      END
