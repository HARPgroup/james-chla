C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CONGRAD (ISTL) 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C **  SUBROUTINE CONGRAD SOLVES THE EXTERNAL MODE BY A CONJUGATE
C **  GRADIENT SCHEME
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C     DIMENSION CGS(LCM),CGW(LCM),CGE(LCM),CGN(LCM),
C    &          PCG(LCM),RCG(LCM),APCG(LCM)
C
C**********************************************************************C
C
C     IF(ISTL.EQ.2) THEN
C      DO L=2,LA
C      CIT=1./CCC(L)
C      CGS(L)=CIT*CCS(L)
C      CGW(L)=CIT*CCW(L)
C      CGE(L)=CIT*CCE(L)
C      CGN(L)=CIT*CCN(L) 
C      END DO
C     ELSE
C      DO L=2,LA
C      CIT=1./CC(L)
C      CGS(L)=CIT*CS(L)
C      CGW(L)=CIT*CW(L)
C      CGE(L)=CIT*CE(L)
C      CGN(L)=CIT*CN(L)
C      END DO
C     END IF
C
C
C     IF(ISTL.EQ.2) THEN
C      DO L=2,LA
C      CG(L)=CCC(L)
C      CGS(L)=CCS(L)
C      CGW(L)=CCW(L)
C      CGE(L)=CCE(L)
C      CGN(L)=CCN(L)
C      END DO
C     ELSE
C      DO L=2,LA
C      CG(L)=CC(L)
C      CGS(L)=CS(L)
C      CGW(L)=CW(L)
C      CGE(L)=CE(L)
C      CGN(L)=CN(L)
C      END DO
C     END IF
C
C     DO L=2,LA
C     LN=LNC(L)
C     LS=LSC(L)
C     FPTMP(L)=FP(L)-CG(L)*PAM(L)-CGS(L)*PAM(LS)-CGW(L)*PAM(L-1)
C    $        -CGE(L)*PAM(L+1)-CGN(L)*PAM(LN)
C     END DO
C
C**********************************************************************C
C
      IF(ISCRAY.EQ.0) THEN
        TTMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
!$OMP PARALLEL
!$OMP DO PRIVATE(LN,LS)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     $        +CCN(L)*P(LN)-FPTMP(L)
      RCG(L)=-RSD
COLD     PCG(L)=-RSD
      PCG(L)=-RSD*CCCI(L)
      END DO
!$OMP END DO
!$OMP END PARALLEL
C
      RPCG=0.
!$OMP PARALLEL
!$OMP DO REDUCTION(+:RPCG)      
      DO L=2,LA
      RPCG=RPCG+RCG(L)*RCG(L)*CCCI(L)
      END DO
!$OMP END DO      
!$OMP END PARALLEL
C
      ITER=0
C
C**********************************************************************C
C
  100 CONTINUE
C
      ITER=ITER+1
C
!$OMP PARALLEL
!$OMP DO PRIVATE(LN,LS)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      APCG(L)=CCC(L)*PCG(L)+CCS(L)*PCG(LS)+CCW(L)*PCG(L-1)
     $       +CCE(L)*PCG(L+1)+CCN(L)*PCG(LN)
      END DO
C
!$OMP ENDDO
!$OMP END PARALLEL
      PAPCG=0.
COLD     RPCG=0.
C
!$OMP PARALLEL
!$OMP DO REDUCTION(+:PAPCG)
      DO L=2,LA
      PAPCG=PAPCG+APCG(L)*PCG(L)
COLD     RPCG=RPCG+RCG(L)*PCG(L)
      END DO
!$OMP END DO
!$OMP END PARALLEL
C
      ALPHA=RPCG/PAPCG
C
!$OMP PARALLEL
!$OMP DO 
      DO L=2,LA
      P(L)=P(L)+ALPHA*PCG(L)
      END DO
!$OMP END DO
!$OMP END PARALLEL      
C
C**********************************************************************C
C
      RSQ=0.
C
      RPCGN=0.
!$OMP PARALLEL
!$OMP DO REDUCTION(+:RPCGN)      
      DO L=2,LA
      RCG(L)=RCG(L)-ALPHA*APCG(L)
      RPCGN=RPCGN+RCG(L)*RCG(L)*CCCI(L)
      END DO
!$OMP END DO
!$OMP END PARALLEL
C
!$OMP PARALLEL
!$OMP DO REDUCTION(+:RSQ)
      DO L=2,LA
      RSQ=RSQ+RCG(L)*RCG(L)
      END DO
!$OMP ENDDO
!$OMP END PARALLEL      
C
      IF (RSQ .LE. RSQM) GO TO 200
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C
COLD     BETA=0.
C
COLD     DO L=2,LA
COLD     BETA=BETA+RCG(L)*APCG(L)
COLD     END DO
C
COLD     BETA=-BETA/PAPCG
      BETA=RPCGN/RPCG
      RPCG=RPCGN
C
!$OMP PARALLEL
!$OMP DO
      DO L=2,LA
COLD     PCG(L)=RCG(L)+BETA*PCG(L)
      PCG(L)=CCCI(L)*RCG(L)+BETA*PCG(L)
      END DO
!$OMP END DO
!$OMP END PARALLEL
C
      GO TO 100
C
  600 FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
C**********************************************************************C
C
C ** CALCULATE FINAL RESIDUAL
C
  200 CONTINUE
C
      RSQ=0.
C
!$OMP DO PRIVATE(LN,LS,RSD) REDUCTION(+:RSQ)
      DO L=2,LA
      LN=LNC(L)
      LS=LSC(L)
      RSD=CCC(L)*P(L)+CCS(L)*P(LS)+CCW(L)*P(L-1)+CCE(L)*P(L+1)
     $        +CCN(L)*P(LN)-FPTMP(L)
      RSD=RSD*CCCI(L)
      RSQ=RSQ+RSD*RSD
      END DO
!$OMP END DO      
C
      IF(ISCRAY.EQ.0) THEN
        TCONG=TCONG+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCONG=TCONG+T2TMP-T1TMP
        WTCONG=WTCONG+(WT2TMP-WT1TMP)*0.001
      END IF
C
      RETURN
      END
