C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CGRAD2 (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CGRAD2 SOLVES THE EXTERNAL MODE BY A CONJUGATE 
C **  GRADIENT SCHEME DOMAIN DECOMPOSITION PARALLELISM    
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      DIMENSION RSDTMP(16),RPCTMP(16),PAPTMP(16)
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
      DO ND=1,NDM
       RSDTMP(ND)=0.
       RPCTMP(ND)=0.
       PAPTMP(ND)=0.
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3S(L)=P(LSC(L))
        TVAR3W(L)=P(L-1   )
        TVAR3E(L)=P(L+1   )
        TVAR3N(L)=P(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       RCG(L)=-CCC(L)*P(L)-CCS(L)*TVAR3S(L)-CCW(L)*TVAR3W(L)
     $               -CCE(L)*TVAR3E(L)-CCN(L)*TVAR3N(L)+FPTMP(L)
       END DO
      END DO
C
      RPCG=0.
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        PCG(L)=RCG(L)*CCCI(L)
        RPCTMP(ND)=RPCTMP(ND)+RCG(L)*RCG(L)*CCCI(L)
       END DO
      END DO
C
      DO ND=1,NDM
       RPCG=RPCG+RPCTMP(ND)
      END DO
C
      ITER=0
C
C**********************************************************************C
C
  100 CONTINUE
C
      ITER=ITER+1
C
      DO ND=1,NDM
       RSDTMP(ND)=0.
       RPCTMP(ND)=0.
       PAPTMP(ND)=0.
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3S(L)=PCG(LSC(L))
        TVAR3W(L)=PCG(L-1   )
        TVAR3E(L)=PCG(L+1   )
        TVAR3N(L)=PCG(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       APCG(L)=CCC(L)*PCG(L)+CCS(L)*TVAR3S(L)+CCW(L)*TVAR3W(L)
     $       +CCE(L)*TVAR3E(L)+CCN(L)*TVAR3N(L)
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       PAPTMP(ND)=PAPTMP(ND)+APCG(L)*PCG(L)
       END DO
      END DO
C
      PAPCG=0.
      DO ND=1,NDM
       PAPCG=PAPCG+PAPTMP(ND)
      END DO
C
      ALPHA=RPCG/PAPCG
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       P(L)=P(L)+ALPHA*PCG(L)
       END DO
      END DO
C
C**********************************************************************C
C
      DO ND=1,NDM
       RSDTMP(ND)=0.
       RPCTMP(ND)=0.
       PAPTMP(ND)=0.
      END DO
C  
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       RCG(L)=RCG(L)-ALPHA*APCG(L)
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       RPCTMP(ND)=RPCTMP(ND)+RCG(L)*RCG(L)*CCCI(L)
       RSDTMP(ND)=RSDTMP(ND)+RCG(L)*RCG(L)
       END DO
      END DO
C
      RPCGN=0.
      RSQ=0. 
      DO ND=1,NDM
       RPCGN=RPCGN+RPCTMP(ND)
       RSQ=RSQ+RSDTMP(ND)
      END DO
C
      IF (RSQ .LE. RSQM) GO TO 200
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C  
      BETA=RPCGN/RPCG
      RPCG=RPCGN
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       PCG(L)=CCCI(L)*RCG(L)+BETA*PCG(L)
       END DO
      END DO
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
      DO ND=1,NDM
       RSDTMP(ND)=0.
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        TVAR3S(L)=P(LSC(L))
        TVAR3W(L)=P(L-1   )
        TVAR3E(L)=P(L+1   )
        TVAR3N(L)=P(LNC(L))
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
       PCG(L)=( CCC(L)*P(L)+CCS(L)*TVAR3S(L)+CCW(L)*TVAR3W(L)
     $                     +CCE(L)*TVAR3E(L)
     $                     +CCN(L)*TVAR3N(L)-FPTMP(L) )*CCCI(L)
       END DO
      END DO
C
      DO ND=1,NDM
       LF=2+(ND-1)*LDM
       LL=LF+LDM-1
       DO L=LF,LL
        RSDTMP(ND)=RSDTMP(ND)+PCG(L)*PCG(L)
       END DO
      END DO
C
      RSQ=0.
      DO ND=1,NDM
       RSQ=RSQ+RSDTMP(ND)
      END DO      
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
