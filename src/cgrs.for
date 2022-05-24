C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CGRS (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CGRS SOLVES THE EXTERNAL MODE BY A RED-BLACK
C **  ORDERED REDUCED SYSTEM CONJUGATEC GRADIENT SCHEME
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C     DIMENSION CGSR(LCM),CGWR(LCM),CGER(LCM),CGNR(LCM),
C    &          CGSB(LCM),CGWB(LCM),CGEB(LCM),CGNB(LCM),
C    &          PCG(LCM),RCG(LCM),APCG(LCM),PBTMP(LCM),FPRT(LCM)
CDHP  DIMENSION PSOUT(LCM),PWEST(LCM),PEAST(LCM),PNORT(LCM),RSDT(LCM)
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
      IF (ISHP.GE.1) GO TO 1000
C
C**********************************************************************C
C
      IF(ISTL.EQ.2) THEN
       DO L=1,LC
       CGSR(L)=CCSR(L)
       CGWR(L)=CCWR(L)
       CGER(L)=CCER(L)
       CGNR(L)=CCNR(L)
       CGSB(L)=CCSB(L)
       CGWB(L)=CCWB(L)
       CGEB(L)=CCEB(L)
       CGNB(L)=CCNB(L)
       END DO
      ELSE
       DO L=1,LC
       CGSR(L)=CSR(L)
       CGWR(L)=CWR(L)
       CGER(L)=CER(L)
       CGNR(L)=CNR(L)
       CGSB(L)=CSB(L)
       CGWB(L)=CWB(L)
       CGEB(L)=CEB(L)
       CGNB(L)=CNB(L)
       END DO
      END IF
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      FPRT(L)=FPR(L)-CGSR(L)*FPB(LS)-CGWR(L)*FPB(LW)
     $              -CGER(L)*FPB(LE)-CGNR(L)*FPB(LN)
      END DO
C
C**********************************************************************C
C
C **  LOAD P INTO PRED AND PBLK
C
      DO LR=1,NRC
      L=LRC(LR)
      PRED(LR)=P(L)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      PBLK(LB)=P(L)
      END DO
C
C**********************************************************************C
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBTMP(L)=CGSB(L)*PRED(LS)+CGWB(L)*PRED(LW)
     $        +CGEB(L)*PRED(LE)+CGNB(L)*PRED(LN)
      END DO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)-CGSR(L)*PBTMP(LS)-CGWR(L)*PBTMP(LW)
     $           -CGER(L)*PBTMP(LE)-CGNR(L)*PBTMP(LN)-FPRT(L)
      RCG(L)=-RSD
      PCG(L)=-RSD
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
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBTMP(L)=CGSB(L)*PCG(LS)+CGWB(L)*PCG(LW)
     $        +CGEB(L)*PCG(LE)+CGNB(L)*PCG(LN)
      END DO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      APCG(L)=PCG(L)-CGSR(L)*PBTMP(LS)-CGWR(L)*PBTMP(LW)
     $              -CGER(L)*PBTMP(LE)-CGNR(L)*PBTMP(LN)
      END DO
C
      PAPCG=0.
      RPCG=0.
C
      DO L=1,NRC
      PAPCG=PAPCG+APCG(L)*PCG(L)
      RPCG=RPCG+RCG(L)*PCG(L)
      END DO
C
      ALPHA=RPCG/PAPCG
C
      DO L=1,NRC
      PRED(L)=PRED(L)+ALPHA*PCG(L)
      END DO
C
C**********************************************************************C
C
      RSQ=0.
C  
      DO L=1,NRC
      RCG(L)=RCG(L)-ALPHA*APCG(L)
      END DO
C
      DO L=1,NRC
      RSQ=RSQ+RCG(L)*RCG(L)
      END DO
C
      IF (RSQ .LE. RSQM) GO TO 200
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C  
      BETA=0.
C
      DO L=1,NRC
      BETA=BETA+RCG(L)*APCG(L)
      END DO
C
      BETA=-BETA/PAPCG
C
      DO L=1,NRC
      PCG(L)=RCG(L)+BETA*PCG(L)
      END DO
C
      GO TO 100
C
C**********************************************************************C
C
C **  CALCULATE PBLK AND FINAL RESIDUAL
C
  200 CONTINUE
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      PBLK(L)=FPB(L)-CGSB(L)*PRED(LS)-CGWB(L)*PRED(LW)
     $              -CGEB(L)*PRED(LE)-CGNB(L)*PRED(LN)
      END DO
C
      RSQ=0.
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CGSB(L)*PRED(LS)+CGWB(L)*PRED(LW)
     $           +CGEB(L)*PRED(LE)+CGNB(L)*PRED(LN)-FPB(L)
      RSQ=RSQ+RSD*RSD
      END DO
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CGSR(L)*PBLK(LS)+CGWR(L)*PBLK(LW)
     $           +CGER(L)*PBLK(LE)+CGNR(L)*PBLK(LN)-FPR(L)
      RSQ=RSQ+RSD*RSD
      END DO
C
C
C**********************************************************************C
C
C **  LOAD PRED AND PBLK INTO P
C
      DO LR=1,NRC
      L=LRC(LR)
      P(L)=PRED(LR)
      END DO
C
      DO LB=1,NBC
      L=LBC(LB)
      P(L)=PBLK(LB)
      END DO
C
      IF(ISCRAY.EQ.0) THEN
        TCGRS=TCGRS+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCGRS=TCGRS+T2TMP-T1TMP
        WTCGRS=WTCGRS+(WT2TMP-WT1TMP)*0.001
      END IF
C
      RETURN
C
C**********************************************************************C
C
C **  HP 9000 S700 OPTIMIZED VERSION OF SUBOURTINE
C
 1000 CONTINUE
C
C**********************************************************************C
C
      IF(ISTL.EQ.2) THEN
C       DO L=1,LC
C       CGSR(L)=CCSR(L)
C       CGWR(L)=CCWR(L)
C       CGER(L)=CCER(L)
C       CGNR(L)=CCNR(L)
C       CGSB(L)=CCSB(L)
C       CGWB(L)=CCWB(L)
C       CGEB(L)=CCEB(L)
C       CGNB(L)=CCNB(L)
C       END DO
CDHP   CALL vec_$copy(CCSR(1),CGSR(1),LC)
CDHP   CALL vec_$copy(CCWR(1),CGWR(1),LC)
CDHP   CALL vec_$copy(CCER(1),CGER(1),LC)
CDHP   CALL vec_$copy(CCNR(1),CGNR(1),LC)
CDHP   CALL vec_$copy(CCSB(1),CGSB(1),LC)
CDHP   CALL vec_$copy(CCWB(1),CGWB(1),LC)
CDHP   CALL vec_$copy(CCEB(1),CGEB(1),LC)
CDHP   CALL vec_$copy(CCNB(1),CGNB(1),LC)
      ELSE
C       DO L=1,LC
C       CGSR(L)=CSR(L)
C       CGWR(L)=CWR(L)
C       CGER(L)=CER(L)
C       CGNR(L)=CNR(L)
C       CGSB(L)=CSB(L)
C       CGWB(L)=CWB(L)
C       CGEB(L)=CEB(L)
C       CGNB(L)=CNB(L)
C       END DO
CDHP   CALL vec_$copy(CSR(1),CGSR(1),LC)
CDHP   CALL vec_$copy(CWR(1),CGWR(1),LC)
CDHP   CALL vec_$copy(CER(1),CGER(1),LC)
CDHP   CALL vec_$copy(CNR(1),CGNR(1),LC)
CDHP   CALL vec_$copy(CSB(1),CGSB(1),LC)
CDHP   CALL vec_$copy(CWB(1),CGWB(1),LC)
CDHP   CALL vec_$copy(CEB(1),CGEB(1),LC)
CDHP   CALL vec_$copy(CNB(1),CGNB(1),LC)
      END IF
C
C     DO L=1,NRC
C     FPRT(L)=FPR(L)-CGSR(L)*FPB(LBSRC(L))-CGWR(L)*FPB(LBWRC(L))
C    $              -CGER(L)*FPB(LBERC(L))-CGNR(L)*FPB(LBNRC(L))
C     END DO
CDHP  CALL vec_$gather(FPB(1),LBSRC(1),NRC,PSOUT(1))
CDHP  CALL vec_$gather(FPB(1),LBWRC(1),NRC,PWEST(1))
CDHP  CALL vec_$gather(FPB(1),LBERC(1),NRC,PEAST(1))
CDHP  CALL vec_$gather(FPB(1),LBNRC(1),NRC,PNORT(1))
C     DO L=1,NRC
C     FPRT(L)=FPR(L)-CGSR(L)*PSOUT(L)-CGWR(L)*PWEST(L)
C    $              -CGER(L)*PEAST(L)-CGNR(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$mult_rSUB_vector(CGSR(1),PSOUT(1),FPR(1),NRC,FPRT(1))
CDHP  CALL vec_$mult_rSUB_vector(CGWR(1),PWEST(1),FPRT(1),NRC,FPRT(1))
CDHP  CALL vec_$mult_rSUB_vector(CGER(1),PEAST(1),FPRT(1),NRC,FPRT(1))
CDHP  CALL vec_$mult_rSUB_vector(CGNR(1),PNORT(1),FPRT(1),NRC,FPRT(1))
C    
C**********************************************************************C
C
C **  LOAD P INTO PRED AND PBLK
C
C     DO LR=1,NRC
C     PRED(LR)=P(LRC(LR))
C     END DO
CDHP  CALL vec_$gather(P(1),LRC(1),NRC,PRED(1))
C
C     DO LB=1,NBC
C     PBLK(LB)=P(LBC(LB))
C     END DO
CDHP  CALL vec_$gather(P(1),LBC(1),NBC,PBLK(1))
C
C**********************************************************************C
C
C     DO L=1,NBC
C     PBTMP(L)=CGSB(L)*PRED(LRSBC(L))+CGWB(L)*PRED(LRWBC(L))
C    $        +CGEB(L)*PRED(LREBC(L))+CGNB(L)*PRED(LRNBC(L))
C     END DO
CDHP  CALL vec_$gather(PRED(1),LRSBC(1),NBC,PSOUT(1))
CDHP  CALL vec_$gather(PRED(1),LRWBC(1),NBC,PWEST(1))
CDHP  CALL vec_$gather(PRED(1),LREBC(1),NBC,PEAST(1))
CDHP  CALL vec_$gather(PRED(1),LRNBC(1),NBC,PNORT(1))
C     DO L=1,NBC
C     PBTMP(L)=CGSB(L)*PSOUT(L)+CGWB(L)*PWEST(L)
C    $        +CGEB(L)*PEAST(L)+CGNB(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$mult_vector(CGSB(1),PSOUT(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGWB(1),PWEST(1),PBTMP(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGEB(1),PEAST(1),PBTMP(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGNB(1),PNORT(1),PBTMP(1),NBC,PBTMP(1))
C    
C     DO L=1,NRC
C     RSDT(L)=PRED(L)-CGSR(L)*PBTMP(LBSRC(L))-CGWR(L)*PBTMP(LBWRC(L))
C    $               -CGER(L)*PBTMP(LBERC(L))-CGNR(L)*PBTMP(LBNRC(L))
C    $               -FPRT(L)
C     RCG(L)=-RSDT(L)
C     PCG(L)=-RSDT(L)
C     END DO
CDHP  CALL vec_$gather(PBTMP(1),LBSRC(1),NRC,PSOUT(1))
CDHP  CALL vec_$gather(PBTMP(1),LBWRC(1),NRC,PWEST(1))
CDHP  CALL vec_$gather(PBTMP(1),LBERC(1),NRC,PEAST(1))
CDHP  CALL vec_$gather(PBTMP(1),LBNRC(1),NRC,PNORT(1))
C     DO L=1,NRC
C     RCG(L)=FPRT(L)-PRED(L)+CGSR(L)*PSOUT(L)+CGWR(L)*PWEST(L)
C    $                      +CGER(L)*PEAST(L)+CGNR(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$SUB_vector(FPRT(1),PRED(1),NRC,RCG(1))
CDHP  CALL vec_$mult_add_vector(CGSR(1),PSOUT(1),RCG(1),NRC,RCG(1))
CDHP  CALL vec_$mult_add_vector(CGWR(1),PWEST(1),RCG(1),NRC,RCG(1))
CDHP  CALL vec_$mult_add_vector(CGER(1),PEAST(1),RCG(1),NRC,RCG(1))
CDHP  CALL vec_$mult_add_vector(CGNR(1),PNORT(1),RCG(1),NRC,RCG(1))
CDHP  CALL vec_$copy(RCG(1),PCG(1),NRC)
C
      ITER=0
C
C**********************************************************************C
C
  101 CONTINUE
C
      ITER=ITER+1
C
C     DO L=1,NBC
C     PBTMP(L)=CGSB(L)*PCG(LRSBC(L))+CGWB(L)*PCG(LRWBC(L))
C    $        +CGEB(L)*PCG(LREBC(L))+CGNB(L)*PCG(LRNBC(L))
C     END DO
CDHP  CALL vec_$gather(PCG(1),LRSBC(1),NBC,PSOUT(1))
CDHP  CALL vec_$gather(PCG(1),LRWBC(1),NBC,PWEST(1))
CDHP  CALL vec_$gather(PCG(1),LREBC(1),NBC,PEAST(1))
CDHP  CALL vec_$gather(PCG(1),LRNBC(1),NBC,PNORT(1))
C     DO L=1,NBC
C     PBTMP(L)=CGSB(L)*PSOUT(L)+CGWB(L)*PWEST(L)
C    $        +CGEB(L)*PEAST(L)+CGNB(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$mult_vector(CGSB(1),PSOUT(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGWB(1),PWEST(1),PBTMP(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGEB(1),PEAST(1),PBTMP(1),NBC,PBTMP(1))
CDHP  CALL vec_$mult_add_vector(CGNB(1),PNORT(1),PBTMP(1),NBC,PBTMP(1))
C
C     DO L=1,NRC
C     APCG(L)=PCG(L)-CGSR(L)*PBTMP(LBSRC(L))-CGWR(L)*PBTMP(LBWRC(L))
C    $              -CGER(L)*PBTMP(LBERC(L))-CGNR(L)*PBTMP(LBNRC(L))
C     END DO
CDHP  CALL vec_$gather(PBTMP(1),LBSRC(1),NRC,PSOUT(1))
CDHP  CALL vec_$gather(PBTMP(1),LBWRC(1),NRC,PWEST(1))
CDHP  CALL vec_$gather(PBTMP(1),LBERC(1),NRC,PEAST(1))
CDHP  CALL vec_$gather(PBTMP(1),LBNRC(1),NRC,PNORT(1))
C     DO L=1,NRC
C     APCG(L)=PCG(L)-CGSR(L)*PSOUT(L)-CGWR(L)*PWEST(L)
C    $              -CGER(L)*PEAST(L)-CGNR(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$mult_rSUB_vector(CGSR(1),PSOUT(1),PCG(1),NRC,APCG(1))
CDHP  CALL vec_$mult_rSUB_vector(CGWR(1),PWEST(1),APCG(1),NRC,APCG(1))
CDHP  CALL vec_$mult_rSUB_vector(CGER(1),PEAST(1),APCG(1),NRC,APCG(1))
CDHP  CALL vec_$mult_rSUB_vector(CGNR(1),PNORT(1),APCG(1),NRC,APCG(1))
C
C     PAPCG=0.
C     RPCG=0.
C     DO L=1,NRC
C     PAPCG=PAPCG+APCG(L)*PCG(L)
C     RPCG=RPCG+RCG(L)*PCG(L)
C     END DO
CDHP  PAPCG=vec_$dot( APCG(1),PCG(1),NRC)
CDHP  RPCG=vec_$dot( RCG(1),PCG(1),NRC)
C
CDHP  ALPHA=RPCG/PAPCG
C
C     DO L=1,NRC
C     PRED(L)=PRED(L)+ALPHA*PCG(L)
C     END DO
CDHP  CALL vec_$mult_add(PRED(1),PCG(1),NRC,ALPHA,PRED(1))
C
C**********************************************************************C
C
C     RSQ=0.
C  
C     DO L=1,NRC
C     RCG(L)=RCG(L)-ALPHA*APCG(L)
C     END DO
CDHP  CALL vec_$mult_add(RCG(1),APCG(1),NRC,-ALPHA,RCG(1))
C
C     DO L=1,NRC
C     RSQ=RSQ+RCG(L)*RCG(L)
C     END DO
CDHP  RSQ=vec_$dot(RCG(1),RCG(1),NRC)
C
      IF (RSQ .LE. RSQM) GO TO 201
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C  
C     BETA=0.
C     DO L=1,NRC
C     BETA=BETA+RCG(L)*APCG(L)
C     END DO
CDHP  BETA=vec_$dot(RCG(1),APCG(1),NRC)
C
CDHP  BETA=-BETA/PAPCG
C
C     DO L=1,NRC
C     PCG(L)=RCG(L)+BETA*PCG(L)
C     END DO
CDHP  CALL vec_$mult_add(RCG(1),PCG(1),NRC,BETA,PCG(1))
C
      GO TO 101
C
C**********************************************************************C
C
C **  CALCULATE PBLK AND FINAL RESIDUAL
C
  201 CONTINUE
C
C     DO L=1,NBC
C     PBLK(L)=FPB(L)-CGSB(L)*PRED(LRSBC(L))-CGWB(L)*PRED(LRWBC(L))
C    $              -CGEB(L)*PRED(LREBC(L))-CGNB(L)*PRED(LRNBC(L))
C     END DO
CDHP  CALL vec_$gather(PRED(1),LRSBC(1),NBC,PSOUT(1))
CDHP  CALL vec_$gather(PRED(1),LRWBC(1),NBC,PWEST(1))
CDHP  CALL vec_$gather(PRED(1),LREBC(1),NBC,PEAST(1))
CDHP  CALL vec_$gather(PRED(1),LRNBC(1),NBC,PNORT(1))
C     DO L=1,NBC
C     PBLK(L)=FPB(L)-CGSB(L)*PSOUT(L)-CGWB(L)*PWEST(L)
C    $              -CGEB(L)*PEAST(L)-CGNB(L)*PNORT(L)
C     END DO
CDHP  CALL vec_$mult_rSUB_vector(CGSB(1),PSOUT(1),FPB(1),NBC,PBLK(1))
CDHP  CALL vec_$mult_rSUB_vector(CGWB(1),PWEST(1),PBLK(1),NBC,PBLK(1))
CDHP  CALL vec_$mult_rSUB_vector(CGEB(1),PEAST(1),PBLK(1),NBC,PBLK(1))
CDHP  CALL vec_$mult_rSUB_vector(CGNB(1),PNORT(1),PBLK(1),NBC,PBLK(1))
C
C     DO L=1,NBC
C     RSDT(L)=PBLK(L)+CGSB(L)*PSOUT(L)+CGWB(L)*PWEST(L)
C    $               +CGEB(L)*PEAST(L)+CGNB(L)*PNORT(L)
C    $               -FPB(L)
C     RSQ=RSQ+RSDT(L)*RSDT(L)
C     END DO
CDHP  CALL vec_$SUB_vector(PBLK(1),FPB(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGSB(1),PSOUT(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGWB(1),PWEST(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGEB(1),PEAST(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGNB(1),PNORT(1),RSDT(1),NBC,RSDT(1))
CDHP  RSQ=vec_$dot(RSDT(1),RSDT(1),NBC)
C
CDHP  CALL vec_$gather(PBLK(1),LBSRC(1),NRC,PSOUT(1))
CDHP  CALL vec_$gather(PBLK(1),LBWRC(1),NRC,PWEST(1))
CDHP  CALL vec_$gather(PBLK(1),LBERC(1),NRC,PEAST(1))
CDHP  CALL vec_$gather(PBLK(1),LBNRC(1),NRC,PNORT(1))
C     DO L=1,NRC
C     RSDT(L)=PRED(L)+CGSR(L)*PSOUT(L)+CGWR(L)*PWEST(L)
C    $               +CGER(L)*PEAST(L)+CGNR(L)*PNORT(L)
C    $               -FPR(L)
C     RSQ=RSQ+RSDT(L)*RSDT(L)
C     END DO
CDHP  CALL vec_$SUB_vector(PRED(1),FPR(1),NRC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGSR(1),PSOUT(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGWR(1),PWEST(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGER(1),PEAST(1),RSDT(1),NBC,RSDT(1))
CDHP  CALL vec_$mult_add_vector(CGNR(1),PNORT(1),RSDT(1),NBC,RSDT(1))
CDHP  RSQ=RSQ+vec_$dot(RSDT(1),RSDT(1),NRC)
C
C**********************************************************************C
C
C **  LOAD PRED AND PBLK INTO P
C
C     DO LR=1,NRC
C     L=LRC(LR)
C     P(L)=PRED(LR)
C     END DO
CDHP  CALL vec_$scatter(PRED(1),LRC(1),NRC,P(1))
C
C     DO LB=1,NBC
C     L=LBC(LB)
C     P(L)=PBLK(LB)
C     END DO
CDHP  CALL vec_$scatter(PBLK(1),LBC(1),NBC,P(1))
C
      IF(ISCRAY.EQ.0) THEN
        TCGRS=TCGRS+SECNDS(TTMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TCGRS=TCGRS+T2TMP-T1TMP
        WTCGRS=WTCGRS+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
  600 FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
      RETURN
      END
