C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RELAX (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RELAX SOLVES THE FINITE DIFFERENCE FORM
C **  OF A PSEUDO HEMHOLTZ EQUATION
C **
C **              CS(L)*P(LS)+CW(L)*P(L-1)                
C **              +CC(L)*P(L)+CE(L)*P(L+1)                
C **              +CN(L)*P(LN) = FP(L)                    
C **                                                      
C **  BY SUCCESSIVE OVER RELAXATION USING A RED-BLACK ORDERING 
C **  WITH CONVERGENCE MEASURED BY A GLOBAL SQUARE ERROR RSQ.  
C **  NON-CONVERGENCE IS SIGNALED WHEN THE ITERATIONS EXCEED A 
C **  MAXIMUM.                                                 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
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
C **  SELECT TIME LEVEL
C
      IF (ISTL.EQ.3) GO TO 250
C
C**********************************************************************C
C
C **  BEGIN TWO TIME LEVEL (CRANK-NICHOLSON) SOR PROCEDURE
C
      ITER=1
C
  200 CONTINUE
      RSQ=0.
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CCSR(L)*PBLK(LS)+CCWR(L)*PBLK(LW)+CCER(L)*PBLK(LE)
     $   +CCNR(L)*PBLK(LN)-FPR(L)
      PRED(L)=PRED(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      END DO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CCSB(L)*PRED(LS)+CCWB(L)*PRED(LW)+CCEB(L)*PRED(LE)
     $   +CCNB(L)*PRED(LN)-FPB(L)
      PBLK(L)=PBLK(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      END DO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      IF (RSQ .LE. RSQM) GO TO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C  
      ITER=ITER+1
      GO TO 200
C     
C**********************************************************************C
C**********************************************************************C
C
  250 CONTINUE
C
C**********************************************************************C
C**********************************************************************C
C
C **  BEGIN THREE TIME LEVEL (LEAP-FROG) SOR PROCEDURE
C
      ITER=1
C
  300 CONTINUE
      RSQ=0.
C
C**********************************************************************C
C
C **  RED CELL LOOP
C
      DO L=1,NRC
      LN=LBNRC(L)
      LS=LBSRC(L)
      LE=LBERC(L)
      LW=LBWRC(L)
      RSD=PRED(L)+CSR(L)*PBLK(LS)+CWR(L)*PBLK(LW)+CER(L)*PBLK(LE)
     $   +CNR(L)*PBLK(LN)-FPR(L)
      PRED(L)=PRED(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      END DO
C
C**********************************************************************C
C
C **  BLACK CELL LOOP
C
      DO L=1,NBC
      LN=LRNBC(L)
      LS=LRSBC(L)
      LE=LREBC(L)
      LW=LRWBC(L)
      RSD=PBLK(L)+CSB(L)*PRED(LS)+CWB(L)*PRED(LW)+CEB(L)*PRED(LE)
     $   +CNB(L)*PRED(LN)-FPB(L)
      PBLK(L)=PBLK(L)-RP*RSD
      RSQ=RSQ+RSD*RSD
      END DO
C
C**********************************************************************C
C
C **  CHECK SQUARED RESIDUAL CONVERGENCE CRITERIA
C
      IF (RSQ .LE. RSQM) GO TO 400
C
C **  CHECK MAXIMUM ITERATION CRITERIA
C
      IF (ITER .GE. ITERM) THEN
       WRITE(6,600)
       STOP
      END IF
C  
      ITER=ITER+1
      GO TO 300
C     
C**********************************************************************C
C
C **  LOAD PRED AND PBLK INTO P
C
  400 CONTINUE
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
C**********************************************************************C
C
  600 FORMAT(1X,'MAXIMUM ITERATIONS EXCEEDED IN EXTERNAL SOLUTION')
C
C**********************************************************************C
C
      RETURN
      END
