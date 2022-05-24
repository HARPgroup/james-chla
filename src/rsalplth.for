C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSALPLTH_ac(ICON,CONC)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RSALPLTH WRITES FILES FOR RESIDUAL SCALAR FIELD 
C **  CONTOURING IN HORIZONTAL PLANES       
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      DIMENSION CONC(LCM,KCM)
      REAL DBS(10)
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF (JSRSPH(ICON).NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
      LINES=LA-1
      LEVELS=2
      LEVELSS=3
      DBS(1)=0.
      DBS(2)=99.
      DBS(3)=-99.
C
      IF (ISTRAN(1).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL SALINITY CONTOURS'
        LUN=11
        OPEN(LUN,FILE='rsalcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsalcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(2).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL TEMPERATURE CONTOURS'
        LUN=12
        OPEN(LUN,FILE='rtemcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rtemcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(3).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL DYE CONC CONTOURS'
        LUN=13
        OPEN(LUN,FILE='rdyecnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rdyecnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(6).GE.1) THEN
        TITLE='RESIDUAL HORIZ COHESIVE SEDIMENT CONC CONTOURS'
        LUN=14
        OPEN(LUN,FILE='rsedcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsedcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(7).GE.1) THEN
        TITLE='RESIDUAL HORIZ NONCOH SEDIMENT CONC CONTOURS'
        LUN=15
        OPEN(LUN,FILE='rsndcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsndcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(5).GE.1) THEN
        TITLE='RESIDUAL HORIZ TOXIC CONTAM CONC CONTOURS'
        LUN=16
        OPEN(LUN,FILE='rtoxcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rtoxcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
        TITLE='RESIDUAL HORIZ TOXIC PART FRAC CONTOURS'
        LUNF=26
        OPEN(LUNF,FILE='rtxpcnh.out',STATUS='UNKNOWN')
        CLOSE(LUNF,STATUS='DELETE')
        OPEN(LUNF,FILE='rtxpcnh.out',STATUS='UNKNOWN')
        WRITE (LUNF,99) TITLE
        WRITE (LUNF,101)LINES,LEVELSS
        WRITE (LUNF,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUNF)
      END IF
C
      IF (ISTRAN(4).GE.1) THEN
        TITLE='RESIDUAL HORIZONTAL SFL CONC CONTOURS'
        LUN=17
        OPEN(LUN,FILE='rsflcnh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='rsflcnh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
      DO ITMP=1,7
       JSRSPH(ITMP)=0
      END DO 
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      IF(ICON.EQ.1) THEN
        LUN=11
        OPEN(LUN,FILE='rsalcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.2) THEN
        LUN=12
        OPEN(LUN,FILE='rtemcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.3) THEN
        LUN=13
        OPEN(LUN,FILE='rdyecnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.6) THEN
        LUN=14
        OPEN(LUN,FILE='rsedcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.7) THEN
        LUN=15
        OPEN(LUN,FILE='rsndcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.5) THEN
        LUN=16
        OPEN(LUN,FILE='rtoxcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
        LUNF=26
        OPEN(LUNF,FILE='rtxpcnh.out',
     $                  ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.4) THEN
        LUN=17
        OPEN(LUN,FILE='rsflcnh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      WRITE (LUN,100)N,TIME
      IF(ICON.EQ.5) THEN
        WRITE (LUNF,100)N,TIME
      END IF
C
      IF(ISPHXY(ICON).EQ.0) THEN
      IF(ICON.LE.3) THEN
       DO L=2,LA
       WRITE(LUN,400)CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.GE.8) THEN
       DO L=2,LA
       WRITE(LUN,400)CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.EQ.6) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L)
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               SEDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L))
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               SNDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
       TOXBT=1000.*TOXBLPF(L,KBT(L),1)
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               TOXBT
       END DO
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       WRITE(LUNF,400)TXPFLPF(L,KC,1,1),
     $               TXPFLPF(L,1,1,1),TOXBLPF(L,KBT(L),1)
       END DO
       CLOSE(LUNF)
      END IF
      IF(ICON.EQ.4) THEN
       DO L=2,LA
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               SFLSBOT(L)
       END DO
      END IF
      END IF
C
C
      IF(ISPHXY(ICON).EQ.1) THEN
      IF(ICON.LE.3) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.GE.8) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.EQ.6) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L)
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SEDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L))
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SNDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
       TOXBT=1000.*TOXBLPF(L,KBT(L),1)
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               TOXBT
       END DO
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       WRITE(LUNF,200)IL(L),JL(L),TXPFLPF(L,KC,1,1),
     $               TXPFLPF(L,1,1,1),TOXBLPF(L,KBT(L),1)
       END DO
       CLOSE(LUNF)
      END IF
      IF(ICON.EQ.4) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SFLSBOT(L)
       END DO
      END IF
      END IF
C
C
      IF(ISPHXY(ICON).EQ.2) THEN
      IF(ICON.LE.3) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.GE.8) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1)
       END DO
      END IF
      IF(ICON.EQ.6) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L)
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               SEDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L))
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               SNDBTLPF(L,KBT(L))
       END DO
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
       TOXBT=1000.*TOXBLPF(L,KBT(L),1)
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               TOXBT
       END DO
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       WRITE(LUNF,200)IL(L),JL(L),DLON(L),DLAT(L),TXPFLPF(L,KC,1,1),
     $               TXPFLPF(L,1,1,1),TOXBLPF(L,KBT(L),1)
       END DO
       CLOSE(LUNF)
      END IF
      IF(ICON.EQ.4) THEN
       DO L=2,LA
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               SFLSBOT(L)
       END DO
      END IF
      END IF
C
      CLOSE(LUN)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
  400 FORMAT(1X,6E14.6)
  420 FORMAT(1X,13E11.3)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
