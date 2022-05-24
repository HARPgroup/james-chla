C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SALPLTH (ICON,CONC)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE SALPLTH WRITES FILES FOR INSTANTANEOUS SCALAR FIELD 
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
      REAL DBS(10), DBSB(0:NSTM)
      CHARACTER*80 TITLE
C
C**********************************************************************C
C
      IF (JSSPH(ICON).NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
      LINES=LA-1
      LEVELS=2
      LEVELSS=3
      DBS(1)=0.
      DBS(2)=99.
      DBS(3)=-99.
      LSEDCL=NSED+NSND
      DO L=0,LSEDCL
       DBSB(L)=FLOAT(L)
      END DO
C
      IF (ISTRAN(1).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZONTAL SALINITY CONTOURS'
        LUN=11
        OPEN(LUN,FILE='salconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='salconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(2).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZONTAL TEMPERATURE CONTOURS'
        LUN=12
        OPEN(LUN,FILE='temconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='temconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(3).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZONTAL DYE CONC CONTOURS'
        LUN=13
        OPEN(LUN,FILE='dyeconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='dyeconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELS
        WRITE (LUN,250)(DBS(L),L=1,LEVELS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(6).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZ COHESIVE SEDIMENT CONC CONTOURS'
        LUN=14
        OPEN(LUN,FILE='sedconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='sedconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(6).GE.1) THEN
        TITLE='INSTANTANEOUS BED SED DEPOSITED CONTOURS gm/m**2'
        LUN=15
        OPEN(LUN,FILE='sbdconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='sbdconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LSEDCL
        WRITE (LUN,250)(DBSB(L),L=0,LSEDCL)
        CLOSE(LUN)
      END IF
C
      IF (ISTRAN(7).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZ NONCOH SEDIMENT CONC CONTOURS'
        LUN=15
        OPEN(LUN,FILE='sndconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='sndconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
C     IF (ISTRAN(7).GE.1) THEN
C       TITLE='INSTANTANEOUS BED SED DEPOSITED CONTOURS gm/m**2'
C       LUN=15
C       OPEN(LUN,FILE='sbdconh.out',STATUS='UNKNOWN')
C       CLOSE(LUN,STATUS='DELETE')
C       OPEN(LUN,FILE='sbdconh.out',STATUS='UNKNOWN')
C       WRITE (LUN,99) TITLE
C       WRITE (LUN,101)LINES,LSEDCL
C       WRITE (LUN,250)(DBSB(L),L=0,LSEDCL)
C       CLOSE(LUN)
C     END IF
C
      IF (ISTRAN(5).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZ TOXIC CONTAM. CONC CONTOURS'
        LUN=16
        OPEN(LUN,FILE='toxconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='toxconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
        TITLE='INSTANTANEOUS HORIZ TOXIC PART FRAC CONTOURS'
        LUNF=26
        OPEN(LUNF,FILE='txpconh.out',STATUS='UNKNOWN')
        CLOSE(LUNF,STATUS='DELETE')
        OPEN(LUNF,FILE='txpconh.out',STATUS='UNKNOWN')
        WRITE (LUNF,99) TITLE
        WRITE (LUNF,101)LINES,LEVELSS
        WRITE (LUNF,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUNF)
      END IF
C
      IF (ISTRAN(4).GE.1) THEN
        TITLE='INSTANTANEOUS HORIZONTAL SFL CONC CONTOURS'
        LUN=17
        OPEN(LUN,FILE='sflconh.out',STATUS='UNKNOWN')
        CLOSE(LUN,STATUS='DELETE')
        OPEN(LUN,FILE='sflconh.out',STATUS='UNKNOWN')
        WRITE (LUN,99) TITLE
        WRITE (LUN,101)LINES,LEVELSS
        WRITE (LUN,250)(DBS(L),L=1,LEVELSS)
        CLOSE(LUN)
      END IF
C
      DO ITMP=1,7
       JSSPH(ITMP)=0
      END DO
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      IF(ICON.EQ.1) THEN
        LUN=11
        OPEN(LUN,FILE='salconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.2) THEN
        LUN=12
        OPEN(LUN,FILE='temconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.3) THEN
        LUN=13
        OPEN(LUN,FILE='dyeconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.6) THEN
        LUN=14
        OPEN(LUN,FILE='sedconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.7) THEN
        LUN=15
        OPEN(LUN,FILE='sndconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.5) THEN
        LUN=16
        OPEN(LUN,FILE='toxconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
        LUNF=26
        OPEN(LUNF,FILE='txpconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.4) THEN
        LUN=17
        OPEN(LUN,FILE='sflconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF(ICON.EQ.6) THEN
        LUB=18
        OPEN(LUB,FILE='sbdconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C     IF(ICON.EQ.7) THEN
C       LUB=18
C       OPEN(LUB,FILE='sbdconh.out',ACCESS='APPEND',STATUS='UNKNOWN')
C     END IF
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      WRITE (LUN,100)N,TIME
      IF(ICON.EQ.5) THEN
        WRITE (LUNF,100)N,TIME
      END IF
      IF(ICON.EQ.6) THEN
        WRITE (LUB,100)N,TIME
      END IF
c     IF(ICON.EQ.7) THEN
c       WRITE (LUB,100)N,TIME
c     END IF
C
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
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               SEDBT(L,KBT(L))
       WRITE(LUB,420)VOLBW2(L,KBT(L)),
     $      (SEDB(L,KBT(L),NX),NX=1,NSED),(SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
       CLOSE(LUB)
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L),1)
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               SNDBT(L,KBT(L))
C      WRITE(LUB,420)VOLBW2(L,KBT(L)),
C    $      (SEDB(L,KBT(L),NX),NX=1,NSED),(SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
C      CLOSE(LUB)
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
       TOXBT=1000.*TOXB(L,KBT(L),1)
       WRITE(LUN,400)CONC(L,KC),CONC(L,1),
     $               TOXBT
       END DO
       DO L=2,LA
       TOXBT=1000.*TOXB(L,KBT(L),1)
       WRITE(LUNF,400)TOXPFTW(L,KC,1),
     $               TOXPFTW(L,1,1),TOXBT
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
       IF(KC.GT.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1)
       IF(KC.EQ.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,1)
       END DO
      END IF
      IF(ICON.GE.8) THEN
       DO L=2,LA
       IF(KC.GT.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1)
       IF(KC.EQ.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,1)
       END DO
      END IF
      IF(ICON.EQ.6) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       IF(KC.GT.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SEDBT(L,KBT(L))
       IF(KC.EQ.1) WRITE(LUN,200)IL(L),JL(L),SEDT(L,1),
     $               SEDBT(L,KBT(L))
c       WRITE(LUB,220)IL(L),JL(L),VOLBW2(L,KBT(L)),
c     $        (SEDB(L,NX),NX=1,NSED),(*SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
       CLOSE(LUB)
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L),1)
       IF(KC.GT.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               SNDBT(L,KBT(L))
       IF(KC.EQ.1) WRITE(LUN,200)IL(L),JL(L),SNDT(L,1),
     $               SNDBT(L,KBT(L))
C      WRITE(LUB,220)IL(L),JL(L),VOLBW2(L,KBT(L)),
C    $        (SEDB(L,NX),NX=1,NSED),(SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
C      CLOSE(LUB)
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
C       TOXBT=1000.*TOXB(L,1)
       IF(KC.GT.1) WRITE(LUN,200)IL(L),JL(L),CONC(L,KC),CONC(L,1),
     $               TOXBT
C HARD WIRE FOR TWO TOXICS
       TOXB1TMP=0.
       TOXB2TMP=0.
       IF(VOLBW2(L,KBT(L)).GT.1.E-12) 
     $       TOXB1TMP=TOXB(L,KBT(L),1)/VOLBW2(L,KBT(L))
c       IF(VOLBW2(L,KBT(L)).GT.1.E-12) 
c     $       TOXB2TMP=TOXB(L,KBT(L),2)/VOLBW2(L,KBT(L))
       IF(KC.EQ.1.AND.NTOX.EQ.1) WRITE(LUN,200)IL(L),JL(L),TOX(L,1,1),
     $             TOXB1TMP
       IF(KC.EQ.1.AND.NTOX.GE.1) WRITE(LUN,200)IL(L),JL(L),TOX(L,1,1),
     $             TOXB1TMP
c       IF(KC.EQ.1.AND.NTOX.GE.2) WRITE(LUN,200)IL(L),JL(L),TOX(L,1,1),
c     $             TOXB1TMP,TOX(L,1,2),TOXB2TMP
C END HARDWIRE
       END DO
       DO L=2,LA
C       TOXBT=1000.*TOXB(L,1)
       IF(KC.GT.1) WRITE(LUNF,200)IL(L),JL(L),TOXPFTW(L,KC,1),
     $               TOXPFTW(L,1,1),TOXBT
C HARD WIRE FOR TWO TOXICS
       IF(KC.EQ.1.AND.NTOX.EQ.1) WRITE(LUNF,200)IL(L),JL(L),
     $       TOXPFTW(L,1,1),TOXPFTB(L,KBT(L),1)
       IF(KC.EQ.1.AND.NTOX.GE.1) WRITE(LUNF,200)IL(L),JL(L),
     $       TOXPFTW(L,1,1),TOXPFTB(L,KBT(L),1)
C       IF(KC.EQ.1.AND.NTOX.GE.2)  WRITE(LUNF,200)IL(L),JL(L),
C     $       TOXPFTW(L,1,1),TOXPFTB(L,1),TOXPFTW(L,1,2),TOXPFTB(L,2)
C END HARDWIRE
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
C      SEDBT=1000*SSG*SEDB(L,KBT(L),1)
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               SEDBT(L,KBT(L))
       WRITE(LUB,220)IL(L),JL(L),DLON(L),DLAT(L),VOLBW2(L,KBT(L)),
     $      (SEDB(L,KBT(L),NX),NX=1,NSED),(SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
       CLOSE(LUB)
      END IF
      IF(ICON.EQ.7) THEN
       DO L=2,LA
C      SEDBT=1000*SSG*SNDB(L,KBT(L),1)
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               SNDBT(L,KBT(L))
C      WRITE(LUB,220)IL(L),JL(L),DLON(L),DLAT(L),VOLBW2(L,KBT(L)),
C    $      (SEDB(L,KBT(L),NX),NX=1,NSED),(SNDB(L,KBT(L),NY),NY=1,NSND)
       END DO
C      CLOSE(LUB)
      END IF
      IF(ICON.EQ.5) THEN
       DO L=2,LA
       TOXBT=1000.*TOXB(L,KBT(L),1)
       WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),CONC(L,KC),CONC(L,1),
     $               TOXBT
       END DO
       DO L=2,LA
       TOXBT=1000.*TOXB(L,KBT(L),1)
       WRITE(LUNF,200)IL(L),JL(L),DLON(L),DLAT(L),TOXPFTW(L,KC,1),
     $               TOXPFTW(L,1,1),TOXBT
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
  200 FORMAT(2I5,1X,8E14.6)
  220 FORMAT(2I5,1X,13E11.3)
  400 FORMAT(1X,6E14.6)
  420 FORMAT(1X,13E11.3)
  250 FORMAT(12E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  220 format(2i5,1x,1p,13e12.4) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
