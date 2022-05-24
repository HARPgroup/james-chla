C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE VELPLTH
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE VELPLTH WRITES A HORIZONTAL INSTANTANEOUS VELOCITY
C **  VECTOR FILE 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      REAL DBS(10)
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7
C
C**********************************************************************C
C
      IF (JSVPH.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      TITLE1='INSTANTANEOUS HORIZ VELOCITY cm/s '
      TITLE2='INSTANTANEOUS BOTTOM STRESS cm2/s2'
      TITLE3='SED TRANS (SURF OR BOT LAYER) kg/s'
      TITLE4='DEPTH INTEGRAED SED TRANS kg/s'
      TITLE5='EFFECTIVE BOTTOM ROUGHNESS cm'
      TITLE6='CURRENT SHEAR VELOCITY cm/s'
      TITLE7='WAVE-CURRENT SHEAR VELOCITY cm/s'
      IF(ISVPH.EQ.1) LINES1=LA-1
      IF(ISVPH.EQ.2) LINES1=NRC
      IF(ISVPH.EQ.3) LINES1=NBC
      LEVELS=2
      LEVELT=1
      DBS(1)=0.
      DBS(2)=99.
C
      OPEN(10,FILE='velvech.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='velvech.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE1
      WRITE (10,101)LINES1,LEVELS
      WRITE (10,250)(DBS(L),L=1,LEVELS)
      CLOSE(10)
C
      OPEN(10,FILE='tauvech.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='tauvech.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE2
      WRITE (10,101)LINES1,LEVELT
      WRITE (10,250)(DBS(L),L=1,LEVELT)
      CLOSE(10)
C
      OPEN(10,FILE='stvech.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='stvech.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE3
      WRITE (10,101)LINES1,LEVELS
      WRITE (10,250)(DBS(L),L=1,LEVELS)
      CLOSE(10)
C
      OPEN(10,FILE='tstvech.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='tstvech.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE4
      WRITE (10,101)LINES1,LEVELT
      WRITE (10,250)(DBS(L),L=1,LEVELT)
      CLOSE(10)
C
      OPEN(10,FILE='zbreffh.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='zbreffh.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE5
      WRITE (10,101)LINES1,LEVELT
      WRITE (10,250)(DBS(L),L=1,LEVELT)
      CLOSE(10)
C
      OPEN(10,FILE='ccustrh.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='ccustrh.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE6
      WRITE (10,101)LINES1,LEVELT
      WRITE (10,250)(DBS(L),L=1,LEVELT)
      CLOSE(10)
C
      OPEN(10,FILE='wcustrh.out',STATUS='UNKNOWN')
      CLOSE(10,STATUS='DELETE')
      OPEN(10,FILE='wcustrh.out',STATUS='UNKNOWN')
      WRITE (10,99) TITLE7
      WRITE (10,101)LINES1,LEVELT
      WRITE (10,250)(DBS(L),L=1,LEVELT)
      CLOSE(10)
C
      JSVPH=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE 
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON 
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014    
C
      OPEN(10,FILE='velvech.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (10,100)N,TIME
      OPEN(11,FILE='tauvech.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (11,100)N,TIME
      OPEN(12,FILE='stvech.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (12,100)N,TIME
      OPEN(13,FILE='tstvech.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (13,100)N,TIME
      OPEN(14,FILE='zbreffh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (14,100)N,TIME
      OPEN(15,FILE='ccustrh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (15,100)N,TIME
      OPEN(16,FILE='wcustrh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (16,100)N,TIME
C
      QBOTTMP=100./CTURB3
      IF(ISVPH.EQ.1) THEN
       DO L=2,LA
       LN=LNC(L)
       UTMP=50.*STCUV(L)*(U(L+1,KC)+U(L,KC))
       VTMP=50.*STCUV(L)*(V(LN,KC)+V(L,KC))
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(U(L+1,1)+U(L,1))
       VTMP=50.*STCUV(L)*(V(LN,1)+V(L,1))
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       IF(KC.GT.1)WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),
     $              VELEKC,VELNKC,VELEKB,VELNKB
       IF(KC.GE.1)WRITE(10,200)IL(L),JL(L),VELEKB,VELNKB
       IF(ISTRAN(4).GE.1) THEN
         UTMP=5000.*STCUV(L)*(TBX(L+1)+TBX(L))
         VTMP=5000.*STCUV(L)*(TBY(LN )+TBY(L))
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,KC)
C    $                        +UHDY(L  ,KC))*SED(L,KC)*DZC(KC)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,KC)
C    $                        +VHDX(L  ,KC))*SED(L,KC)*DZC(KC)
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,1)+UHDY(L,1))*SED(L,1)*DZC(1)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,1)+VHDX(L,1))*SED(L,1)*DZC(1)
         VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC,
     $                VELEKB,VELNKB
         VTMP=0.0
         UTMP=0.0
         DO K=1,KC
C          UTMP=UTMP+0.0005*STCUV(L)*(UHDY(L+1,K)+UHDY(L,K))
C     $                             *SED(L,K)*DZC(K)
C          VTMP=VTMP+0.0005*STCUV(L)*(VHDX(LN ,K)+VHDX(L,K))
C     $                             *SED(L,K)*DZC(K)
         END DO
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
       END IF
       IF (ISWAVE.GE.2) THEN
         ZBREFF=100.*ZBRE(L)
         WRITE(14,200)IL(L),JL(L),DLON(L),DLAT(L),ZBREFF
         QTURBC=QBOTTMP*SQRT(QQ(L,0))
         WRITE(15,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
         QTURBC=QBOTTMP*QQWV2(L)
         WRITE(16,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
       END IF
       END DO
      END IF
C
      IF(ISVPH.EQ.2) THEN 
       DO LR=1,NRC
       L=LRC(LR)
       LN=LNC(L)
       UTMP=50.*STCUV(L)*(U(L+1,KC)+U(L,KC))
       VTMP=50.*STCUV(L)*(V(LN,KC)+V(L,KC))
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(U(L+1,1)+U(L,1))
       VTMP=50.*STCUV(L)*(V(LN,1)+V(L,1))
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC,
     $              VELEKB,VELNKB
       IF(ISTRAN(4).GE.1) THEN
         UTMP=5000.*STCUV(L)*(TBX(L+1)+TBX(L))
         VTMP=5000.*STCUV(L)*(TBY(LN )+TBY(L))
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,KC)
C    $                        +UHDY(L  ,KC))*SED(L,KC)*DZC(KC)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,KC)
C    $                        +VHDX(L  ,KC))*SED(L,KC)*DZC(KC)
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,1)+UHDY(L,1))*SED(L,1)*DZC(1)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,1)+VHDX(L,1))*SED(L,1)*DZC(1)
         VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC,
     $                VELEKB,VELNKB
         VTMP=0.0
         UTMP=0.0
         DO K=1,KC
C          UTMP=UTMP+0.0005*STCUV(L)*(UHDY(L+1,K)+UHDY(L,K))
C     $                             *SED(L,K)*DZC(K)
C          VTMP=VTMP+0.0005*STCUV(L)*(VHDX(LN ,K)+VHDX(L,K))
C     $                             *SED(L,K)*DZC(K)
         END DO
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
       END IF
       IF (ISWAVE.GE.2) THEN
         ZBREFF=100.*ZBRE(L)
         WRITE(14,200)IL(L),JL(L),DLON(L),DLAT(L),ZBREFF
         QTURBC=QBOTTMP*SQRT(QQ(L,0))
         WRITE(15,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
         QTURBC=QBOTTMP*QQWV2(L)
         WRITE(16,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
       END IF
       END DO
      END IF
C
      IF(ISVPH.EQ.3) THEN
       DO LB=1,NBC
       L=LBC(LB)
       LN=LNC(L)
       UTMP=50.*STCUV(L)*(U(L+1,KC)+U(L,KC))
       VTMP=50.*STCUV(L)*(V(LN,KC)+V(L,KC))
       VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
       UTMP=50.*STCUV(L)*(U(L+1,1)+U(L,1))
       VTMP=50.*STCUV(L)*(V(LN,1)+V(L,1))
       VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
       VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
       WRITE(10,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC,
     $              VELEKB,VELNKB
       IF(ISTRAN(4).GE.1) THEN
         UTMP=5000.*STCUV(L)*(TBX(L+1)+TBX(L))
         VTMP=5000.*STCUV(L)*(TBY(LN )+TBY(L))
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,KC)
C    $                        +UHDY(L  ,KC))*SED(L,KC)*DZC(KC)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,KC)
C    $                        +VHDX(L  ,KC))*SED(L,KC)*DZC(KC)
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
C        UTMP=0.0005*STCUV(L)*(UHDY(L+1,1)+UHDY(L,1))*SED(L,1)*DZC(1)
C        VTMP=0.0005*STCUV(L)*(VHDX(LN ,1)+VHDX(L,1))*SED(L,1)*DZC(1)
         VELEKB=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKB=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC,
     $                VELEKB,VELNKB
         VTMP=0.0
         UTMP=0.0
         DO K=1,KC
C          UTMP=UTMP+0.0005*STCUV(L)*(UHDY(L+1,K)+UHDY(L,K))
C     $                             *SED(L,K)*DZC(K)
C          VTMP=VTMP+0.0005*STCUV(L)*(VHDX(LN ,K)+VHDX(L,K))
C     $                             *SED(L,K)*DZC(K)
         END DO
         VELEKC=CUE(L)*UTMP+CVE(L)*VTMP
         VELNKC=CUN(L)*UTMP+CVN(L)*VTMP
         WRITE(13,200)IL(L),JL(L),DLON(L),DLAT(L),VELEKC,VELNKC
       END IF
       IF (ISWAVE.GE.2) THEN
         ZBREFF=100.*ZBRE(L)
         WRITE(14,200)IL(L),JL(L),DLON(L),DLAT(L),ZBREFF
         QTURBC=QBOTTMP*SQRT(QQ(L,0))
         WRITE(15,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
         QTURBC=QBOTTMP*QQWV2(L)
         WRITE(16,200)IL(L),JL(L),DLON(L),DLAT(L),QTURBC
       END IF
       END DO
      END IF
C
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
