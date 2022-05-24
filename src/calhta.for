C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALHTA
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE CALHTA PERFORMS A HARMONIC ANALYSIS FOR THE M2 TIDE
C **  OVER TWO TIDAL CYCLES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE11,TITLE12
C
C**********************************************************************C
C
C **  INITIALIZE ON FIRST ENTRY FOR CURRENT ANALYSIS INTERVAL
C
      IF (NHAR.GT.1) GO TO 1000
C
      DO L=2,LA
      AMCP(L)=0.
      AMSP(L)=0.
      AMCUE(L)=0.
      AMSUE(L)=0.
      AMCVE(L)=0.
      AMSVE(L)=0.
      END DO
C
      DO K=1,KC
      DO L=2,LA
      AMCU(L,K)=0.
      AMSU(L,K)=0.
      AMCV(L,K)=0.
      AMSV(L,K)=0.
      END DO
      END DO
C
c     DO K=1,KS
c     DO L=2,LA
c     AMCW(L,K)=0.
c     AMSW(L,K)=0.
c     AMCAB(L,K)=0.
c     AMSAB(L,K)=0.
c     AMC2AB(L,K)=0.
c     AMS2AB(L,K)=0.
c     END DO
c     END DO
C
      OPEN(1,FILE='surfamp.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='surfamp.out',STATUS='UNKNOWN')
      TITLE1='SURFACE DISPLACEMENT AMPLITUDE CONTOURS'
      OPEN(2,FILE='surfpha.out',STATUS='UNKNOWN')
      CLOSE(2,STATUS='DELETE')
      OPEN(2,FILE='surfpha.out',STATUS='UNKNOWN')
      TITLE2='SURFACE DISPLACEMENT PHASE CONTOURS'
      OPEN(3,FILE='majaxis.out',STATUS='UNKNOWN')
      CLOSE(3,STATUS='DELETE')
      OPEN(3,FILE='majaxis.out',STATUS='UNKNOWN')
      TITLE3='MAJOR AXES OF TIDAL ELLIPSES'
      OPEN(4,FILE='majapha.out',STATUS='UNKNOWN')
      CLOSE(4,STATUS='DELETE')
      OPEN(4,FILE='majapha.out',STATUS='UNKNOWN')
      TITLE4='TIDAL ELLIPSES PHASE CONTOURS'
      OPEN(11,FILE='tidelkc.out',STATUS='UNKNOWN')
      CLOSE(11,STATUS='DELETE')
      OPEN(11,FILE='tidelkc.out',STATUS='UNKNOWN')
      TITLE11='SURFACE TIDAL ELLIPSES'
      OPEN(12,FILE='tidelkb.out',STATUS='UNKNOWN')
      CLOSE(12,STATUS='DELETE')
      OPEN(12,FILE='tidelkb.out',STATUS='UNKNOWN')
      TITLE12='BOTTOM TIDAL ELLIPSES'
C
      LINES=LA-1
      LEVELS=1
      DBS=0.
C
      WRITE (1,99) TITLE1
      WRITE (1,101)LINES,LEVELS
      WRITE (1,250)DBS
      CLOSE(1)
      WRITE (2,99) TITLE2
      WRITE (2,101)LINES,LEVELS
      WRITE (2,250)DBS
      CLOSE(2)
      WRITE (11,99) TITLE11
      WRITE (11,101)LINES,LEVELS
      WRITE (11,250)DBS
      CLOSE(11)
      DBS=99.
      WRITE (12,99) TITLE12
      WRITE (12,101)LINES,LEVELS
      WRITE (12,250)DBS
      CLOSE(12)
C
      LEVELS=2
      DBS1=0.
      DBS2=99.
C
      WRITE (3,99) TITLE3
      WRITE (3,101)LINES,LEVELS
      WRITE (3,250)DBS1,DBS2
      CLOSE(3)
      WRITE (4,99) TITLE4
      WRITE (4,101)LINES,LEVELS
      WRITE (4,250)DBS1,DBS2
      CLOSE(4)
C
C**********************************************************************C
C
C **  ACCUMULATE HARMONIC ANALYSIS
C
 1000 CONTINUE
C
      DO L=2,LA
      LN=LNC(L)
      AMCP(L)=P(L)*WC(NHAR)+AMCP(L)
      AMSP(L)=P(L)*WS(NHAR)+AMSP(L)
      UTMP1=0.5*(UHDYE(L+1)+UHDYE(L))*HPI(L)/DYP(L)
      VTMP1=0.5*(VHDXE(LN)+VHDXE(L))*HPI(L)/DXP(L)
      UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
      VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
      AMCUE(L)=UTMP*WC(NHAR)+AMCUE(L)
      AMSUE(L)=UTMP*WS(NHAR)+AMSUE(L)
      AMCVE(L)=VTMP*WC(NHAR)+AMCVE(L)
      AMSVE(L)=VTMP*WS(NHAR)+AMSVE(L)
      END DO
C
      DO K=1,KC
      DO L=2,LA
      LN=LNC(L)
      UTMP1=0.5*(U(L+1,K)+U(L,K))
      VTMP1=0.5*(V(LN,K)+V(L,K))
      UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
      VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
      AMCU(L,K)=UTMP*WC(NHAR)+AMCU(L,K)
      AMSU(L,K)=UTMP*WS(NHAR)+AMSU(L,K)
      AMCV(L,K)=VTMP*WC(NHAR)+AMCV(L,K)
      AMSV(L,K)=VTMP*WS(NHAR)+AMSV(L,K)
      END DO
      END DO
C
c     DO K=1,KS
c     DO L=2,LA
c     AMCW(L,K)=W(L,K)*WC(NHAR)+AMCW(L,K)
c     AMSW(L,K)=W(L,K)*WS(NHAR)+AMSW(L,K)
c     AMCAB(L,K)=AB(L,K)*WC(NHAR)+AMCAB(L,K)
c     AMSAB(L,K)=AB(L,K)*WS(NHAR)+AMSAB(L,K)
c     AMC2AB(L,K)=AB(L,K)*WC2(NHAR)+AMC2AB(L,K)
c     AMS2AB(L,K)=AB(L,K)*WS2(NHAR)+AMS2AB(L,K)
c     END DO
c     END DO
C
C**********************************************************************C
C
C **  CHECK FOR END OF ANALYSIS INTERVAL
C
      IF (NHAR.LT.NTSPTC2) GO TO 2000
C
C**********************************************************************C
C
C **  COMPLETE HARMONIC ANALYSIS
C
      DO L=2,LA
      AMC=AS*AMCP(L)-ACS*AMSP(L)
      AMS=-ACS*AMCP(L)+AC*AMSP(L)
      AMCP(L)=AMC
      AMSP(L)=AMS
      AMC=AS*AMCUE(L)-ACS*AMSUE(L)
      AMS=-ACS*AMCUE(L)+AC*AMSUE(L)
      AMCUE(L)=AMC
      AMSUE(L)=AMS
      AMC=AS*AMCVE(L)-ACS*AMSVE(L)
      AMS=-ACS*AMCVE(L)+AC*AMSVE(L)
      AMCVE(L)=AMC
      AMSVE(L)=AMS
      END DO
C
      DO K=1,KC
      DO L=2,LA
      AMC=AS*AMCU(L,K)-ACS*AMSU(L,K)
      AMS=-ACS*AMCU(L,K)+AC*AMSU(L,K)
      AMCU(L,K)=AMC
      AMSU(L,K)=AMS
      AMC=AS*AMCV(L,K)-ACS*AMSV(L,K)
      AMS=-ACS*AMCV(L,K)+AC*AMSV(L,K)
      AMCV(L,K)=AMC
      AMSV(L,K)=AMS
      END DO
      END DO
C
c     DO K=1,KS
c     DO L=2,LA
c     AMC=AS*AMCW(L,K)-ACS*AMSW(L,K)
c     AMS=-ACS*AMCW(L,K)+AC*AMSW(L,K)
c     AMCW(L,K)=AMC
c     AMSW(L,K)=AMS
c     AMC=AS*AMCAB(L,K)-ACS*AMSAB(L,K)
c     AMS=-ACS*AMCAB(L,K)+AC*AMSAB(L,K)
c     AMCAB(L,K)=AMC
c     AMSAB(L,K)=AMS
c     AMC=AS2*AMC2AB(L,K)-ACS2*AMS2AB(L,K)
c     AMS=-ACS2*AMC2AB(L,K)+AC2*AMS2AB(L,K)
c     AMC2AB(L,K)=AMC
c     AMS2AB(L,K)=AMS
c     END DO
c     END DO
C
      NHAR=0
C
      OPEN(1,FILE='surfamp.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (1,100)N
      OPEN(2,FILE='surfpha.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (2,100)N
C
      DO L=2,LA
      SSURFAMP=GI*SQRT(AMCP(L)*AMCP(L)+AMSP(L)*AMSP(L))
       IF (AMCP(L).EQ.0.0.AND.AMSP(L).EQ.0.0) THEN
         PHI=999999.
        ELSE
         PHI=ATAN2(AMSP(L),AMCP(L))
       END IF
      SSURFPHS=TIDALP*PHI/(3600.*PI2)
      WRITE(1,200)IL(L),JL(L),DLON(L),DLAT(L),SSURFAMP
      WRITE(2,200)IL(L),JL(L),DLON(L),DLAT(L),SSURFPHS
      END DO
C
      CLOSE(1)
      CLOSE(2)
C
      OPEN(3,FILE='majaxis.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (3,100)N
      OPEN(4,FILE='majapha.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (4,100)N
      OPEN(11,FILE='tidelkc.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (11,100)N
      OPEN(12,FILE='tidelkb.out',ACCESS='APPEND',STATUS='UNKNOWN')
      WRITE (12,100)N
C
      DO L=2,LA
      TERM1=AMCU(L,KC)+AMSV(L,KC)
      TERM2=AMCV(L,KC)-AMSU(L,KC)
      TERM3=AMCU(L,KC)-AMSV(L,KC)
      TERM4=AMCV(L,KC)+AMSU(L,KC)
      RPLUS=0.5*SQRT(TERM1*TERM1+TERM2*TERM2)
      RMINS=0.5*SQRT(TERM3*TERM3+TERM4*TERM4)
C     APLUS=ATAN2(TERM2,TERM1)
C     AMINS=ATAN2(TERM4,TERM3)
      IF(TERM1.EQ.0.0.AND.TERM2.EQ.0.0) THEN
        APLUS=999999.
       ELSE
        APLUS=ATAN2(TERM2,TERM1)
      END IF
      IF(TERM3.EQ.0.0.AND.TERM4.EQ.0.0) THEN
        AMINS=999999.
       ELSE
        AMINS=ATAN2(TERM4,TERM3)
      END IF
      RRMAJ=RPLUS+RMINS
      RRMIN=ABS(RPLUS-RMINS)
      AACCWX=0.5*(APLUS+AMINS)
      RMAJUKC=RRMAJ*COS(AACCWX)
      RMAJVKC=RRMAJ*SIN(AACCWX)
      IF(RMAJUKC.LT.0.0) THEN
        RMAJUKC=-RMAJUKC
        RMAJVKC=-RMAJVKC
      END IF
      PHASEKC=(0.25/PI)*TIDALP*(AMINS-APLUS)/3600.
      WRITE(11,200)IL(L),JL(L),DLON(L),DLAT(L),AACCWX,RRMAJ,RRMIN
      TERM1=AMCU(L,1)+AMSV(L,1)
      TERM2=AMCV(L,1)-AMSU(L,1)
      TERM3=AMCU(L,1)-AMSV(L,1)
      TERM4=AMCV(L,1)+AMSU(L,1)
      RPLUS=0.5*SQRT(TERM1*TERM1+TERM2*TERM2)
      RMINS=0.5*SQRT(TERM3*TERM3+TERM4*TERM4)
C     APLUS=ATAN2(TERM2,TERM1)
C     AMINS=ATAN2(TERM4,TERM3)
      IF(TERM1.EQ.0.0.AND.TERM2.EQ.0.0) THEN
        APLUS=999999.
       ELSE
        APLUS=ATAN2(TERM2,TERM1)
      END IF
      IF(TERM3.EQ.0.0.AND.TERM4.EQ.0.0) THEN
        AMINS=999999.
       ELSE
        AMINS=ATAN2(TERM4,TERM3)
      END IF
      RRMAJ=RPLUS+RMINS
      RRMIN=ABS(RPLUS-RMINS)
      AACCWX=0.5*(APLUS+AMINS)
      RMAJUKB=RRMAJ*COS(AACCWX)
      RMAJVKB=RRMAJ*SIN(AACCWX)
      IF(RMAJUKB.LT.0.0) THEN
        RMAJUKB=-RMAJUKB
        RMAJVKB=-RMAJVKB
      END IF
      PHASEKB=(0.25/PI)*TIDALP*(AMINS-APLUS)/3600.
      WRITE(3,200)IL(L),JL(L),DLON(L),DLAT(L),RMAJUKC,RMAJVKC,
     $              RMAJUKB,RMAJVKB
      WRITE(4,200)IL(L),JL(L),DLON(L),DLAT(L),PHASEKC,PHASEKB
      WRITE(12,200)IL(L),JL(L),DLON(L),DLAT(L),AACCWX,RRMAJ,RRMIN
      END DO
C
      CLOSE(3)
      CLOSE(4)
      CLOSE(11)
      CLOSE(12)
C
C**********************************************************************C
C
 2000 CONTINUE
C
      NHAR=NHAR+1
C
C**********************************************************************C
C
   99 FORMAT(A80)
  100 FORMAT(I10)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E13.5)
  250 FORMAT(12E12.4)
C
      RETURN
      END
