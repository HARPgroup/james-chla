C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WAVEBL
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      character*9 FNWAVE
      character*1 CFNWAVE(0:9)
C
C**********************************************************************C
C
C **  INITIALIZE AND INPUT WAVE INFORMATION
C
      IF (JSWAVE.EQ.1) GO TO 100
C
      JSWRPH=1
C
      DO L=1,LC
      HMPW(L)=0.
      HMCW(L)=0.
      HMUW(L)=0.
      HMVW(L)=0.
      WVWHA(L)=0.
C     WVACOS(L)=0.
C     WVASIN(L)=0.
      WVKHP(L)=0.
      WVKHC(L)=0.
      WVKHU(L)=0.
      WVKHV(L)=0.
      WVTMP1(L)=0.
      WVTMP2(L)=0.
      WVTMP3(L)=0.
      WVTMP4(L)=0.
      UWVMAG(L)=0.
      VWVMAG(L)=0.
      WVENEP(L)=0.
      UWVSQ(L)=0.
      QQWC(L)=1.E-12
      QQWCR(L)=1.E-12
      QQWV1(L)=1.E-12
      QQWV2(L)=1.E-12
      QQWV3(L)=1.E-12
      WACCWE(L)=0.
      END DO
C
      DO K=1,KC
      DO L=1,LC
      WVHUU(L,K)=0.
      WVHVV(L,K)=0.
      WVHUV(L,K)=0.
      WVPP(L,K)=0.
      WVPU(L,K)=0.
      WVPV(L,K)=0.
      WVDISP(L,K)=0.
      UWVRE(L,K)=0.
      UWVIM(L,K)=0.
      VWVRE(L,K)=0.
      VWVIM(L,K)=0.
      FXWAVE(L,K)=0.
      FYWAVE(L,K)=0.
      END DO
      END DO
C
      OPEN(1,FILE='wavebl.inp',STATUS='UNKNOWN')
C
      DO NSKIP=1,11
      READ(1,1,IOSTAT=ISO)
      IF(ISO.GT.0) GO TO 1081
      END DO
C
      READ(1,*,IOSTAT=ISO)NWVDAT,CVTWHA,ISWCBL,NWUPDT,NTSWV,ISDZBR
      IF(ISO.GT.0) GO TO 1082
C
      CLOSE(1)
C
      FNWAVE='wv001.inp'
C
      OPEN(1,FILE=FNWAVE)
C
      DO L=2,LA
      READ(1,*,IOSTAT=ISO)IWVH,IWDIR,IWPRD1,IWPRD2
      IF(ISO.GT.0) GO TO 1083
      HMPW(L)=HMP(L)
      WVWHA(L)=0.001*FLOAT(IWVH)
      WACCWE(L)=FLOAT(IWDIR)/57.29578
      IF(IWPRD2.GT.0) WVFRQL(L)=0.01*FLOAT(IWPRD2)
      IF(IWPRD2.LE.0) WVFRQL(L)=1.
      WVFRQL(L)=2.*PI/WVFRQL(L)
      END DO
C
      CLOSE(1)
C
      OPEN(1,FILE='wavebl.dia',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
C

      NTMP=0
      WRITE(6,666)NTMP,FNWAVE
      WRITE(8,666)NTMP,FNWAVE
C
      JSWAVE=1
      ITWCBL1=1
      ITWCBL2=0
      ITWCBL3=0
      NWCUNT=1
C
      GO TO 200
C
C**********************************************************************C
C
C ** UPDATE WAVE FIELD AT NWUPDT INTERVALS
C
  100 CONTINUE
c
      JSWRPH=0
C
      IF(JSWAVE.EQ.1) THEN
        IF(NWCUNT.EQ.NWUPDT)THEN
          NWCUNT=1
          ITWCBL1=ITWCBL1+1
          IF(ITWCBL1.GT.9) THEN
           ITWCBL1=0
           ITWCBL2=ITWCBL2+1
          END IF
          IF(ITWCBL2.GT.9) THEN
           ITWCBL2=0
           ITWCBL3=ITWCBL3+1
          END IF
          CFNWAVE(0)='0'
          CFNWAVE(1)='1'
          CFNWAVE(2)='2'
          CFNWAVE(3)='3'
          CFNWAVE(4)='4'
          CFNWAVE(5)='5'
          CFNWAVE(6)='6'
          CFNWAVE(7)='7'
          CFNWAVE(8)='8'
          CFNWAVE(9)='9'
          FNWAVE='wv' // CFNWAVE(ITWCBL3) // CFNWAVE(ITWCBL2)
     $                // CFNWAVE(ITWCBL1)// '.inp'
        ELSE
          NWCUNT=NWCUNT+1
          RETURN
        END IF
      END IF
C
      OPEN(1,FILE=FNWAVE)
C
      DO L=2,LA
      READ(1,*,IOSTAT=ISO)IWVH,IWDIR,IWPRD1,IWPRDP
      IF(ISO.GT.0) GO TO 1083
      HMPW(L)=HMP(L)
      WVWHA(L)=0.5*(0.001*FLOAT(IWVH)+WVWHA(L))
      WACCWE(L)=FLOAT(IWDIR)/57.29578
      IF(IWPRDP.GT.0) WVFRQL(L)=0.01*FLOAT(IWPRDP)
      IF(IWPRDP.LE.0) WVFRQL(L)=1.
      WVFRQL(L)=2.*PI/WVFRQL(L)
      END DO
C
      CLOSE(1)
C
      WRITE(6,666)N,FNWAVE
      WRITE(8,666)N,FNWAVE
C
  666 FORMAT(' UPDATED WAVE FIELD N,FNWAVE = ',I12,A12)
C
C**********************************************************************C
C
C **  GENERATE WAVE TABLE
C
  200 CONTINUE
C
      HMXTMP=0.
      WVFRQM=0.
      DO L=2,LA
       HMXTMP=MAX(HMXTMP,HMPW(L))
       IF(WVWHA(L).GT.0.) WVFRQM=MAX(WVFRQM,WVFRQL(L))
      END DO
C
      FKHMAX=1.5*GI*WVFRQM*WVFRQM*HMXTMP
      RKHTMP=0.001
   10 CONTINUE
      FKHTMP=RKHTMP*TANH(RKHTMP)
      IF (FKHTMP.LT.FKHMAX) THEN
        RKHTMP=2.*RKHTMP
        GO TO 10
       ELSE
        DKH=RKHTMP/1000.
      END IF
C
      RKHTAB(1)=0.
      FUNKH(1)=0.
      DO NKH=2,1001
       RKHTAB(NKH)=RKHTAB(NKH-1)+DKH
       FUNKH(NKH)=RKHTAB(NKH)*TANH(RKHTAB(NKH))
      END DO
C
      DO L=2,LA
       HFFDG=GI*WVFRQL(L)*WVFRQL(L)*HMPW(L)
       WVKHP(L)=1.
       IF(WVWHA(L).GT.0.) WVKHP(L)=VALKH(HFFDG)
      END DO
C
      IF(JSWRPH.EQ.1)THEN
        OPEN(1,FILE='wvtab.out',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='wvtab.out',STATUS='UNKNOWN')
        DO NKH=1,1001
         WRITE(1,111)RKHTAB(NKH),FUNKH(NKH)
        END DO
        CLOSE(1)
      END IF
C
      GO TO 400
C
 1081 WRITE(6,1091)
      STOP
 1082 WRITE(6,1092)
      STOP
c 1083 WRITE(6,1093) NWV ! bug, NWV used before assigned a value, Ji, 10/10/00
 1083 WRITE(6,1093) IWVH
      STOP
 1084 WRITE(6,1094) NWV
      STOP
C
    1 FORMAT(120X)
 1091 FORMAT('  READ ERROR ON FILE wave.inp , HEADER')
 1092 FORMAT('  READ ERROR ON FILE wave.inp , 1ST DATA')
 1093 FORMAT('  READ ERROR ON FILE wave.inp , 2ND DATA, NWV = ',I5)
 1094 FORMAT('  READ ERROR ON FILE wave.inp , 3RD DATA, NWV = ',I5)
  111 FORMAT(2E14.4)
C
C**********************************************************************C
C
  400 CONTINUE
C
      DO L=2,LA
       IF(HMP(L).LT.0.55) WVWHA(L)=0. ! Ji, Hardwired by John
       IF(MVEGL(L).NE.MVEGOW) WVWHA(L)=0. !Ji, Hardwired by John,
cJi, It makes only deep water and mvegl(l)=1 has wind wave, 1/29/00
       IWVRDC=0
       IF(MVEGL(L).EQ.MVEGOW) THEN
         IF(MVEGL(L-1).NE.MVEGOW) IWVRDC=1
         IF(MVEGL(L+1).NE.MVEGOW) IWVRDC=1
         IF(MVEGL(LSC(L)).NE.MVEGOW) IWVRDC=1
         IF(MVEGL(LNC(L)).NE.MVEGOW) IWVRDC=1
       END IF
       IF(IWVRDC.GT.0) WVWHA(L)=0.5*WVWHA(L)
      END DO
C
C**********************************************************************C
C
C **  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING
C **  THE WAVE TURBULENT INTENSITY, QQWV
C **  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE
C
       OPEN(1,FILE='wavebl.dia')
       CLOSE(1,STATUS='DELETE')
       OPEN(1,FILE='wavebl.dia')
C
      DO L=2,LA
       AEXTMP=0.5*WVWHA(L)/SINH(WVKHP(L))
       UWORBIT=AEXTMP*WVFRQL(L)
       UWVSQ(L)=UWORBIT*UWORBIT
       VISMUD=1.E-6
       IF(ISMUD.GE.1) VISMUD=CSEDVIS(SED(L,1,1))
       DELWAVE=SQRT(VISMUD/WVFRQL(L))
       REYWAVE=UWORBIT*AEXTMP/VISMUD
       QQWV1(L)=0.
       QQWV2(L)=0.
C
C ** LAMINAR WAVE BOUNDARY LAYER
C
       IF(REYWAVE.LT.5.E5) THEN
         CDTMP=0.
         IF(REYWAVE.GT.0.) CDTMP=1./SQRT(REYWAVE)
         QQWV2(L)=CDTMP*UWORBIT*UWORBIT
       END IF
C
C ** TURBULENT WAVE BOUNDARY LAYER
C
       IF(REYWAVE.GE.5.E5) THEN
C
C ** TURBULENT SMOOTH WAVE BOUNDARY LAYER
C
         IF(ZBR(L).LE.0.) THEN
           CDTMP=0.012/(REYWAVE**0.123)
           QQWV1(L)=CDTMP*UWORBIT*UWORBIT
C
C ** TURBULENT ROUGH WAVE BOUNDARY LAYER
C
          ELSE
	exp=0.0 ! bug, Ji, 10/10/00
           TMPVAL=30.*ZBR(L)/AEXTMP
           TMPVAL=5.5*(TMPVAL**0.2)-6.3
           CDTMP=0.5*EXP*TMPVAL   !Bug, Ji, exp is no defined!, 3/5/00
           QQWV1(L)=CDTMP*UWORBIT*UWORBIT
           ZBRE(L)=ZBR(L)
c          IF(QQ(L,O).GT.0.) THEN  ! Bug, Ji, o->0, 3/5/00
           IF(QQ(L,0).GT.0.) THEN
             TMPVAL=UWORBIT*SQRT( AEXTMP/(30.*ZBR(L)) )
             USTARC=SQRT(QQ(L,0)/CTURB2)
             TMPVAL=TMPVAL/USTARC
             ZBRE(L)=ZBR(L)*(1.+0.19*TMPVAL)
           END IF
         END IF
       END IF
C
c       TVAL1=SINH(WVKHP(L))
c       TVAL2=CSEDVIS(SED(L,1,1))
        WRITE(1,600)L,IL(L),JL(L),WVWHA(L),WVFRQL(L),AEXTMP,UWORBIT,
     $             VISMUD,REYWAVE,CDTMP,QQWV1(L),QQWV2(L)
      END DO
C
       CLOSE(1)
C
  600 FORMAT(3I5,9E12.4)
C
c      DO L=2,LA
c      AEXTMP=0.5*WVWHA(L)/SINH(WVKHP(L))
c      UWVSQ(L)=AEXTMP*AEXTMP*WVFRQ*WVFRQ
c      CDRGTMP=(30.*ZBR(L)/AEXTMP)**0.2
c      CDRGTMP=5.57*CDRGTMP-6.13
c      CDRGTMP=EXP(CDRGTMP)
c      CDRGTMP=MIN(CDRGTMP,0.22)
c      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
c      QQWV1(L)=CTURB2*TMPVAL
c      TAUTMP=TMPVAL/TAUR(NSED+1)
c      CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
c      ZBRE(L)=CORZBR*ZBR(L)
c      CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2
c      CDRGTMP=5.57*CDRGTMP-6.13
c      CDRGTMP=EXP(CDRGTMP)
c      CDRGTMP=MIN(CDRGTMP,0.22)
c      TMPVAL=0.5*CDRGTMP*UWVSQ(L)
c      QQWV2(L)=CTURB2*TMPVAL
c      IF(ISTRAN(7).EQ.0) QQWV2(L)=0.
c      END DO
C
C      IF (ISRESTI.NE.0) THEN
C       OPEN(1,FILE='wvqwcp.inp',STATUS='UNKNOWN')
C       DO L=2,LA
C       READ(1,*)IDUM,JDUM,QQWV1(L),QQWV2(L),QQWV2(L),QQWC(L),QQWCR(L)
C       END DO
C      END IF
C
C**********************************************************************C
C
C      IF (ISWAVE.GE.1) CALL WRSPLTH
C
C**********************************************************************C
C
      RETURN
      END
