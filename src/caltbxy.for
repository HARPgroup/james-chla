C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALTBXY
C
C**********************************************************************C
C
C **  SUBROUTINE CALTBXY CALCULATES BOTTOM FRICTION OR DRAG
C **  COEFFICIENTS IN QUADRATIC LAW FORM REFERENCED TO NEAR
C **  BOTTOM OR DEPTH AVERAGED HORIZONTAL VELOCITIES
C **  FOR VEGETATION RESISTANCE IN DEPTH INTEGRATED FLOW
C **  THE COEFFICIENT REPRESENTS BOTTOM AND WATER COLUMN VEGETATION
C **  RESISTANCE
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C **  IF WAVE-CURRENT BBL MODEL IS ACTIVE, GO TO WAVE CURRENT BBL
C
      IF(ISWCBL.GE.1) GO TO 1947
C
C**********************************************************************C
C
C **  INITIALIZE IMPLICIT BOTTOM FRICTION AND SET DIAGNOSTIC FILES
C **  ON FIRST CALL
C
      IF (JSTBXY.EQ.1) GO TO 100
C
      IF(ISITB.GE.1) THEN
        IF(ISITB.EQ.1) THEN
          RITB1=0.5
          RITB=0.5
          CDLIMIT=1.
         ELSE
          RITB1=0.0
          RITB=1.0
          CDLIMIT=10.
        END IF
       ELSE
        RITB1=1.0
        RITB=0.0
        CDLIMIT=0.5
      END IF
C
C       CDMAXU=CDLIMIT*H1U(L)/( DT*UMAGTMP )
C       CDMAXV=CDLIMIT*H1V(L)/( DT*VMAGTMP )
C
      IF(ISVEG.GE.2) THEN
       OPEN(1,FILE='cbot.log',STATUS='UNKNOWN')
       CLOSE(1,STATUS='DELETE')
      END IF
C
      DO L=2,LA
       STBXO(L)=STBX(L)
       STBYO(L)=STBY(L)
      END DO
C
      N=-2
      JSTBXY=1
C
  100 CONTINUE
C
      IF(ISITB.GE.1) THEN
        IF(ISITB.EQ.1) THEN
          CDLIMIT=1.
         ELSE
          CDLIMIT=10.
        END IF
       ELSE
        CDLIMIT=0.5
      END IF
C
C**********************************************************************C
C
C **  NORMAL ENTRY INTO STANDARD AND VEGE RESISTANCE CALCULATION
C
      IF(ISVEG.GE.2) THEN
       OPEN(1,FILE='cbot.log',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      CDTOTUM=0.
      CDTOTVM=0.
      CDMAXUM=0.
      CDMAXVM=0.
      IF(ISVEG.EQ.0) UVEGSCL=1.E-12
C
      IF (KC.GT.1) GO TO 200
C
C ** For Single Layer Model setup:
C
       DO L=2,LA
C      STBXO(L)=STBX(L)
C      STBYO(L)=STBY(L)
C      STBX(L)=STBX(L)*.16/((LOG(HMU(L)/ZBR(L))-1.)**2)
C      STBY(L)=STBY(L)*.16/((LOG(HMV(L)/ZBR(L))-1.)**2)
       UMAGTMP=SQRT( U(L,1)*U(L,1)+VU(L)*VU(L) )
       VMAGTMP=SQRT( UV(L)*UV(L)+V(L,1)*V(L,1) )
       IF(N.EQ.-2) THEN
        UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L) )
        VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1) )
       END IF
       UMAGTMP=MAX(UMAGTMP,UVEGSCL)
       VMAGTMP=MAX(VMAGTMP,UVEGSCL)
       UMAGTMP=MAX(UMAGTMP,1.E-12)
       VMAGTMP=MAX(VMAGTMP,1.E-12)
c       CDMAXU=H1U(L)/( 4.*DT*UMAGTMP )
c       CDMAXV=H1V(L)/( 4.*DT*VMAGTMP )
       CDMAXU=CDLIMIT*H1U(L)/( DT*UMAGTMP )
       CDMAXV=CDLIMIT*H1V(L)/( DT*VMAGTMP )
       CDMAXUM=MAX(CDMAXUM,CDMAXU)
       CDMAXVM=MAX(CDMAXVM,CDMAXV)
       CDVEGU=0.
       CDVEGV=0.
       CDBOTU=0.
       CDBOTV=0.
C      IF(ISRESTI.EQ.0) GO TO 7777
       IF(ISRESTI.EQ.0.AND.N.LE.1) GO TO 7777
       IF(ISVEG.GE.1) THEN
         LS=LSC(L)
         M=MVEGL(L)
         MW=MVEGL(L-1)
         MS=MVEGL(LS)
C        IF(M.EQ.MVEGOW) GO TO 7777
         RVEGUL=0.
         RVEGVL=0.
         RVEGUM=1.
         RVEGVM=1.
         CPVEGU=0.5
         IF(ISVEGL.EQ.1) CPVEGU=CPVEGU + 10.E-6/(
     $                   (BPVEG(MW)+BPVEG(M))*UMAGTMP )
         IF (CPVEGU.GT.1.0) THEN
C          CALCULATE R FOR LAMINAR FLOW
           CPVEGU=CPVEGU-0.5
           RVEGUL=SQRT( (2.5*SCVEG(M)*HU(L)*HU(L)*RDLPSQ(M)/PVEGZ(M))
     $        +(2.5*SCVEG(MW)*HU(L-1)*HU(L-1)*RDLPSQ(MW)/PVEGZ(MW)) )
           IF(N.EQ.-2) THEN
           RVEGUL=SQRT( (2.5*SCVEG(M)*H1U(L)*H1U(L)*RDLPSQ(M)/PVEGZ(M))
     $        +(2.5*SCVEG(MW)*H1U(L-1)*H1U(L-1)*RDLPSQ(MW)/PVEGZ(MW)) )
           END IF
           RVEGUM=0.
         END IF
         CPVEGU=SCVEG(M)*CPVEGU
         CPVEGV=0.5
         IF(ISVEGL.EQ.1) CPVEGV=CPVEGV + 10.E-6/(
     $                    (BPVEG(MS)+BPVEG(M))*VMAGTMP )
         IF (CPVEGV.GT.1.0) THEN
C          CALCULATE R FOR LAMINAR FLOW
           CPVEGV=CPVEGV-0.5
           RVEGVL=SQRT( (2.5*SCVEG(M)*HV(L)*HV(L)*RDLPSQ(M)/PVEGZ(M))
     $          +(2.5*SCVEG(MS)*HV(LS)*HV(LS)*RDLPSQ(MS)/PVEGZ(MS)) )
           IF(N.EQ.-2) THEN
           RVEGVL=SQRT( (2.5*SCVEG(M)*H1V(L)*H1V(L)*RDLPSQ(M)/PVEGZ(M))
     $          +(2.5*SCVEG(MS)*H1V(LS)*H1V(LS)*RDLPSQ(MS)/PVEGZ(MS)) )
           END IF
           RVEGVM=0.
         END IF
         CPVEGV=SCVEG(M)*CPVEGV
         CPTMPU=0.5*CPVEGU*( (BDLPSQ(M)*HU(L)/PVEGZ(M))
     $                      +(BDLPSQ(MW)*HU(L-1)/PVEGZ(MW)) )
         CPTMPV=0.5*CPVEGV*( (BDLPSQ(M)*HV(L)/PVEGZ(M))
     $                      +(BDLPSQ(MS)*HV(LS)/PVEGZ(MS)) )
         RVEGU=1.41*( ( 0.5*HU(L)/BPVEG(M)
     $                 +0.5*HU(L-1)/BPVEG(MW) )**.6667 )
         RVEGV=1.41*( ( 0.5*HV(L)/BPVEG(M)
     $                 +0.5*HV(LS)/BPVEG(MS) )**.6667 )
         IF(N.EQ.-2) THEN
          CPTMPU=0.5*CPVEGU*( (BDLPSQ(M)*H1U(L)/PVEGZ(M))
     $                      +(BDLPSQ(MW)*H1U(L-1)/PVEGZ(MW)) )
          CPTMPV=0.5*CPVEGV*( (BDLPSQ(M)*H1V(L)/PVEGZ(M))
     $                      +(BDLPSQ(MS)*H1V(LS)/PVEGZ(MS)) )
          RVEGU=1.41*( ( 0.5*H1U(L)/BPVEG(M)
     $                 +0.5*H1U(L-1)/BPVEG(MW) )**.6667 )
          RVEGV=1.41*( ( 0.5*H1V(L)/BPVEG(M)
     $                 +0.5*H1V(LS)/BPVEG(MS) )**.6667 )
         END IF
         RVEGU=RVEGU*( CPTMPU**.3333)
         RVEGU=RVEGUM*RVEGU+RVEGUL
         RVEGV=RVEGV*( CPTMPV**.3333)
         RVEGV=RVEGVM*RVEGV+RVEGVL
         FRVEGU=RVEGU/( RVEGU-TANH(RVEGU) )
         FRVEGV=RVEGV/( RVEGV-TANH(RVEGV) )
         CDVEGU=CPTMPU*FRVEGU
         CDVEGV=CPTMPV*FRVEGV
         GO TO 7778
       END IF
 7777  CONTINUE
       HUDZBR=H1U(L)/ZBR(L)
       IF(HUDZBR.LT.7.5) HUDZBR=7.5
       HVDZBR=H1V(L)/ZBR(L)
       IF(HVDZBR.LT.7.5) HVDZBR=7.5
       CDBOTU=.16/( (LOG( HUDZBR ) -1.)**2)
       CDBOTV=.16/( (LOG( HVDZBR ) -1.)**2)
 7778  CONTINUE
       CDTOTU=CDBOTU+CDVEGU
       CDTOTV=CDBOTV+CDVEGV
       IF(ISVEG.EQ.2) THEN
         IF(CDTOTU.GT.CDMAXU) THEN
           IF(RVEGUM.EQ.1.) THEN
             WRITE(1,1717)N,IL(L),JL(L),CDTOTU,CDMAXU
            ELSE
             WRITE(1,1727)N,IL(L),JL(L),CDTOTU,CDMAXU
           END IF
         END IF
         IF(CDTOTV.GT.CDMAXV) THEN
           IF(RVEGVM.EQ.1.) THEN
             WRITE(1,1718)N,IL(L),JL(L),CDTOTV,CDMAXV
            ELSE
             WRITE(1,1728)N,IL(L),JL(L),CDTOTV,CDMAXV
           END IF
         END IF
       END IF
       CDTOTUM=MAX(CDTOTU,CDTOTUM)
       CDTOTVM=MAX(CDTOTV,CDTOTVM)
c mrm  begin special case for Christina project:
c mrm  13 Dec 1998 M.R. Morton
c mrm       CDTOTU=MIN(CDTOTU,CDMAXU)
c mrm       CDTOTV=MIN(CDTOTV,CDMAXV)
c mrm  If zbr rougness is gt 1.00, zbr value will be directly assigned
c mrm  to the Cd drag coefficient.  This was done to get subcritical
c mrm  flows in steep channels.  This change will have no impact on
c mrm  the tidal portion of the model as long as zbr is < 1.00.
       if (zbr(L) .le. 1.00) then
         CDTOTU=MIN(CDTOTU,CDMAXU)
         CDTOTV=MIN(CDTOTV,CDMAXV)
       else
         cdtotu = zbr(L) - 1.0
         cdtotv = zbr(L) - 1.0
       end if
c mrm  end special case
       STBX(L)=STBXO(L)*CDTOTU
       STBY(L)=STBYO(L)*CDTOTV
       END DO
C
       IF(ISVEG.GE.2) WRITE(1,1719)N,CDTOTUM,CDTOTVM
       IF(ISVEG.GE.2) WRITE(1,1729)N,CDMAXUM,CDMAXVM
C
       GO TO 300
C
C **  MULTIPLE LAYERS
C
  200 CONTINUE
C
       DO L=2,LA
       UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )
       VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )
C       CDMAXU=STBXO(L)*H1U(L)/( 4.*DT*UMAGTMP )
C       CDMAXV=STBYO(L)*H1V(L)/( 4.*DT*VMAGTMP )
       CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DT*UMAGTMP )
       CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DT*VMAGTMP )
       HURTMP=MAX(ZBR(L),H1U(L))
       HVRTMP=MAX(ZBR(L),H1V(L))
       DZHUDZBR=1.+0.5*DZC(1)*HURTMP/ZBR(L)
       DZHVDZBR=1.+0.5*DZC(1)*HVRTMP/ZBR(L)
c       DZHUDZBR=0.5*DZC(1)*H1U(L)/ZBR(L)
c       DZHVDZBR=0.5*DZC(1)*H1V(L)/ZBR(L)
c       DZHUDZBR=MAX(DZHUDZBR,1.35)
c       DZHVDZBR=MAX(DZHVDZBR,1.35)
       STBX(L)=AVCON*STBXO(L)*.16/((LOG(DZHUDZBR))**2)
       STBY(L)=AVCON*STBYO(L)*.16/((LOG(DZHVDZBR))**2)
       STBX(L)=MIN(CDMAXU,STBX(L))
       STBY(L)=MIN(CDMAXV,STBY(L))
       END DO
C
       IF(N.EQ.-2) THEN
       DO L=2,LA
C      STBXO(L)=STBX(L)
C      STBYO(L)=STBY(L)
C      STBX(L)=AVCON*STBX(L)*.16/((LOG(0.5*DZC(1)*HMU(L)/ZBR(L)))**2)
C      STBY(L)=AVCON*STBY(L)*.16/((LOG(0.5*DZC(1)*HMV(L)/ZBR(L)))**2)
       UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )
       VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )
C       CDMAXU=STBXO(L)*H1U(L)/( 4.*DT*UMAGTMP )
C       CDMAXV=STBYO(L)*H1V(L)/( 4.*DT*VMAGTMP )
       CDMAXU=CDLIMIT*STBXO(L)*H1U(L)/( DT*UMAGTMP )
       CDMAXV=CDLIMIT*STBYO(L)*H1V(L)/( DT*VMAGTMP )
       HURTMP=MAX(ZBR(L),H1U(L))
       HVRTMP=MAX(ZBR(L),H1V(L))
       DZHUDZBR=1.+0.5*DZC(1)*HURTMP/ZBR(L)
       DZHVDZBR=1.+0.5*DZC(1)*HVRTMP/ZBR(L)
c       DZHUDZBR=0.5*DZC(1)*H1U(L)/ZBR(L)
c       DZHVDZBR=0.5*DZC(1)*H1V(L)/ZBR(L)
c       DZHUDZBR=MAX(DZHUDZBR,1.35)
c       DZHVDZBR=MAX(DZHVDZBR,1.35)
       STBX(L)=AVCON*STBXO(L)*.16/((LOG(DZHUDZBR))**2)
       STBY(L)=AVCON*STBYO(L)*.16/((LOG(DZHVDZBR))**2)
       STBX(L)=MIN(CDMAXU,STBX(L))
       STBY(L)=MIN(CDMAXV,STBY(L))
       END DO
       END IF
C
  300 CONTINUE
C
      IF(ISVEG.GE.2) CLOSE(1)
C
      GO TO 1948
C
C**********************************************************************C
C
C **  ENTER HERE FOR WAVE-CURRENT BOUNDARY LAYER
C
 1947 CONTINUE
C
      IF(JSTBXY.EQ.0) THEN
        DO L=2,LA
        STBXO(L)=STBX(L)
        STBYO(L)=STBY(L)
        END DO
        N=0
        JSTBXY=1
        IF(ISDZBR.GE.1) THEN
         OPEN(1,FILE='zbremx.out',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
        END IF
      END IF
C
C      COSWC=1.
c      N=0                             !comm out in second call
c      IF (N.LT.NTSWV) THEN            !comm out in second call
c        WVFACT=FLOAT(N)/FLOAT(NTSWV)  !comm out in second call
c       ELSE                           !comm out in second call
c        WVFACT=1.0                    !comm out in second call
c      END IF                          !comm out in second call
C
       IF(ISDZBR.EQ.N) THEN
         OPEN(1,FILE='cddiag.out',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE='cddiag.out',STATUS='UNKNOWN')
       END IF
C
       NTMP=MAX(N,1)
       IF (NTMP.LT.NTSWV) THEN
         TMPVALW=FLOAT(NTMP)/FLOAT(NTSWV)
         WVFACT=0.5-0.5*COS(PI*TMPVALW)
        ELSE
         WVFACT=1.0
       END IF
C
       DO L=2,LA
       QQWCTMP=SQRT( QQWV2(L)*QQWV2(L)+QQ(L,0)*QQ(L,0) )
       TWCTMP=QQWCTMP/CTURB2
       TAUTMP=TWCTMP/TAUR(NSED+1)
       CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
       ZBRE(L)=CORZBR*ZBR(L)
       AEXTMP=WVWHA(L)/SINH(WVKHP(L))
       CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2
       CDRGTMP=5.57*CDRGTMP-6.13
       CDRGTMP=EXP(CDRGTMP)
       CDRGTMP=MIN(CDRGTMP,0.22)
       TAUTMP=0.5*CDRGTMP*UWVSQ(L)
       QQWV2(L)=CTURB2*TAUTMP
       QQWC(L)=SQRT( QQWV2(L)*QQWV2(L)+QQ(L,0)*QQ(L,0) )
       TWCTMP=QQWC(L)/CTURB2
       TAUB=QQWV1(L)/CTURB2
       TAUE=TWCTMP/TAUN(NSED+1)
       RIPAMP=0.
       RIPSTP=0.
       IF(TAUB.GT.TAUN(NSED+1).AND.TAUB.LE.TAUD(NSED+1)) THEN
         RIPAMP=0.22/(TAUE**0.16)
         RIPSTP=0.16/(TAUE**0.04)
       END IF
       IF(TAUB.GT.TAUD(NSED+1)) THEN
         RIPAMP=0.78/(TAUE**1.5)
         RIPSTP=0.41/TAUE
       END IF
       RIPAMP=RIPAMP*WVWHA(L)/SINH(WVKHP(L))
       TMPVAL=0.
       IF (RIPAMP.GT.0.) TMPVAL=LOG(RIPAMP/ZBRE(L))-1.
       TMPVAL=MAX(TMPVAL,0.)
       RIPFAC=1.+3.125*TMPVAL*TMPVAL*RIPSTP
       QQWV3(L)=RIPFAC*QQWV2(L)
       QQWCR(L)=SQRT( QQWV3(L)*QQWV3(L)+QQ(L,0)*QQ(L,0) )
       END DO
C
       ZBRMAX=-(1.E+12)*ZBRADJ
       ZBRMIN=(1.E+12)*ZBRADJ
       CDRGMAX=-1.E+12
       CDRGMIN=1.E+12
       WVDTMP=0.4/(WVFRQ*CTURB3)
       RKZTURB=0.4/CTURB3
C
       DO L=2,LA
       LS=LSC(L)
       LN=LNC(L)
       UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
       VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
       CURANG=ATAN2(VTMP,UTMP)
       COSWC=COS(CURANG-WACCWE(L))
       UMAGTMP=SQRT( U1(L,1)*U1(L,1)+V1U(L)*V1U(L)+1.E-12 )
       VMAGTMP=SQRT( U1V(L)*U1V(L)+V1(L,1)*V1(L,1)+1.E-12 )
       CDMAXU=STBXO(L)*H1U(L)/( 4.*DT*UMAGTMP )
       CDMAXV=STBYO(L)*H1V(L)/( 4.*DT*VMAGTMP )
       CDTMPU=-1.
       CDTMPV=-1.
       QWCTMPU=0.5*( QQWV2(L)+QQWV2(L+1) )
       QWCTMPV=0.5*( QQWV2(L)+QQWV2(LS ) )
       IF(ISWCBL.EQ.2) THEN
         QWCTMPU=0.5*( QQWC(L)+QQWC(L+1) )
         QWCTMPV=0.5*( QQWC(L)+QQWC(LS ) )
       END IF
       WVDELU=WVDTMP*SQRT(QWCTMPU)
       WVDELV=WVDTMP*SQRT(QWCTMPV)
       QWCTMPU=0.5*( QQWCR(L)+QQWCR(L+1) )
       QWCTMPV=0.5*( QQWCR(L)+QQWCR(LS ) )
       QWCTMPU=SQRT(QWCTMPU)
       QWCTMPV=SQRT(QWCTMPV)
       QCTMPU=0.5*( QQ(L,0)+QQ(L+1,0) )
       QCTMPV=0.5*( QQ(L,0)+QQ(LS ,0) )
       QWDQCU=QWCTMPU/SQRT(QCTMPU)
       QWDQCV=QWCTMPV/SQRT(QCTMPV)
       HZREFU=DZC(1)*H1U(L)
       HZREFV=DZC(1)*H1V(L)
       ZBREU=0.5*(ZBRE(L)+ZBRE(L+1))
       ZBREV=0.5*(ZBRE(L)+ZBRE(LS ))
       ZDHZRU=ZBREU/HZREFU
       ZDHZRV=ZBREV/HZREFV
       HZRUDZ=1./ZDHZRU
       HZRVDZ=1./ZDHZRV
       DWUD2Z=0.5*WVDELU/ZBREU
       DWVD2Z=0.5*WVDELV/ZBREV
       DWUDZ=2.*DWUD2Z
       DWVDZ=2.*DWVD2Z
       DWUDHR=WVDELU/HZREFU
       DWVDHR=WVDELV/HZREFV
       CDTMPUX=RKZTURB*QWCTMPU
       CDTMPVY=RKZTURB*QWCTMPV
       JWCBLU=0
       JWCBLV=0
       IF( HZRUDZ.LE.DWUD2Z) THEN
          CDTMPU=CDTMPUX/( (1.+ZDHZRU)*LOG(1.+HZRUDZ)-1. )
          JWCBLU=1
       END IF
       IF( HZRVDZ.LE.DWVD2Z) THEN
          CDTMPV=CDTMPVY/( (1.+ZDHZRV)*LOG(1.+HZRVDZ)-1. )
          JWCBLV=1
       END IF
       IF( HZRUDZ.GT.DWUD2Z.AND.HZRUDZ.LE.DWUDZ) THEN
          BOTTMP=(1.+ZDHZRU)*LOG(1.+DWUD2Z)-0.5*DWUDHR
     $          +0.5*HZRUDZ*(1.-0.5*DWUDHR)*(1.-0.5*DWUDHR)/(1.+DWUD2Z)
          CDTMPU=CDTMPUX/BOTTMP
          JWCBLU=2
       END IF
       IF( HZRVDZ.GT.DWVD2Z.AND.HZRVDZ.LE.DWVDZ) THEN
          BOTTMP=(1.+ZDHZRV)*LOG(1.+DWVD2Z)-0.5*DWVDHR
     $          +0.5*HZRVDZ*(1.-0.5*DWVDHR)*(1.-0.5*DWVDHR)/(1.+DWVD2Z)
          CDTMPV=CDTMPVY/BOTTMP
          JWCBLV=2
       END IF
       IF( HZRUDZ.GT.DWUDZ) THEN
          BOTTMP=QWDQCU*( (1.+ZDHZRU)*(LOG(1.+HZRUDZ)-LOG(1.+DWUDZ))
     $          +DWUDHR-1. )
          BOTTMP=BOTTMP+(1.+ZDHZRU)*LOG(1.+DWUD2Z)
     $          +DWUD2Z*(1.-1.25*DWUDHR-ZDHZRU)/(1.+DWUD2Z)
          CDTMPU=CDTMPUX/BOTTMP
          JWCBLU=3
       END IF
       IF( HZRVDZ.GT.DWVDZ) THEN
          BOTTMP=QWDQCV*( (1.+ZDHZRV)*(LOG(1.+HZRVDZ)-LOG(1.+DWVDZ))
     $          +DWVDHR-1. )
          BOTTMP=BOTTMP+(1.+ZDHZRV)*LOG(1.+DWVD2Z)
     $          +DWVD2Z*(1.-1.25*DWVDHR-ZDHZRV)/(1.+DWVD2Z)
          CDTMPV=CDTMPVY/BOTTMP
          JWCBLV=3
       END IF
       CDTMPU=CDTMPU/UMAGTMP
       CDTMPV=CDTMPV/VMAGTMP
C TMP DIAG
       IF(ISDZBR.EQ.N) THEN
       WRITE(1,1779) IL(L),JL(L),JWCBLU,JWCBLV
       WRITE(1,1780)
       WRITE(1,1781) ZBREU,WVDELU,HZREFU,CDTMPU,CDMAXU
       WRITE(1,1782)
       WRITE(1,1781) ZBREV,WVDELV,HZREFV,CDTMPV,CDMAXV
       END IF
C TMP DIAG
       IF(CDTMPU.LE.0.) CDTMPU=CDMAXU
       IF(CDTMPV.LE.0.) CDTMPV=CDMAXV
       STBX(L)=AVCON*STBXO(L)*CDTMPU
       STBY(L)=AVCON*STBYO(L)*CDTMPV
       STBX(L)=MIN(CDMAXU,STBX(L),0.11)
       STBY(L)=MIN(CDMAXV,STBY(L),0.11)
       END DO
C
       IF(ISDZBR.EQ.N) CLOSE(1)
C
       IF(ISDZBR.GE.1) THEN
C
       DO L=2,LA
       IF(ZBRE(L).GT.ZBRMAX) THEN
         ZBRMAX=ZBRE(L)
         LZBMAX=L
       END IF
       IF(ZBRE(L).LT.ZBRMIN) THEN
         ZBRMIN=ZBRE(L)
         LZBMIN=L
       END IF
       IF(STBX(L).GT.CDRGMAX) THEN
         CDRGMAX=STBX(L)
         LCDMAX=L
       END IF
       IF(STBX(L).LT.CDRGMIN) THEN
         CDRGMIN=STBX(L)
         LCDMIN=L
       END IF
       IF(STBY(L).GT.CDRGMAX) THEN
         CDRGMAX=STBY(L)
         LCDMAX=L
       END IF
       IF(STBY(L).LT.CDRGMIN) THEN
         CDRGMIN=STBY(L)
         LCDMIN=L
       END IF
       END DO
C
      OPEN(1,FILE='zbremx.out',STATUS='UNKNOWN',ACCESS='APPEND')
      HOTLYMX=DZC(1)*H1P(LZBMAX)
      HOTLYMN=DZC(1)*H1P(LZBMIN)
      WRITE(1,1739)N,IL(LZBMAX),JL(LZBMAX),ZBRMAX,HOTLYMX
      WRITE(1,1749)N,IL(LZBMIN),JL(LZBMIN),ZBRMIN,HOTLYMN
      WRITE(1,1759)N,IL(LCDMAX),JL(LCDMAX),CDRGMAX,STBX(LCDMAX),
     $             STBY(LCDMAX)
      WRITE(1,1769)N,IL(LCDMIN),JL(LCDMIN),CDRGMIN,STBX(LCDMIN),
     $             STBY(LCDMIN)
      CLOSE(1)
C
      END IF
C
 1948 CONTINUE
C
C**********************************************************************C
C
 1717 FORMAT(' N,I,J = ',I10,2I5,'   CDTOTU,CDMAXU = ',2F15.10)
 1718 FORMAT(' N,I,J = ',I10,2I5,'   CDTOTV,CDMAXV = ',2F15.10)
 1727 FORMAT(' N,I,J = ',I10,2I5,'   LAM CDTOTU,CDMAXU = ',2F15.10)
 1728 FORMAT(' N,I,J = ',I10,2I5,'   LAM CDTOTV,CDMAXV = ',2F15.10)
 1719 FORMAT(' N = ',I10,'  CDTOTUM,CDTOTVM = ',2F15.10)
 1729 FORMAT(' N = ',I10,'  CDMAXUM,CDMAXVM = ',2F15.10)
 1739 FORMAT(' N,I,J = ',I10,2I5,'  ZBRMAX,HBTLYMX = ',2E14.6)
 1749 FORMAT(' N,I,J = ',I10,2I5,'  ZBRMIN,HBTLYMN = ',2E14.6)
 1759 FORMAT(' N,I,J = ',I10,2I5,'  CDRGMAX,STBX,STBY = ',3E14.6)
 1769 FORMAT(' N,I,J = ',I10,2I5,'  CDRGMIN,STBX,STBY = ',3E14.6)
 1779 FORMAT(' I, J, JWCBLU, JWCBLV = ',4I8)
 1780 FORMAT('    ZBREU        WVDELU        HZREFU        CDTMPU    ',
     $     1X,'  CDMAXU')
 1781 FORMAT(5E12.4)
 1782 FORMAT('    ZBREV        WVDELV        HZREFV        CDTMPV    ',
     $     1X,'  CDMAXV')
C
C**********************************************************************C
C
      RETURN
      END
