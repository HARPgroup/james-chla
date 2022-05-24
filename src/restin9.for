C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTIN9
c
c   modified based on restout9.for for binary restart, ji, 10/31/00
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE RESTIN1 READS A RESTART FILE
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
c     DIMENSION TDUMMY(KCM)
C
C**********************************************************************C
C
      OPEN(1,FILE='restart.inp',form='binary')
C
C**********************************************************************C
C
      ISBELVC=0
C
      READ(1) NREST,time
      write(6,*) "Hot start from :  ", Nrest,time
C
      read(1) HP,H1P,HWQ,H2WQ,BELV
     *,UHDYE,UHDY1E,VHDXE,VHDX1E
     *,U,U1,V,V1,QQ,QQ1,QQL,QQL1,DML
      IF(ISCO(1).EQ.1) read(1) SAL,SAL1
      IF(ISCO(2).EQ.1) read(1) TEM,TEM1,TEMB
      IF(ISCO(3).EQ.1) read(1) DYE,DYE1
      IF(ISCO(4).EQ.1) read(1) SFLSBOT,SFL,SFLSBOT,SFL2
      IF(ISCO(5).EQ.1) read(1) TOXB,TOX,TOXB1,TOX1
      IF(ISCO(6).EQ.1) read(1) SEDB,SED,SEDB1,SED1
      IF(ISCO(7).EQ.1) read(1) SNDB,SND,SNDB1,SND1
      IF(ISCI(6).EQ.1.OR.ISCI(7).EQ.1) read(1)
     +    HBED,HBED1,VDRBED,VDRBED1
C
      DO M=1,4
      IF(ISCI(M).EQ.1)THEN
C
       DO LL=1,NCBS
       read(1)(NLOS(LL,K,M),K=1,KC)
       read(1)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       read(1)(NLOW(LL,K,M),K=1,KC)
       read(1)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       read(1)(NLOE(LL,K,M),K=1,KC)
       read(1)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       read(1)(NLON(LL,K,M),K=1,KC)
       read(1)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END IF
      END DO
C
      IF(ISCI(5).EQ.1)THEN
      DO NT=1,NTOX
      M=MSVTOX(NT)
C
       DO LL=1,NCBS
       read(1)(NLOS(LL,K,M),K=1,KC)
       read(1)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       read(1)(NLOW(LL,K,M),K=1,KC)
       read(1)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       read(1)(NLOE(LL,K,M),K=1,KC)
       read(1)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       read(1)(NLON(LL,K,M),K=1,KC)
       read(1)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
      IF(ISCI(6).EQ.1)THEN
      DO NT=1,NSED
      M=MSVSED(NT)
C
       DO LL=1,NCBS
       read(1)(NLOS(LL,K,M),K=1,KC)
       read(1)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       read(1)(NLOW(LL,K,M),K=1,KC)
       read(1)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       read(1)(NLOE(LL,K,M),K=1,KC)
       read(1)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       read(1)(NLON(LL,K,M),K=1,KC)
       read(1)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
      IF(ISCI(7).EQ.1)THEN
      DO NT=1,NSND
      M=MSVSND(NT)
C
       DO LL=1,NCBS
       read(1)(NLOS(LL,K,M),K=1,KC)
       read(1)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       read(1)(NLOW(LL,K,M),K=1,KC)
       read(1)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       read(1)(NLOE(LL,K,M),K=1,KC)
       read(1)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       read(1)(NLON(LL,K,M),K=1,KC)
       read(1)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
        read(1) QSUME,QSUM
        read(1) NXX
c check
      if(Nxx.ne.nrest) then
      write(6,*) " Restart Wrong : ", NREST,NXX
 !     stop
      else
      write(6,*) " Read restart.inp successfully "
      endif
C
c     IF (MDCHH.GE.1) THEN
c       DO NMD=1,MDCHH
c       READ(1,*,ERR=1032)ITMP1,JTMP1,ITMP2,JTMP2,
c    $                    ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
C       WRITE(6,6667)NMD,ITMP1,JTMP1,ITMP2,JTMP2,
C    $                   ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
c       END DO
c     END IF
C
c     IF (ISGWIE.GE.1) THEN
c       DO L=2,LA
c       READ(1,*,ERR=1033)AGWELV(L),AGWELV1(L)
c       END DO
c     END IF
C
      CLOSE(1)
C
 6666 FORMAT(3I10,F12.6)
 6667 FORMAT(7I5,2X,E12.4,2X,E12.4)
C
      DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      SED(1,K,1)=0.
      SND(1,K,1)=0.
      TOX(1,K,1)=0.
      SFL(1,K)=0.
      CWQ(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SED1(1,K,1)=0.
      SND1(1,K,1)=0.
      TOX1(1,K,1)=0.
      SFL2(1,K)=0.
      CWQ2(1,K)=0.
      VHDX1(1,K)=0.
      UHDY1(1,K)=0.
      VHDXWQ(1,K)=0.
      UHDYWQ(1,K)=0.
      SAL(LC,K)=0.
      TEM(LC,K)=0.
      DYE(LC,K)=0.
      SED(LC,K,1)=0.
      SND(LC,K,1)=0.
      TOX(LC,K,1)=0.
      SFL(LC,K)=0.
      CWQ(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SED1(LC,K,1)=0.
      SND1(LC,K,1)=0.
      TOX1(LC,K,1)=0.
      SFL2(LC,K)=0.
      CWQ2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXWQ(LC,K)=0.
      UHDYWQ(LC,K)=0.
      END DO
C
      DO L=2,LA
      UHDYE(L)=SUB(L)*UHDYE(L)
      UHDY1E(L)=SUB(L)*UHDY1E(L)
      VHDXE(L)=SVB(L)*VHDXE(L)
      VHDX1E(L)=SVB(L)*VHDX1E(L)
      END DO
C
      DO K=1,KC
      DO L=2,LA
      U(L,K)=SUB(L)*U(L,K)
      U1(L,K)=SUB(L)*U1(L,K)
      V(L,K)=SVB(L)*V(L,K)
      V1(L,K)=SVB(L)*V1(L,K)
      END DO
      END DO
C
C**********************************************************************C
C
      DO L=2,LA
      LS=LSC(L)
      H1U(L)=0.5*(H1P(L)+H1P(L-1))
      H1V(L)=0.5*(H1P(L)+H1P(LS))
      P1(L)=G*(H1P(L)+BELV(L))
      HU(L)=0.5*(HP(L)+HP(L-1))
      HV(L)=0.5*(HP(L)+HP(LS))
      P(L)=G*(HP(L)+BELV(L))
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
      END DO
C
      H1U(1)=H1U(2)
      H1V(1)=H1V(2)
      P1(1)=P1(2)
      HU(1)=HU(2)
      HV(1)=HV(2)
      P(1)=P(2)
      HPI(1)=1./HP(2)
      HUI(1)=1./HU(2)
      HVI(1)=1./HV(2)
      H1UI(1)=1./H1U(2)
      H1VI(1)=1./H1V(2)
C
      H1U(LC)=H1U(LA)
      H1V(LC)=H1V(LA)
      P1(LC)=P1(LA)
      HU(LC)=HU(LA)
      HV(LC)=HV(LA)
      P(LC)=P(LA)
      HPI(LC)=1./HP(LA)
      HUI(LC)=1./HU(LA)
      HVI(LC)=1./HV(LA)
      H1UI(LC)=1./H1U(LA)
      H1VI(LC)=1./H1V(LA)
C
      DO K=1,KC
      DO L=2,LA
      UHDY1(L,K)=DYU(L)*H1U(L)*U1(L,K)
      VHDX1(L,K)=DXV(L)*H1V(L)*V1(L,K)
      UHDY(L,K)=DYU(L)*HU(L)*U(L,K)
      VHDX(L,K)=DXV(L)*HV(L)*V(L,K)
      SAL(L,K)=MAX(SAL(L,K),0.)
      SAL1(L,K)=MAX(SAL1(L,K),0.)
      END DO
      END DO
C
C**********************************************************************C
C
C **  CORRECT FOR CHANGED BOTTOM ELEV
C
c     IF (ISRESTI.EQ.-1.AND.ISBELVC.EQ.1)THEN
C
c       DO L=2,LA
c       UHE(L)=0.
c       VHE(L)=0.
c       END DO
C
c       DO K=1,KC
c       DO L=2,LA
c       UHE(L)=UHE(L)+UHDY1(L,K)
c       VHE(L)=VHE(L)+VHDX1(L,K)
c       END DO
c       END DO
C
c       DO L=2,LA
c       IF(UHE(L).NE.0.)THEN
c         RTMP=UHDY1E(L)/UHE(L)
c         DO K=1,KC
c         U1(L,K)=RTMP*U1(L,K)
c         END DO
c       END IF
c       IF(VHE(L).NE.0.)THEN
c         RTMP=VHDX1E(L)/VHE(L)
c         DO K=1,KC
c         V1(L,K)=RTMP*V1(L,K)
c         END DO
c       END IF
c       END DO
C
c       DO L=2,LA
c       UHE(L)=0.
c       VHE(L)=0.
c       END DO
C
c       DO K=1,KC
c       DO L=2,LA
c       UHE(L)=UHE(L)+UHDY(L,K)
c       VHE(L)=VHE(L)+VHDX(L,K)
c       END DO
c       END DO
C
c       DO L=2,LA
c       IF(UHE(L).NE.0.)THEN
c         RTMP=UHDYE(L)/UHE(L)
c         DO K=1,KC
c         U(L,K)=RTMP*U(L,K)
c         END DO
c       END IF
c       IF(VHE(L).NE.0.)THEN
c         RTMP=VHDXE(L)/VHE(L)
c         DO K=1,KC
c         V(L,K)=RTMP*V(L,K)
c         END DO
c       END IF
c       END DO
C
c       DO K=1,KC
c       DO L=2,LA
c       UHDY1(L,K)=DYU(L)*H1U(L)*U1(L,K)
c       VHDX1(L,K)=DXV(L)*H1V(L)*V1(L,K)
c       UHDY(L,K)=DYU(L)*HU(L)*U(L,K)
c       VHDX(L,K)=DXV(L)*HV(L)*V(L,K)
c       END DO
c       END DO
C
c     END IF
C
C**********************************************************************C
C
      N=0
      IF (ISDRY.EQ.0) THEN
        CALL CALTSXY
        CALL CALQVS (2)
      END IF
C
      DO K=1,KS
      RDZC=DZC(K)
      DO L=2,LA
      LN=LNC(L)
      W(L,K)=SWB(L)*(W(L,K-1)
     $       -RDZC*(UHDY(L+1,K)-UHDY(L,K)-UHDYE(L+1)+UHDYE(L)
     $       +VHDX(LN,K)-VHDX(L,K)-VHDXE(LN)+VHDXE(L))*DXYIP(L))
     $        +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      W1(L,K)=SWB(L)*(W1(L,K-1)
     $       -RDZC*(UHDY1(L+1,K)-UHDY1(L,K)-UHDY1E(L+1)+UHDY1E(L)
     $       +VHDX1(LN,K)-VHDX1(L,K)-VHDX1E(LN)+VHDX1E(L))*DXYIP(L))
     $        +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      END DO
      END DO
C
      OPEN(1,FILE='newsalt.out',STATUS='UNKNOWN')   ! Ji, 10/31/00, for checking only
      IONE=1
      WRITE(1,101)IONE
      DO L=2,LA
      WRITE(1,102)L,IL(L),JL(L),(SAL(L,K),K=1,KC)
      END DO
      CLOSE(1)
c
C ** READ BED TEMPERATURE  HTBED1 HTBED2
C
c     IF (ISCO(2).GE.1.) THEN               ! Ji, modified, 10/31/00
c      OPEN(1,FILE='temp.rst',STATUS='UNKNOWN')
c      DO L=2,LA
c       READ(1,*)LDUM,IDUM,JDUM,(TDUMMY(K),K=1,KC),TEMB(L)
c      END DO
c      CLOSE(1)
c     END IF
C
C
  101 FORMAT(I5)
  102 FORMAT(3I5,12F8.2)
C
C**********************************************************************C
C
C **  WRITE READ ERRORS ON RESTART
C
      GO TO 3000
 1000 WRITE(6,2000)
      STOP
 1001 WRITE(6,2001)L
      STOP
 1002 WRITE(6,2002)L
      STOP
 1003 WRITE(6,2003)L
      STOP
 1004 WRITE(6,2004)L
      STOP
 1005 WRITE(6,2005)L
      STOP
 1006 WRITE(6,2006)L
      STOP
 1007 WRITE(6,2007)L
      STOP
 1008 WRITE(6,2008)L
      STOP
 1009 WRITE(6,2009)L
      STOP
 1010 WRITE(6,2010)L
      STOP
 1011 WRITE(6,2011)L
      STOP
 1012 WRITE(6,2012)L
      STOP
 1013 WRITE(6,2013)L
      STOP
 1014 WRITE(6,2014)L
      STOP
 1015 WRITE(6,2015)L
      STOP
 1016 WRITE(6,2016)L
      STOP
 1017 WRITE(6,2017)L
      STOP
 1018 WRITE(6,2018)L
      STOP
 1019 WRITE(6,2019)L
      STOP
 1020 WRITE(6,2020)L
      STOP
 1021 WRITE(6,2021)L
      STOP
 1022 WRITE(6,2022)L
      STOP
 1023 WRITE(6,2023)L
      STOP
 1024 WRITE(6,2024)L
      STOP
 1025 WRITE(6,2025)L
      STOP
 1026 WRITE(6,2026)L
      STOP
 1027 WRITE(6,2027)L
      STOP
 1028 WRITE(6,2028)L
      STOP
 1029 WRITE(6,2029)L
      STOP
 1030 WRITE(6,2030)L
      STOP
 1031 WRITE(6,2031)L
      STOP
 1032 WRITE(6,2032)NMD
      STOP
 1033 WRITE(6,2033)L
      STOP
 3000 CONTINUE
C
C**********************************************************************C
C
  600 FORMAT(2X,'I,J,BELVOLD,BELVNEW',2I5,2F12.2)
  906 FORMAT(5E15.7)
  907 FORMAT(12E12.4)
  908 FORMAT(12I10)
 2000 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1000')
 2001 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1001 L =',I6)
 2002 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1002 L =',I6)
 2003 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1003 L =',I6)
 2004 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1004 L =',I6)
 2005 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1005 L =',I6)
 2006 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1006 L =',I6)
 2007 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1007 L =',I6)
 2008 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1008 L =',I6)
 2009 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1009 L =',I6)
 2010 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1010 L =',I6)
 2011 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1011 L =',I6)
 2012 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1012 L =',I6)
 2013 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1013 L =',I6)
 2014 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1014 L =',I6)
 2015 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1015 L =',I6)
 2016 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1016 L =',I6)
 2017 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1017 L =',I6)
 2018 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1018 L =',I6)
 2019 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1019 L =',I6)
 2020 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1020 L =',I6)
 2021 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1021 L =',I6)
 2022 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1022 L =',I6)
 2023 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1023 L =',I6)
 2024 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1024 L =',I6)
 2025 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1025 L =',I6)
 2026 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1026 L =',I6)
 2027 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1027 L =',I6)
 2028 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1028 L =',I6)
 2029 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1029 L =',I6)
 2030 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1030 L =',I6)
 2031 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1031 L =',I6)
 2032 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1032 NMD =',I6)
 2033 FORMAT(1X,'READ ERROR ON FILE restart.inp ERR 1033 L =',I6)
C
C**********************************************************************C
C
      RETURN
      END
