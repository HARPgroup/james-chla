C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTOUT(IRSTYP)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RESTOUT WRITES A RESTART FILE
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      IF(IRSTYP.EQ.0) THEN
        OPEN(99,FILE='restart.out',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='restart.out',STATUS='UNKNOWN')
      END IF
C
      IF(IRSTYP.EQ.1) THEN
        OPEN(99,FILE='crashst.out',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='crashst.out',STATUS='UNKNOWN')
      END IF
C
C**********************************************************************C
C
      IF (ISRESTO.EQ.-11) THEN
C
      DO L=1,LC
C     HP(L)=HMP(L)
C     H1P(L)=HMP(L)
C     HWQ(L)=HMP(L)
C     H2WQ(L)=HMP(L)
      HP(L)=-BELV(L)
      H1P(L)=-BELV(L)
      HWQ(L)=-BELV(L)
      H2WQ(L)=-BELV(L)
      UHDYE(L)=0.
      UHDY1E(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      END DO
C
      DO K=0,KC
      DO L=1,LC
      QQ(L,K)=QQMIN
      QQ1(L,K)=QQMIN
      QQL(L,K)=QQLMIN
      QQL1(L,K)=QQLMIN
      DML(L,K)=DMLMIN
      END DO
      END DO
C
      DO K=1,KC
      DO L=1,LC
      U(L,K)=0.
      U1(L,K)=0.
      V(L,K)=0.
      V1(L,K)=0.
      END DO
      END DO
C
      DO K=1,KC
      DO L=1,LC
      SAL(L,K)=MAX(SAL(L,K),0.)
      SAL1(L,K)=MAX(SAL1(L,K),0.)
      END DO
      END DO
C
      END IF
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014
      WRITE(99,909)N,TIME    ! they are NOT used in restin1.for, Ji, 2/15/99
C
      DO L=2,LA
c     WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L)
      WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
      WRITE(99,907)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      WRITE(99,907)(U(L,K),K=1,KC)
      WRITE(99,907)(U1(L,K),K=1,KC)
      WRITE(99,907)(V(L,K),K=1,KC)
      WRITE(99,907)(V1(L,K),K=1,KC)
      WRITE(99,907)(QQ(L,K),K=0,KC)
      WRITE(99,907)(QQ1(L,K),K=0,KC)
      WRITE(99,907)(QQL(L,K),K=0,KC)
      WRITE(99,907)(QQL1(L,K),K=0,KC)
      WRITE(99,907)(DML(L,K),K=0,KC)
      IF(ISCO(1).EQ.1)THEN
       WRITE(99,907)(SAL(L,K),K=1,KC)
       WRITE(99,907)(SAL1(L,K),K=1,KC)
      END IF
      IF(ISCO(2).EQ.1)THEN
       WRITE(99,907)(TEM(L,K),K=1,KC)
       WRITE(99,907)(TEM1(L,K),K=1,KC)
       WRITE(99,907)TEMB(L)              ! Ji, 10/31/00
      END IF
      IF(ISCO(3).EQ.1)THEN
       WRITE(99,907)(DYE(L,K),K=1,KC)
       WRITE(99,907)(DYE1(L,K),K=1,KC)
      END IF
      IF(ISCO(4).EQ.1)THEN
       WRITE(99,907)SFLSBOT(L),(SFL(L,K),K=1,KC)
       WRITE(99,907)SFLSBOT(L),(SFL2(L,K),K=1,KC)
      END IF
      IF(ISCO(5).EQ.1)THEN
       DO NT=1,NTOX
       WRITE(99,907)(TOXB(L,K,NT),K=1,KB)
       WRITE(99,907)(TOX(L,K,NT),K=1,KC)
       WRITE(99,907)(TOXB1(L,K,NT),K=1,KB)
       WRITE(99,907)(TOX1(L,K,NT),K=1,KC)
       END DO
      END IF
      IF(ISCO(6).EQ.1)THEN
       DO NS=1,NSED
       WRITE(99,907)(SEDB(L,K,NS),K=1,KB)
c      WRITE(99,907)(SED1(L,K,NS),K=1,KC)  ! bug, Ji, 10/31/00
       WRITE(99,907)(SED (L,K,NS),K=1,KC)
       WRITE(99,907)(SEDB1(L,K,NS),K=1,KB)
       WRITE(99,907)(SED1(L,K,NS),K=1,KC)
       END DO
      END IF
      IF(ISCO(7).EQ.1)THEN
       DO NS=1,NSND
       WRITE(99,907)(SNDB(L,K,NS),K=1,KB)
       WRITE(99,907)(SND(L,K,NS),K=1,KC)
       WRITE(99,907)(SNDB1(L,K,NS),K=1,KB)
       WRITE(99,907)(SND1(L,K,NS),K=1,KC)
       END DO
      END IF
      IF(ISCI(6).EQ.1.OR.ISCI(7).EQ.1)THEN
       WRITE(99,907)(HBED(L,K),K=1,KB)
       WRITE(99,907)(HBED1(L,K),K=1,KB)
       WRITE(99,907)(VDRBED(L,K),K=1,KB)
       WRITE(99,907)(VDRBED1(L,K),K=1,KB)
      END IF
      END DO
C
      DO M=1,4
      IF(ISCO(M).EQ.1)THEN
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END IF
      END DO
C
      IF(ISCO(5).EQ.1)THEN
      DO NT=1,NTOX
       M=MSVTOX(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
      IF(ISCO(6).EQ.1)THEN
      DO NT=1,NSED
       M=MSVSED(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
      IF(ISCO(7).EQ.1)THEN
      DO NT=1,NSND
       M=MSVSND(NT)
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
       WRITE(99,907)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99,908)(NLON(LL,K,M),K=1,KC)
       WRITE(99,907)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
      DO L=2,LA
        WRITE(99,907)QSUME(L),(QSUM(L,K),K=1,KC)
      END DO
C
      IF (MDCHH.GE.1) THEN
        DO NMD=1,MDCHH
        WRITE(99,910)IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),
     $               IMDCHV(NMD),JMDCHV(NMD),QCHANU(NMD),QCHANV(NMD)
        END DO
      END IF
C
      IF (ISGWIE.GE.1) THEN
        DO L=2,LA
        WRITE(99,907)AGWELV(L),AGWELV1(L)
        END DO
      END IF
C
C
      IF (ISWAVE.GE.1) THEN
       OPEN(1,FILE='wvqwcp.out',STATUS='UNKNOWN')
       CLOSE(1, STATUS='DELETE')
       OPEN(1,FILE='wvqwcp.out',STATUS='UNKNOWN')
       DO L=2,LA
       WRITE(1,911)IL(L),JL(L),QQWV1(L),QQWV2(L),QQWV3(L),QQWC(L),
     $                         QQWCR(L),QQ(L,0)
       END DO
       CLOSE(1)
      END IF
C
c     IF (ISCO(2).GE.1.) THEN                        ! Ji, 10/31/00
c      OPEN(1,FILE='temp.rst',STATUS='UNKNOWN')      ! not needed, temb is save with tem
c      CLOSE(1, STATUS='DELETE')
c      OPEN(1,FILE='temp.rst',STATUS='UNKNOWN')
c      DO L=2,LA
c       WRITE(1,912)L,IL(L),JL(L),(TEM(L,K),K=1,KC),TEMB(L)
c      END DO
c      CLOSE(1)
c     END IF
C
C**********************************************************************C
C
      CLOSE(99)
C
C     IF(ISRESTO.EQ.-2) CALL OUT3D
C
C**********************************************************************C
C
  906 FORMAT(5E17.8)
  907 FORMAT(13E17.8)
cmrm      906 FORMAT(1p,5E17.8)
cmrm      907 FORMAT(1p,13E17.8)
  908 FORMAT(12I10)
  909 FORMAT(I20,4X,F12.4)
  910 FORMAT(6I5,2X,E17.8,2X,E17.8)
  911 FORMAT(2I5,2X,6E13.4)
cmrm      910 FORMAT(6I5,2X,1p,E17.8,2X,E17.8)
cmrm      911 FORMAT(2I5,2X,1p,6E13.4)
  912 FORMAT(3I5,12F7.3)
C
C**********************************************************************C
C
      RETURN
      END
