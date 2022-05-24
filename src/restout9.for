C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTOUT9
c
c    Modified based on restout.for for binary saving, Ji, 10/31/00
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
c     IF(IRSTYP.EQ.0) THEN
        OPEN(99,FILE='restart.bin',STATUS='UNKNOWN',form='binary')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='restart.bin',STATUS='new',form='binary')
c     END IF
C
C**********************************************************************C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
      WRITE(99) N,TIME    ! they are NOT used in restin1.for, Ji, 2/15/99
C
      WRITE(99) HP,H1P,HWQ,H2WQ,BELV
     *,UHDYE,UHDY1E,VHDXE,VHDX1E
     *,U,U1,V,V1,QQ,QQ1,QQL,QQL1,DML
      IF(ISCO(1).EQ.1) WRITE(99) SAL,SAL1
      IF(ISCO(2).EQ.1) WRITE(99) TEM,TEM1,TEMB
      IF(ISCO(3).EQ.1) WRITE(99) DYE,DYE1
      IF(ISCO(4).EQ.1) WRITE(99) SFLSBOT,SFL,SFLSBOT,SFL2
      IF(ISCO(5).EQ.1) WRITE(99) TOXB,TOX,TOXB1,TOX1
      IF(ISCO(6).EQ.1) WRITE(99) SEDB,SED,SEDB1,SED1
      IF(ISCO(7).EQ.1) WRITE(99) SNDB,SND,SNDB1,SND1
      IF(ISCI(6).EQ.1.OR.ISCI(7).EQ.1) WRITE(99)
     +    HBED,HBED1,VDRBED,VDRBED1
C
      DO M=1,4
      IF(ISCO(M).EQ.1)THEN
C
       DO LL=1,NCBS
       DO K=1,KC
       NLOS(LL,K,M)=NLOS(LL,K,M)-N
       END DO
       WRITE(99)(NLOS(LL,K,M),K=1,KC)
       WRITE(99)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99)(NLOW(LL,K,M),K=1,KC)
       WRITE(99)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99)(NLOE(LL,K,M),K=1,KC)
       WRITE(99)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99)(NLON(LL,K,M),K=1,KC)
       WRITE(99)(CLON(LL,K,M),K=1,KC)
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
       WRITE(99)(NLOS(LL,K,M),K=1,KC)
       WRITE(99)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99)(NLOW(LL,K,M),K=1,KC)
       WRITE(99)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99)(NLOE(LL,K,M),K=1,KC)
       WRITE(99)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99)(NLON(LL,K,M),K=1,KC)
       WRITE(99)(CLON(LL,K,M),K=1,KC)
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
       WRITE(99)(NLOS(LL,K,M),K=1,KC)
       WRITE(99)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99)(NLOW(LL,K,M),K=1,KC)
       WRITE(99)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99)(NLOE(LL,K,M),K=1,KC)
       WRITE(99)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99)(NLON(LL,K,M),K=1,KC)
       WRITE(99)(CLON(LL,K,M),K=1,KC)
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
       WRITE(99)(NLOS(LL,K,M),K=1,KC)
       WRITE(99)(CLOS(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBW
       DO K=1,KC
       NLOW(LL,K,M)=NLOW(LL,K,M)-N
       END DO
       WRITE(99)(NLOW(LL,K,M),K=1,KC)
       WRITE(99)(CLOW(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBE
       DO K=1,KC
       NLOE(LL,K,M)=NLOE(LL,K,M)-N
       END DO
       WRITE(99)(NLOE(LL,K,M),K=1,KC)
       WRITE(99)(CLOE(LL,K,M),K=1,KC)
       END DO
C
       DO LL=1,NCBN
       DO K=1,KC
       NLON(LL,K,M)=NLON(LL,K,M)-N
       END DO
       WRITE(99)(NLON(LL,K,M),K=1,KC)
       WRITE(99)(CLON(LL,K,M),K=1,KC)
       END DO
C
      END DO
      END IF
C
        WRITE(99) QSUME,QSUM
C
        write(99) N  ! write again for checking
c
c     IF (MDCHH.GE.1) THEN
c       DO NMD=1,MDCHH
c       WRITE(99,910)IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),
c    $               IMDCHV(NMD),JMDCHV(NMD),QCHANU(NMD),QCHANV(NMD)
c       END DO
c     END IF
C
c     IF (ISGWIE.GE.1) THEN
c       DO L=2,LA
c       WRITE(99,907)AGWELV(L),AGWELV1(L)
c       END DO
c     END IF
C
C
c     IF (ISWAVE.GE.1) THEN
c      OPEN(1,FILE='wvqwcp.out',STATUS='UNKNOWN')
c      CLOSE(1, STATUS='DELETE')
c      OPEN(1,FILE='wvqwcp.out',STATUS='UNKNOWN')
c      DO L=2,LA
c      WRITE(1,911)IL(L),JL(L),QQWV1(L),QQWV2(L),QQWV3(L),QQWC(L),
c    $                         QQWCR(L),QQ(L,0)
c      END DO
c      CLOSE(1)
c     END IF
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
