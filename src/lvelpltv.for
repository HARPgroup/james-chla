C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE LVELPLTV
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE VELPLTV WRITES A FILE FOR VERTICAL PLANE CONTOURING
C **  OF VELOCITY NORMAL TO AN ARBITARY SEQUENCE OF (I,J) POINTS AND
C **  AND VERTICAL PLANE TANGENTIAL-VERTICAL VELOCITY VECTORS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      REAL RVELN (KCM,100)
      REAL RVELT (KCM,100)
      REAL RW (KCM,100),HLRAVG(LCM)
C
      CHARACTER*80 TITLE10,TITLE20
      CHARACTER*80 TITLE40,TITLE50
      CHARACTER*80 TITLE70,TITLE80
C
C**********************************************************************C
C
C **  WRITE HEADINGS
C
      TITLE10='NORMAL LAG MEAN VEL CONTOURS'
      TITLE20='NORMAL AVG LAG MEAN VEL CONTOURS'
      TITLE40='TANGENTIAL LAG MEAN VEL CONTOURS'
      TITLE50='TANGENTIAL AVG LAG MEAN VEL CONTOURS'
      TITLE70='TANGENTIAL LAG MEAN VEL VECTORS'
      TITLE80='TANGENTIAL AVG LAG MEAN VEL VECTORS'
C
      LEVELS=KC
C
      IF (ISECVPV.GE.1) THEN
        OPEN(11,FILE='lmvcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='alvcnv1.out',STATUS='UNKNOWN')
        OPEN(41,FILE='lmvcvt1.out',STATUS='UNKNOWN')
        OPEN(51,FILE='alvcvt1.out',STATUS='UNKNOWN')
        OPEN(71,FILE='lmvvcv1.out',STATUS='UNKNOWN')
        OPEN(81,FILE='alvvcv1.out',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        CLOSE(41,STATUS='DELETE')
        CLOSE(51,STATUS='DELETE')
        CLOSE(71,STATUS='DELETE')
        CLOSE(81,STATUS='DELETE')
        OPEN(11,FILE='lmvcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='alvcnv1.out',STATUS='UNKNOWN')
        OPEN(41,FILE='lmvcvt1.out',STATUS='UNKNOWN')
        OPEN(51,FILE='alvcvt1.out',STATUS='UNKNOWN')
        OPEN(71,FILE='lmvvcv1.out',STATUS='UNKNOWN')
        OPEN(81,FILE='alvvcv1.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.2) THEN
        OPEN(12,FILE='lmvcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='alvcnv2.out',STATUS='UNKNOWN')
        OPEN(42,FILE='lmvcvt2.out',STATUS='UNKNOWN')
        OPEN(52,FILE='alvcvt2.out',STATUS='UNKNOWN')
        OPEN(72,FILE='lmvvcv2.out',STATUS='UNKNOWN')
        OPEN(82,FILE='alvvcv2.out',STATUS='UNKNOWN')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        CLOSE(42,STATUS='DELETE')
        CLOSE(52,STATUS='DELETE')
        CLOSE(72,STATUS='DELETE')
        CLOSE(82,STATUS='DELETE')
        OPEN(12,FILE='lmvcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='alvcnv2.out',STATUS='UNKNOWN')
        OPEN(42,FILE='lmvcvt2.out',STATUS='UNKNOWN')
        OPEN(52,FILE='alvcvt2.out',STATUS='UNKNOWN')
        OPEN(72,FILE='lmvvcv2.out',STATUS='UNKNOWN')
        OPEN(82,FILE='alvvcv2.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.3) THEN
        OPEN(13,FILE='lmvcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='alvcnv3.out',STATUS='UNKNOWN')
        OPEN(43,FILE='lmvcvt3.out',STATUS='UNKNOWN')
        OPEN(53,FILE='alvcvt3.out',STATUS='UNKNOWN')
        OPEN(73,FILE='lmvvcv3.out',STATUS='UNKNOWN')
        OPEN(83,FILE='alvvcv3.out',STATUS='UNKNOWN')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        CLOSE(43,STATUS='DELETE')
        CLOSE(53,STATUS='DELETE')
        CLOSE(73,STATUS='DELETE')
        CLOSE(83,STATUS='DELETE')
        OPEN(13,FILE='lmvcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='alvcnv3.out',STATUS='UNKNOWN')
        OPEN(43,FILE='lmvcvt3.out',STATUS='UNKNOWN')
        OPEN(53,FILE='alvcvt3.out',STATUS='UNKNOWN')
        OPEN(73,FILE='lmvvcv3.out',STATUS='UNKNOWN')
        OPEN(83,FILE='alvvcv3.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.4) THEN
        OPEN(14,FILE='lmvcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='alvcnv4.out',STATUS='UNKNOWN')
        OPEN(44,FILE='lmvcvt4.out',STATUS='UNKNOWN')
        OPEN(54,FILE='alvcvt4.out',STATUS='UNKNOWN')
        OPEN(74,FILE='lmvvcv4.out',STATUS='UNKNOWN')
        OPEN(84,FILE='alvvcv4.out',STATUS='UNKNOWN')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        CLOSE(44,STATUS='DELETE')
        CLOSE(54,STATUS='DELETE')
        CLOSE(74,STATUS='DELETE')
        CLOSE(84,STATUS='DELETE')
        OPEN(14,FILE='lmvcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='alvcnv4.out',STATUS='UNKNOWN')
        OPEN(44,FILE='lmvcvt4.out',STATUS='UNKNOWN')
        OPEN(54,FILE='alvcvt4.out',STATUS='UNKNOWN')
        OPEN(74,FILE='lmvvcv4.out',STATUS='UNKNOWN')
        OPEN(84,FILE='alvvcv4.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.5) THEN
        OPEN(15,FILE='lmvcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='alvcnv5.out',STATUS='UNKNOWN')
        OPEN(45,FILE='lmvcvt5.out',STATUS='UNKNOWN')
        OPEN(55,FILE='alvcvt5.out',STATUS='UNKNOWN')
        OPEN(75,FILE='lmvvcv5.out',STATUS='UNKNOWN')
        OPEN(85,FILE='alvvcv5.out',STATUS='UNKNOWN')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        CLOSE(45,STATUS='DELETE')
        CLOSE(55,STATUS='DELETE')
        CLOSE(75,STATUS='DELETE')
        CLOSE(85,STATUS='DELETE')
        OPEN(15,FILE='lmvcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='alvcnv5.out',STATUS='UNKNOWN')
        OPEN(45,FILE='lmvcvt5.out',STATUS='UNKNOWN')
        OPEN(55,FILE='alvcvt5.out',STATUS='UNKNOWN')
        OPEN(75,FILE='lmvvcv5.out',STATUS='UNKNOWN')
        OPEN(85,FILE='alvvcv5.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.6) THEN
        OPEN(16,FILE='lmvcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='alvcnv6.out',STATUS='UNKNOWN')
        OPEN(46,FILE='lmvcvt6.out',STATUS='UNKNOWN')
        OPEN(56,FILE='alvcvt6.out',STATUS='UNKNOWN')
        OPEN(76,FILE='lmvvcv6.out',STATUS='UNKNOWN')
        OPEN(86,FILE='alvvcv6.out',STATUS='UNKNOWN')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        CLOSE(46,STATUS='DELETE')
        CLOSE(56,STATUS='DELETE')
        CLOSE(76,STATUS='DELETE')
        CLOSE(86,STATUS='DELETE')
        OPEN(16,FILE='lmvcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='alvcnv6.out',STATUS='UNKNOWN')
        OPEN(46,FILE='lmvcvt6.out',STATUS='UNKNOWN')
        OPEN(56,FILE='alvcvt6.out',STATUS='UNKNOWN')
        OPEN(76,FILE='lmvvcv6.out',STATUS='UNKNOWN')
        OPEN(86,FILE='alvvcv6.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.7) THEN
        OPEN(17,FILE='lmvcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='alvcnv7.out',STATUS='UNKNOWN')
        OPEN(47,FILE='lmvcvt7.out',STATUS='UNKNOWN')
        OPEN(57,FILE='alvcvt7.out',STATUS='UNKNOWN')
        OPEN(77,FILE='lmvvcv7.out',STATUS='UNKNOWN')
        OPEN(87,FILE='alvvcv7.out',STATUS='UNKNOWN')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        CLOSE(47,STATUS='DELETE')
        CLOSE(57,STATUS='DELETE')
        CLOSE(77,STATUS='DELETE')
        CLOSE(87,STATUS='DELETE')
        OPEN(17,FILE='lmvcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='alvcnv7.out',STATUS='UNKNOWN')
        OPEN(47,FILE='lmvcvt7.out',STATUS='UNKNOWN')
        OPEN(57,FILE='alvcvt7.out',STATUS='UNKNOWN')
        OPEN(77,FILE='lmvvcv7.out',STATUS='UNKNOWN')
        OPEN(87,FILE='alvvcv7.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.8) THEN
        OPEN(18,FILE='lmvcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='alvcnv8.out',STATUS='UNKNOWN')
        OPEN(48,FILE='lmvcvt8.out',STATUS='UNKNOWN')
        OPEN(58,FILE='alvcvt8.out',STATUS='UNKNOWN')
        OPEN(78,FILE='lmvvcv8.out',STATUS='UNKNOWN')
        OPEN(88,FILE='alvvcv8.out',STATUS='UNKNOWN')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        CLOSE(48,STATUS='DELETE')
        CLOSE(58,STATUS='DELETE')
        CLOSE(78,STATUS='DELETE')
        CLOSE(88,STATUS='DELETE')
        OPEN(18,FILE='lmvcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='alvcnv8.out',STATUS='UNKNOWN')
        OPEN(48,FILE='lmvcvt8.out',STATUS='UNKNOWN')
        OPEN(58,FILE='alvcvt8.out',STATUS='UNKNOWN')
        OPEN(78,FILE='lmvvcv8.out',STATUS='UNKNOWN')
        OPEN(88,FILE='alvvcv8.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.9) THEN
        OPEN(19,FILE='lmvcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='alvcnv9.out',STATUS='UNKNOWN')
        OPEN(49,FILE='lmvcvt9.out',STATUS='UNKNOWN')
        OPEN(59,FILE='alvcvt9.out',STATUS='UNKNOWN')
        OPEN(79,FILE='lmvvcv9.out',STATUS='UNKNOWN')
        OPEN(89,FILE='alvvcv9.out',STATUS='UNKNOWN')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        CLOSE(49,STATUS='DELETE')
        CLOSE(59,STATUS='DELETE')
        CLOSE(79,STATUS='DELETE')
        CLOSE(89,STATUS='DELETE')
        OPEN(19,FILE='lmvcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='alvcnv9.out',STATUS='UNKNOWN')
        OPEN(49,FILE='lmvcvt9.out',STATUS='UNKNOWN')
        OPEN(59,FILE='alvcvt9.out',STATUS='UNKNOWN')
        OPEN(79,FILE='lmvvcv9.out',STATUS='UNKNOWN')
        OPEN(89,FILE='alvvcv9.out',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LUN4=40+IS
      LUN5=50+IS
      LUN7=70+IS
      LUN8=80+IS
      LINES=NIJVPV(IS)
      WRITE (LUN1,99) TITLE10,CVTITLE(LUN1)
      WRITE (LUN2,99) TITLE20,CVTITLE(LUN2)
      WRITE (LUN4,99) TITLE40,CVTITLE(LUN4)
      WRITE (LUN5,99) TITLE50,CVTITLE(LUN5)
      WRITE (LUN7,99) TITLE70,CVTITLE(LUN7)
      WRITE (LUN8,99) TITLE80,CVTITLE(LUN8)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN4,101)LINES,LEVELS
      WRITE (LUN5,101)LINES,LEVELS
      WRITE (LUN7,101)LINES,LEVELS
      WRITE (LUN8,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      WRITE (LUN4,250)(ZZ(K),K=1,KC)
      WRITE (LUN5,250)(ZZ(K),K=1,KC)
      WRITE (LUN7,250)(ZZ(K),K=1,KC)
      WRITE (LUN8,250)(ZZ(K),K=1,KC)
      END DO
C
C----------------------------------------------------------------------C
C
      DO M=1,MLRPDRT
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN4=40+IS
      LUN7=70+IS
      WRITE (LUN1,100)NLRPDRT(M)
      WRITE (LUN4,100)NLRPDRT(M)
      WRITE (LUN7,100)NLRPDRT(M)
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
        DO K=1,KC
        RVELN(K,NN)=100.*(XLRPD(NLRPDL(L),K,M)*COSC
     $            +YLRPD(NLRPDL(L),K,M)*SINC)
        RVELT(K,NN)=-100.*XLRPD(NLRPDL(L),K,M)*SINC
     $            +100.*YLRPD(NLRPDL(L),K,M)*COSC
        RW(K,NN)=100.*ZLRPD(NLRPDL(L),K,M)
        END DO
       END DO
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
c      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(RVELN(K,NN),K=1,KC)
       WRITE(LUN4,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RW(K,NN),K=1,KC)
       END DO
      END DO
C
      END DO
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN4=40+IS
      LUN7=70+IS
      CLOSE(LUN1)
      CLOSE(LUN4)
      CLOSE(LUN7)
      END DO
C
C----------------------------------------------------------------------C
C     
      DO NLR=1,NLRPD
      HLRAVG(NLR)=0.
      END DO
      DO M=1,MLRPDRT
      DO NLR=1,NLRPD
      HLRAVG(NLR)=HLRAVG(NLR)+HLRPD(NLR,M)
      END DO
      END DO
      DO NLR=1,NLRPD
      HLRAVG(NLR)=HLRAVG(NLR)/FLOAT(MLRPDRT)
      END DO
C     
      M=MLRAVG
C
      DO IS=1,ISECVPV
      LUN2=20+IS
      LUN5=50+IS
      LUN8=80+IS
      WRITE (LUN2,100)N
      WRITE (LUN5,100)N
      WRITE (LUN8,100)N
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
        DO K=1,KC
        RVELN(K,NN)=100.*(XLRPD(NLRPDL(L),K,M)*COSC
     $            +YLRPD(NLRPDL(L),K,M)*SINC)
        RVELT(K,NN)=-100.*XLRPD(NLRPDL(L),K,M)*SINC
     $            +100.*YLRPD(NLRPDL(L),K,M)*COSC
        RW(K,NN)=100.*ZLRPD(NLRPDL(L),K,M)
        END DO
       END DO
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
c      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,250)(RVELN(K,NN),K=1,KC)
       WRITE(LUN5,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN8,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN8,250)(RW(K,NN),K=1,KC)
       END DO
      CLOSE(LUN2)
      CLOSE(LUN5)
      CLOSE(LUN8)
      END DO
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10)
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
