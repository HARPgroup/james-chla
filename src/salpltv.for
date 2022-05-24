C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SALPLTV(ITMP)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE SALPLTV WRITES A FILE FOR VERTICAL PLANE CONTOURING
C **  OF SALINITY, DYE CONCENTRATION, AND SEDIMENT CONCENTRATION 
C **  ALONG AN ARBITARY SEQUENCE OF (I,J) POINTS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7
C     DIMENSION QQCTR(KCM)
C
C**********************************************************************C
C
      IF (JSSPV(ITMP).NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      TITLE1='INSTANTANEOUS SALINITY CONTOURS'
      TITLE2='INSTANTANEOUS TEMPERATURE CONTOURS'
      TITLE3='INSTANTANEOUS DYE CONC CONTOURS'
      TITLE4='INSTANT COHESIVE SED CONC CONTOURS'
      TITLE5='INSTANT NONCOH SED CONC CONTOURS'
      TITLE6='INSTANT TOXIC CONTAM CONC CONTOURS'
      TITLE7='INSTANT SHELLFISH LARVAE CONTOURS'
C
      IF(ITMP.EQ.1) THEN
      IF (ISTRAN(1).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(11,FILE='salcnv1.out',STATUS='UNKNOWN')
          CLOSE(11,STATUS='DELETE')
          OPEN(11,FILE='salcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(12,FILE='salcnv2.out',STATUS='UNKNOWN')
          CLOSE(12,STATUS='DELETE')
          OPEN(12,FILE='salcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(13,FILE='salcnv3.out',STATUS='UNKNOWN')
          CLOSE(13,STATUS='DELETE')
          OPEN(13,FILE='salcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(14,FILE='salcnv4.out',STATUS='UNKNOWN')
          CLOSE(14,STATUS='DELETE')
          OPEN(14,FILE='salcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(15,FILE='salcnv5.out',STATUS='UNKNOWN')
          CLOSE(15,STATUS='DELETE')
          OPEN(15,FILE='salcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(16,FILE='salcnv6.out',STATUS='UNKNOWN')
          CLOSE(16,STATUS='DELETE')
          OPEN(16,FILE='salcnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(17,FILE='salcnv7.out',STATUS='UNKNOWN')
          CLOSE(17,STATUS='DELETE')
          OPEN(17,FILE='salcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(18,FILE='salcnv8.out',STATUS='UNKNOWN')
          CLOSE(18,STATUS='DELETE')
          OPEN(18,FILE='salcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(19,FILE='salcnv9.out',STATUS='UNKNOWN')
          CLOSE(19,STATUS='DELETE')
          OPEN(19,FILE='salcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=10+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE1,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.2) THEN
      IF (ISTRAN(2).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(21,FILE='temcnv1.out',STATUS='UNKNOWN')
          CLOSE(21 ,STATUS='DELETE')
          OPEN(21,FILE='temcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(22,FILE='temcnv2.out',STATUS='UNKNOWN')
          CLOSE(22,STATUS='DELETE')
          OPEN(22,FILE='temcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(23,FILE='temcnv3.out',STATUS='UNKNOWN')
          CLOSE(23,STATUS='DELETE')
          OPEN(23,FILE='temcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(24,FILE='temcnv4.out',STATUS='UNKNOWN')
          CLOSE(24,STATUS='DELETE')
          OPEN(24,FILE='temcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(25,FILE='temcnv5.out',STATUS='UNKNOWN')
          CLOSE(25,STATUS='DELETE')
          OPEN(25,FILE='temcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(26,FILE='temcnv6.out',STATUS='UNKNOWN')
          CLOSE(26,STATUS='DELETE')
          OPEN(25,FILE='temcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(27,FILE='temcnv7.out',STATUS='UNKNOWN')
          CLOSE(27,STATUS='DELETE')
          OPEN(27,FILE='temcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(28,FILE='temcnv8.out',STATUS='UNKNOWN')
          CLOSE(28,STATUS='DELETE')
          OPEN(28,FILE='temcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(29,FILE='temcnv9.out',STATUS='UNKNOWN')
          CLOSE(29,STATUS='DELETE')
          OPEN(29,FILE='temcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=20+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE2,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.3) THEN
      IF (ISTRAN(3).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(31,FILE='dyecnv1.out',STATUS='UNKNOWN')
          CLOSE(31 ,STATUS='DELETE')
          OPEN(31,FILE='dyecnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(32,FILE='dyecnv2.out',STATUS='UNKNOWN')
          CLOSE(32,STATUS='DELETE')
          OPEN(32,FILE='dyecnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(33,FILE='dyecnv3.out',STATUS='UNKNOWN')
          CLOSE(33,STATUS='DELETE')
          OPEN(33,FILE='dyecnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(34,FILE='dyecnv4.out',STATUS='UNKNOWN')
          CLOSE(34,STATUS='DELETE')
          OPEN(34,FILE='dyecnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(35,FILE='dyecnv5.out',STATUS='UNKNOWN')
          CLOSE(35,STATUS='DELETE')
          OPEN(35,FILE='dyecnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(36,FILE='dyecnv6.out',STATUS='UNKNOWN')
          CLOSE(36,STATUS='DELETE')
          OPEN(36,FILE='dyecnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(37,FILE='dyecnv7.out',STATUS='UNKNOWN')
          CLOSE(37,STATUS='DELETE')
          OPEN(37,FILE='dyecnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(38,FILE='dyecnv8.out',STATUS='UNKNOWN')
          CLOSE(38,STATUS='DELETE')
          OPEN(38,FILE='dyecnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(39,FILE='dyecnv9.out',STATUS='UNKNOWN')
          CLOSE(39,STATUS='DELETE')
          OPEN(39,FILE='dyecnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=30+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE3,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.6) THEN
      IF (ISTRAN(6).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(41,FILE='sedcnv1.out',STATUS='UNKNOWN')
          CLOSE(41,STATUS='DELETE')
          OPEN(41,FILE='sedcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(42,FILE='sedcnv2.out',STATUS='UNKNOWN')
          CLOSE(42,STATUS='DELETE')
          OPEN(42,FILE='sedcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(43,FILE='sedcnv3.out',STATUS='UNKNOWN')
          CLOSE(43,STATUS='DELETE')
          OPEN(43,FILE='sedcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(44,FILE='sedcnv4.out',STATUS='UNKNOWN')
          CLOSE(44,STATUS='DELETE')
          OPEN(44,FILE='sedcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(45,FILE='sedcnv5.out',STATUS='UNKNOWN')
          CLOSE(45,STATUS='DELETE')
          OPEN(45,FILE='sedcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(46,FILE='sedcnv6.out',STATUS='UNKNOWN')
          CLOSE(46,STATUS='DELETE')
          OPEN(46,FILE='sedcnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(47,FILE='sedcnv7.out',STATUS='UNKNOWN')
          CLOSE(47,STATUS='DELETE')
          OPEN(47,FILE='sedcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(48,FILE='sedcnv8.out',STATUS='UNKNOWN')
          CLOSE(48,STATUS='DELETE')
          OPEN(48,FILE='sedcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(49,FILE='sedcnv9.out',STATUS='UNKNOWN')
          CLOSE(49,STATUS='DELETE')
          OPEN(49,FILE='sedcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=40+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE4,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.7) THEN
      IF (ISTRAN(7).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(51,FILE='sndcnv1.out',STATUS='UNKNOWN')
          CLOSE(51,STATUS='DELETE')
          OPEN(51,FILE='sndcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(52,FILE='sndcnv2.out',STATUS='UNKNOWN')
          CLOSE(52,STATUS='DELETE')
          OPEN(52,FILE='sndcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(53,FILE='sndcnv3.out',STATUS='UNKNOWN')
          CLOSE(53,STATUS='DELETE')
          OPEN(53,FILE='sndcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(54,FILE='sndcnv4.out',STATUS='UNKNOWN')
          CLOSE(54,STATUS='DELETE')
          OPEN(54,FILE='sndcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(55,FILE='sndcnv5.out',STATUS='UNKNOWN')
          CLOSE(55,STATUS='DELETE')
          OPEN(55,FILE='sndcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(56,FILE='sndcnv6.out',STATUS='UNKNOWN')
          CLOSE(56,STATUS='DELETE')
          OPEN(56,FILE='sndcnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(57,FILE='sndcnv7.out',STATUS='UNKNOWN')
          CLOSE(57,STATUS='DELETE')
          OPEN(57,FILE='sndcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(58,FILE='sndcnv8.out',STATUS='UNKNOWN')
          CLOSE(58,STATUS='DELETE')
          OPEN(58,FILE='sndcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(59,FILE='sndcnv9.out',STATUS='UNKNOWN')
          CLOSE(59,STATUS='DELETE')
          OPEN(59,FILE='sndcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=50+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE5,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.5) THEN
      IF (ISTRAN(5).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(61,FILE='toxcnv1.out',STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE='toxcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(62,FILE='toxcnv2.out',STATUS='UNKNOWN')
          CLOSE(62,STATUS='DELETE')
          OPEN(62,FILE='toxcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(63,FILE='toxcnv3.out',STATUS='UNKNOWN')
          CLOSE(63,STATUS='DELETE')
          OPEN(63,FILE='toxcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(64,FILE='toxcnv4.out',STATUS='UNKNOWN')
          CLOSE(64,STATUS='DELETE')
          OPEN(64,FILE='toxcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(65,FILE='toxcnv5.out',STATUS='UNKNOWN')
          CLOSE(65,STATUS='DELETE')
          OPEN(65,FILE='toxcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(66,FILE='toxcnv6.out',STATUS='UNKNOWN')
          CLOSE(66,STATUS='DELETE')
          OPEN(66,FILE='toxcnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(67,FILE='toxcnv7.out',STATUS='UNKNOWN')
          CLOSE(67,STATUS='DELETE')
          OPEN(67,FILE='toxcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(68,FILE='toxcnv8.out',STATUS='UNKNOWN')
          CLOSE(68,STATUS='DELETE')
          OPEN(68,FILE='toxcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(69,FILE='toxcnv9.out',STATUS='UNKNOWN')
          CLOSE(69,STATUS='DELETE')
          OPEN(69,FILE='toxcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=60+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE6,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      IF(ITMP.EQ.4) THEN
      IF (ISTRAN(4).GE.1) THEN
        IF (ISECSPV.GE.1) THEN
          OPEN(71,FILE='sflcnv1.out',STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE='sflcnv1.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.2) THEN
          OPEN(72,FILE='sflcnv2.out',STATUS='UNKNOWN')
          CLOSE(72,STATUS='DELETE')
          OPEN(72,FILE='sflcnv2.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.3) THEN
          OPEN(73,FILE='sflcnv3.out',STATUS='UNKNOWN')
          CLOSE(73,STATUS='DELETE')
          OPEN(73,FILE='sflcnv3.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.4) THEN
          OPEN(74,FILE='sflcnv4.out',STATUS='UNKNOWN')
          CLOSE(74,STATUS='DELETE')
          OPEN(74,FILE='sflcnv4.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.5) THEN
          OPEN(75,FILE='sflcnv5.out',STATUS='UNKNOWN')
          CLOSE(75,STATUS='DELETE')
          OPEN(75,FILE='sflcnv5.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.6) THEN
          OPEN(76,FILE='sflcnv6.out',STATUS='UNKNOWN')
          CLOSE(76,STATUS='DELETE')
          OPEN(76,FILE='sflcnv6.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.7) THEN
          OPEN(77,FILE='sflcnv7.out',STATUS='UNKNOWN')
          CLOSE(77,STATUS='DELETE')
          OPEN(77,FILE='sflcnv7.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.8) THEN
          OPEN(78,FILE='sflcnv8.out',STATUS='UNKNOWN')
          CLOSE(78,STATUS='DELETE')
          OPEN(78,FILE='sflcnv8.out',STATUS='UNKNOWN')
        END IF
        IF (ISECSPV.GE.9) THEN
          OPEN(79,FILE='sflcnv9.out',STATUS='UNKNOWN')
          CLOSE(79,STATUS='DELETE')
          OPEN(79,FILE='sflcnv9.out',STATUS='UNKNOWN')
        END IF
         DO IS=1,ISECSPV
         LUN=70+IS 
         LINES=NIJSPV(IS)
         LEVELS=KC
         WRITE (LUN,99) TITLE7,CCTITLE(LUN)
         WRITE (LUN,101)LINES,LEVELS
         WRITE (LUN,250)(ZZ(K),K=1,KC)
         CLOSE(LUN)
         END DO
      END IF
      END IF
C
      JSSPV(ITMP)=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      IF(ITMP.EQ.1.AND.ISTRAN(1).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(11,FILE='salcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(12,FILE='salcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(13,FILE='salcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(14,FILE='salcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(15,FILE='salcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(16,FILE='salcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(17,FILE='salcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(18,FILE='salcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(19,FILE='salcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=10+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         WRITE(LUN,250)(SAL(L,K),K=1,KC)
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.2.AND.ISTRAN(2).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(21,FILE='temcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(22,FILE='temcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(23,FILE='temcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(24,FILE='temcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
         IF (ISECSPV.GE.5)
     $    OPEN(25,FILE='temcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(26,FILE='temcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(27,FILE='temcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(28,FILE='temcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(29,FILE='temcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=20+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         WRITE(LUN,250)(TEM(L,K),K=1,KC)
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.3.AND.ISTRAN(3).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(31,FILE='dyecnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(32,FILE='dyecnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(33,FILE='dyecnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(34,FILE='dyecnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(35,FILE='dyecnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(36,FILE='dyecnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(37,FILE='dyecnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(38,FILE='dyecnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(39,FILE='dyecnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=30+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         WRITE(LUN,250)(DYE(L,K),K=1,KC)
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.6.AND.ISTRAN(6).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(41,FILE='sedcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(42,FILE='sedcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(43,FILE='sedcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(44,FILE='sedcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(45,FILE='sedcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(46,FILE='sedcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(47,FILE='sedcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(48,FILE='sedcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(49,FILE='sedcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=40+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         DO NSC=1,NSED
           WRITE(LUN,250)(SEDT(L,K),K=1,KC)
         END DO
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.7.AND.ISTRAN(7).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(51,FILE='sndcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(52,FILE='sndcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(53,FILE='sndcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(54,FILE='sndcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(55,FILE='sndcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(56,FILE='sndcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(57,FILE='sndcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(58,FILE='sndcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(59,FILE='sndcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=50+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         DO NSN=1,NSND
           WRITE(LUN,250)(SNDT(L,K),K=1,KC)
         END DO
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.5.AND.ISTRAN(5).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(61,FILE='toxcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(62,FILE='toxcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(63,FILE='toxcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(64,FILE='toxcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(65,FILE='toxcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(66,FILE='toxcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(67,FILE='toxcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(68,FILE='toxcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(69,FILE='toxcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=60+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         DO NT=1,NTOX
           WRITE(LUN,250)(TOX(L,K,1),K=1,KC)
         END DO
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
      IF(ITMP.EQ.4.AND.ISTRAN(4).GE.1) THEN
        IF (ISECSPV.GE.1)
     $    OPEN(71,FILE='sflcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.2)
     $    OPEN(72,FILE='sflcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.3)
     $    OPEN(73,FILE='sflcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.4)
     $    OPEN(74,FILE='sflcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.5)
     $    OPEN(75,FILE='sflcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.6)
     $    OPEN(76,FILE='sflcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.7)
     $    OPEN(77,FILE='sflcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.8)
     $    OPEN(78,FILE='sflcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        IF (ISECSPV.GE.9)
     $    OPEN(79,FILE='sflcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        DO IS=1,ISECSPV
        LUN=70+IS
        WRITE (LUN,100)N,TIME
         DO NN=1,NIJSPV(IS)
         I=ISPV(NN,IS)
         J=JSPV(NN,IS)
         L=LIJ(I,J)
         ZETA=P(L)*GI-SBPLTV(ITMP)*(HMP(L)+BELV(L))
         HBTMP=HMP(L)
C        HBTMP=SHPLTV*HMP(L)+SBPLTV(ITMP)*BELV(L)
C        WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
         WRITE(LUN,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
         WRITE(LUN,250)(SFL(L,K),K=1,KC)
         END DO
        CLOSE(LUN)
       END DO
      END IF
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(30E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
