C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RVELPLTV
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
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
      REAL RW (KCM,100)
      REAL PVELN (KCM,100)
      REAL PVELT (KCM,100)
      REAL PW (KCM,100)
      REAL RLVELN (KCM,100)
      REAL RLVELT (KCM,100)
      REAL RLW (KCM,100)
      CHARACTER*80 TITLE10,TITLE20,TITLE30
      CHARACTER*80 TITLE40,TITLE50,TITLE60
      CHARACTER*80 TITLE70,TITLE80,TITLE90
C
C**********************************************************************C
C
      IF (JSRVPV.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      TITLE10='NORMAL ERT VELOCITY CONTOURS'
      TITLE20='NORMAL VPT VELOCITY CONTOURS'
      TITLE30='NORMAL MMT VELOCITY CONTOURS'
      TITLE40='TANGENTIAL ERT VELOCITY CONTOURS'
      TITLE50='TANGENTIAL VPT VELOCITY CONTOURS'
      TITLE60='TANGENTIAL MMT VELOCITY CONTOURS'
      TITLE70='TANGENTIAL ERT VELOCITY VECTORS'
      TITLE80='TANGENTIAL VPT VELOCITY VECTORS'
      TITLE90='TANGENTIAL MMT VELOCITY VECTORS'
C
      LEVELS=KC
C
      IF (ISECVPV.GE.1) THEN
        OPEN(11,FILE='rvlcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='pvlcnv1.out',STATUS='UNKNOWN')
        OPEN(31,FILE='mvlcnv1.out',STATUS='UNKNOWN')
        OPEN(41,FILE='rvlcvt1.out',STATUS='UNKNOWN')
        OPEN(51,FILE='pvlcvt1.out',STATUS='UNKNOWN')
        OPEN(61,FILE='mvlcvt1.out',STATUS='UNKNOWN')
        OPEN(71,FILE='rvlvcv1.out',STATUS='UNKNOWN')
        OPEN(81,FILE='pvlvcv1.out',STATUS='UNKNOWN')
        OPEN(91,FILE='mvlvcv1.out',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        CLOSE(31,STATUS='DELETE')
        CLOSE(41,STATUS='DELETE')
        CLOSE(51,STATUS='DELETE')
        CLOSE(61,STATUS='DELETE')
        CLOSE(71,STATUS='DELETE')
        CLOSE(81,STATUS='DELETE')
        CLOSE(91,STATUS='DELETE')
        OPEN(11,FILE='rvlcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='pvlcnv1.out',STATUS='UNKNOWN')
        OPEN(31,FILE='mvlcnv1.out',STATUS='UNKNOWN')
        OPEN(41,FILE='rvlcvt1.out',STATUS='UNKNOWN')
        OPEN(51,FILE='pvlcvt1.out',STATUS='UNKNOWN')
        OPEN(61,FILE='mvlcvt1.out',STATUS='UNKNOWN')
        OPEN(71,FILE='rvlvcv1.out',STATUS='UNKNOWN')
        OPEN(81,FILE='pvlvcv1.out',STATUS='UNKNOWN')
        OPEN(91,FILE='mvlvcv1.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.2) THEN
        OPEN(12,FILE='rvlcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='pvlcnv2.out',STATUS='UNKNOWN')
        OPEN(32,FILE='mvlcnv2.out',STATUS='UNKNOWN')
        OPEN(42,FILE='rvlcvt2.out',STATUS='UNKNOWN')
        OPEN(52,FILE='pvlcvt2.out',STATUS='UNKNOWN')
        OPEN(62,FILE='mvlcvt2.out',STATUS='UNKNOWN')
        OPEN(72,FILE='rvlvcv2.out',STATUS='UNKNOWN')
        OPEN(82,FILE='pvlvcv2.out',STATUS='UNKNOWN')
        OPEN(92,FILE='mvlvcv2.out',STATUS='UNKNOWN')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        CLOSE(32,STATUS='DELETE')
        CLOSE(42,STATUS='DELETE')
        CLOSE(52,STATUS='DELETE')
        CLOSE(62,STATUS='DELETE')
        CLOSE(72,STATUS='DELETE')
        CLOSE(82,STATUS='DELETE')
        CLOSE(92,STATUS='DELETE')
        OPEN(12,FILE='rvlcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='pvlcnv2.out',STATUS='UNKNOWN')
        OPEN(32,FILE='mvlcnv2.out',STATUS='UNKNOWN')
        OPEN(42,FILE='rvlcvt2.out',STATUS='UNKNOWN')
        OPEN(52,FILE='pvlcvt2.out',STATUS='UNKNOWN')
        OPEN(62,FILE='mvlcvt2.out',STATUS='UNKNOWN')
        OPEN(72,FILE='rvlvcv2.out',STATUS='UNKNOWN')
        OPEN(82,FILE='pvlvcv2.out',STATUS='UNKNOWN')
        OPEN(92,FILE='mvlvcv2.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.3) THEN
        OPEN(13,FILE='rvlcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='pvlcnv3.out',STATUS='UNKNOWN')
        OPEN(33,FILE='mvlcnv3.out',STATUS='UNKNOWN')
        OPEN(43,FILE='rvlcvt3.out',STATUS='UNKNOWN')
        OPEN(53,FILE='pvlcvt3.out',STATUS='UNKNOWN')
        OPEN(63,FILE='mvlcvt3.out',STATUS='UNKNOWN')
        OPEN(73,FILE='rvlvcv3.out',STATUS='UNKNOWN')
        OPEN(83,FILE='pvlvcv3.out',STATUS='UNKNOWN')
        OPEN(93,FILE='mvlvcv3.out',STATUS='UNKNOWN')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        CLOSE(33,STATUS='DELETE')
        CLOSE(43,STATUS='DELETE')
        CLOSE(53,STATUS='DELETE')
        CLOSE(63,STATUS='DELETE')
        CLOSE(73,STATUS='DELETE')
        CLOSE(83,STATUS='DELETE')
        CLOSE(93,STATUS='DELETE')
        OPEN(13,FILE='rvlcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='pvlcnv3.out',STATUS='UNKNOWN')
        OPEN(33,FILE='mvlcnv3.out',STATUS='UNKNOWN')
        OPEN(43,FILE='rvlcvt3.out',STATUS='UNKNOWN')
        OPEN(53,FILE='pvlcvt3.out',STATUS='UNKNOWN')
        OPEN(63,FILE='mvlcvt3.out',STATUS='UNKNOWN')
        OPEN(73,FILE='rvlvcv3.out',STATUS='UNKNOWN')
        OPEN(83,FILE='pvlvcv3.out',STATUS='UNKNOWN')
        OPEN(93,FILE='mvlvcv3.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.4) THEN
        OPEN(14,FILE='rvlcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='pvlcnv4.out',STATUS='UNKNOWN')
        OPEN(34,FILE='mvlcnv4.out',STATUS='UNKNOWN')
        OPEN(44,FILE='rvlcvt4.out',STATUS='UNKNOWN')
        OPEN(54,FILE='pvlcvt4.out',STATUS='UNKNOWN')
        OPEN(64,FILE='mvlcvt4.out',STATUS='UNKNOWN')
        OPEN(74,FILE='rvlvcv4.out',STATUS='UNKNOWN')
        OPEN(84,FILE='pvlvcv4.out',STATUS='UNKNOWN')
        OPEN(94,FILE='mvlvcv4.out',STATUS='UNKNOWN')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        CLOSE(34,STATUS='DELETE')
        CLOSE(44,STATUS='DELETE')
        CLOSE(54,STATUS='DELETE')
        CLOSE(64,STATUS='DELETE')
        CLOSE(74,STATUS='DELETE')
        CLOSE(84,STATUS='DELETE')
        CLOSE(94,STATUS='DELETE')
        OPEN(14,FILE='rvlcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='pvlcnv4.out',STATUS='UNKNOWN')
        OPEN(34,FILE='mvlcnv4.out',STATUS='UNKNOWN')
        OPEN(44,FILE='rvlcvt4.out',STATUS='UNKNOWN')
        OPEN(54,FILE='pvlcvt4.out',STATUS='UNKNOWN')
        OPEN(64,FILE='mvlcvt4.out',STATUS='UNKNOWN')
        OPEN(74,FILE='rvlvcv4.out',STATUS='UNKNOWN')
        OPEN(84,FILE='pvlvcv4.out',STATUS='UNKNOWN')
        OPEN(94,FILE='mvlvcv4.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.5) THEN
        OPEN(15,FILE='rvlcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='pvlcnv5.out',STATUS='UNKNOWN')
        OPEN(35,FILE='mvlcnv5.out',STATUS='UNKNOWN')
        OPEN(45,FILE='rvlcvt5.out',STATUS='UNKNOWN')
        OPEN(55,FILE='pvlcvt5.out',STATUS='UNKNOWN')
        OPEN(65,FILE='mvlcvt5.out',STATUS='UNKNOWN')
        OPEN(75,FILE='rvlvcv5.out',STATUS='UNKNOWN')
        OPEN(85,FILE='pvlvcv5.out',STATUS='UNKNOWN')
        OPEN(95,FILE='mvlvcv5.out',STATUS='UNKNOWN')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        CLOSE(35,STATUS='DELETE')
        CLOSE(45,STATUS='DELETE')
        CLOSE(55,STATUS='DELETE')
        CLOSE(65,STATUS='DELETE')
        CLOSE(75,STATUS='DELETE')
        CLOSE(85,STATUS='DELETE')
        CLOSE(95,STATUS='DELETE')
        OPEN(15,FILE='rvlcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='pvlcnv5.out',STATUS='UNKNOWN')
        OPEN(35,FILE='mvlcnv5.out',STATUS='UNKNOWN')
        OPEN(45,FILE='rvlcvt5.out',STATUS='UNKNOWN')
        OPEN(55,FILE='pvlcvt5.out',STATUS='UNKNOWN')
        OPEN(65,FILE='mvlcvt5.out',STATUS='UNKNOWN')
        OPEN(75,FILE='rvlvcv5.out',STATUS='UNKNOWN')
        OPEN(85,FILE='pvlvcv5.out',STATUS='UNKNOWN')
        OPEN(95,FILE='mvlvcv5.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.6) THEN
        OPEN(16,FILE='rvlcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='pvlcnv6.out',STATUS='UNKNOWN')
        OPEN(36,FILE='mvlcnv6.out',STATUS='UNKNOWN')
        OPEN(46,FILE='rvlcvt6.out',STATUS='UNKNOWN')
        OPEN(56,FILE='pvlcvt6.out',STATUS='UNKNOWN')
        OPEN(66,FILE='mvlcvt6.out',STATUS='UNKNOWN')
        OPEN(76,FILE='rvlvcv6.out',STATUS='UNKNOWN')
        OPEN(86,FILE='pvlvcv6.out',STATUS='UNKNOWN')
        OPEN(96,FILE='mvlvcv6.out',STATUS='UNKNOWN')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        CLOSE(36,STATUS='DELETE')
        CLOSE(46,STATUS='DELETE')
        CLOSE(56,STATUS='DELETE')
        CLOSE(66,STATUS='DELETE')
        CLOSE(76,STATUS='DELETE')
        CLOSE(86,STATUS='DELETE')
        CLOSE(96,STATUS='DELETE')
        OPEN(16,FILE='rvlcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='pvlcnv6.out',STATUS='UNKNOWN')
        OPEN(36,FILE='mvlcnv6.out',STATUS='UNKNOWN')
        OPEN(46,FILE='rvlcvt6.out',STATUS='UNKNOWN')
        OPEN(56,FILE='pvlcvt6.out',STATUS='UNKNOWN')
        OPEN(66,FILE='mvlcvt6.out',STATUS='UNKNOWN')
        OPEN(76,FILE='rvlvcv6.out',STATUS='UNKNOWN')
        OPEN(86,FILE='pvlvcv6.out',STATUS='UNKNOWN')
        OPEN(96,FILE='mvlvcv6.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.7) THEN
        OPEN(17,FILE='rvlcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='pvlcnv7.out',STATUS='UNKNOWN')
        OPEN(37,FILE='mvlcnv7.out',STATUS='UNKNOWN')
        OPEN(47,FILE='rvlcvt7.out',STATUS='UNKNOWN')
        OPEN(57,FILE='pvlcvt7.out',STATUS='UNKNOWN')
        OPEN(67,FILE='mvlcvt7.out',STATUS='UNKNOWN')
        OPEN(77,FILE='rvlvcv7.out',STATUS='UNKNOWN')
        OPEN(87,FILE='pvlvcv7.out',STATUS='UNKNOWN')
        OPEN(97,FILE='mvlvcv7.out',STATUS='UNKNOWN')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        CLOSE(37,STATUS='DELETE')
        CLOSE(47,STATUS='DELETE')
        CLOSE(57,STATUS='DELETE')
        CLOSE(67,STATUS='DELETE')
        CLOSE(77,STATUS='DELETE')
        CLOSE(87,STATUS='DELETE')
        CLOSE(97,STATUS='DELETE')
        OPEN(17,FILE='rvlcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='pvlcnv7.out',STATUS='UNKNOWN')
        OPEN(37,FILE='mvlcnv7.out',STATUS='UNKNOWN')
        OPEN(47,FILE='rvlcvt7.out',STATUS='UNKNOWN')
        OPEN(57,FILE='pvlcvt7.out',STATUS='UNKNOWN')
        OPEN(67,FILE='mvlcvt7.out',STATUS='UNKNOWN')
        OPEN(77,FILE='rvlvcv7.out',STATUS='UNKNOWN')
        OPEN(87,FILE='pvlvcv7.out',STATUS='UNKNOWN')
        OPEN(97,FILE='mvlvcv7.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.8) THEN
        OPEN(18,FILE='rvlcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='pvlcnv8.out',STATUS='UNKNOWN')
        OPEN(38,FILE='mvlcnv8.out',STATUS='UNKNOWN')
        OPEN(48,FILE='rvlcvt8.out',STATUS='UNKNOWN')
        OPEN(58,FILE='pvlcvt8.out',STATUS='UNKNOWN')
        OPEN(68,FILE='mvlcvt8.out',STATUS='UNKNOWN')
        OPEN(78,FILE='rvlvcv8.out',STATUS='UNKNOWN')
        OPEN(88,FILE='pvlvcv8.out',STATUS='UNKNOWN')
        OPEN(98,FILE='mvlvcv8.out',STATUS='UNKNOWN')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        CLOSE(38,STATUS='DELETE')
        CLOSE(48,STATUS='DELETE')
        CLOSE(58,STATUS='DELETE')
        CLOSE(68,STATUS='DELETE')
        CLOSE(78,STATUS='DELETE')
        CLOSE(88,STATUS='DELETE')
        CLOSE(98,STATUS='DELETE')
        OPEN(18,FILE='rvlcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='pvlcnv8.out',STATUS='UNKNOWN')
        OPEN(38,FILE='mvlcnv8.out',STATUS='UNKNOWN')
        OPEN(48,FILE='rvlcvt8.out',STATUS='UNKNOWN')
        OPEN(58,FILE='pvlcvt8.out',STATUS='UNKNOWN')
        OPEN(68,FILE='mvlcvt8.out',STATUS='UNKNOWN')
        OPEN(78,FILE='rvlvcv8.out',STATUS='UNKNOWN')
        OPEN(88,FILE='pvlvcv8.out',STATUS='UNKNOWN')
        OPEN(98,FILE='mvlvcv8.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.9) THEN
        OPEN(19,FILE='rvlcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='pvlcnv9.out',STATUS='UNKNOWN')
        OPEN(39,FILE='mvlcnv9.out',STATUS='UNKNOWN')
        OPEN(49,FILE='rvlcvt9.out',STATUS='UNKNOWN')
        OPEN(59,FILE='pvlcvt9.out',STATUS='UNKNOWN')
        OPEN(69,FILE='mvlcvt9.out',STATUS='UNKNOWN')
        OPEN(79,FILE='rvlvcv9.out',STATUS='UNKNOWN')
        OPEN(89,FILE='pvlvcv9.out',STATUS='UNKNOWN')
        OPEN(99,FILE='mvlvcv9.out',STATUS='UNKNOWN')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        CLOSE(39,STATUS='DELETE')
        CLOSE(49,STATUS='DELETE')
        CLOSE(59,STATUS='DELETE')
        CLOSE(69,STATUS='DELETE')
        CLOSE(79,STATUS='DELETE')
        CLOSE(89,STATUS='DELETE')
        CLOSE(99,STATUS='DELETE')
        OPEN(19,FILE='rvlcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='pvlcnv9.out',STATUS='UNKNOWN')
        OPEN(39,FILE='mvlcnv9.out',STATUS='UNKNOWN')
        OPEN(49,FILE='rvlcvt9.out',STATUS='UNKNOWN')
        OPEN(59,FILE='pvlcvt9.out',STATUS='UNKNOWN')
        OPEN(69,FILE='mvlcvt9.out',STATUS='UNKNOWN')
        OPEN(79,FILE='rvlvcv9.out',STATUS='UNKNOWN')
        OPEN(89,FILE='pvlvcv9.out',STATUS='UNKNOWN')
        OPEN(99,FILE='mvlvcv9.out',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LUN4=40+IS
      LUN5=50+IS
      LUN6=60+IS
      LUN7=70+IS
      LUN8=80+IS
      LUN9=90+IS
      LINES=NIJVPV(IS)
      WRITE (LUN1,99) TITLE10,CVTITLE(LUN1)
      WRITE (LUN2,99) TITLE20,CVTITLE(LUN2)
      WRITE (LUN3,99) TITLE30,CVTITLE(LUN3)
      WRITE (LUN4,99) TITLE40,CVTITLE(LUN4)
      WRITE (LUN5,99) TITLE50,CVTITLE(LUN5)
      WRITE (LUN6,99) TITLE60,CVTITLE(LUN6)
      WRITE (LUN7,99) TITLE70,CVTITLE(LUN7)
      WRITE (LUN8,99) TITLE80,CVTITLE(LUN8)
      WRITE (LUN9,99) TITLE90,CVTITLE(LUN9)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN3,101)LINES,LEVELS
      WRITE (LUN4,101)LINES,LEVELS
      WRITE (LUN5,101)LINES,LEVELS
      WRITE (LUN6,101)LINES,LEVELS
      WRITE (LUN7,101)LINES,LEVELS
      WRITE (LUN8,101)LINES,LEVELS
      WRITE (LUN9,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      WRITE (LUN3,250)(ZZ(K),K=1,KC)
      WRITE (LUN4,250)(ZZ(K),K=1,KC)
      WRITE (LUN5,250)(ZZ(K),K=1,KC)
      WRITE (LUN6,250)(ZZ(K),K=1,KC)
      WRITE (LUN7,250)(ZZ(K),K=1,KC)
      WRITE (LUN8,250)(ZZ(K),K=1,KC)
      WRITE (LUN9,250)(ZZ(K),K=1,KC)
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      CLOSE(LUN4)
      CLOSE(LUN5)
      CLOSE(LUN6)
      CLOSE(LUN7)
      CLOSE(LUN8)
      CLOSE(LUN9)
      END DO
C
      JSRVPV=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      IF (ISECVPV.GE.1) THEN
        OPEN(11,FILE='rvlcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(21,FILE='pvlcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(31,FILE='mvlcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(41,FILE='rvlcvt1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(51,FILE='pvlcvt1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(61,FILE='mvlcvt1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(71,FILE='rvlvcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(81,FILE='pvlvcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(91,FILE='mvlvcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.2) THEN
        OPEN(12,FILE='rvlcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(22,FILE='pvlcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(32,FILE='mvlcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(42,FILE='rvlcvt2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(52,FILE='pvlcvt2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(62,FILE='mvlcvt2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(72,FILE='rvlvcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(82,FILE='pvlvcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(92,FILE='mvlvcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.3) THEN
        OPEN(13,FILE='rvlcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(23,FILE='pvlcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(33,FILE='mvlcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(43,FILE='rvlcvt3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(53,FILE='pvlcvt3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(63,FILE='mvlcvt3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(73,FILE='rvlvcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(83,FILE='pvlvcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(93,FILE='mvlvcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.4) THEN
        OPEN(14,FILE='rvlcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(24,FILE='pvlcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(34,FILE='mvlcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(44,FILE='rvlcvt4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(54,FILE='pvlcvt4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(64,FILE='mvlcvt4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(74,FILE='rvlvcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(84,FILE='pvlvcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(94,FILE='mvlvcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.5) THEN
        OPEN(15,FILE='rvlcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(25,FILE='pvlcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(35,FILE='mvlcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(45,FILE='rvlcvt5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(55,FILE='pvlcvt5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(65,FILE='mvlcvt5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(75,FILE='rvlvcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(85,FILE='pvlvcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(95,FILE='mvlvcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.6) THEN
        OPEN(16,FILE='rvlcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(26,FILE='pvlcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(36,FILE='mvlcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(46,FILE='rvlcvt6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(56,FILE='pvlcvt6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(66,FILE='mvlcvt6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(76,FILE='rvlvcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(86,FILE='pvlvcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(96,FILE='mvlvcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.7) THEN
        OPEN(17,FILE='rvlcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(27,FILE='pvlcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(37,FILE='mvlcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(47,FILE='rvlcvt7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(57,FILE='pvlcvt7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(67,FILE='mvlcvt7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(77,FILE='rvlvcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(87,FILE='pvlvcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(97,FILE='mvlvcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.8) THEN
        OPEN(18,FILE='rvlcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(28,FILE='pvlcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(38,FILE='mvlcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(48,FILE='rvlcvt8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(58,FILE='pvlcvt8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(68,FILE='mvlcvt8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(78,FILE='rvlvcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(88,FILE='pvlvcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(98,FILE='mvlvcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.9) THEN
        OPEN(19,FILE='rvlcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(29,FILE='pvlcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(39,FILE='mvlcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(49,FILE='rvlcvt9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(59,FILE='pvlcvt9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(69,FILE='mvlcvt9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(79,FILE='rvlvcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(89,FILE='pvlvcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(99,FILE='mvlvcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LUN4=40+IS
      LUN5=50+IS
      LUN6=60+IS
      LUN7=70+IS
      LUN8=80+IS
      LUN9=90+IS
      WRITE (LUN1,100)N,TIME
      WRITE (LUN2,100)N,TIME
      WRITE (LUN3,100)N,TIME
      WRITE (LUN4,100)N,TIME
      WRITE (LUN5,100)N,TIME
      WRITE (LUN6,100)N,TIME
      WRITE (LUN7,100)N,TIME
      WRITE (LUN8,100)N,TIME
      WRITE (LUN9,100)N,TIME
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       LN=LNC(L)
       LS=LSC(L)
        DO K=1,KC
        RVELN(K,NN)=50.*((UHLPF(L+1,K)+UHLPF(L,K))*COSC
     $            +(VHLPF(LN,K)+VHLPF(L,K))*SINC)/HLPF(L)
        RVELT(K,NN)=-50.*((UHLPF(L+1,K)+UHLPF(L,K))*SINC
     $            -(VHLPF(LN,K)+VHLPF(L,K))*COSC)/HLPF(L)
        RW(K,NN)=50.*(WLPF(L,K)+WLPF(L,K-1))
        PVELN(K,NN)=50.*((UVPT(L+1,K)+UVPT(L,K))*COSC
     $            +(VVPT(LN,K)+VVPT(L,K))*SINC)/HLPF(L)
        PVELT(K,NN)=-50.*((UVPT(L+1,K)+UVPT(L,K))*SINC
     $            -(VVPT(LN,K)+VVPT(L,K))*COSC)/HLPF(L)
        PW(K,NN)=50.*(WVPT(L,K)+WVPT(L,K-1))
        RLVELN(K,NN)=RVELN(K,NN)+PVELN(K,NN)
        RLVELT(K,NN)=RVELT(K,NN)+PVELT(K,NN)
        RLW(K,NN)=RW(K,NN)+PW(K,NN)
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
c      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN6,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
c      WRITE(LUN9,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN4,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN5,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN6,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN7,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN8,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN9,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(RVELN(K,NN),K=1,KC)
       WRITE(LUN2,250)(PVELN(K,NN),K=1,KC)
       WRITE(LUN3,250)(RLVELN(K,NN),K=1,KC)
       WRITE(LUN4,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN5,250)(PVELT(K,NN),K=1,KC)
       WRITE(LUN6,250)(RLVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RVELT(K,NN),K=1,KC)
       WRITE(LUN8,250)(PVELT(K,NN),K=1,KC)
       WRITE(LUN9,250)(RLVELT(K,NN),K=1,KC)
       WRITE(LUN7,250)(RW(K,NN),K=1,KC)
       WRITE(LUN8,250)(PW(K,NN),K=1,KC)
       WRITE(LUN9,250)(RLW(K,NN),K=1,KC)
       END DO
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      CLOSE(LUN4)
      CLOSE(LUN5)
      CLOSE(LUN6)
      CLOSE(LUN7)
      CLOSE(LUN8)
      CLOSE(LUN9)
      END DO
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
cmrm  200 FORMAT(2I5,1X,6E14.6)
cmrm  250 FORMAT(12E12.4)
  200 FORMAT(2I5,1X,1p,6E13.5)
  250 FORMAT(1p,20E11.3)
C
C**********************************************************************C
C
      RETURN
      END
