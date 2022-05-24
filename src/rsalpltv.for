C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RSALPLTV(ITMP)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE RSALPLTV WRITES A FILE FOR VERTICAL PLANE CONTOURING
C **  OF RESIDUAL SALINITY AND VERTICAL DIFFUSIVITY ALONG AN ARBITARY
C **  SEQUENCE OF (I,J) POINTS 
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE5
      DIMENSION ABTMP(KCM)
C
C**********************************************************************C
C
      IF(ITMP.EQ.2) RETURN
      IF(ITMP.EQ.3) RETURN
      IF(ITMP.EQ.4) RETURN
      IF(ITMP.GE.5) GO TO 1000
      IF (JSRSPV(ITMP).NE.1) GO TO 300
C
C**********************************************************************C
C
      TITLE1='RESIDUAL SALINITY CONTOURS'
      TITLE2='RESIDUAL VERTICAL DIFFUSIVITY CONTOURS'
      TITLE3='FLUX AVG RESID VERT DIFF CONTOURS'
C
      IF (ISECSPV.GE.1) THEN
        OPEN(11,FILE='rsalcv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='rviscv1.out',STATUS='UNKNOWN')
        OPEN(31,FILE='rvefcv1.out',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        CLOSE(31,STATUS='DELETE')
        OPEN(11,FILE='rsalcv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='rviscv1.out',STATUS='UNKNOWN')
        OPEN(31,FILE='rvefcv1.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.2) THEN
        OPEN(12,FILE='rsalcv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='rviscv2.out',STATUS='UNKNOWN')
        OPEN(32,FILE='rvefcv2.out',STATUS='UNKNOWN')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        CLOSE(32,STATUS='DELETE')
        OPEN(12,FILE='rsalcv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='rviscv2.out',STATUS='UNKNOWN')
        OPEN(32,FILE='rvefcv2.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        OPEN(13,FILE='rsalcv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='rviscv3.out',STATUS='UNKNOWN')
        OPEN(33,FILE='rvefcv3.out',STATUS='UNKNOWN')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        CLOSE(33,STATUS='DELETE')
        OPEN(13,FILE='rsalcv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='rviscv3.out',STATUS='UNKNOWN')
        OPEN(33,FILE='rvefcv3.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        OPEN(14,FILE='rsalcv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='rviscv4.out',STATUS='UNKNOWN')
        OPEN(34,FILE='rvefcv4.out',STATUS='UNKNOWN')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        CLOSE(34,STATUS='DELETE')
        OPEN(14,FILE='rsalcv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='rviscv4.out',STATUS='UNKNOWN')
        OPEN(34,FILE='rvefcv4.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.5) THEN
        OPEN(15,FILE='rsalcv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='rviscv5.out',STATUS='UNKNOWN')
        OPEN(35,FILE='rvefcv5.out',STATUS='UNKNOWN')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        CLOSE(35,STATUS='DELETE')
        OPEN(15,FILE='rsalcv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='rviscv5.out',STATUS='UNKNOWN')
        OPEN(35,FILE='rvefcv5.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        OPEN(16,FILE='rsalcv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='rviscv6.out',STATUS='UNKNOWN')
        OPEN(36,FILE='rvefcv6.out',STATUS='UNKNOWN')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        CLOSE(36,STATUS='DELETE')
        OPEN(16,FILE='rsalcv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='rviscv6.out',STATUS='UNKNOWN')
        OPEN(36,FILE='rvefcv6.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        OPEN(17,FILE='rsalcv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='rviscv7.out',STATUS='UNKNOWN')
        OPEN(37,FILE='rvefcv7.out',STATUS='UNKNOWN')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        CLOSE(37,STATUS='DELETE')
        OPEN(17,FILE='rsalcv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='rviscv7.out',STATUS='UNKNOWN')
        OPEN(37,FILE='rvefcv7.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.8) THEN
        OPEN(18,FILE='rsalcv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='rviscv8.out',STATUS='UNKNOWN')
        OPEN(38,FILE='rvefcv8.out',STATUS='UNKNOWN')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        CLOSE(38,STATUS='DELETE')
        OPEN(18,FILE='rsalcv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='rviscv8.out',STATUS='UNKNOWN')
        OPEN(38,FILE='rvefcv8.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        OPEN(19,FILE='rsalcv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='rviscv9.out',STATUS='UNKNOWN')
        OPEN(39,FILE='rvefcv9.out',STATUS='UNKNOWN')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        CLOSE(39,STATUS='DELETE')
        OPEN(19,FILE='rsalcv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='rviscv9.out',STATUS='UNKNOWN')
        OPEN(39,FILE='rvefcv9.out',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LINES=NIJSPV(IS)
      LEVELS=KC
      WRITE (LUN1,99) TITLE1,CCTITLE(LUN1)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,99) TITLE2,CCTITLE(LUN2)
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      WRITE (LUN3,99) TITLE3,CCTITLE(LUN2)
      WRITE (LUN3,101)LINES,LEVELS
      WRITE (LUN3,250)(ZZ(K),K=1,KC)
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
CREMOVE LATER      CLOSE(LUN5)
      END DO
C
      JSRSPV(ITMP)=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
      IF (ISECSPV.GE.1) THEN
        OPEN(11,FILE='rsalcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(21,FILE='rviscv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(31,FILE='rvefcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.2) THEN
        OPEN(12,FILE='rsalcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(22,FILE='rviscv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(32,FILE='rvefcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        OPEN(13,FILE='rsalcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(23,FILE='rviscv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(33,FILE='rvefcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        OPEN(14,FILE='rsalcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(24,FILE='rviscv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(34,FILE='rvefcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
       END IF 
      IF (ISECSPV.GE.5) THEN
        OPEN(15,FILE='rsalcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(25,FILE='rviscv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(35,FILE='rvefcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        OPEN(16,FILE='rsalcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(26,FILE='rviscv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(36,FILE='rvefcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        OPEN(17,FILE='rsalcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(27,FILE='rviscv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(37,FILE='rvefcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.8) THEN
        OPEN(18,FILE='rsalcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(28,FILE='rviscv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(38,FILE='rvefcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        OPEN(19,FILE='rsalcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(29,FILE='rviscv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(39,FILE='rvefcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      WRITE (LUN1,100)N,TIME
      WRITE (LUN2,100)N,TIME
      WRITE (LUN3,100)N,TIME
       DO NN=1,NIJSPV(IS)
       I=ISPV(NN,IS)
       J=JSPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(SALLPF(L,K),K=1,KC)
C      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
C      WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       ABTMP(1)=5000.*ABLPF(L,1)*HLPF(L)
       ABTMP(KC)=5000.*ABLPF(L,KS)*HLPF(L)
        DO K=2,KS
        ABTMP(K)=5000.*(ABLPF(L,K-1)+ABLPF(L,K))*HLPF(L)
        END DO
       WRITE(LUN2,250)(ABTMP(K),K=1,KC)
C      ABTMP(1)=5000.*ABEFF(L,1)*HLPF(L)
C      ABTMP(KC)=5000.*ABEFF(L,KS)*HLPF(L)
C       DO K=2,KS
C       ABTMP(K)=5000.*(ABEFF(L,K-1)+ABEFF(L,K))*HLPF(L)
C       END DO
       ABTMP(1)=-50.*ABEFF(L,1)
       ABTMP(KC)=50.*ABEFF(L,KS)
        DO K=2,KS
        ABTMP(K)=-50.*(ABEFF(L,K-1)+ABEFF(L,K))
        END DO
       WRITE(LUN3,250)(ABTMP(K),K=1,KC)
       END DO
      CLOSE(LUN1)
      CLOSE(LUN2)
      CLOSE(LUN3)
      END DO
C
      GO TO 2000
C
C**********************************************************************C
C
 1000 CONTINUE
C
      IF (JSRSPV(ITMP).NE.1) GO TO 1300

      TITLE1='RESIDUAL TOXIC CONTAMIANT CONTOURS'
      TITLE2='RESIDUAL COHESIVE SED CONTOURS'
      TITLE2='RESIDUAL NONCOHESIVE SED CONTOURS'
C
      IF (ISECSPV.GE.1) THEN
        IF(ITMP.EQ.5) OPEN(11,FILE='rtoxcv1.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(21,FILE='rsedcv1.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(31,FILE='rsndcv1.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(11,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(21,STATUS='DELETE')
        IF(ITMP.EQ.8) CLOSE(31,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(11,FILE='rtoxcv1.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(21,FILE='rsedcv1.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.8) OPEN(31,FILE='rsndcv1.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.2) THEN
        IF(ITMP.EQ.5) OPEN(12,FILE='rtoxcv2.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(22,FILE='rsedcv2.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(32,FILE='rsndcv2.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(12,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(22,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(32,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(12,FILE='rtoxcv2.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(22,FILE='rsedcv2.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(32,FILE='rsndcv2.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        IF(ITMP.EQ.5) OPEN(13,FILE='rtoxcv3.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(23,FILE='rsedcv3.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(33,FILE='rsndcv3.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5)  CLOSE(13,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(23,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(33,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(13,FILE='rtoxcv3.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(23,FILE='rsedcv3.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(33,FILE='rsndcv3.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        IF(ITMP.EQ.5) OPEN(14,FILE='rtoxcv4.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(24,FILE='rsedcv4.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(34,FILE='rsndcv4.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(14,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(24,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(34,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(14,FILE='rtoxcv4.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(24,FILE='rsedcv4.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(34,FILE='rsndcv4.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.5) THEN
        IF(ITMP.EQ.5) OPEN(15,FILE='rtoxcv5.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(25,FILE='rsedcv5.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(35,FILE='rsndcv5.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(15,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(25,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(35,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(15,FILE='rtoxcv5.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(25,FILE='rsedcv5.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(35,FILE='rsndcv5.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        IF(ITMP.EQ.5) OPEN(16,FILE='rtoxcv6.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(26,FILE='rsedcv6.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(36,FILE='rsndcv6.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(16,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(26,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(36,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(16,FILE='rtoxcv6.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(26,FILE='rsedcv6.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(36,FILE='rsndcv6.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        IF(ITMP.EQ.5) OPEN(17,FILE='rtoxcv7.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(27,FILE='rsedcv7.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(37,FILE='rsndcv7.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(17,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(27,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(37,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(17,FILE='rtoxcv7.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(27,FILE='rsedcv7.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(37,FILE='rsndcv7.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.8) THEN
        IF(ITMP.EQ.5) OPEN(18,FILE='rtoxcv8.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(28,FILE='rsedcv8.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(38,FILE='rsndcv8.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(18,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(28,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(38,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(18,FILE='rtoxcv8.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(28,FILE='rsedcv8.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(38,FILE='rsndcv8.out',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        IF(ITMP.EQ.5) OPEN(19,FILE='rtoxcv9.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(29,FILE='rsedcv9.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(39,FILE='rsndcv9.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.5) CLOSE(19,STATUS='DELETE')
        IF(ITMP.EQ.6) CLOSE(29,STATUS='DELETE')
        IF(ITMP.EQ.7) CLOSE(39,STATUS='DELETE')
        IF(ITMP.EQ.5) OPEN(19,FILE='rtoxcv9.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.6) OPEN(29,FILE='rsedcv9.out',STATUS='UNKNOWN')
        IF(ITMP.EQ.7) OPEN(39,FILE='rsndcv9.out',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN1=10+IS
      LUN2=20+IS
      LUN3=30+IS
      LINES=NIJSPV(IS)
      LEVELS=KC
      IF(ITMP.EQ.5) WRITE (LUN1,99) TITLE1,CCTITLE(LUN1)
      IF(ITMP.EQ.5) WRITE (LUN1,101)LINES,LEVELS
      IF(ITMP.EQ.5) WRITE (LUN1,250)(ZZ(K),K=1,KC)
      IF(ITMP.EQ.6) WRITE (LUN2,99) TITLE2,CCTITLE(LUN2)
      IF(ITMP.EQ.6) WRITE (LUN2,101)LINES,LEVELS
      IF(ITMP.EQ.6) WRITE (LUN2,250)(ZZ(K),K=1,KC)
      IF(ITMP.EQ.7) WRITE (LUN3,99) TITLE3,CCTITLE(LUN2)
      IF(ITMP.EQ.7) WRITE (LUN3,101)LINES,LEVELS
      IF(ITMP.EQ.7) WRITE (LUN3,250)(ZZ(K),K=1,KC)
      IF(ITMP.EQ.5) CLOSE(LUN1)
      IF(ITMP.EQ.6) CLOSE(LUN2)
      IF(ITMP.EQ.7) CLOSE(LUN3)
      END DO
C
      JSRSPV(ITMP)=0
C
C----------------------------------------------------------------------C
C
 1300 CONTINUE
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON    
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
C ** TOXICS
C
      IF (ITMP.EQ.5) THEN
C  
      IF (ISECSPV.GE.1) THEN
        OPEN(11,FILE='rtoxcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.2) THEN
        OPEN(12,FILE='rtoxcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        OPEN(13,FILE='rtoxcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        OPEN(14,FILE='rtoxcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.5) THEN
        OPEN(15,FILE='rtoxcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        OPEN(16,FILE='rtoxcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        OPEN(17,FILE='rtoxcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.8) THEN
        OPEN(18,FILE='rtoxcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        OPEN(19,FILE='rtoxcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN1=10+IS
      WRITE (LUN1,100)N,TIME
       DO NN=1,NIJSPV(IS)
       I=ISPV(NN,IS)
       J=JSPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(TOXLPF(L,K,1),K=1,KC)
       END DO
      CLOSE(LUN1)
      END DO
C
      END IF
C
C **  COHESIVE SEDIMENT
C
      IF (ITMP.EQ.6) THEN
C
      IF (ISECSPV.GE.1) THEN
        OPEN(21,FILE='rsedcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.2) THEN
        OPEN(22,FILE='rsedcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        OPEN(23,FILE='rsedcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        OPEN(24,FILE='rsedcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.5) THEN
        OPEN(25,FILE='rsedcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        OPEN(26,FILE='rsedcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        OPEN(27,FILE='rsedcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.8) THEN
        OPEN(28,FILE='rsedcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        OPEN(29,FILE='rsedcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN2=20+IS
      WRITE (LUN2,100)N,TIME
       DO NN=1,NIJSPV(IS)
       I=ISPV(NN,IS)
       J=JSPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,250)(SEDTLPF(L,K),K=1,KC)
       END DO
      CLOSE(LUN2)
      END DO
C
      END IF
C
C ** NONCHOESIVE SEDIMENT
C
      IF (ITMP.EQ.7) THEN
C
      IF (ISECSPV.GE.1) THEN
        OPEN(31,FILE='rsndcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.2) THEN
        OPEN(32,FILE='rsndcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.3) THEN
        OPEN(33,FILE='rsndcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.4) THEN
        OPEN(34,FILE='rsndcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.5) THEN
        OPEN(35,FILE='rsndcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.6) THEN
        OPEN(36,FILE='rsndcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.7) THEN
        OPEN(37,FILE='rsndcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECSPV.GE.8) THEN
        OPEN(38,FILE='rsndcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECSPV.GE.9) THEN
        OPEN(39,FILE='rsndcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECSPV
      LUN3=30+IS
      WRITE (LUN3,100)N,TIME
       DO NN=1,NIJSPV(IS)
       I=ISPV(NN,IS)
       J=JSPV(NN,IS)
       L=LIJ(I,J)
       ZETA=HLPF(L)-HMP(L)
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN3,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN3,250)(SNDTLPF(L,K),K=1,KC)
       END DO
      CLOSE(LUN3)
      END DO
C
      END IF
C
C**********************************************************************C
C
 2000 CONTINUE
C
C**********************************************************************C
C
   99 FORMAT(A40,2X,A20)
  100 FORMAT(I10,F12.4)
  101 FORMAT(2I10)
C 200 FORMAT(2I4,1X,10F12.6)
C 250 FORMAT(12F10.6)
  200 FORMAT(2I5,1X,6E14.6)
  250 FORMAT(12E12.4)
cmrm  200 FORMAT(2I5,1X,1p,6E13.5) 
cmrm  250 FORMAT(1p,12E11.3)
C
C**********************************************************************C
C
      RETURN
      END
