C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE VELPLTV
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE VELPLTV WRITES A FIL FOR VERTICAL PLANE CONTOURING
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
      REAL VELN (KCM,100)
      REAL VELT (KCM,100)
      REAL WZ (KCM,100)
      CHARACTER*80 TITLE1,TITLE2
C
C**********************************************************************C
C
      IF (JSVPV.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
C **  WRITE HEADINGS
C
      LEVELS=KC
      TITLE1='INSTANTANEOUS NORMAL VELOCITY CONTOURS'
      TITLE2='INSTANTANEOUS TANGENTIAL VELOCITY VECTORS'
C
      IF (ISECVPV.GE.1) THEN
        OPEN(11,FILE='velcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='velvcv1.out',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        CLOSE(21,STATUS='DELETE')
        OPEN(11,FILE='velcnv1.out',STATUS='UNKNOWN')
        OPEN(21,FILE='velvcv1.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.2) THEN
        OPEN(12,FILE='velcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='velvcv2.out',STATUS='UNKNOWN')
        CLOSE(12,STATUS='DELETE')
        CLOSE(22,STATUS='DELETE')
        OPEN(12,FILE='velcnv2.out',STATUS='UNKNOWN')
        OPEN(22,FILE='velvcv2.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.3) THEN
        OPEN(13,FILE='velcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='velvcv3.out',STATUS='UNKNOWN')
        CLOSE(13,STATUS='DELETE')
        CLOSE(23,STATUS='DELETE')
        OPEN(13,FILE='velcnv3.out',STATUS='UNKNOWN')
        OPEN(23,FILE='velvcv3.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.4) THEN
        OPEN(14,FILE='velcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='velvcv4.out',STATUS='UNKNOWN')
        CLOSE(14,STATUS='DELETE')
        CLOSE(24,STATUS='DELETE')
        OPEN(14,FILE='velcnv4.out',STATUS='UNKNOWN')
        OPEN(24,FILE='velvcv4.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.5) THEN
        OPEN(15,FILE='velcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='velvcv5.out',STATUS='UNKNOWN')
        CLOSE(15,STATUS='DELETE')
        CLOSE(25,STATUS='DELETE')
        OPEN(15,FILE='velcnv5.out',STATUS='UNKNOWN')
        OPEN(25,FILE='velvcv5.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.6) THEN
        OPEN(16,FILE='velcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='velvcv6.out',STATUS='UNKNOWN')
        CLOSE(16,STATUS='DELETE')
        CLOSE(26,STATUS='DELETE')
        OPEN(16,FILE='velcnv6.out',STATUS='UNKNOWN')
        OPEN(26,FILE='velvcv6.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.7) THEN
        OPEN(17,FILE='velcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='velvcv7.out',STATUS='UNKNOWN')
        CLOSE(17,STATUS='DELETE')
        CLOSE(27,STATUS='DELETE')
        OPEN(17,FILE='velcnv7.out',STATUS='UNKNOWN')
        OPEN(27,FILE='velvcv7.out',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.8) THEN
        OPEN(18,FILE='velcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='velvcv8.out',STATUS='UNKNOWN')
        CLOSE(18,STATUS='DELETE')
        CLOSE(28,STATUS='DELETE')
        OPEN(18,FILE='velcnv8.out',STATUS='UNKNOWN')
        OPEN(28,FILE='velvcv8.out',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.9) THEN
        OPEN(19,FILE='velcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='velvcv9.out',STATUS='UNKNOWN')
        CLOSE(19,STATUS='DELETE')
        CLOSE(29,STATUS='DELETE')
        OPEN(19,FILE='velcnv9.out',STATUS='UNKNOWN')
        OPEN(29,FILE='velvcv9.out',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      LINES=NIJVPV(IS)
      WRITE (LUN1,99)TITLE1,CVTITLE(LUN1)
      WRITE (LUN2,99)TITLE2,CVTITLE(LUN2)
      WRITE (LUN1,101)LINES,LEVELS
      WRITE (LUN2,101)LINES,LEVELS
      WRITE (LUN1,250)(ZZ(K),K=1,KC)
      WRITE (LUN2,250)(ZZ(K),K=1,KC)
      CLOSE(LUN1)
      CLOSE(LUN2)
      END DO
C
      JSVPV=0
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
        OPEN(11,FILE='velcnv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(21,FILE='velvcv1.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.2) THEN
        OPEN(12,FILE='velcnv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(22,FILE='velvcv2.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.3) THEN
        OPEN(13,FILE='velcnv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(23,FILE='velvcv3.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.4) THEN
        OPEN(14,FILE='velcnv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(24,FILE='velvcv4.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.5) THEN
        OPEN(15,FILE='velcnv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(25,FILE='velvcv5.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.6) THEN
        OPEN(16,FILE='velcnv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(26,FILE='velvcv6.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.7) THEN
        OPEN(17,FILE='velcnv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(27,FILE='velvcv7.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF 
      IF (ISECVPV.GE.8) THEN
        OPEN(18,FILE='velcnv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(28,FILE='velvcv8.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
      IF (ISECVPV.GE.9) THEN
        OPEN(19,FILE='velcnv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(29,FILE='velvcv9.out',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
      DO IS=1,ISECVPV
      LUN1=10+IS
      LUN2=20+IS
      WRITE (LUN1,100)N,TIME
      WRITE (LUN2,100)N,TIME
      COSC=COS(PI*ANGVPV(IS)/180.)
      SINC=SIN(PI*ANGVPV(IS)/180.)
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       LN=LNC(L)
       LS=LSC(L)
        DO K=1,KC
        VELN(K,NN)=50.*((U(L+1,K)+U(L,K))*COSC+(V(LN,K)
     $                                         +V(L,K))*SINC)
        VELT(K,NN)=-50.*((U(L+1,K)+U(L,K))*SINC-(V(LN,K)
     $                                          +V(L,K))*COSC)
        WZ(K,NN)=50.*(W(L,K)+W(L,K-1))+GI*ZZ(K)*(DTI*(P(L)-P1(L))
     $         +50.*(U(L+1,K)*(P(L+1)-P(L))*DXIU(L+1)
     $              +U(L,K)*(P(L)-P(L-1))*DXIU(L)
     $              +V(LN,K)*(P(LN)-P(L))*DYIV(LN)
     $              +V(L,K)*(P(L)-P(LS))*DYIV(L)))
     $         +50.*(1.-ZZ(K))*(U(L+1,K)*(BELV(L+1)-BELV(L))*DXIU(L+1)
     $                         +U(L,K)*(BELV(L)-BELV(L-1))*DXIU(L)
     $                         +V(LN,K)*(BELV(LN)-BELV(L))*DYIV(LN)
     $                         +V(L,K)*(BELV(L)-BELV(LS))*DYIV(L))
c    $         -50.*(1.-ZZ(K))*(U(L+1,K)*(HMP(L+1)-HMP(L))*DXIU(L+1)
c    $                         +U(L,K)*(HMP(L)-HMP(L-1))*DXIU(L)
c    $                         +V(LN,K)*(HMP(LN)-HMP(L))*DYIV(LN)
c    $                         +V(L,K)*(HMP(L)-HMP(LS))*DYIV(L))
c
        END DO
       END DO
       DO NN=1,NIJVPV(IS)
       I=IVPV(NN,IS)
       J=JVPV(NN,IS)
       L=LIJ(I,J)
       ZETA=P(L)*GI-SBPLTV(1)*(HMP(L)+BELV(L))
       HBTMP=HMP(L)
C      HBTMP=SHPLTV*HMP(L)+SBPLTV*BELV(L)
C      WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
C      WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HMP(L)
       WRITE(LUN1,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN2,200)IL(L),JL(L),DLON(L),DLAT(L),ZETA,HBTMP
       WRITE(LUN1,250)(VELN(K,NN),K=1,KC)
       WRITE(LUN2,250)(VELT(K,NN),K=1,KC)
       WRITE(LUN2,250)(WZ(K,NN),K=1,KC)
       END DO
      CLOSE(LUN1)
      CLOSE(LUN2)
      END DO
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
