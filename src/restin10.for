C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE RESTIN10
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE RESTINP READS A RESTART FILE GENERATED BY A 
C **  PRE SEPTEMBER 8, 1992 VERSION OF EFDC.FOR
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      OPEN(1,FILE='restart.inp',STATUS='UNKNOWN')
C
C**********************************************************************C
C
      READ(1,*,ERR=1000)NREST
C
      DO L=2,LA
      READ(1,*,ERR=1000)P(L),P1(L),UHDYE(L),UHDY1E(L),VHDXE(L),
     $                         VHDX1E(L)
      READ(1,*,ERR=1000)(U(L,K),K=1,KC)
      READ(1,*,ERR=1000)(U1(L,K),K=1,KC)
      READ(1,*,ERR=1000)(V(L,K),K=1,KC)
      READ(1,*,ERR=1000)(V1(L,K),K=1,KC)
      READ(1,*,ERR=1000)(W(L,K),K=1,KS)
      READ(1,*,ERR=1000)(W1(L,K),K=1,KS)
      READ(1,*,ERR=1000)(QQ(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQ1(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQL(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQL1(L,K),K=0,KC)
      READ(1,*,ERR=1000)(DML(L,K),K=0,KC)
      IF(ISCI(1).EQ.1)THEN
       READ(1,*,ERR=1000)(SAL(L,K),K=1,KC)
       READ(1,*,ERR=1000)(SAL1(L,K),K=1,KC)
      END IF
      IF(ISCI(2).EQ.1)THEN
       READ(1,*,ERR=1000)(TEM(L,K),K=1,KC)
       READ(1,*,ERR=1000)(TEM1(L,K),K=1,KC)
      END IF
      IF(ISCI(3).EQ.1)THEN
       READ(1,*,ERR=1000)(DYE(L,K),K=1,KC)
       READ(1,*,ERR=1000)(DYE1(L,K),K=1,KC)
      END IF
      IF(ISCI(4).EQ.1)THEN
       READ(1,*,ERR=1000) SEDB(L,1,1),(SED(L,K,1),K=1,KC)
       READ(1,*,ERR=1000) SEDB1(L,1,1),(SED1(L,K,1),K=1,KC)
      END IF
      IF(ISCI(5).EQ.1)THEN
       READ(1,*,ERR=1000)(SFL(L,K),K=1,KC)
       READ(1,*,ERR=1000)(SFL2(L,K),K=1,KC)
      END IF
      END DO
C
      IF(ISCI(1).EQ.1)THEN
       DO LL=1,NCBS
       L=LCBS(LL)
       READ(1,*,ERR=1000)(NLOS(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(CLOS(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBW
       L=LCBW(LL)      
       READ(1,*,ERR=1000)(NLOW(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(CLOW(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBE
       L=LCBE(LL)      
       READ(1,*,ERR=1000)(NLOE(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(CLOE(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBN
       L=LCBN(LL)
       READ(1,*,ERR=1000)(NLON(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(CLON(LL,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
       END DO
      END IF
C
      IF(ISCI(2).EQ.1)THEN
       DO LL=1,NCBS
       L=LCBS(LL)
       READ(1,*,ERR=1000)(NLOS(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(CLOS(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBW
       L=LCBW(LL)      
       READ(1,*,ERR=1000)(NLOW(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(CLOW(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBE
       L=LCBE(LL)      
       READ(1,*,ERR=1000)(NLOE(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(CLOE(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBN
       L=LCBN(LL)
       READ(1,*,ERR=1000)(NLON(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(CLON(LL,K,2),K=1,KC)
       READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
       END DO
      END IF
C
      IF(ISCI(3).EQ.1)THEN
       DO LL=1,NCBS
       L=LCBS(LL)
       READ(1,*,ERR=1000)(NLOS(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(CLOS(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBW
       L=LCBW(LL)      
       READ(1,*,ERR=1000)(NLOW(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(CLOW(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBE
       L=LCBE(LL)      
       READ(1,*,ERR=1000)(NLOE(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(CLOE(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBN
       L=LCBN(LL)
       READ(1,*,ERR=1000)(NLON(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(CLON(LL,K,3),K=1,KC)
       READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
       END DO
      END IF
C
      IF(ISCI(4).EQ.1)THEN
       DO LL=1,NCBS
       L=LCBS(LL)
       READ(1,*,ERR=1000)(NLOS(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(CLOS(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
       END DO
       DO LL=1,NCBW
       L=LCBW(LL)      
       READ(1,*,ERR=1000)(NLOW(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(CLOW(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
       END DO
       DO LL=1,NCBE
       L=LCBE(LL)      
       READ(1,*,ERR=1000)(NLOE(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(CLOE(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBN
       L=LCBN(LL)
       READ(1,*,ERR=1000)(NLON(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(CLON(LL,K,4),K=1,KC)
       READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
       READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
       END DO
      END IF
C
      IF(ISCI(5).EQ.1)THEN
       DO LL=1,NCBS
       L=LCBS(LL)
       READ(1,*,ERR=1000)(NLOS(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(CLOS(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBW
       L=LCBW(LL)      
       READ(1,*,ERR=1000)(NLOW(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(CLOW(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBE
       L=LCBE(LL)      
       READ(1,*,ERR=1000)(NLOE(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(CLOE(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDYWQ(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(UHDYWQ(LC,K),K=1,KC)
       END DO
       DO LL=1,NCBN
       L=LCBN(LL)
       READ(1,*,ERR=1000)(NLON(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(CLON(LL,K,5),K=1,KC)
       READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDXWQ(LC,K),K=1,KC)
       READ(1,*,ERR=1000)(VHDXWQ(LC,K),K=1,KC)
       END DO
      END IF
C
      CLOSE(1)
C
      DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      SED(1,K,1)=0.
      SFL(1,K)=0.
      CWQ(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SED1(1,K,1)=0.
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
      SFL(LC,K)=0.
      CWQ(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SED1(LC,K,1)=0.
      SFL2(LC,K)=0.
      CWQ2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXWQ(LC,K)=0.
      UHDYWQ(LC,K)=0.
      END DO
C
C**********************************************************************C
C
      DO L=2,LA
      LS=LSC(L)
      H1U(L)=0.5*GI*(P1(L)+P1(L-1))-0.5*(BELV(L)+BELV(L-1))
      H1V(L)=0.5*GI*(P1(L)+P1(LS))-0.5*(BELV(L)+BELV(LS))
      H1P(L)=GI*P1(L)-BELV(L)
      HU(L)=0.5*GI*(P(L)+P(L-1))-0.5*(BELV(L)+BELV(L-1))
      HV(L)=0.5*GI*(P(L)+P(LS))-0.5*(BELV(L)+BELV(LS))
      HP(L)=GI*P(L)-BELV(L)
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
      H2WQ(L)=HP(L)
      END DO
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
C **  WRITE READ ERRORS ON RESTART
C
      GO TO 1002
 1000 WRITE(6,1001)
 1001 FORMAT(1X,'READ ERROR ON FILE restart.inp ')
      STOP
 1002 CONTINUE
C
C**********************************************************************C
C
  907 FORMAT(12E12.4)
  908 FORMAT(12I10)
C
C**********************************************************************C
C
      RETURN
      END
