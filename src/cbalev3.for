C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CBALEV3 
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINES CBALEV CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, 
C **  AND ENERGY BALANCES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
C     DIMENSION CONT(LCM,KCM)
C
C**********************************************************************C
C
C **  ACCUMULATE INTERNAL SOURCES AND SINKS
C
C----------------------------------------------------------------------C
C
      DO L=2,LA
      VOLOUTE=VOLOUTE-QSUME(L)
      END DO
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
      DO LL=1,NQSIJ
      L=LQS(LL)
C
      PPEOUTE=PPEOUTE-QSS(K,LL)*G*( 0.5*(BELV(L)+BELV(L-1))
     $       +0.125*(HP(L)+H2P(L)+HP(L-1)+H2P(L-1))*(Z(K)+Z(K-1)) )
C
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
      IF (ISTRAN(1).GE.1) THEN
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=SAL(L,K)
       END DO 
      END DO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
       DO K=1,KC
       SALOUTE=SALOUTE
     $           -MAX(QSS(K,NS),0.)*CQS(K,NS,1)
     $           -MIN(QSS(K,NS),0.)*SAL(L,K) 
     $           -MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,1)
     $           -MIN(QSERT(K,NQSTMP),0.)*SAL(L,K)
       END DO
      END DO
C
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF (ID.EQ.0.AND.JD.EQ.0) THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      END IF
       DO K=1,KC
       SALOUTE=SALOUTE+QCTLT(K,NCTL)*CONT(LU,K)
     &         -RQWD*QCTLT(K,NCTL)*CONT(LU,K)
       END DO
      END DO
C
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       SALOUTE=SALOUTE+
     $   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC) THEN
         SALOUTE=SALOUTE-
     $     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))
     $      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
       END IF
      END DO
C
C
      END IF
C
C----------------------------------------------------------------------C
C
      IF (ISTRAN(3).GE.1) THEN
C
      DO K=1,KC
       DO L=2,LC
       CONT(L,K)=DYE(L,K)
       END DO 
      END DO
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
       DO K=1,KC
       DYEOUTE=DYEOUTE
     $           -MAX(QSS(K,NS),0.)*CQS(K,NS,3)
     $           -MIN(QSS(K,NS),0.)*DYE(L,K) 
     $           -MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,3)
     $           -MIN(QSERT(K,NQSTMP),0.)*DYE(L,K)
       END DO
      END DO
C
      DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF (ID.EQ.0.AND.JD.EQ.0) THEN
        LD=LC
        RQWD=0.
       ELSE
        LD=LIJ(ID,JD)
      END IF
       DO K=1,KC
      DYEOUTE=DYEOUTE+QCTLT(K,NCTL)*CONT(LU,K)
     &        -RQWD*QCTLT(K,NCTL)*CONT(LU,K)
       END DO
      END DO
C
      DO NWR=1,NQWR
       IU=IQWRU(NWR)
       JU=JQWRU(NWR)
       KU=KQWRU(NWR)
       ID=IQWRD(NWR)
       JD=JQWRD(NWR)
       KD=KQWRD(NWR)
       LU=LIJ(IU,JU)
       LD=LIJ(ID,JD)
       NQSTMP=NQWRSERQ(NWR)
       NCSTMP=NQWRSERQ(NWR)
       DYEOUTE=DYEOUTE+
     $   ( (QWR(NWR)+QWRSERT(NQSTMP))*CONT(LU,KU) )
       IF(LD.NE.1.OR.LD.NE.LC) THEN
         DYEOUTE=DYEOUTE-
     $     ( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))
     $      +QSERT(K,NQSTMP)*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
       END IF
      END DO
C
      END IF
C
C**********************************************************************C
C
      RETURN
      END
