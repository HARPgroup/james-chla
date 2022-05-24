C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALFQC (ISTL,MVAR,M,CON,CON1)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
C
C **  SUBROUTINE CALFQC CALCULATES MASS SOURCES AND SINKS ASSOCIATED
C **  WITH CONSTANT AND TIME SERIES INFLOWS AND OUTFLOWS; CONTROL
C **  STRUCTURE INFLOWS AND OUTLOWS; WITHDRAWAL AND RETURN STRUCTURE
C **  OUTFLOWS; AND  EMBEDED CHANNEL INFLOWS AND OUTFLOWS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION CON(LCM,KCM),CON1(LCM,KCM)
C
C**********************************************************************C
C
C **  VOLUMETRIC SOURCE-SINK FLUXES
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC
       FQC(L,K)=0.
       END DO
      END DO
C
      IF(MVAR.EQ.4) GO TO 1000
      IF(MVAR.EQ.8) GO TO 1500
C
C**********************************************************************C
C
      IF(ISTL.EQ.2) THEN
C     CALCULATE FOR TWO TIME LEVELS
C
C----------------------------------------------------------------------C
C
C **  STANDARD VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       FQC(L,K)=FQC(L,K)                    ! Discharge C out? when Q<0, Ji, 8/3/99
     $         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     $         +MIN(QSS(K,NS),0.)*CON1(L,K) !CON1 is from call statement                                      
     $         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     $         +MIN(QSERT(K,NQSTMP),0.)*CON1(L,K)
C    $         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
C    $         +MIN(QSS(K,NS),0.)*0.5*(CON(L,K)+CON1(L,K))
C    $         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
C    $         +MIN(QSERT(K,NQSTMP),0.)*0.5*(CON(L,K)+CON1(L,K))
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (2TL)
C
      IF(NQJPIJ.GT.0) THEN
      DO NJP=1,NQJPIJ
      IF(ICALJP(NJP).GT.0) THEN
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
       QVJPTMP=0.
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       END DO
       FQC(LJP,KTMP)=FQC(LJP,KTMP)
     $         +QQCJP(NJP)*CQCJP(1,NJP,M)
     $         +QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)
      END IF
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  CONTROL STRUCTURES (2TL)
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
       FQC(LU,K)=FQC(LU,K)
     $          -QCTLT(K,NCTL)*CON1(LU,K)
       FQC(LD,K)=FQC(LD,K)
     $          +RQWD*QCTLT(K,NCTL)*CON1(LU,K)
C      FQC(LU,K)=FQC(LU,K)
C    $          -QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K))
C      FQC(LD,K)=FQC(LD,K)
C    $          +RQWD*QCTLT(K,NCTL)*0.5*(CON(LU,K)+CON1(LU,K))
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  WITHDRAWAL CONCENTRATION RISE RETURN (2TL)
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
       FQC(LU,KU)=FQC(LU,KU)
     $          -(QWR(NWR)+QWRSERT(NQSTMP))*CON1(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)
     $          +QWR(NWR)*(CON1(LU,KU)+CQWR(NWR,M))
     $          +QWRSERT(NQSTMP)*(CON1(LU,KU)+CQWRSERT(NCSTMP,M))
C      FQC(LU,KU)=FQC(LU,K)
C    $          -(QWR(K,NWR)+QSERT(K,NQSTMP))*0.5*(CON(LU,K)
C    $                                            +CON1(LU,K))
C      FQC(LD,KD)=FQC(LD,K)
C    $          +QWR(K,NWR)*(0.5*(CON(LU,K)+CON1(LU,K))+CQWR(K,NWR,M))
C    $          +QSERT(K,NQSTMP)*(0.5*(CON(LU,K)+CON1(LU,K))
C    $                                +CSERT(K,NCSTMP,M))
      END DO
C
C----------------------------------------------------------------------C
C
C **  SUBGRID SCALE CHANNEL EXCHANGE (2TL)
C
      IF (MDCHH.GE.1) THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF (MDCHTYP(NMD).EQ.1) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        END IF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
     $               +MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               +MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
     $               +MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               +MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
     $               -MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               -MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
     $               -MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               -MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
C    $             +0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
C    $             +0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
C    $             +0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
C    $             +0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
C       FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
C    $             -0.5*MAX(QUKTMP,0.)*(CON1(LMDCHUT,K)+CON(LMDCHUT,K))
C    $             -0.5*MIN(QUKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
C       FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
C    $             -0.5*MAX(QVKTMP,0.)*(CON1(LMDCHVT,K)+CON(LMDCHVT,K))
C    $             -0.5*MIN(QVKTMP,0.)*(CON1(LMDCHHT,K)+CON(LMDCHHT,K))
        END DO
        END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  GROUNDWATER, EVAP, RAINFALL (2TL)
C
      IF (ISGWIE.NE.0) THEN
        DO L=2,LA
         FQC(L,1)=FQC(L,1)-RIFTR(L)*CON1(L,1)
        END DO
      END IF
C
      IF(M.EQ.2) THEN
        IF(ISTOPT(2).EQ.0.OR.ISTOPT(2).EQ.3) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TEMO*DXYP(L)
          END DO
        END IF
        IF(ISTOPT(2).EQ.1.OR.ISTOPT(2).EQ.2) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
          END DO
        END IF
      END IF
C
      IF(M.EQ.2) THEN
        IF(ISTOPT(2).EQ.0) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)-EVAPSW(L)*CON1(L,KC)
          END DO
        END IF
      END IF
C
C**********************************************************************C
C
      ELSE
C     CALCULATE FOR THREE TIME LEVELS
C
C----------------------------------------------------------------------C
C
C **  STANDARD VOLUMETRICS SOURCE SINK LOCATIONS (3TL)
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       FQC(L,K)=FQC(L,K)
     $         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     $         +MIN(QSS(K,NS),0.)*CON1(L,K)
     $         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     $         +MIN(QSERT(K,NQSTMP),0.)*CON1(L,K)
C    $         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
C    $         +MIN(QSS(K,NS),0.)*CON(L,K)
C    $         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
C    $         +MIN(QSERT(K,NQSTMP),0.)*CON(L,K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  JET-PLUME VOLUMETRICS SOURCE SINK LOCATIONS (3TL)
C
      IF(NQJPIJ.GT.0) THEN
      DO NJP=1,NQJPIJ
      IF(ICALJP(NJP).GT.0) THEN
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
       QVJPTMP=0.
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       END DO
       FQC(LJP,KTMP)=FQC(LJP,KTMP)
     $         +QQCJP(NJP)*CQCJP(1,NJP,M)
     $         +QVJPTMP*CSERT(1,NCSERJP(NJP,M),M)
      END IF
      END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  CONTROL STRUCTURES (3TL)
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
       FQC(LU,K)=FQC(LU,K)
     $          -QCTLT(K,NCTL)*CON1(LU,K)
       FQC(LD,K)=FQC(LD,K)
     $          +RQWD*QCTLT(K,NCTL)*CON1(LU,K)
C      FQC(LU,K)=FQC(LU,K)
C    $          -QCTLT(K,NCTL)*CON(LU,K)
C      FQC(LD,K)=FQC(LD,K)
C    $          +RQWD*QCTLT(K,NCTL)*CON(LU,K)
       END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  WITHDRAWAL CONCENTRATION RISE RETURN (3TL)
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
       FQC(LU,KU)=FQC(LU,KU)
     $          -(QWR(NWR)+QWRSERT(NQSTMP))*CON1(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)
     $          +QWR(NWR)*(CON1(LU,KU)+CQWR(NWR,M))
     $          +QWRSERT(NQSTMP)*(CON1(LU,KU)+CQWRSERT(NCSTMP,M))
C      FQC(LU,K)=FQC(LU,K)
C    $          -(QWR(K,NWR)+QSERT(K,NQSTMP))*CON(LU,K)
C      FQC(LD,K)=FQC(LD,K)
C    $          +QWR(K,NWR)*(CON(LU,K)+CQWR(K,NWR,M))
C    $          +QSERT(K,NQSTMP)*(CON(LU,K)+CSERT(K,NCSTMP,M))
      END DO
C
C----------------------------------------------------------------------C
C
C **  SUBGRID SCALE CHANNEL EXCHANGE (3TL)
C
      IF (MDCHH.GE.1) THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF (MDCHTYP(NMD).EQ.1) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        END IF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
     $               +MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               +MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
     $               +MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               +MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
     $               -MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               -MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
     $               -MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               -MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
C       FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
C    $               +MAX(QUKTMP,0.)*CON(LMDCHUT,K)
C    $               +MIN(QUKTMP,0.)*CON(LMDCHHT,K)
C    $               +MAX(QVKTMP,0.)*CON(LMDCHVT,K)
C    $               +MIN(QVKTMP,0.)*CON(LMDCHHT,K)
C       FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
C    $               -MAX(QUKTMP,0.)*CON(LMDCHUT,K)
C    $               -MIN(QUKTMP,0.)*CON(LMDCHHT,K)
C       FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
C    $               -MAX(QVKTMP,0.)*CON(LMDCHVT,K)
C    $               -MIN(QVKTMP,0.)*CON(LMDCHHT,K)
        END DO
        END DO
      END IF
C
C----------------------------------------------------------------------C
C
C **  GROUNDWATER, EVAP, RAINFALL (3TL)
C
      IF (ISGWIE.NE.0) THEN
        DO L=2,LA
         FQC(L,1)=FQC(L,1)-RIFTR(L)*CON1(L,1)
        END DO
      END IF
C
      IF(M.EQ.2) THEN
        IF(ISTOPT(2).EQ.0.OR.ISTOPT(2).EQ.3) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TEMO*DXYP(L)
          END DO
        END IF
        IF(ISTOPT(2).EQ.1.OR.ISTOPT(2).EQ.2) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)+RAINT(L)*TATMT(L)*DXYP(L)
          END DO
        END IF
      END IF
C
      IF(M.EQ.2) THEN
        IF(ISTOPT(2).EQ.0) THEN
          DO L=2,LA
           FQC(L,KC)=FQC(L,KC)-EVAPSW(L)*CON1(L,KC)
          END DO
        END IF
      END IF
C
      END IF
      GO TO 2000
C
C**********************************************************************C
C
 1000 CONTINUE
C
      DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,M)
       DO K=1,KC
       FQC(L,K)=FQC(L,K)
     $         +MAX(QSS(K,NS),0.)*CQS(K,NS,M)
     $         +MIN(QSS(K,NS),0.)*CON1(L,K)
     $         +MAX(QSERT(K,NQSTMP),0.)*CSERT(K,NCSTMP,M)
     $         +MIN(QSERT(K,NQSTMP),0.)*CON1(L,K)
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
       FQC(LU,K)=FQC(LU,K)
     $          -QCTLT(K,NCTL)*CON1(LU,K)
       FQC(LD,K)=FQC(LD,K)
     $          +RQWD*QCTLT(K,NCTL)*CON1(LU,K)
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
       FQC(LU,KU)=FQC(LU,KU)
     $          -(QWR(NWR)+QWRSERT(NQSTMP))*CON1(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)
     $          +QWR(NWR)*(SFLKILL*CON1(LU,KU)+CQWR(NWR,M))
     $          +QWRSERT(NQSTMP)*(SFLKILL*CON1(LU,KU)
     $                                   +CQWRSERT(NCSTMP,M))
C      FQC(LU,K)=FQC(LU,K)
C    $          -(QWR(K,NWR)+QSERT(K,NQSTMP))*CON(LU,K)
C      FQC(LD,K)=FQC(LD,K)
C    $          +QWR(K,NWR)*(CON(LU,K)+CQWR(K,NWR,M))
C    $          +QSERT(K,NQSTMP)*(CON(LU,K)+CSERT(K,NCSTMP,M))
      END DO
C
      IF (MDCHH.GE.1) THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF (MDCHTYP(NMD).EQ.1) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        END IF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
     $               +MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               +MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
     $               +MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               +MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
     $               -MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               -MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
     $               -MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               -MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        END DO
        END DO
      END IF
C
      GO TO 2000
C
C**********************************************************************C
C
 1500 CONTINUE
c
cBug, Ji, 12/30/99, This section of the code is for WQ model,  but
c the point source load (NQSIJ) is not included. Point source loads of
c WQ might(?) be added by wqpsl.inp. Therefore, (-q) point source
c should not work here. We have to use withdraw to get WQ variables out
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
       FQC(LU,K)=FQC(LU,K)
     $          -QCTLT(K,NCTL)*CON1(LU,K)
       FQC(LD,K)=FQC(LD,K)
     $          +RQWD*QCTLT(K,NCTL)*CON1(LU,K)
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
        FQC(LU,KU)=FQC(LU,KU)
     $          -(QWR(NWR)+QWRSERT(NQSTMP))*CON1(LU,KU)
       FQC(LD,KD)=FQC(LD,KD)
     $          +QWR(NWR)*(CON1(LU,KU)+CQWR(NWR,M))
     $          +QWRSERT(NQSTMP)*(CON1(LU,KU)+CQWRSERT(NCSTMP,M))
      END DO
C
      IF (MDCHH.GE.1) THEN
        DO K=1,KC
        DO NMD=1,MDCHH
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF (MDCHTYP(NMD).EQ.1) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.2) THEN
          QVKTMP=QCHANV(NMD)*DZC(K)
          QUKTMP=0.
        END IF
        IF (MDCHTYP(NMD).EQ.3) THEN
          QUKTMP=QCHANU(NMD)*DZC(K)
          QVKTMP=QCHANV(NMD)*DZC(K)
        END IF
        FQC(LMDCHHT,K)=FQC(LMDCHHT,K)
     $               +MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               +MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
     $               +MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               +MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHUT,K)=FQC(LMDCHUT,K)
     $               -MAX(QUKTMP,0.)*CON1(LMDCHUT,K)
     $               -MIN(QUKTMP,0.)*CON1(LMDCHHT,K)
        FQC(LMDCHVT,K)=FQC(LMDCHVT,K)
     $               -MAX(QVKTMP,0.)*CON1(LMDCHVT,K)
     $               -MIN(QVKTMP,0.)*CON1(LMDCHHT,K)
        END DO
        END DO
      END IF
C
C**********************************************************************C
C
 2000 CONTINUE
C
      DO K=1,KC
      RTMP=DZIC(K)
       DO L=1,LC
       FQC(L,K)=RTMP*FQC(L,K)
       END DO
      END DO
C
C**********************************************************************C
C
      RETURN
      END
