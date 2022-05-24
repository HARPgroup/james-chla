C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE CALQVS (ISTL)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C ** SUBROUTINE CALQVS UPDATES TIME VARIABLE VOLUME SOURCES
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
      DELT=DT2
      IF (ISTL.EQ.2) DELT=DT
      NTOX1=NTOX
      IF(NPCB.GT.0)NTOX1=NPCB      
C
C**********************************************************************C
C
C **  INITIALIZE NULL (0) FLOW SERIES
C
      QWRSERT(0)=0.
      DO K=1,KC
      QSERT(K,0)=0.
      QCTLT(K,0)=0.
      QCTLTO(K,0)=0.
      END DO
C
C **  INITIALIZE TOTAL FLOW SERIES
C
      DO L=1,LC
      QSUME(L)=0.
      END DO
C
      DO K=1,KC
       DO L=1,LC
       QSUM(L,K)=0.
       END DO
      END DO
C
C**********************************************************************C
C
C **  VOLUME SOURCE/SINK INTERPOLATION
C
       DO NS=1,NQSER
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
        ELSE
         TIME=DT*FLOAT(N-1)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
       END IF
 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
        ELSE
         TIME=(SECDLAST-DT)/TCQSER(NS)+TBEGIN*(TCON/TCQSER(NS))
        END IF        
       ENDIF
C
       M1=MQTLAST(NS)
  100  CONTINUE
       M2=M1+1
c     write(652,*) M2,NS,time,tqser(m2,ns)    ! check qser.inp
       IF (TIME.GT.TQSER(M2,NS)) THEN
         M1=M2
         GO TO 100
        ELSE
         MQTLAST(NS)=M1
       END IF
C
       TDIFF=TQSER(M2,NS)-TQSER(M1,NS)
       WTM1=(TQSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TQSER(M1,NS))/TDIFF
       DO K=1,KC
       QSERT(K,NS)=WTM1*QSER(M1,K,NS)+WTM2*QSER(M2,K,NS)
       END DO
!       if(NS.EQ.8) then
!       write(991,*)QSERT(KC,NS),QSER(M2,KC,NS)
!       endif
C
       END DO
c
      IF(N.EQ.1) THEN
      DO LL=1,NQSIJ
      L=LQS(LL)
      ITYP=LCT(L)
      IF (ITYP.EQ.0.OR.ITYP.GE.8) THEN
        WRITE(6,6111)LL,IQS(LL),JQS(LL)
        WRITE(8,6111)LL,IQS(LL),JQS(LL)
      END IF
      END DO
      END IF
C
      DO LL=1,NQSIJ
      NS=NQSERQ(LL)
      L=LQS(LL)
       DO K=1,KC
       QSUM(L,K)=QSUM(L,K)+RQSMUL(LL)*(QSS(K,LL)+QSERT(K,NS))
       END DO
      END DO
C
C**********************************************************************C
C
C **  CONTROL STRUCTURES AND TIDAL INLETS
C
      IF(ISCRAY.EQ.0) THEN
        T1TMP=SECNDS(0.0)
       ELSE
        T1TMP=SECOND( )
        CALL TIMEF(WT1TMP)
      END IF
C
      DO NCTL=1,NQCTL
      IF(NQCTYP(NCTL).LE.1)THEN
      NCTLT=NQCTLQ(NCTL)
      RQDW=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF (ID.EQ.0.AND.JD.EQ.0) THEN
        LD=LC
        HDW=0.
        RQDW=0.
       ELSE
        LD=LIJ(ID,JD)
        HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
      END IF
      DELH=HCTLUM(NCTLT)*HUP-HCTLDM(NCTLT)*HDW
      IF (DELH.LE.0.) THEN
        DO K=1,KC
        QCTLT(K,NCTL)=0.
        END DO
       ELSE
        IF(NQCTYP(NCTL).EQ.1)DELH=SQRT(DELH)
        M1=0
        M2=1
 500    M1=M1+1
        M2=M2+1
        IF (M2.GT.MQCTL(NCTLT)) THEN
           WRITE(6,6666)
           WRITE(6,6667)NCTL,NCTLT,IU,JU,ID,JD
           WRITE(6,6668)HUP,HP(LU),HDW,HP(LD)
           WRITE(8,6666)
           WRITE(8,6667)NCTL,NCTLT,IU,JU,ID,JD
           WRITE(8,6668)HUP,HP(LU),HDW,HP(LD)
           STOP
        END IF
        IF(DELH.GE.HDIFCTL(M1,NCTLT).AND.DELH.LE.HDIFCTL(M2,NCTLT))THEN
          TDIFF=HDIFCTL(M2,NCTLT)-HDIFCTL(M1,NCTLT)
          WTM1=(HDIFCTL(M2,NCTLT)-DELH)/TDIFF
          WTM2=(DELH-HDIFCTL(M1,NCTLT))/TDIFF
           DO K=1,KC
           QCTLT(K,NCTL)=WTM1*QCTL(M1,1,K,NCTLT)
     $                  +WTM2*QCTL(M2,1,K,NCTLT)
           END DO
         ELSE
          GO TO 500
        END IF
      END IF
      IF(NQCTYP(NCTL).EQ.1) THEN
      IF(ISTL.EQ.3) THEN
        DO K=1,KC
         QCTLST(K,NCTL)=QCTLT(K,NCTL)
         RTMP=QCTLTO(K,NCTL)
     $       +DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
         QCTLT(K,NCTL)=RTMP/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
         QCTLTO(K,NCTL)=QCTLT(K,NCTL)
         QCTLSTO(K,NCTL)=QCTLST(K,NCTL)
        END DO
       ELSE
        DO K=1,KC
         QCTLST(K,NCTL)=QCTLT(K,NCTL)
         RTMP=QCTLTO(K,NCTL)
     $       +DT*AQCTL(NCTLT)*QCTLST(K,NCTL)*QCTLST(K,NCTL)
         QCTLT(K,NCTL)=RTMP/(1.+DT*AQCTL(NCTLT)*QCTLTO(K,NCTL))
         QCTLT(K,NCTL)=0.5*(QCTLT(K,NCTL)+QCTLTO(K,NCTL))
        END DO
      END IF
      END IF
      DO K=1,KC
       QSUM(LU,K)=QSUM(LU,K)-RQCMUL(NCTL)*QCTLT(K,NCTL)
       QSUM(LD,K)=QSUM(LD,K)+RQCMUL(NCTL)*RQDW*QCTLT(K,NCTL)
       if(mod(N,NTSPTC/48).eq.0) then
       write(800,*)time,RQCMUL(NCTL)*QCTLT(K,NCTL),k,NCTL
       endif
      END DO
      END IF
      END DO
C
      DO NCTL=1,NQCTL
      IF(NQCTYP(NCTL).EQ.2)THEN
      NCTLT=NQCTLQ(NCTL)
      RQDW=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
      IF(HUP.LT.HDIFCTL(1,NCTLT)) THEN
        DO K=1,KC
         QCTLT(K,NCTL)=0.
        END DO
        GO TO 560
      END IF
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      LD=LIJ(ID,JD)
      HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
      HTMPD=HDIFCTD(1,NCTLT)+0.001
      HDW=MAX(HDW,HTMPD)
      MU1=0
      MU2=1
      MD1=0
      MD2=1
 555  MU1=MU1+1
      MU2=MU1+1
      IF (MU2.GT.MQCTL(NCTLT)) THEN
        WRITE(6,6676)
        WRITE(6,6677)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(6,6678)HUP,HP(LU),HDW,HP(LD)
        WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     $               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
        WRITE(8,6676)
        WRITE(8,6677)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(8,6678)HUP,HP(LU),HDW,HP(LD)
        WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     $               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
        STOP
      END IF  
      IF(HUP.GE.HDIFCTL(MU1,NCTLT).AND.HUP.LE.HDIFCTL(MU2,NCTLT))THEN
        TDIFFU=HDIFCTL(MU2,NCTLT)-HDIFCTL(MU1,NCTLT)
        WTM1U=(HDIFCTL(MU2,NCTLT)-HUP)/TDIFFU
        WTM2U=(HUP-HDIFCTL(MU1,NCTLT))/TDIFFU
       ELSE
        GO TO 555
      END IF
 556  MD1=MD1+1
      MD2=MD1+1
      IF (MD2.GT.MQCTL(NCTLT)) THEN
        WRITE(6,6686)
        WRITE(6,6687)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(6,6688)HUP,HP(LU),HDW,HP(LD)
        WRITE(6,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     $               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
        WRITE(8,6686)
        WRITE(8,6687)NCTL,NCTLT,IU,JU,ID,JD
        WRITE(8,6688)HUP,HP(LU),HDW,HP(LD)
        WRITE(8,6679)HDIFCTL(1,NCTLT),HDIFCTL(MQCTL(NCTLT),NCTLT),
     $               HDIFCTD(1,NCTLT),HDIFCTD(MQCTL(NCTLT),NCTLT)
        STOP
      END IF  
      IF(HDW.GE.HDIFCTD(MD1,NCTLT).AND.HDW.LE.HDIFCTD(MD2,NCTLT))THEN
        TDIFFD=HDIFCTD(MD2,NCTLT)-HDIFCTD(MD1,NCTLT)
        WTM1D=(HDIFCTD(MD2,NCTLT)-HDW)/TDIFFD
        WTM2D=(HDW-HDIFCTD(MD1,NCTLT))/TDIFFD
       ELSE
        GO TO 556
      END IF
      DO K=1,KC
      QCTLT(K,NCTL)=WTM1U*( WTM1D*QCTL(MU1,MD1,K,NCTLT)
     $                    +WTM2D*QCTL(MU1,MD2,K,NCTLT) )
     $             +WTM2U*( WTM1D*QCTL(MU2,MD1,K,NCTLT)
     $                    +WTM2D*QCTL(MU2,MD2,K,NCTLT) )
      END DO
  560 CONTINUE
      DO K=1,KC
       QSUM(LU,K)=QSUM(LU,K)-RQCMUL(NCTL)*QCTLT(K,NCTL)
       QSUM(LD,K)=QSUM(LD,K)+RQCMUL(NCTL)*RQDW*QCTLT(K,NCTL)
      END DO
      END IF
      END DO
C
      IF(ISCRAY.EQ.0) THEN
        TQCLT=TQCLT+SECNDS(T1TMP)
       ELSE
        T2TMP=SECOND( )
        CALL TIMEF(WT2TMP)
        TQCLT=TQCLT+T2TMP-T1TMP
        WTQCLT=WTQCLT+(WT2TMP-WT1TMP)*0.001
      END IF
C
C**********************************************************************C
C
C **  FLOW WITHDRAWAL AND RETURN
C
      NTMP=4+NSED+NSND+NTOX1
C
      DO NC=1,NTMP
        CQWRSERT(0,NC)=0.
      END DO
C
      DO NS=1,NQWRSR
C
       IF (ISTL.EQ.2) THEN
         TIME=DT*(FLOAT(N)-0.5)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
        ELSE
         TIME=DT*FLOAT(N-1)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
       END IF
 
! %change timestep
!      
        IF(NCSTEP.GT.0) THEN   ! J.S. 12/24/2010
        IF (ISTL.EQ.2) THEN
         TIME=(SECDLAST-0.5*DT)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
        ELSE
         TIME=(SECDLAST-DT)/TCQWRSR(NS)+TBEGIN*(TCON/TCQWRSR(NS))
        END IF        
       ENDIF
C
       M1=MQWRTLST(NS)
  200  CONTINUE
       M2=M1+1
       IF (TIME.GT.TQWRSER(M2,NS)) THEN
         M1=M2
         GO TO 200
        ELSE
         MQWRTLST(NS)=M1
       END IF      
C
       TDIFF=TQWRSER(M2,NS)-TQWRSER(M1,NS)
       WTM1=(TQWRSER(M2,NS)-TIME)/TDIFF
       WTM2=(TIME-TQWRSER(M1,NS))/TDIFF
       QWRSERT(NS)=WTM1*QWRSER(M1,NS)+WTM2*QWRSER(M2,NS)
       DO NC=1,NTMP
        CQWRSERT(NS,NC)=WTM1*CQWRSER(M1,NS,NC)+WTM2*CQWRSER(M2,NS,NC)
       END DO
C
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
      NS=NQWRSERQ(NWR)
      QSUM(LU,KU)=QSUM(LU,KU)-QWR(NWR)-QWRSERT(NS)
      QSUM(LD,KD)=QSUM(LD,KD)+QWR(NWR)+QWRSERT(NS)
      END DO
C
C**********************************************************************C
C
C **  CALL JPEFDC AND PLACE JET-PLUME VOLUMES SOURCES
C
!      IF(NQJPIJ.GT.0.AND.N.EQ.1) CALL JPEFDC
C
      IF(NQJPIJ.GT.0.AND.ISTL.EQ.3) THEN   
        IF(NUDJPC(1).EQ.NUDJP(1)) THEN
!          CALL JPEFDC
          NUDJPC(1)=1
         ELSE
          NUDJPC(1)=NUDJPC(1)+1
        END IF
      END IF
C
C **  PLACE JET-PLUME VOLUMES SOURCES
C
      IF(NQJPIJ.GT.0) THEN
      DO NJP=1,NQJPIJ
      IF(ICALJP(NJP).GT.0) THEN
       LJP=LIJ(IQJP(NJP),JQJP(NJP))
       KTMP=KEFFJP(NJP)
       QVJPTMP=QQCJP(NJP)
       DO K=1,KC
        QVJPTMP=QVJPTMP+QSERT(K,NQSERJP(NJP))
       END DO
       QSUM(LJP,KTMP)=QSUM(LJP,KTMP)+QVJPTMP
      END IF
      END DO
      END IF
C
C**********************************************************************C
C
C **  GROUND WATER INTERACTION, EVAPORATION AND RAINFALL
C
      IF (ISGWIE.EQ.0) THEN
        IF(EVAPCVT.LT.0.) THEN
          DO L=2,LA
           SVPW=(10.**((0.7859+0.03477*TEM(L,KC))/
     $              (1.+0.00412*TEM(L,KC))))
          EVAPT(L)=CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW-VPA(L))/PATMT(L) 
     & *abs(evapcvt)  ! added to adjust evaporation rate for TestWQ's Chloride modeling, Ji, 9/22/02      
           QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*(RAINT(L)-EVAPT(L))
          END DO
         ELSE
          DO L=2,LA
           QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*(RAINT(L)-EVAPT(L))
          END DO
        END IF
       ELSE
        DO L=2,LA
        QSUM(L,KC)=QSUM(L,KC)+DXYP(L)*RAINT(L)
        END DO
      END IF
C
C**********************************************************************C
C
C **  DETERMINE NET EXTERNAL VOLUME SOURCE/SINK
C
      DO K=1,KC
       DO L=1,LC
       QSUME(L)=QSUME(L)+QSUM(L,K)
       END DO
      END DO
C
C**********************************************************************C
C
C **  UPDATE ZERO DIMENSION VOLUME BALANCE
C 
C----------------------------------------------------------------------C
C
c     IF (ISDRY.GE.1.AND.ISTL.EQ.3) THEN
c       VOLADD=0.
c       DO L=2,LA
c       IF (SPB(L).NE.0) THEN
c         VOLADD=VOLADD+QSUME(L)
c       END IF
c       END DO
c       VOLADD=VOLADD*DT
c       VOLZERD=VOLZERD+VOLADD
c       VETZERD=VETZERD+VOLADD+DT*QETTMP
c     END IF
C
c5303 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
C                         
C**********************************************************************C
C
C **  WRITE DIAGNOSTIC FILE FOR VOLUME SOURCES,SINKS, ETC
C
      IF (ISDIQ.EQ.2.AND.ISTL.EQ.2) ITMPD=1
      IF (ISDIQ.EQ.1) ITMPD=1
C 
      NTT=4+NTOX1+NSED+NSND
C 
      IF(ITMPD.EQ.1) THEN
       IF(N.EQ.NTSPTC.OR.N.EQ.1) THEN
         OPEN(1,FILE='qdiag.out',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE='qdiag1.out',STATUS='UNKNOWN')
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE='qdiag1.out',STATUS='UNKNOWN')
        ELSE
         OPEN(1,FILE='qdiag.out',ACCESS='APPEND',STATUS='UNKNOWN')
       END IF
       WRITE(1,101)N
C
        DO LL=1,NQSIJ
        NQSTMP=NQSERQ(LL)
        NCSTMP=NCSERQ(LL,1)
        L=LQS(LL)
        I=IL(L)
        J=JL(L)
        WRITE(1,102)I,J
        WRITE(1,216)LL,L,(QSS(K,LL),K=1,KC)
        DO NT=1,NTT
          WRITE(1,217)LL,NT,(CQS(K,LL,NT),K=1,KC)
        END DO
        WRITE(1,104)
        WRITE(1,105)I,J
        WRITE(1,206)LL,L,(QSERT(K,NQSTMP),K=1,KC)
        DO NT=1,NTT
         NCSTMP=NCSERQ(LL,NT)
         WRITE(1,207)LL,NT,NCSTMP,(CSERT(K,NCSTMP,NT),K=1,KC)
         END DO
        WRITE(1,104)   
        END DO
C
        DO NCTL=1,NQCTL
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF(IU.EQ.0.AND.JU.EQ.0) THEN
          LU=0
          HUP=0.
         ELSE
          LU=LIJ(IU,JU)
          HUP=HP(LU)+BELV(LU)+HCTLUA(NCTLT)
        END IF
        IF(ID.EQ.0.AND.JD.EQ.0) THEN
          LD=0
          HDW=0.
         ELSE
          LD=LIJ(ID,JD)
          HDW=HP(LD)+BELV(LD)+HCTLDA(NCTLT)
        END IF
        WRITE(1,107)IU,JU,LU,NCTLT,HUP
         DO K=1,KC
         WRITE(1,108)K,QCTLT(K,NCTL)
         END DO
        WRITE(1,104)   
        WRITE(1,109)ID,JD,LD,NCTLT,HDW
         DO K=1,KC
         WRITE(1,108)K,QCTLT(K,NCTL)
         END DO
        WRITE(1,104) 
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
        WRITE(1,110)IU,JU
        WRITE(1,111)KU,QWR(NWR),CQWR(NWR,1),CQWR(NWR,2)
        WRITE(1,104)
        WRITE(1,112)ID,JD
        WRITE(1,111)KD,QWR(NWR),CQWR(NWR,1),CQWR(NWR,2)
        WRITE(1,104)   
        WRITE(1,113)IU,JU
        WRITE(1,114)KU,QWRSERT(NQSTMP),CQWRSERT(NQSTMP,1),
     $                 CQWRSERT(NQSTMP,2)
        WRITE(1,104)
        WRITE(1,115)ID,JD
        WRITE(1,114)KD,QWRSERT(NQSTMP),CQWRSERT(NQSTMP,1),
     $                 CQWRSERT(NQSTMP,2)
        WRITE(1,104)   
        END DO
C
        CLOSE(1)
      END IF
C
  101 FORMAT(1X,'SOURCE/SINK DIAGNOSTICS AT TIME STEP =',I8,//)
  102 FORMAT(3X,'CONST NQSIJ SOURCE/SINK FLOW AT I =',I5,' J =',I5,/)
  103 FORMAT(5X,'K =',I5,5X,'QSS(K) = ',E12.4,5X,'CQS(K,1) = ',E12.4,
     $       5X,'CQS(K,5) = ',E12.4)
  203 FORMAT(5X,'K =',I5,5X,'QSS(K) = ',E12.4,5X,'CQS(K, ) = ',
     $       5X, 12E12.4)
  104 FORMAT(/)
  105 FORMAT(3X,'TIME VAR NQSIJ SOURCE/SINK FLOW AT I =',I5,' J=',I5,/)
  106 FORMAT(5X,'K =',I5,5X,'QSERT(K) = ',E12.4,
     $       5X,'CSERT(K,1) = ',E12.4,5X,'CSERT(K,5) = ',E12.4)
  206 FORMAT(5X,'NQ,LQ     =',2I4,7X,'QSERT() = ',12E12.4)
  207 FORMAT(5X,'NQ,NT,NCQ =',3I4,3X,'CSERT() = ',12E12.4)
  216 FORMAT(5X,'NQ,LQ =',2I4,3X,'QSS() = ',12E12.4)
  217 FORMAT(5X,'NQ,NT =',2I4,3X,'CQS() = ',12E12.4)
  107 FORMAT(3X,'UPSTRM CONTROLED SINK FLOW AT I =',I5,' J =',I5,
     $          ' L =',I5,'  NQCTLT =',I5,'  HUP = ',E12.4/)
  108 FORMAT(5X,'K =',I5,5X,'QCTL(K) = ',2E12.4) 
  109 FORMAT(3X,'DWNSTRM CONTROLED SOURCE FLOW AT I =',I5,' J =',I5,
     $          ' L =',I5,'  NQCTLT =',I5,'  HDW = ',E12.4/)
  110 FORMAT(3X,'UPSTRM CONST WITHDRW SINK FLOW AT I =',I5,' J =',I5,/)
  111 FORMAT(5X,'K =',I5,5X,'QWR(K) = ',E12.4,
     $       5X,'CQWR(1) = ',E12.4,5X,'CQWR(2) = ',E12.4) 
  112 FORMAT(3X,'DWNSTRM CONST RETN SOURCE FLOW AT I =',I5,' J =',I5,/)
  113 FORMAT(3X,'UPSTRM VAR WITHDRW SINK FLOW AT I =',I5,' J =',I5,/)
  114 FORMAT(5X,'K =',I5,5X,'QSERT(K) = ',E12.4,
     $       5X,'CSERT(K,1) = ',E12.4,5X,'CSERT(K,2) = ',E12.4) 
  115 FORMAT(3X,'DWNSTRM VAR RETN SOURCE FLOW AT I =',I5,' J =',I5,/)
 6666 FORMAT(' SINGLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS ')
 6667 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6668 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6676 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, UP ')
 6677 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6678 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6679 FORMAT(' HUF,HUL,HDF,HDL = ',4(2X,E12.4))
 6686 FORMAT(' DOUBLE VAL CONTROL STRUCTURE TABLE OUT OF BOUNDS, DW ')
 6687 FORMAT(' NCTL,NCTLT,IU,JU,ID,JD = ',6I5)
 6688 FORMAT(' SELU,HU,SELD,HD = ',4(2X,E12.4))
 6111 FORMAT(' INVALID NQSIJ LOCATION, NQSIJ,I,J = ',3I5)
C  
C**********************************************************************C
C
      RETURN
      END
