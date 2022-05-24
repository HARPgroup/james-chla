C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
       SUBROUTINE DUMP
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE DUMP WRITES FULL FIELD DUMPS OF MODEL VARIABLES
C **  AT SPECIFIED TIME INTERVALS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
C**********************************************************************C
C
c      DIMENSION P(LCM),SEDBT(LCM),SNDBT(LCM),VOLBW2(LCM),
c     $          TOXB(LCM,NTXM),TOXPFTB(LCM,NTXM),U(LCM,KCM),          
c     $          V(LCM,KCM),W(LCM,KCM),SAL(LCM,KCM),TEM(LCM,KCM),
c     $          DYE(LCM,KCM),SEDT(LCM,KCM),
c     $          SNDT(LCM,KCM),TOX(LCM,KCM,NTXM),
c     $          TOXPFTW(LCM,KCM,NTXM),LNC(LCM),ISTRAN(7)
C
C      CHARACTER*13 FNDSEL,FNDUUU,FNDVVV,FNDWWW,FNDSAL,FNDTEM,FNDDYE,
C     $          FNDSDW,FNDSDB,FNDSNW,FNDSNB,FNDBDH,FNDTXW(NTXM),
C     $          FNDTXB(NTXM),FNDTPW(NTXM),FNDTPB(NTXM)
C          
C**********************************************************************C
C
      DIMENSION TXWMAX(NTXM),TXWMIN(NTXM),TXBMAX(NTXM),TXBMIN(NTXM)
      DIMENSION DMPVALL(LCM-2),DMPVAL(LCM-2,KCM)
      DIMENSION IDMPVALL(LCM-2),IDMPVAL(LCM-2,KCM)
      INTEGER*4 IB08VALL(LCM-2),IB08VAL(LCM-2,KCM)
      INTEGER*4 IB16VALL(LCM-2),IB16VAL(LCM-2,KCM)
C      DIMENSION IB08VALL(LCM-2),IB08VAL(LCM-2,KCM)
C      DIMENSION IB16VALL(LCM-2),IB16VAL(LCM-2,KCM)
      CHARACTER*1 CZTT(0:9)
      CHARACTER*1 CCHTMF,CCHTMS
      CHARACTER*2 CNTTOX(NTXM)
C
C**********************************************************************C
C
C **  INFORMATION FOR TESTING AS STAND ALONE PROGRAM
C
C      LC=514
C      LA=LC-1
C      KC=10
C      NTOX=20
C
C      DO L=1,LC
C       P(L)=FLOAT(L-1)
C       SEDBT(L)=FLOAT(L-1)
C       SNDBT(L)=FLOAT(L-1)
C      END DO
C      DO K=1,KC
C       DO L=1,LC
C        U(L,K)=FLOAT(L-1)
C        V(L,K)=FLOAT(L-1)
C        W(L,K)=FLOAT(L-1)
C        SAL(L,K)=FLOAT(L-1)
C        TEM(L,K)=FLOAT(L-1)
C        DYE(L,K)=FLOAT(L-1)
C        SEDT(L,K)=FLOAT(L-1)
C        SNDT(L,K)=FLOAT(L-1)
C       END DO
C      END DO
C
C      DO NT=1,NTOX
C      DO L=1,LC
C       TOXB(L,KB,NT)=FLOAT(L-1)
C       TOXPFTB(L,NT)=FLOAT(L-1)/512.
C      END DO
C      DO K=1,KC
C       DO L=1,LC
C        TOX(L,K,NT)=FLOAT(L-1)
C        TOXPFTW(L,K,NT)=FLOAT(L-1)/512.
C       END DO
C      END DO
C      END DO
C
C      JSDUMP=1
C      ISDUMP=2 
C      ISADMP=0
C      NSDUMP=1 
C      TSDUMP=1. 
C      TEDUMP=1. 
C      ISDMPP=1 
C      ISDMPU=1 
C      ISDMPW=1 
C      ISDMPT=1 
C      IADJDMP=-32768
C      NTOX=20
C
C      DO N=1,7
C       ISTRAN(N)=1
C      END DO
C
C      GI=1./9.8
C
C      DT=1.
C      TCON=1
C      TBEGIN=0.
C      N=1
C
C**********************************************************************C
C
      IF (JSDUMP.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C
      CZTT(0)='0'
      CZTT(1)='1'
      CZTT(2)='2'
      CZTT(3)='3'
      CZTT(4)='4'
      CZTT(5)='5'
      CZTT(6)='6'
      CZTT(7)='7'
      CZTT(8)='8'
      CZTT(9)='9'
C
      DO MLTM=1,NTOX
      MSDIG=MOD(MLTM,10)
      MTMP=MLTM-MSDIG
      MFDIG=MTMP/10
      CCHTMF=CZTT(MFDIG)
      CCHTMS=CZTT(MSDIG)
      CNTTOX(MLTM)= CCHTMF // CCHTMS
      END DO
C
C  ISDUMP=1, ASCII INTERGER OUTPUT
C
      IF(ISDUMP.EQ.1) THEN
      FNDSEL='seldmpi.asc'
      FNDUUU='uuudmpi.asc'
      FNDVVV='vvvdmpi.asc'
      FNDWWW='wwwdmpi.asc'
      FNDSAL='saldmpi.asc'
      FNDTEM='temdmpi.asc'
      FNDDYE='dyedmpi.asc'
      FNDSDW='sdwdmpi.asc'
      FNDSDB='sdbdmpi.asc'
      FNDSNW='snwdmpi.asc'
      FNDSNB='snbdmpi.asc'
      FNDBDH='bdhdmpi.asc'
      DO NT=1,NTOX
       FNDTXW(NT)='tw'// CNTTOX(NT) // 'dpi.asc'
       FNDTXB(NT)='tb'// CNTTOX(NT) // 'dpi.asc'
       FNDTPW(NT)='fw'// CNTTOX(NT) // 'dpi.asc'
       FNDTPB(NT)='fb'// CNTTOX(NT) // 'dpi.asc'
      END DO
      END IF
C
C  ISDUMP=2, 16/8 BIT BINARY INTERGER OUTPUT
C
      IF(ISDUMP.EQ.2) THEN
      FNDSEL='seldmpi.bin'
      FNDUUU='uuudmpi.bin'
      FNDVVV='vvvdmpi.bin'
      FNDWWW='wwwdmpi.bin'
      FNDSAL='saldmpi.bin'
      FNDTEM='temdmpi.bin'
      FNDDYE='dyedmpi.bin'
      FNDSDW='sdwdmpi.bin'
      FNDSDB='sdbdmpi.bin'
      FNDSNW='snwdmpi.bin'
      FNDSNB='snbdmpi.bin'
      FNDBDH='bdhdmpi.bin'
      DO NT=1,NTOX
       FNDTXW(NT)='tw'// CNTTOX(NT) // 'dpi.bin'
       FNDTXB(NT)='tb'// CNTTOX(NT) // 'dpi.bin'
       FNDTPW(NT)='fw'// CNTTOX(NT) // 'dpi.bin'
       FNDTPB(NT)='fb'// CNTTOX(NT) // 'dpi.bin'
      END DO
      END IF
C
C  ISDUMP=3, ASCII FLOATING POINT OUTPUT
C
      IF(ISDUMP.EQ.3) THEN
      FNDSEL='seldmpf.asc'
      FNDUUU='uuudmpf.asc'
      FNDVVV='vvvdmpf.asc'
      FNDWWW='wwwdmpf.asc'
      FNDSAL='saldmpf.asc'
      FNDTEM='temdmpf.asc'
      FNDDYE='dyedmpf.asc'
      FNDSDW='sdwdmpf.asc'
      FNDSDB='sdbdmpf.asc'
      FNDSNW='snwdmpf.asc'
      FNDSNB='snbdmpf.asc'
      FNDBDH='bdhdmpf.asc'
      DO NT=1,NTOX
       FNDTXW(NT)='tw'// CNTTOX(NT) // 'dpf.asc'
       FNDTXB(NT)='tb'// CNTTOX(NT) // 'dpf.asc'
       FNDTPW(NT)='fw'// CNTTOX(NT) // 'dpf.asc'
       FNDTPB(NT)='fb'// CNTTOX(NT) // 'dpf.asc'
      END DO
      END IF
C
C  ISDUMP=4, 32/64 BIT BINARY FLOATING POINT OUTPUT
C
      IF(ISDUMP.EQ.4) THEN
      FNDSEL='seldmpf.bin'
      FNDUUU='uuudmpf.bin'
      FNDVVV='vvvdmpf.bin'
      FNDWWW='wwwdmpf.bin'
      FNDSAL='saldmpf.bin'
      FNDTEM='temdmpf.bin'
      FNDDYE='dyedmpf.bin'
      FNDSDW='sdwdmpf.bin'
      FNDSDB='sdbdmpf.bin'
      FNDSNW='snwdmpf.bin'
      FNDSNB='snbdmpf.bin'
      FNDBDH='bdhdmpf.bin'
      DO NT=1,NTOX
       FNDTXW(NT)='tw'// CNTTOX(NT) // 'dpf.bin'
       FNDTXB(NT)='tb'// CNTTOX(NT) // 'dpf.bin'
       FNDTPW(NT)='fw'// CNTTOX(NT) // 'dpf.bin'
       FNDTPB(NT)='fb'// CNTTOX(NT) // 'dpf.bin'
      END DO
      END IF
C
      IF(ISADMP.EQ.0) THEN
C
        OPEN(1,FILE=FNDSEL)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDUUU)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDVVV)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDWWW)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDSAL)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDTEM)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDDYE)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDSDW)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDSDB)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDSNW)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDSNB)
        CLOSE(1,STATUS='DELETE')
C
        OPEN(1,FILE=FNDBDH)
        CLOSE(1,STATUS='DELETE')
C
        DO NT=1,NTOX
         OPEN(1,FILE=FNDTXW(NT))
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE=FNDTXB(NT))
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE=FNDTPW(NT))
         CLOSE(1,STATUS='DELETE')
         OPEN(1,FILE=FNDTPB(NT))
         CLOSE(1,STATUS='DELETE')
        END DO
C
      END IF
C
      JSDUMP=0
C
C**********************************************************************C
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      DO K=1,KC
       DO L=1,LC-2
        DMPVAL(L,K)=0.
        IDMPVAL(L,K)=0
        IB08VAL(L,K)=0
        IB16VAL(L,K)=0
       END DO
      END DO
C
      DO L=1,LC-2
       DMPVALL(L)=0.
       IDMPVALL(L)=0
       IB08VALL(L)=0
       IB16VALL(L)=0
      END DO
C
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/86400.
      IF(NCSTEP.GT.0) TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
      R1=1.
      R0=0.
C
C**********************************************************************C
C
C **  IF (ISDUMP EQUAL 1 OR 2, SCALE VARIABLES AND WRITE INTEGER 
C **  DUMP FILES
C
      IF (ISDUMP.LE.2) THEN
C
C **  SCALE VARIABLES
C
      SELMAX=-1.E12
      SELMIN=1.E12
      UUUMAX=-1.E12
      UUUMIN=1.E12
      VVVMAX=-1.E12
      VVVMIN=1.E12
      WWWMAX=-1.E12
      WWWMIN=1.E12
      SALMAX=-1.E12
      SALMIN=1.E12
      TEMMAX=-1.E12
      TEMMIN=1.E12
      DYEMAX=-1.E12
      DYEMIN=1.E12
      SDWMAX=-1.E12
      SDWMIN=1.E12
      SDBMAX=-1.E12
      SDBMIN=1.E12
      SNWMAX=-1.E12
      SNWMIN=1.E12
      SNWMAX=-1.E12
      SNWMIN=1.E12
      SNBMAX=-1.E12
      SNBMIN=1.E12
      BDHMAX=-1.E12
      BDHMIN=1.E12
      DO NT=1,NTOX
       TXWMAX(NT)=-1.E12
       TXWMIN(NT)=1.E12
       TXBMAX(NT)=-1.E12
       TXBMIN(NT)=1.E12
      END DO       
C
      IF(ISDMPP.GE.1) THEN
        DO L=2,LA
         SELMAX=MAX(SELMAX,P(L))
         SELMIN=MIN(SELMIN,P(L))
        END DO
      END IF
      SELMAX=GI*SELMAX
      SELMIN=GI*SELMIN
C
      IF(ISDMPU.GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         UTMP=0.5*(U(L,K)+U(L+1,K))
         VTMP=0.5*(V(L,K)+V(LNC(L),K))
         UUUMAX=MAX(UUUMAX,UTMP)
         UUUMIN=MIN(UUUMIN,UTMP)
         VVVMAX=MAX(VVVMAX,UTMP)
         VVVMIN=MIN(VVVMIN,VTMP)
        END DO
        END DO
      END IF
C
      IF(ISDMPW.GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         WTMP=0.5*(W(L,K)+W(L,K-1))
         WWWMAX=MAX(WWWMAX,WTMP)
         WWWMIN=MIN(WWWMIN,WTMP)
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         SALMAX=MAX(SALMAX,SAL(L,K))
         SALMIN=MIN(SALMIN,SAL(L,K))
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         TEMMAX=MAX(TEMMAX,TEM(L,K))
         TEMMIN=MIN(TEMMIN,TEM(L,K))
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         DYEMAX=MAX(DYEMAX,DYE(L,K))
         DYEMIN=MIN(DYEMIN,DYE(L,K))
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         SDWMAX=MAX(SDWMAX,SEDT(L,K))
         SDWMIN=MIN(SDWMIN,SEDT(L,K))
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        DO K=1,KC
        DO L=2,LA
         SNWMAX=MAX(SNWMAX,SNDT(L,K))
         SNWMIN=MIN(SNWMIN,SNDT(L,K))
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
        DO NT=1,NTOX
        DO K=1,KC
        DO L=2,LA
         TXWMAX(NT)=MAX(TXWMAX(NT),TOX(L,K,NT))
         TXWMIN(NT)=MIN(TXWMIN(NT),TOX(L,K,NT))
        END DO
        END DO
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        DO L=2,LA
         SDBMAX=MAX(SDBMAX,SEDBT(L,KBT(L)))
         SDBMIN=MIN(SDBMIN,SEDBT(L,KBT(L)))
        END DO
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        DO L=2,LA
         SNBMAX=MAX(SNBMAX,SNDBT(L,KBT(L)))
         SNBMIN=MIN(SNBMIN,SNDBT(L,KBT(L)))
        END DO
      END IF
C
      IF(ISDMPT.GE.1) THEN
      IF(ISTRAN(7).GE.1.OR.ISTRAN(6).GE.1) THEN
        DO L=2,LA
         BDHMAX=MAX(BDHMAX,VOLBW2(L,KBT(L)))
         BDHMIN=MIN(BDHMIN,VOLBW2(L,KBT(L)))
        END DO
      END IF
      END IF
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
        DO NT=1,NTOX
        DO L=2,LA
         TXBMAX(NT)=MAX(TXBMAX(NT),TOXB(L,KBT(L),NT))
         TXBMIN(NT)=MIN(TXBMIN(NT),TOXB(L,KBT(L),NT))
        END DO
        END DO
      END IF
C
C **  WRITE ARRAYS
C
      IF(ISDUMP.EQ.1) RSCALE=65535.
      IF(ISDUMP.EQ.2) RSCALE=65535.
C
C **  WATER SURFACE ELEVATION
C
      IF(ISDMPP.GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSEL,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSEL,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SELMAX-SELMIN)
        DO L=2,LA
         DMPVALL(L-1)=SCALE*(GI*P(L)-SELMIN)
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SELMAX,SELMIN
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,SELMAX,SELMIN
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END IF
C
C **  U VELOCITY COMPONENT
C
      IF(ISDMPU.GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDUUU,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDUUU,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(UUUMAX-UUUMIN)
        DO K=1,KC
        DO L=2,LA
         UUUTMP=0.5*(U(L,K)+U(L+1,K))
         DMPVAL(L-1,K)=SCALE*(UUUTMP-UUUMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,UUUMAX,UUUMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,UUUMAX,UUUMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  V VELOCITY COMPONENT
C
      IF(ISDMPU.GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDVVV,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDVVV,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(VVVMAX-VVVMIN)
        DO K=1,KC
        DO L=2,LA
         VVVTMP=0.5*(V(L,K)+V(LNC(L),K))
         DMPVAL(L-1,K)=SCALE*(VVVTMP-VVVMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,VVVMAX,VVVMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,VVVMAX,VVVMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  W VELOCITY COMPONENT
C
      IF(ISDMPW.GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDWWW,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDWWW,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(WWWMAX-WWWMIN)
        DO K=1,KC
        DO L=2,LA
         WWWTMP=0.5*(W(L,K)+W(L,K-1))
         DMPVAL(L-1,K)=SCALE*(WWWTMP-WWWMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,WWWMAX,WWWMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,WWWMAX,WWWMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  SALINITY
C
      IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSAL,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSAL,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SALMAX-SALMIN)
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(SAL(L,K)-SALMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SALMAX,SALMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,SALMAX,SALMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  TEMPATURE
C
      IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTEM,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDTEM,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(TEMMAX-TEMMIN)
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(TEM(L,K)-TEMMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,TEMMAX,TEMMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,TEMMAX,TEMMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  DYE
C
      IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDDYE,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDDYE,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(DYEMAX-DYEMIN)
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(DYE(L,K)-DYEMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,DYEMAX,DYEMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,DYEMAX,DYEMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL COHESIVE SEDIMENT WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSDW,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSDW,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SDWMAX-SDWMIN)
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(SEDT(L,K)-SDWMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SDWMAX,SDWMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,SDWMAX,SDWMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL NONCOHESIVE SEDIMENT IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSNW,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSNW,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SNWMAX-SNWMIN)
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(SNDT(L,K)-SNWMIN)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SNWMAX,SNWMIN
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,SNWMAX,SNWMIN
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL TOXIC CONTAMINANTS IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTXW(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDTXW(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(TXWMAX(NT)-TXWMIN(NT))
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*(TOX(L,K,NT)-TXWMIN(NT))
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,TXWMAX(NT),TXWMIN(NT)
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,TXWMAX(NT),TXWMIN(NT)
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOTAL TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTPW(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDTPW(NT),ACCESS='APPEND',FORM='UNFORMATTED')
C        SCALE=100.
        SCALE=RSCALE
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SCALE*TOXPFTW(L,K,NT)
         IDMPVAL(L-1,K)=NINT(DMPVAL(L-1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,R1,R0
C          WRITE(1,102)IDMPVAL
          WRITE(1,101)IDMPVAL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO K=1,KC
          DO L=2,LA
C           IB08VAL(L-1,K)=IDMPVAL(L-1,K)
           IB16VAL(L-1,K)=IDMPVAL(L-1,K)+IADJDMP
          END DO
          END DO
          WRITE(1)TIME,R1,R0
C          WRITE(1)IB08VAL
          WRITE(1)IB16VAL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOTAL COHESIVE SEDIMENT IN BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSDB,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSDB,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SDBMAX-SDBMIN)
        DO L=2,LA
         DMPVALL(L-1)=SCALE*(SEDBT(L,KBT(L))-SDBMIN)
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SDBMAX,SDBMIN
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,SDBMAX,SDBMIN
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL NONCOHESIVE SEDIMENT IN BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDSNB,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDSNB,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(SNBMAX-SNBMIN)
        DO L=2,LA
         DMPVALL(L-1)=SCALE*(SNDBT(L,KBT(L))-SNBMIN)
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,SNBMAX,SNBMIN
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,SNBMAX,SNBMIN
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END IF
C
C **  THICKNESS OF SEDIMENT BED
C
      IF(ISDMPT.GE.1) THEN
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDBDH,ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDBDH,ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(BDHMAX-BDHMIN)
        DO L=2,LA
         DMPVALL(L-1)=SCALE*(VOLBW2(L,KBT(L))-BDHMIN)
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,BDHMAX,BDHMIN
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,BDHMAX,BDHMIN
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END IF
      END IF
C
C **  TOTAL TOXIC CONTAMINANTS IN SEDIMENT BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTXB(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDTXB(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        SCALE=RSCALE/(TXBMAX(NT)-TXBMIN(NT))
        DO L=2,LA
         DMPVALL(L-1)=SCALE*(TOXB(L,KBT(L),NT)-TXBMIN(NT))
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,TXBMAX(NT),TXBMIN(NT)
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,TXBMAX(NT),TXBMIN(NT)
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOXIC PARTICULATE FRACTION IS SEDIMENT BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.1) OPEN(1,FILE=FNDTPB(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.2) 
     $     OPEN(1,FILE=FNDTPB(NT),ACCESS='APPEND',FORM='UNFORMATTED')
C        SCALE=100.
        SCALE=RSCALE
        DO L=2,LA
         DMPVALL(L-1)=SCALE*TOXPFTB(L,KB,NT)
         IDMPVALL(L-1)=NINT(DMPVALL(L-1))
        END DO
        IF(ISDUMP.EQ.1) THEN
          WRITE(1,*)TIME,R1,R0
C          WRITE(1,102)IDMPVALL
          WRITE(1,101)IDMPVALL
        END IF
        IF(ISDUMP.EQ.2) THEN
          DO L=2,LA
C           IB08VALL(L-1)=IDMPVALL(L-1)
           IB16VALL(L-1)=IDMPVALL(L-1)+IADJDMP
          END DO
          WRITE(1)TIME,R1,R0
C          WRITE(1)IB08VALL
          WRITE(1)IB16VALL
        END IF
        CLOSE(1)
      END DO
      END IF
C
      END IF
C
C**********************************************************************C
C
C **  IF (ISDUMP EQUAL 3 OR 4, WRITE FLOATING POINT
C **  DUMP FILES
C
      IF (ISDUMP.GE.3) THEN
C
C **  WATER SURFACE ELEVATION
C
      IF(ISDMPP.GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSEL,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSEL,ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=GI*P(L)
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVALL
          WRITE(1,111)(DMPVALL(L),L=1,LA-1)
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END IF
C
C **  U VELOCITY COMPONENT
C
      IF(ISDMPU.GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDUUU,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDUUU,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=0.5*(U(L,K)+U(L+1,K))
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  V VELOCITY COMPONENT
C
      IF(ISDMPU.GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDVVV,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDVVV,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=0.5*(V(L,K)+V(LNC(L),K))
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  W VELOCITY COMPONENT
C
      IF(ISDMPW.GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDWWW,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDWWW,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=0.5*(W(L,K)+W(L,K-1))
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  SALINITY
C
      IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSAL,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSAL,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SAL(L,K)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  TEMPERATURE
C
      IF(ISDMPT.GE.1.AND.ISTRAN(2).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDTEM,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDTEM,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=TEM(L,K)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  DYE
C
      IF(ISDMPT.GE.1.AND.ISTRAN(3).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDDYE,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDDYE,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=DYE(L,K)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL COHESIVE SEDIMENT IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSDW,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSDW,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SEDT(L,K)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL NONCOHESIVE SEDIMENT IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSNW,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSNW,ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=SNDT(L,K)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL TOXIC CONTAMINANTS IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDTXW(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDTXW(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=TOX(L,K,NT)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOTAL TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDTPW(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDTPW(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        DO K=1,KC
        DO L=2,LA
         DMPVAL(L-1,K)=TOXPFTW(L,K,NT)
        END DO
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVAL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOTAL COHESIVE SEDIMENT IN BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(6).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSDB,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSDB,ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=SEDBT(L,KBT(L))
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END IF
C
C **  TOTAL NONCOHESIVE SEDIMENT IN BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDSNB,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDSNB,ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=SNDBT(L,KBT(L))
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END IF
C
C **  THICKNESS OF SEDIMENT BED
C
      IF(ISDMPT.GE.1) THEN
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDBDH,ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDBDH,ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=VOLBW2(L,KB)
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END IF
      END IF
C
C **  TOTAL TOXIC CONTAMINANTS IN SEDIMENT BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDTXB(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDTXB(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=TOXB(L,KB,NT)
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END DO
      END IF
C
C **  TOXIC PARTICULATE FRACTION IS SEDIMENT BED
C
      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
        IF(ISDUMP.EQ.3) OPEN(1,FILE=FNDTPB(NT),ACCESS='APPEND')
        IF(ISDUMP.EQ.4) 
     $     OPEN(1,FILE=FNDTPB(NT),ACCESS='APPEND',FORM='UNFORMATTED')
        DO L=2,LA
         DMPVALL(L-1)=SCALE*TOXPFTB(L,KB,NT)
        END DO
        IF(ISDUMP.EQ.3) THEN
          WRITE(1,*)TIME
C          WRITE(1,111)DMPVAL
          DO L=1,LA-1
            WRITE(1,111)(DMPVAL(L,K), K=1,KC)
          END DO
        END IF
        IF(ISDUMP.EQ.4) THEN
          WRITE(1)TIME
          WRITE(1)DMPVALL
        END IF
        CLOSE(1)
      END DO
      END IF
C
      END IF
C
C**********************************************************************C
C
C **  CHECK BY READING BINARY FILES
C
C      IF(ISDUMP.EQ.2) THEN
C
C      OPEN(2,FILE='dump.chk')
C      CLOSE(2,STATUS='DELETE')
C      OPEN(2,FILE='dump.chk')
C
C      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
C        OPEN(1,FILE=FNDTPB(1),FORM='UNFORMATTED')
C        READ(1)TIME,RMAX,RMIN
C        READ(1)IB08VALL
C        CLOSE(1)
C        DO L=1,LC-2
C         DMPVALL(L)=FLOAT(IB08VALL(L))/100.
C        END DO
C        WRITE(2,201)
C        WRITE(2,205)TIME,RMAX,RMIN
C        WRITE(2,205)DMPVALL
C      END IF
C
C      IF(ISDMPT.GE.1.AND.ISTRAN(5).GE.1) THEN
C        OPEN(1,FILE=FNDTPW(1),FORM='UNFORMATTED')
C        READ(1)TIME,RMAX,RMIN
C        READ(1)IB08VAL
C        CLOSE(1)
C        DO K=1,KC
C        DO L=1,LC-2
C         DMPVAL(L,K)=FLOAT(IB08VAL(L,K))/100.
C        END DO
C        END DO
C        WRITE(2,202)
C        WRITE(2,205)TIME,RMAX,RMIN
C        WRITE(2,205)DMPVAL
C      END IF
C
C      IF(ISDMPP.GE.1) THEN
C        OPEN(1,FILE=FNDSEL,FORM='UNFORMATTED')
C        READ(1)TIME,SELMAX,SELMIN
C        SELMAX=SELMAX/GI
C        SELMIN=SELMIN/GI
C        READ(1)IB16VALL
C        CLOSE(1)
C        DO L=1,LC-2
C         DMPVALL(L)=FLOAT(IB16VALL(L))
C         DMPVALL(L)=DMPVALL(L)-FLOAT(IADJDMP)
C        END DO
C        TMPVAL=(SELMAX-SELMIN)/RSCALE
C        DO L=1,LC-2
C         DMPVALL(L)=TMPVAL*DMPVALL(L)+SELMIN
C        END DO
C        WRITE(2,203)
C        WRITE(2,205)TIME,SELMAX,SELMIN
C        WRITE(2,205)DMPVALL
C      END IF
C
C      IF(ISDMPT.GE.1.AND.ISTRAN(1).GE.1) THEN
C        OPEN(1,FILE=FNDSAL,FORM='UNFORMATTED')
C        READ(1)TIME,SALMAX,SALMIN
C        READ(1)IB16VAL
C        CLOSE(1)
C        DO K=1,KC
C        DO L=1,LC-2
C         DMPVAL(L,K)=FLOAT(IB16VAL(L,K))
C         DMPVAL(L,K)=DMPVAL(L,K)-FLOAT(IADJDMP)
C        END DO
C        END DO
C        TMPVAL=(SALMAX-SALMIN)/RSCALE
C        DO K=1,KC
C        DO L=1,LC-2
C         DMPVAL(L,K)=TMPVAL*DMPVAL(L,K)+SALMIN
C        END DO
C        END DO
C        WRITE(2,204)
C        WRITE(2,205)TIME,SALMAX,SALMIN
C        WRITE(2,205)DMPVAL
C      END IF
C
C      CLOSE(2)
C
C      END IF
C
C**********************************************************************C
C
  100 FORMAT(A80)
  101 FORMAT(8I6)
  102 FORMAT(8I4)
c  111 FORMAT(8E14.6)
  111 FORMAT(10E12.4)
  201 FORMAT(//,' CHECK 2D  8 BIT VARIABLE',/)
  202 FORMAT(//,' CHECK 3D  8 BIT VARIABLE',/)
  203 FORMAT(//,' CHECK 2D 16 BIT VARIABLE',/)
  204 FORMAT(//,' CHECK 3D 16 BIT VARIABLE',/)
  205 FORMAT(8F8.2)
C
C**********************************************************************C
C
      RETURN
      END
