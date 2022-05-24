C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE TMSR_asc
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE TMSR WRITES TIME SERIES FILES FOR SURFACE ELEVATON,
C **  VELOCITY, CONCENTRATION, AND VOLUME SOURCES AT SPECIFIED
C **  (I,J) POINTS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      INCLUDE 'wq.par' 
      INCLUDE 'wqcom.cmn'      
C
C**********************************************************************C
C
C     DIMENSION ATMP(KCM),BTMP(KCM)
      CHARACTER*80 TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,TITLE6,TITLE7,
     $         TITLE11,TITLE12,TITLE13,TITLE14,TITLE15,TITLE16,TITLE17
      CHARACTER*10 CTUNIT
      CHARACTER*1 CZTT(0:9)
      CHARACTER*1 CCHTMF,CCHTMS
      CHARACTER*2 CNTOX(25)
C 
C**********************************************************************C
C
      NTOX1=NTOX

      IF(NTOX.EQ.0)NTOX1=NPCB
      
      IF (JSTMSR.NE.1) GO TO 300
C
C----------------------------------------------------------------------C
C

      
      IF(MLTMSR.GT.MLTMSRM) THEN
        WRITE (6,600)
        STOP
      END IF
      IF(MLTMSR.GT.99) THEN
        WRITE (6,601)
        STOP
      END IF
C
  600 FORMAT(' NUMBER OF TIME SER LOC, MLTMSR, EXCEEDS DIM, MLTMSRM')
  601 FORMAT(' NUMBER OF TIME SERIES LOCATIONS EXCEED 99')
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
      DO MLTM=1,MLTMSR
      MSDIG=MOD(MLTM,10)
      MTMP=MLTM-MSDIG
      MFDIG=MTMP/10
      CCHTMF=CZTT(MFDIG)
      CCHTMS=CZTT(MSDIG)
      CNTMSR(MLTM)= CCHTMF // CCHTMS
      END DO
C
      IF(TCTMSR.EQ.1.) CTUNIT='SECONDS'
      IF(TCTMSR.EQ.60.) CTUNIT='MINUTES'
      IF(TCTMSR.EQ.3600.) CTUNIT='HOURS'
      IF(TCTMSR.EQ.86400.) CTUNIT='DAYS'
C
       CNTOX( 1)= '01'
       CNTOX( 2)= '02'
       CNTOX( 3)= '03'
       CNTOX( 4)= '04'
       CNTOX( 5)= '05'
       CNTOX( 6)= '06'
       CNTOX( 7)= '07'
       CNTOX( 8)= '08'
       CNTOX( 9)= '09'
       CNTOX(10)= '10'
       CNTOX(11)= '11'
       CNTOX(12)= '12'
       CNTOX(13)= '13'
       CNTOX(14)= '14'
       CNTOX(15)= '15'
       CNTOX(16)= '16'
       CNTOX(17)= '17'
       CNTOX(18)= '18'
       CNTOX(19)= '19'
       CNTOX(20)= '20'
       CNTOX(21)= '21'
       CNTOX(22)= '22'
       CNTOX(23)= '23'
       CNTOX(24)= '24'
       CNTOX(25)= '25'
C
C **  WRITE HEADINGS
C
      TITLE1=' SALINITY (PSU) TIME SERIES, K=1,KC'
      TITLE2=' TEMPERATURE (DEG C) TIME SERIES, K=1,KC'
      TITLE3=' DYE CONC (KG/M**3) TIME SERIES, K=1,KC'
      TITLE4=' SED CONC (MG/LITER) TIME SERIES, K=1,KC'
      TITLE5=' TOXIC CONC (UG/LITER) TIME SERIES K=1,KC'
      TITLE6=' VISCOSITY (CM**2/S) TIME SERIES, K=1,KS'
      TITLE7=' DIFFUSIVITY (CM**2/S) TIME SERIES, K=1,KS'
      TITLE11=' SURFACE ELEVATION & DEPTH (METERS) TIME SERIES'
      TITLE12=' EXT MODE E,N VEL (CM/S), BOT TAU (CM/S)**2 TIME SERIES'
      TITLE13=' EXT MODE U,V TRANSPORT (M**3/S) TIME SERIES'
      TITLE14=' INT MODE EAST VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE15=' INT MODE NORTH VEL (CM/S) TIME SERIES, K=1,KC'
      TITLE16=' EXT MODE VOLUME S/S (M**3/S) TIME SERIES'
      TITLE17=' INT MODE VOL S/S (M**3/S) TIME SERIES, K=1,KC'
C
      IF (ISTMSR.EQ.2) THEN
      DO MLTM=1,MLTMSR
        IF (MTMSRC(MLTM).EQ.1) THEN
          IF(ISTRAN(1).GE.1) THEN
            FNSAL(MLTM)='salts' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(2).GE.1) THEN
            FNTEM(MLTM)='temts' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(3).GE.1) THEN
            FNDYE(MLTM)='dyets' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(4).GE.1) THEN
            FNSFL(MLTM)='sflts' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(6).GE.1) THEN
            FNSED(MLTM)='sedts' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(7).GE.1) THEN
            FNSND(MLTM)='sndts' // CNTMSR(MLTM) // '.out'
          END IF
          IF (ISTRAN(5).GE.1) THEN
            DO NT=1,NTOX1
            FNTOX(MLTM,NT)='t' // CNTOX(NT) // 'ts' //
     $                           CNTMSR(MLTM) // '.out'
            END DO
          END IF
        END IF
        IF (MTMSRA(MLTM).EQ.1) THEN
          FNAVV(MLTM)='avvts' // CNTMSR(MLTM) // '.out'
          FNAVB(MLTM)='avbts' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRP(MLTM).EQ.1) THEN
          FNSEL(MLTM)='selts' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRUE(MLTM).EQ.1) THEN
          FNUVE(MLTM)='uvets' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRUT(MLTM).EQ.1) THEN
          FNUVT(MLTM)='uvtts' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRU(MLTM).EQ.1) THEN
          FNU3D(MLTM)='u3dts' // CNTMSR(MLTM) // '.out'
          FNV3D(MLTM)='v3dts' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRQE(MLTM).EQ.1) THEN
          FNQQE(MLTM)='qqets' // CNTMSR(MLTM) // '.out'
        END IF
        IF (MTMSRQ(MLTM).EQ.1) THEN
          FNQ3D(MLTM)='q3dts' // CNTMSR(MLTM) // '.out'
        END IF
      END DO
      JSTMSR=0
      END IF
C
      IF(JSTMSR.EQ.0) GO TO 300
C
      DO MLTM=1,MLTMSR
        IF (MTMSRC(MLTM).EQ.1) THEN
          IF(ISTRAN(1).GE.1) THEN
            FNSAL(MLTM)='salts' // CNTMSR(MLTM) // '.out'
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            CLOSE(11,STATUS='DELETE')
            OPEN(11,FILE=FNSAL(MLTM),STATUS='UNKNOWN')
            WRITE (11,100) TITLE1
            WRITE (11,101) CLTMSR(MLTM)
            WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (11,102) CTUNIT
            CLOSE(11)
          END IF
          IF (ISTRAN(2).GE.1) THEN
            FNTEM(MLTM)='temts' // CNTMSR(MLTM) // '.out'
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            CLOSE(21,STATUS='DELETE')
            OPEN(21,FILE=FNTEM(MLTM),STATUS='UNKNOWN')
            WRITE (21,100) TITLE2
            WRITE (21,101) CLTMSR(MLTM)
            WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (21,102) CTUNIT
            CLOSE(21)
          END IF
          IF (ISTRAN(3).GE.1) THEN
            FNDYE(MLTM)='dyets' // CNTMSR(MLTM) // '.out'
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNDYE(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          END IF
          IF (ISTRAN(4).GE.1) THEN
            FNDYE(MLTM)='sflts' // CNTMSR(MLTM) // '.out'
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            CLOSE(31,STATUS='DELETE')
            OPEN(31,FILE=FNSFL(MLTM),STATUS='UNKNOWN')
            WRITE (31,100) TITLE3
            WRITE (31,101) CLTMSR(MLTM)
            WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (31,102) CTUNIT
            CLOSE(31)
          END IF
          IF (ISTRAN(6).GE.1) THEN
            FNSED(MLTM)='sedts' // CNTMSR(MLTM) // '.out'
            OPEN(41,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNSED(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
          END IF
          IF (ISTRAN(7).GE.1) THEN
            FNSND(MLTM)='sndts' // CNTMSR(MLTM) // '.out'
            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
            CLOSE(41,STATUS='DELETE')
            OPEN(41,FILE=FNSND(MLTM),STATUS='UNKNOWN')
            WRITE (41,100) TITLE4
            WRITE (41,101) CLTMSR(MLTM)
            WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (41,102) CTUNIT
            CLOSE(41)
          END IF
          IF (ISTRAN(5).GE.1) THEN
            DO NT=1,NTOX1
            FNTOX(MLTM,NT)='t'// CNTOX(NT) // 'ts' //
     $                        CNTMSR(MLTM) // '.out'
            OPEN(51,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
            CLOSE(51,STATUS='DELETE')
            OPEN(51,FILE=FNTOX(MLTM,NT),STATUS='UNKNOWN')
            WRITE (51,100) TITLE5
            WRITE (51,101) CLTMSR(MLTM)
            WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
            WRITE (51,102) CTUNIT
            if(NPCB.GT.0) THEN
            WRITE (51,*)'1 TIME'
            WRITE (51,*)'2 TOXB2'
            WRITE (51,*)'3 TOXB1'
            WRITE (51,*)'4-3+2KC TOX TOXPFTW'
            WRITE (51,*)'22 SEDPOC2'
            WRITE (51,*)'23 SEDPOC1'
            WRITE (51,*)'24 TOXPFTW1'
            WRITE (51,*)'25 TOXPFTW2'
            WRITE (51,*)'26 27 28 POC Alae DOC'
            WRITE (51,*)'29 30 TOXB2*B2 TOXB1*B1 no correction'
            WRITE (51,*)'31 32 pro2 pro1'
            WRITE (51,*)'33 34 dep2 dep1'
            endif
            CLOSE(51)
            END DO
          END IF
        END IF
        IF (MTMSRA(MLTM).EQ.1) THEN
          FNAVV(MLTM)='avvts' // CNTMSR(MLTM) // '.out'
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE=FNAVV(MLTM),STATUS='UNKNOWN')
          WRITE (61,100) TITLE6
          WRITE (61,101) CLTMSR(MLTM)
          WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (61,102) CTUNIT
          CLOSE(61)
          FNAVB(MLTM)='avbts' // CNTMSR(MLTM) // '.out'
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNAVB(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE7
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
        END IF
        IF (MTMSRP(MLTM).EQ.1) THEN
          FNSEL(MLTM)='selts' // CNTMSR(MLTM) // '.out'
          OPEN(11,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          CLOSE(11,STATUS='DELETE')
          OPEN(11,FILE=FNSEL(MLTM),STATUS='UNKNOWN')
          WRITE (11,100) TITLE11
          WRITE (11,101) CLTMSR(MLTM)
          WRITE (11,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (11,102) CTUNIT
          CLOSE(11)
        END IF
        IF (MTMSRUE(MLTM).EQ.1) THEN
          FNUVE(MLTM)='uvets' // CNTMSR(MLTM) // '.out'
          OPEN(21,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          CLOSE(21,STATUS='DELETE')
          OPEN(21,FILE=FNUVE(MLTM),STATUS='UNKNOWN')
          WRITE (21,100) TITLE12
          WRITE (21,101) CLTMSR(MLTM)
          WRITE (21,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (21,102) CTUNIT
          CLOSE(21)
        END IF
        IF (MTMSRUT(MLTM).EQ.1) THEN
          FNUVT(MLTM)='uvtts' // CNTMSR(MLTM) // '.out'
          OPEN(31,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          CLOSE(31,STATUS='DELETE')
          OPEN(31,FILE=FNUVT(MLTM),STATUS='UNKNOWN')
          WRITE (31,100) TITLE13
          WRITE (31,101) CLTMSR(MLTM)
          WRITE (31,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (31,102) CTUNIT
          CLOSE(31)
        END IF
        IF (MTMSRU(MLTM).EQ.1) THEN
          FNU3D(MLTM)='u3dts' // CNTMSR(MLTM) // '.out'
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          CLOSE(41,STATUS='DELETE')
          OPEN(41,FILE=FNU3D(MLTM),STATUS='UNKNOWN')
          WRITE (41,100) TITLE14
          WRITE (41,101) CLTMSR(MLTM)
          WRITE (41,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (41,102) CTUNIT
          CLOSE(41)
          FNV3D(MLTM)='v3dts' // CNTMSR(MLTM) // '.out'
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          CLOSE(51,STATUS='DELETE')
          OPEN(51,FILE=FNV3D(MLTM),STATUS='UNKNOWN')
          WRITE (51,100) TITLE15
          WRITE (51,101) CLTMSR(MLTM)
          WRITE (51,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (51,102) CTUNIT
          CLOSE(51)
        END IF
        IF (MTMSRQE(MLTM).EQ.1) THEN
          FNQQE(MLTM)='qqets' // CNTMSR(MLTM) // '.out'
          OPEN(61,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          CLOSE(61,STATUS='DELETE')
          OPEN(61,FILE=FNQQE(MLTM),STATUS='UNKNOWN')
          WRITE (61,100) TITLE16
          WRITE (61,101) CLTMSR(MLTM)
          WRITE (61,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (61,102) CTUNIT
          CLOSE(61)
        END IF
        IF (MTMSRQ(MLTM).EQ.1) THEN
          FNQ3D(MLTM)='q3dts' // CNTMSR(MLTM) // '.out'
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          CLOSE(71,STATUS='DELETE')
          OPEN(71,FILE=FNQ3D(MLTM),STATUS='UNKNOWN')
          WRITE (71,100) TITLE17
          WRITE (71,101) CLTMSR(MLTM)
          WRITE (71,103)ILTMSR(MLTM),JLTMSR(MLTM)
          WRITE (71,102) CTUNIT
          CLOSE(71)
        END IF
      END DO
C
      JSTMSR=0
C
C----------------------------------------------------------------------C
C
  300 CONTINUE
C
C----------------------------------------------------------------------C
C
      TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
C
      FOURDPI=4./PI
C
C **  STEP CURRENT TIME INTERVALS FOR WRITE SCENARIOS
C
      DO NTSSS=1,NTSSTSP
       DO MTSSS=1,MTSSTSP(NTSSS)
c	write(302,*) NTSSS,MTSSS,N,NTSSTSP,MTSSTSP(NTSSS)
c	write(302,*) time,tsstrt(MTSSS,NTSSS),TSSTOP(MTSSS,NTSSS)
cwrite(301,3011) NTSSS,MTSSS,N,time,TSSTRT(MTSSS,NTSSS),
c    *  TSSTOP(MTSSS,NTSSS)
3011	format(3i6,999e12.4)
        IF(TIME.GE.TSSTRT(MTSSS,NTSSS)) THEN
        IF(TIME.LE.TSSTOP(MTSSS,NTSSS)) THEN
          MTSCUR(NTSSS)=MTSSS
        END IF
        END IF
       END DO
      END DO
C
      DO MLTM=1,MLTMSR
       NTSSS=NTSSSS(MLTM)
       MTSCC=MTSCUR(NTSSS)
c      write(302,*) MLTM,NTSSS,MTSCC
       IF(TIME.GE.TSSTRT(MTSCC,NTSSS)) THEN
       IF(TIME.LE.TSSTOP(MTSCC,NTSSS)) THEN
        I=ILTMSR(MLTM)
        J=JLTMSR(MLTM)
        L=LIJ(I,J)
        LN=LNC(L)
        IF (MTMSRC(MLTM).EQ.1) THEN
          IF(ISTRAN(1).GE.1) THEN
            OPEN(11,FILE=FNSAL(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
            WRITE(11,201)TIME,(SAL(L,K),K=1,KC)
            CLOSE(11)
          END IF
          IF (ISTRAN(2).GE.1) THEN
            OPEN(21,FILE=FNTEM(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
            WRITE(21,201)TIME,(TEM(L,K),K=1,KC)
            CLOSE(21)
          END IF
          IF (ISTRAN(3).GE.1) THEN
            OPEN(31,FILE=FNDYE(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
            WRITE(31,201)TIME,(DYE(L,K),K=1,KC)
            CLOSE(31)
          END IF
          IF (ISTRAN(4).GE.1) THEN
            OPEN(31,FILE=FNSFL(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
            WRITE(31,201)TIME,(SFL(L,K),K=1,KC)
            CLOSE(31)
          END IF
          IF (ISTRAN(6).GE.1) THEN
            OPEN(41,FILE=FNSED(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
c           WRITE(41,201)TIME,SEDBT(L,KBT(L)),(SEDT(L,K),K=1,KC) ! Ji, 10/22/00
            WRITE(41,201)TIME,(SEDBT(L,K),k=1,KB),(SEDT(L,K),K=1,KC)
            CLOSE(41)
          END IF
          IF (ISTRAN(7).GE.1) THEN
            OPEN(41,FILE=FNSND(MLTM),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
            WRITE(41,201)TIME,SNDBT(L,KBT(L)),(SNDT(L,K),K=1,KC)
            CLOSE(41)
          END IF

          IF (ISTRAN(5).GE.1) THEN
            DO NT=1,NTOX1
            OPEN(51,FILE=FNTOX(MLTM,NT),ACCESS='APPEND',
     $                                  STATUS='UNKNOWN')
           k1=1
           a_carb=(WQV(L,K1,4)+WQV(L,K1,5))
           a_agle=(WQV(L,K1,1)+WQV(L,K1,2)+WQV(L,K1,3))
           a_doc=WQV(L,K1,6)
           TOXBTMP=TOXB(L,KBP,NT)
           TOXBTMP1=TOXB(L,KBP-1,NT)
           IF(ITXBDUT(NT).EQ.1) THEN
           TOXBTMP=TOXB(L,KBP,NT)/HBEDP(L,KBP)/
     &     (2500.0*(1-PORBEDP(L,KBP)))
           TOXBTMP1=TOXB(L,KBP-1,NT)/HBEDP(L,KBP-1)
     &      /(2500.0*(1-PORBEDP(L,KBP-1)))
           ENDIF
  !          IF(VOLBW2(L,KB).GT.1.E-12)TOXBTMP=TOXB(L,KB,NT)/VOLBW2(L,KB)
            WRITE(51,201)TIME,TOXBTMP,TOXBTMP1,
     &      (TOX(L,K,NT),K=1,KC),
     &      (TOXPFTW(L,K,NT),K=1,KC),SEDPOC(L,2),SEDPOC(L,1),
     &       TOXPFTB(L,2,NT),TOXPFTB(L,1,NT),
     &       a_carb, a_agle,a_doc,TOXB(L,KBP,NT),TOXB(L,KBP-1,NT),
     &       PORBEDP(L,KBP),PORBEDP(L,KBP-1),HBEDP(L,2),HBEDP(L,1)
            CLOSE(51)
            END DO
          END IF
        END IF
        IF (MTMSRA(MLTM).EQ.1) THEN
          OPEN(61,FILE=FNAVV(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          OPEN(71,FILE=FNAVB(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
           DO K=1,KS
           ATMP(K)=10000.*AV(L,K)*HP(L)
           END DO
          WRITE(61,201)TIME,(ATMP(K),K=1,KS)
           DO K=1,KS
           ATMP(K)=10000.*AB(L,K)*HP(L)
           END DO
          WRITE(71,201)TIME,(ATMP(K),K=1,KS)
          CLOSE(61)
          CLOSE(71)
        END IF
        IF (MTMSRP(MLTM).EQ.1) THEN
          OPEN(11,FILE=FNSEL(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          PPTMP=GI*P(L)
          HHTMP=PPTMP-BELV(L)
          GWELTMP=AGWELV(L)-BELAGW(L)
          WRITE(11,201)TIME,PPTMP,HHTMP,GWELTMP
          CLOSE(11)
        END IF
        IF (MTMSRUE(MLTM).EQ.1) THEN
          OPEN(21,FILE=FNUVE(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          UTMP1=50.*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=50.*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          IF(SPB(L).EQ.0) THEN
            UTMP1=2.*UTMP1
            VTMP1=2.*VTMP1
          END IF
          UTMP=UTMP1 !CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=VTMP1 !CUN(L)*UTMP1+CVN(L)*VTMP1
          UTMP1=5000.*(TBX(L+1)+TBX(L))
          VTMP1=5000.*(TBY(LN)+TBY(L))
          TBEAST=CUE(L)*UTMP1+CVE(L)*VTMP1
          TBNORT=CUN(L)*UTMP1+CVN(L)*VTMP1
          TAUB=1000.*QQ(L,0)/CTURB2
          TAUW1=1000.*QQWV1(L)
          TAUW2=1000.*QQWV2(L)
         ! UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
         ! VTMP=0.5*STCUV(L)*(V(LN,1)+V(L,1))
          CURANG=ATAN2(VTMP,UTMP)
          TAUB2=TAUB*TAUB+0.5*(TAUW2*TAUW2)
     $           +FOURDPI*TAUB*TAUW2*COS(CURANG-WACCWE(L))
          TAUB2=MAX(TAUB2,0.)
          TAUTOT=SQRT(TAUB2)
          WRITE(21,201)TIME,UTMP,VTMP,TBEAST,TBNORT,TAUB,TAUW1,
     $                 TAUW2,TAUTOT
          CLOSE(21)
        END IF
        IF (MTMSRUT(MLTM).EQ.1) THEN
          OPEN(31,FILE=FNUVT(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          WRITE(31,201)TIME,UHDYE(L),VHDXE(L)
          CLOSE(31)
        END IF
        IF (MTMSRU(MLTM).EQ.1) THEN
          OPEN(41,FILE=FNU3D(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          OPEN(51,FILE=FNV3D(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          RUVTMP=50.
          IF(SPB(L).EQ.0) RUVTMP=100.
          DO K=1,KC
           UTMP1=RUVTMP*(U(L+1,K)+U(L,K))
           VTMP1=RUVTMP*(V(LN,K)+V(L,K))
           ATMP(K)=UTMP1 !CUE(L)*UTMP1+CVE(L)*VTMP1
           BTMP(K)=VTMP1 !CUN(L)*UTMP1+CVN(L)*VTMP1
          END DO
          WRITE(41,201)TIME,(ATMP(K),K=1,KC) ! rotated into cartisina direction & in cm/s!
          WRITE(51,201)TIME,(BTMP(K),K=1,KC)
          CLOSE(41)
          CLOSE(51)
        END IF
        IF (MTMSRQE(MLTM).EQ.1) THEN
          OPEN(61,FILE=FNQQE(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          QRNT=DXYP(L)*RAINT(L)
          WRITE(61,201)TIME,QSUME(L),QRNT,EVAPSW(L),EVAPGW(L),RIFTR(L)
          CLOSE(61)
        END IF
        IF (MTMSRQ(MLTM).EQ.1) THEN
          OPEN(71,FILE=FNQ3D(MLTM),ACCESS='APPEND',STATUS='UNKNOWN')
          WRITE (71,201)TIME,(QSUM(L,K),K=1,KC)
          CLOSE(71)
        END IF
       END IF
       END IF
      END DO
C
C**********************************************************************C
C
  100 FORMAT(A80)
  101 FORMAT(1X,'AT LOCATION  ',A20)
  102 FORMAT(1X,'TIME IN FIRST COLUMN HAS UNITS OF ',A10)
  103 FORMAT(1X,'CELL I,J = ',2I5)
  201 FORMAT(F12.5,130E12.4)
C
C**********************************************************************C
C
      RETURN
      END
