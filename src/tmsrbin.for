C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE TMSR
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 APRIL 1998
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
      IF(JSTMSR.EQ.0) GO TO 300
C
      DO ii=511,526
      close(ii)
      ENDDO
      
      DO MLTM=1,1
        IF (MTMSRC(MLTM).EQ.1) THEN
          IF(ISTRAN(1).GE.1) THEN
            FNSAL(MLTM)='salts'//iyear//'.bin'
            OPEN(511,FILE=FNSAL(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(511,STATUS='DELETE')
            OPEN(511,FILE=FNSAL(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (511) KC+2
          END IF
          IF (ISTRAN(2).GE.1) THEN
            FNTEM(MLTM)='temts'//iyear//'.bin'
            OPEN(512,FILE=FNTEM(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(512,STATUS='DELETE')
            OPEN(512,FILE=FNTEM(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (512) KC+2
          END IF
          IF (ISTRAN(3).GE.1) THEN
            FNDYE(MLTM)='dyets'//iyear//'.bin'
            OPEN(513,FILE=FNDYE(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(513,STATUS='DELETE')
            OPEN(513,FILE=FNDYE(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (513) KC+2
           END IF
          IF (ISTRAN(4).GE.1) THEN
            FNSFL(MLTM)='sflts' //iyear//'.bin'
            OPEN(514,FILE=FNSFL(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(514,STATUS='DELETE')
            OPEN(514,FILE=FNSFL(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (514) KC+2
          END IF
          IF (ISTRAN(6).GE.1) THEN
            FNSED(MLTM)='sedts'//iyear//'.bin'
            OPEN(515,FILE=FNSED(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(515,STATUS='DELETE')
            OPEN(515,FILE=FNSED(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (515) KC+3
          END IF
          IF (ISTRAN(7).GE.1) THEN
            FNSND(MLTM)='sndts'//iyear//'.bin'
            OPEN(516,FILE=FNSND(MLTM),STATUS='UNKNOWN',form='unformatted')
            CLOSE(516,STATUS='DELETE')
            OPEN(516,FILE=FNSND(MLTM),STATUS='UNKNOWN',form='unformatted')
            WRITE (516)  KC+3
          END IF
          IF (ISTRAN(5).GE.1) THEN
            FNTOX(MLTM,1)='tox'//iyear//'.bin'
            OPEN(517,FILE=FNTOX(MLTM,1),STATUS='UNKNOWN',form='unformatted')
            CLOSE(517,STATUS='DELETE')
            OPEN(517,FILE=FNTOX(MLTM,1),STATUS='UNKNOWN',form='unformatted')
            WRITE (517) NTOX,KC+4
          END IF
        END IF
        IF (MTMSRA(MLTM).EQ.1) THEN
          FNAVV(MLTM)='avvts'//iyear//'.bin'
          OPEN(518,FILE=FNAVV(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(518,STATUS='DELETE')
          OPEN(518,FILE=FNAVV(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (518)KC+1
          FNAVB(MLTM)='avbts'//iyear//'.bin'
          OPEN(519,FILE=FNAVB(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(519,STATUS='DELETE')
          OPEN(519,FILE=FNAVB(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (519) KC+1
        END IF
        IF (MTMSRP(MLTM).EQ.1) THEN
          FNSEL(MLTM)='selts'//iyear//'.bin'
          OPEN(520,FILE=FNSEL(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(520,STATUS='DELETE')
          OPEN(520,FILE=FNSEL(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (520) 5
        END IF
        IF (MTMSRUE(MLTM).EQ.1) THEN
          FNUVE(MLTM)='uvets' //iyear//'.bin'
          OPEN(521,FILE=FNUVE(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(521,STATUS='DELETE')
          OPEN(521,FILE=FNUVE(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (521) 6
        END IF
        IF (MTMSRUT(MLTM).EQ.1) THEN
          FNUVT(MLTM)='uvtts' //iyear//'.bin'
          OPEN(522,FILE=FNUVT(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(522,STATUS='DELETE')
          OPEN(522,FILE=FNUVT(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE(522)3
        END IF
        IF (MTMSRU(MLTM).EQ.1) THEN
          FNU3D(MLTM)='u3dts' //iyear//'.bin'
          OPEN(523,FILE=FNU3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(523,STATUS='DELETE')
          OPEN(523,FILE=FNU3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (523) KC+3
          FNV3D(MLTM)='v3dts'//iyear//'.bin'
          OPEN(524,FILE=FNV3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(524,STATUS='DELETE')
          OPEN(524,FILE=FNV3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (524) KC+3
        END IF
        IF (MTMSRQE(MLTM).EQ.1) THEN
          FNQQE(MLTM)='qqets' //iyear//'.bin'
          OPEN(525,FILE=FNQQE(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(525,STATUS='DELETE')
          OPEN(525,FILE=FNQQE(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (525) 7
        END IF
        IF (MTMSRQ(MLTM).EQ.1) THEN
          FNQ3D(MLTM)='q3dts' //iyear//'.bin'
          OPEN(526,FILE=FNQ3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          CLOSE(526,STATUS='DELETE')
          OPEN(526,FILE=FNQ3D(MLTM),STATUS='UNKNOWN',form='unformatted')
          WRITE (526) KC+1
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
      IF(NCSTEP.GT.0)TIME=SECDLAST/TCON+TBEGIN  !% J.S. 1/31/2014 
C
C **  STEP CURRENT TIME INTERVALS FOR WRITE SCENARIOS
C
      DO NTSSS=1,NTSSTSP
       DO MTSSS=1,MTSSTSP(NTSSS)
        IF(TIME.GE.TSSTRT(MTSSS,NTSSS)) THEN
        IF(TIME.LE.TSSTOP(MTSSS,NTSSS)) THEN
          MTSCUR(NTSSS)=MTSSS
        END IF
        END IF
       END DO
      END DO
C
      DO MLTM=1,MLTMSR  ! number of stations to output
       NTSSS=NTSSSS(MLTM)
       MTSCC=MTSCUR(NTSSS)
       IF(TIME.GE.TSSTRT(MTSCC,NTSSS)) THEN
       IF(TIME.LE.TSSTOP(MTSCC,NTSSS)) THEN
        I=ILTMSR(MLTM)
        J=JLTMSR(MLTM)
        L=LIJ(I,J)
        LN=LNC(L)
          IF (MTMSRC(1).EQ.1) THEN               ! output if it is turn on at 1st station
            WRITE(511)TIME,MLTM,(SAL(L,K),K=1,KC)
          END IF
          IF (ISTRAN(2).GE.1) THEN
            WRITE(512)TIME,MLTM,(TEM(L,K),K=1,KC)
          END IF
          IF (ISTRAN(3).GE.1) THEN
            WRITE(513)TIME,MLTM,(DYE(L,K),K=1,KC)
          END IF
          IF (ISTRAN(4).GE.1) THEN
            WRITE(514)TIME,MLTM,(SFL(L,K),K=1,KC)
          END IF
          IF (ISTRAN(6).GE.1) THEN
            SEDBT(L,1)=SEDT(L,1)/(1.0e-6+DYE(L,1))/86400
            WRITE(515)TIME,MLTM,SEDBT(L,1),(SEDT(L,K),K=1,KC)
          END IF
          IF (ISTRAN(7).GE.1) THEN
            WRITE(516)TIME,MLTM,SNDBT(L,1),(SNDT(L,K),K=1,KC)
          END IF
          IF (ISTRAN(5).GE.1) THEN
            DO NT=1,NPCB
           K1=KC
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
            WRITE(517)TIME,TOXBTMP,TOXBTMP1,
     &      (TOX(L,K,NT),K=1,KC),
     &      (TOXPFTW(L,K,NT),K=1,KC),SEDPOC(L,2),SEDPOC(L,1),
     &       TOXPFTB(L,2,NT),TOXPFTB(L,1,NT),
     &       a_carb, a_agle,a_doc,TOXB(L,KBP,NT),TOXB(L,KBP-1,NT),
     &       PORBEDP(L,KBP),PORBEDP(L,KBP-1),HBEDP(L,2),HBEDP(L,1)
            END DO
 !           DO NT=1,NTOX
 !           TOXBTMP=0.
 !!           IF(VOLBW2(L).GT.1.E-12)TOXBTMP=TOXB(L,1,NT)/VOLBW2(L)
 !           WRITE(517)TIME,MLTM,NT,TOXBTMP,(TOX(L,K,NT),K=1,KC)
 !           END DO
          END IF
          IF (MTMSRA(1).EQ.1) THEN
           DO K=1,KS
           ATMP(K)=10000.*AV(L,K)*HP(L)
           END DO
           WRITE(518)TIME,MLTM,(ATMP(K),K=1,KS)
           DO K=1,KS
           ATMP(K)=10000.*AB(L,K)*HP(L)
           END DO
           WRITE(519)TIME,MLTM,(ATMP(K),K=1,KS)
         END IF
         IF (MTMSRP(1).EQ.1) THEN
          PPTMP=GI*P(L)
          HHTMP=PPTMP-BELV(L)
          GWELTMP=AGWELV(L)-BELAGW(L)
          WRITE(520)TIME,MLTM,PPTMP,HHTMP,GWELTMP
         END IF
         IF (MTMSRUE(1).EQ.1) THEN
          UTMP1=50.*(UHDYE(L+1)+UHDYE(L))/(DYP(L)*HP(L))
          VTMP1=50.*(VHDXE(LN)+VHDXE(L))/(DXP(L)*HP(L))
          IF(SPB(L).EQ.0) THEN
            UTMP1=2.*UTMP1
            VTMP1=2.*VTMP1
          END IF
          UTMP=CUE(L)*UTMP1+CVE(L)*VTMP1
          VTMP=CUN(L)*UTMP1+CVN(L)*VTMP1
          UTMP1=5000.*(TBX(L+1)+TBX(L))
          VTMP1=5000.*(TBY(LN)+TBY(L))
          TBEAST=CUE(L)*UTMP1+CVE(L)*VTMP1
          TBNORT=CUN(L)*UTMP1+CVN(L)*VTMP1
          WRITE(521)TIME,MLTM,UTMP,VTMP,TBEAST,TBNORT
         END IF
         IF (MTMSRUT(1).EQ.1) THEN
          WRITE(522)TIME,UHDYE(L),VHDXE(L)
         END IF
         IF (MTMSRU(MLTM).EQ.1) THEN
          RUVTMP=50.
          IF(SPB(L).EQ.0) RUVTMP=100.
          DO K=1,KC
           UTMP1=2.0*RUVTMP*U(L,K)  !(U(L+1,K)+U(L,K))
           VTMP1=2.0*RUVTMP*V(L,K) !(V(LN,K)+V(L,K))
       !    ATMP(K)=CUE(L)*UTMP1+CVE(L)*VTMP1
       !    BTMP(K)=CUN(L)*UTMP1+CVN(L)*VTMP1
           ATMP(K)=UTMP1
           BTMP(K)=VTMP1
          END DO
          WRITE(523)TIME,MLTM,(ATMP(K),K=1,KC) !,(STBX(L+1)+STBX(L))*0.5
          WRITE(524)TIME,MLTM,(BTMP(K),K=1,KC) !,(STBY(LN)+STBY(L))*0.5
         END IF
         IF (MTMSRQE(1).EQ.1) THEN
          QRNT=DXYP(L)*RAINT(L)
          WRITE(525)TIME,MLTM,QSUME(L),QRNT,EVAPSW(L),EVAPGW(L),RIFTR(L)
         END IF
         IF (MTMSRQ(1).EQ.1) THEN
          WRITE (526)TIME,(QSUM(L,K),K=1,KC)
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
  201 FORMAT(F12.5, 1p, 100E12.4)
C
C**********************************************************************C
C
      RETURN
      END
