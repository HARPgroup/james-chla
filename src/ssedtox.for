C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE SSEDTOX(ISTL,CORDT)
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C**********************************************************************C
C
C **  SUBROUTINE SSEDTOX CALCULATES SETTLING AND WATER COLUMN-BED
C **  EXCHANGE OF SEDIMENT AND SORBED TOXIC CONTAMINANTS
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
      DIMENSION STRESSS(0:KSM)
      DIMENSION CTMPDRY(LCM)
      DIMENSION WSETA(LCM,0:KSM,NSTM)
      DIMENSION SEDS(LCM,KCM,NSCM),SNDS(LCM,KCM,NSNM),
     $          TOXS(LCM,KCM,NTXM)
      DIMENSION SEDBS(LCM,KBM,NSCM),SNDBS(LCM,KBM,NSNM),
     $           TOXBS(LCM,KBM,NTXM)
      DIMENSION  CDECAYW(LCM,KCM),CDECAYB(LCM,KBM)
      DIMENSION  SEDDIA50(LCM,KBM),SEDBALL(LCM,KBM)
      DIMENSION  ALOW(LCM,KBM+1),BMNN(LCM,KBM+1),CUPP(LCM,KBM+1),
     $         RRHS(LCM,KBM+1),TOXTMP(LCM,KBM+1),GAMTMP(LCM,KBM+1)
      DIMENSION    QWBDTOP(LCM),QSBDTOP(LCM),ZETATOP(LCM),
     $             ZBEDGT(LCM)
      DIMENSION STRSE(LCM,KBM),HYDCN(LCM,KBM),CCOEF(LCM,KBM),
     $          DZBTR(LCM,KBM),STRSEM(LCM,KBM),PRESE(LCM,KBM),
     $          PRESH(LCM,KBM),PREST(LCM,KBM),STRST(LCM,KBM),
     $          DZBTR1(LCM,KBM),ZBEDG(LCM,KBM),ZBEDC(LCM,KBM),
     $          SGSM1(LCM,KBM)
C
C**********************************************************************C
C
      DELT=DT2
      S3TL=1.0
      S2TL=0.0
      ISUD=1
      IF (ISTL.NE.3) THEN
        DELT=DT
        S3TL=0.0
        S2TL=1.0
        ISUD=0
      END IF
      DELTI=1./DELT
C
      SEDMDGM=SQRT(SEDMDMX*SEDMDMN)
cbegmod
      BEDEX=1.
cendmod
cjh      BEDEX=0.
      NVAL=MOD(N,2)
CJH      IF(NVAL.EQ.0) THEN
CJH        IF(ISTL.EQ.3) BEDEX=1.
CJH      END IF
C
      TIME=DT*FLOAT(N)+TCON*TBEGIN
      TIME=TIME/TCON
      FOURDPI=4./PI
C
      DO L=2,LA
       CTMPDRY(L)=1.
       IF(ISCDRY(L).NE.0) CTMPDRY(L)=0.
      END DO
C
C**********************************************************************C
C
      IF(N.EQ.1) THEN
        OPEN(1,FILE='ssedtox.dia',STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE='ssedtox.dia',STATUS='UNKNOWN')
        OPEN(11,FILE='ssedtox.dia1',STATUS='UNKNOWN')
        CLOSE(11,STATUS='DELETE')
        OPEN(11,FILE='ssedtox.dia1',STATUS='UNKNOWN')
        OPEN(21,FILE='ssedtox.dia2',STATUS='UNKNOWN')
        CLOSE(21,STATUS='DELETE')
        OPEN(21,FILE='ssedtox.dia2',STATUS='UNKNOWN')
        OPEN(31,FILE='ssedtox.dia3',STATUS='UNKNOWN')
        CLOSE(31,STATUS='DELETE')
        OPEN(31,FILE='ssedtox.dia3',STATUS='UNKNOWN')
        OPEN(41,FILE='ssedtox.dia4',STATUS='UNKNOWN')
        CLOSE(41,STATUS='DELETE')
        OPEN(41,FILE='ssedtox.dia4',STATUS='UNKNOWN')
       ELSE
        OPEN(1,FILE='ssedtox.dia',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(11,FILE='ssedtox.dia1',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(21,FILE='ssedtox.dia2',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(31,FILE='ssedtox.dia3',ACCESS='APPEND',STATUS='UNKNOWN')
        OPEN(41,FILE='ssedtox.dia4',ACCESS='APPEND',STATUS='UNKNOWN')
      END IF
C
C**********************************************************************C
C
C **  IF N=1 CALCULATE INITIAL SEDIMENT BED THICKNESS
C **  MOVED TO SUBROUTINE BEDINIT  IN APRIL 10, 1999 VERSION
C
C      IF(N.EQ.1) THEN
C
C      DO K=1,KB
C      DO L=2,LA
C       HBED(L,K)=0.
C       SEDBALL(L,K)=0.
C      END DO
C      END DO
C
C      IF(ISTRAN(6).GE.1) THEN
C      DO NS=1,NSED
C       DO K=1,KB
C       DO L=2,LA
C        HBED(L,K)=HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
C        SEDBALL(L,K)=SEDBALL(L,K)+SEDB(L,K,NS)
C       END DO
C       END DO
C      END DO
C      END IF
C
C      IF(ISTRAN(7).GE.1) THEN
C      DO NX=1,NSND
C       NS=NSED+NX
C       DO K=1,KB
C       DO L=2,LA
C        HBED(L,K)=HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)
C        SEDBALL(L,K)=SEDBALL(L,K)+SNDB(L,K,NX)
C       END DO
C       END DO
C      END DO
C      END IF
C
C      TMPVAL=1./(1.-PORBED(L,K))
C      DO K=1,KB
C      DO L=2,LA
C       HBED(L,K)=TMPVAL*HBED(L,K)
c       HBED(L,K)=MAX(1.E-9,HBED(L,K))
C       VOLBW2(L,K)=HBED(L,K)
C      END DO
C      END DO
C
C ** DIAGNOSTICS OF INITIALIZATION
C
C       OPEN(2,file='depbed.dia')
C       CLOSE(2,STATUS='DELETE')
C       OPEN(2,file='depbed.dia')
C       DO L=2,LA
C        WRITE(2,2222)IL(L),JL(L),HP(L),SEDB(L,1,1),SNDB(L,1,1),HBED(L,1)
C       END DO
C       CLOSE(2)
c
C      END IF
C
C**********************************************************************C
C
C **   UPDATE SEDIMENT PROCESSES
C
C----------------------------------------------------------------------C
C
C **  CALCULATE TOTAL SEDIMENT IN THE BED
C
      DO K=1,KB
      DO L=1,LC
       SEDBT(L,K)=0.
       SNDBT(L,K)=0.
       SEDBALL(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)
c     write(183,1831) N,K,L,sedbt(L,K),sedb(L,K,ns)
1831  format(3i5,999e12.4)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NS=1,NSND
       DO K=1,KB
       DO L=2,LA
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      DO K=1,KB
      DO L=1,LC
       SEDBALL(L,K)=SEDBT(L,K)+SNDBT(L,K)
      END DO
      END DO
C
C----------------------------------------------------------------------C
C
C **  SET SEDIMENT VOLUME FRACTIONS
C
      DO K=1,KB
      DO L=2,LA
       BEDLINIT(L,K)=0.
       BEDDINIT(L,K)=0.
      END DO
      END DO
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SEDB(L,K,NS)
        VFRBED1(L,K,NS)=SDEN(NS)*SEDB1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        VFRBED(L,K,NS)=SDEN(NS)*SNDB(L,K,NX)
        VFRBED1(L,K,NS)=SDEN(NS)*SNDB1(L,K,NX)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        BEDLINIT(L,K)=BEDLINIT(L,K)+VFRBED(L,K,NS)
        BEDDINIT(L,K)=BEDDINIT(L,K)+VFRBED1(L,K,NS)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KB
       DO L=2,LA
c       VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
c       VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/(BEDLINIT(L,K)+1.0e-10)   ! to avoid bedLinit=0.0 -> NaN, Ji, 2/3/03
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/(BEDDINIT(L,K)+1.0e-10)
       END DO
       END DO
      END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
      NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
c       VFRBED(L,K,NS)=VFRBED(L,K,NS)/BEDLINIT(L,K)
c       VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/BEDDINIT(L,K)
        VFRBED(L,K,NS)=VFRBED(L,K,NS)/(BEDLINIT(L,K)+1.0e-10)    ! to avoid bedLinit=0.0 -> NaN, Ji, 2/3/03
        VFRBED1(L,K,NS)=VFRBED1(L,K,NS)/(BEDDINIT(L,K)+1.0e-10)
       END DO
       END DO
      END DO
      END IF
C
      DO L=2,LA
       QWBDTOP(L)=0.
       QSBDTOP(L)=0.
      END DO
C
C----------------------------------------------------------------------C
C
C **  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES
C
      IF(ISTRAN(6).GE.1) THEN
      IF(IWRSP(1).EQ.0) THEN         ! #1
        DO K=1,KB
c
        DO L=2,LA
c         TAURS(L,K)=TAUR(1)    !@, Ji, 2/7/00
          TAURS(L,K)=xLTAUR(L,K)       ! surface layer uses spatially varying xltaur
      if(itmap.eq.0.and.k.lt.kb) taurs(L,k)=taurc(k) ! bottom bed layers use uniform taurc on each layer
          TAURB(L,K)=1.E6
c         WRSPS(L,K)=WRSPO(1)   !@
          WRSPS(L,K)=XLWRSPO(L)
          WRSPB(L,K)=0.0
        END DO
        END DO
c  check
c     if(N.eq.1000) then
c     do L=2,LA
c     write(714,704) L, (taurs(L,K),K=1,kbm)
c704   format(i5,9e12.4)
c      enddo
c      endif
c
      END IF                           ! #1
      IF(IWRSP(1).GE.1) THEN
        DO K=1,KB
        DO L=2,LA
          TAURS(L,K)=CSEDTAUS(BDENBED(L,K),IWRSP(1))
          TAURB(L,K)=CSEDTAUB(BDENBED(L,K),IWRSP(1))
          WRSPS(L,K)=CSEDRESS(BDENBED(L,K),IWRSP(1))
	    aatt=0.0
          WRSPB(L,K)=CSEDRESS(aatt,IWRSP(1))
        END DO
        END DO
      END IF
      END IF
C
C**********************************************************************C
C
C **  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN
C **  BED AND REINITILIZE IF NECESSARY
C
      IF(N.EQ.1.AND.ISTRAN(5).GE.1) THEN
      IF(ISRESTI.EQ.0.OR.ISCI(5).EQ.0) THEN
C
C **  CALCULATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        DO K=1,KB
        DO L=2,LA
          TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        TOXPFTB(L,K,NT)=0.
       END DO
       END DO
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        IF(SEDBALL(L,K).GT.0.0) THEN
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     &                   /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
         ELSE
          TOXPFTB(L,K,NT)=1.
        END IF
       END DO
       END DO
      END DO
C
C **  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC
C **  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG
C **  TO TOXB UNITS OF OF MG/M**2
C
      DO NT=1,NTOX
      IF(ITXBDUT(NT).EQ.0) THEN
        DO K=1,KB
        DO L=2,LA
         TOXB(L,K,NT)=HBED(L,K)*TOXB(L,K,NT)
         TOXB1(L,K,NT)=TOXB(L,K,NT)
        END DO
        END DO
      END IF
      IF(ITXBDUT(NT).EQ.1) THEN
        DO K=1,KB
        DO L=2,LA
         TOXB(L,K,NT)=0.001*TOXB(L,K,NT)*(SEDBT(L,K)+SNDBT(L,K))
     $               /TOXPFTB(L,K,NT)
         TOXB1(L,K,NT)=TOXB(L,K,NT)
        END DO
        END DO
      END IF
      END DO
C
C ** DIAGNOSTICS OF INITIALIZATION
C
      OPEN(2,file='toxbed.dia')
      CLOSE(2,STATUS='DELETE')
      OPEN(2,file='toxbed.dia')
      DO L=2,LA
C       TMP1=-999.
C       TMP2=-999.
C       IF(HBED(L).GT.0.)TMP1=TOXB(L,1)/HBED(L)
C       IF(HBED(L).GT.0.)TMP2=TOXB(L,2)/HBED(L)
C       WRITE(2,2222)IL(L),JL(L),HBED(L),TOXB(L,1),TOXB(L,2),TMP1,TMP2
       TMP1=TOXB(L,1,1)/(HBED(L,1)+1.E-12)
       WRITE(2,2222)IL(L),JL(L),TOXPFTB(L,1,1),TOXB(L,1,1),
     $              TMP1,TOX(L,1,1)
      END DO
      CLOSE(2)
C
      END IF
      END IF
c
 2222 FORMAT(2I5,7E13.4)
C
C**********************************************************************C
C
C **  SAVE OLD VALUES
C
      IF(ISTRAN(5).GE.1) THEN
      DO NT=1,NTOX
       DO K=1,KC
        DO L=2,LA
         TOXS(L,K,NT)=TOX(L,K,NT)
        END DO
       END DO
      END DO
       DO NT=1,NTOX
        DO K=1,KB
        DO L=2,LA
         TOXBS(L,K,NT)=TOXB(L,K,NT)
        END DO
        END DO
       END DO
      END IF
C
      IF(ISTRAN(6).GE.1) THEN
      DO NS=1,NSED
       DO K=1,KC
        DO L=2,LA
         SEDS(L,K,NS)=SED(L,K,NS)
        END DO
       END DO
      END DO
       DO NS=1,NSED
        DO K=1,KB
        DO L=2,LA
         SEDBS(L,K,NS)=SEDB(L,K,NS)
        END DO
        END DO
       END DO
      END IF
C
      IF(ISTRAN(7).GE.1) THEN
      DO NX=1,NSND
       DO K=1,KC
        DO L=2,LA
         SNDS(L,K,NX)=SND(L,K,NX)
        END DO
       END DO
      END DO
       DO NX=1,NSND
        DO K=1,KB
        DO L=2,LA
         SNDBS(L,K,NX)=SNDB(L,K,NX)
        END DO
        END DO
       END DO
      END IF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.1) THEN
        DO NS=1,NSED
         DSEDGMM=1./(1.E6*SSG(NS))   ! =sden, efdc.inp, C39
         IF(WSEDO(NS).GE.1.E-15) THEN
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
c
           K=0
C
           IF(ISEDVW.EQ.0) THEN
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
            END IF
C
           IF(ISEDVW.EQ.1) THEN
             DO L=2,LA
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),0.0,ISEDVW)
             END DO
            END IF
C
           IF(ISEDVW.EQ.2) THEN
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
            END IF
C
           IF(ISEDVW.GE.3) THEN
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=0.5*SQRT(TAUB2)+0.5*QQ(L,1)/CTURB2
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
            END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loop
c
           DO L=2,LA
            SEDF(L,1,NS)=0.
            PROBDEP=0.
            WESE=0.
            TAUBC=QQ(L,0)/CTURB2
            UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
            VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
            CURANG=ATAN2(VTMP,UTMP)
            TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $           +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
            TAUB2=MAX(TAUB2,0.)
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUB.GT.TAURB(L,KBT(L))) THEN
c **  mass erosion
               WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
               WESE=MIN(WESE,WESEMX)
              ELSE
               IF(TAUB.GT.TAURS(L,KBT(L))) THEN
c **  surface erosion
                  WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
                  WESE=MIN(WESE,WESEMX)
                  TAUE=(TAUB-TAURS(L,KBT(L)))/TAURS(L,KBT(L))
                  TAUE=MAX(TAUE,0.0)
                  WESE=WESE*(TAUE**TEXP(NS))
                 ELSE
c **  no erosion
                  WESE=0.0
               END IF
            END IF
c **  set probability of deposition
            IF(TAUB.LT.TAUD(NS)) PROBDEP=(TAUD(NS)-TAUB)/TAUD(NS)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
c            IF(SEDBTMP.LT.0.0) THEN
c              SEDF(L,0,NS)=0.
c              SEDBTMP=SEDB1(L,KBT(L),NS)
c              SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
c            END IF
c            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
c     $                        +S2TL*SEDB1(L,KBT(L),NS)
c            SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0) THEN
              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)+(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            END IF
            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
     $                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     $                  +SEDVDRD*MIN(SEDF(L,0,NS),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.2) THEN
        DO NS=1,NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISEDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.EQ.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),0.0,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.EQ.2) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
             K=1
             DO L=2,LA
              LN=LNC(L)
              SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
              END DO
            END IF
C
           IF(ISEDVW.GE.3) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
             K=1
             DO L=2,LA
              LN=LNC(L)
              STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
            END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loops
c
           K=2
           DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
           END DO
C
           DO L=2,LA
            PROBDEP=0.
            WESE=0.
            TAUBC=QQ(L,0)/CTURB2
            UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
            VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
            CURANG=ATAN2(VTMP,UTMP)
            TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $           +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
            TAUB2=MAX(TAUB2,0.)
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUB.GT.TAURB(L,KBT(L))) THEN
c **  mass erosion
               WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
               WESE=MIN(WESE,WESEMX)
              ELSE
               IF(TAUB.GT.TAURS(L,KBT(L))) THEN
c **  surface erosion
                  WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
                  WESE=MIN(WESE,WESEMX)
                  TAUE=(TAUB-TAURS(L,KBT(L)))/TAURS(L,KBT(L))
                  TAUE=MAX(TAUE,0.0)
                  WESE=WESE*(TAUE**TEXP(NS))
                 ELSE
c **  no erosion
                  WESE=0.0
               END IF
            END IF
c **  set probability of deposition
            IF(TAUB.LT.TAUD(NS)) PROBDEP=(TAUD(NS)-TAUB)/TAUD(NS)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
c            IF(SEDBTMP.LT.0.0) THEN
c              SEDF(L,0,NS)=0.
c              SEDBTMP=SEDB1(L,KBT(L),NS)
c              SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
c            END IF
c            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
c     $                        +S2TL*SEDB1(L,KBT(L),NS)
c            SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0) THEN
              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)+(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            END IF
            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
     $                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     $                  +SEDVDRD*MIN(SEDF(L,0,NS),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.EQ.2
C
           IF(ISTOPT(6).EQ.1) THEN
C
           DO L=2,LA
            CRNUM=1.+DELT*WSETA(L,1,NS)*HPI(L)*DZIC(2)
            GRADSED=(SED(L,2,NS)-SED(L,1,NS))/(DZC(2)+DZC(1))
            SEDAVG=0.5*(SED(L,2,NS)+SED(L,1,NS)+1.E-16)
            WSETA(L,1,NS)=-CRNUM*DZC(2)*WSETA(L,1,NS)*GRADSED/SEDAVG
           END DO
C
           DO L=2,LA
            AA11=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
            AA12=-MAX(WSETA(L,1,NS),0.)
            AA21=MIN(WSETA(L,1,NS),0.)
            AA22=DELTI*DZC(2)*HP(L)+MAX(WSETA(L,1,NS),0.)
            BB11=DELTI*DZC(1)*HP(L)*SED(L,1,NS)
            BB22=DELTI*DZC(2)*HP(L)*SED(L,2,NS)
            DETI=1./(AA11*AA22-AA12*AA21)
            SED(L,1,NS)=DETI*( BB11*AA22-BB22*AA12 )
            SED(L,2,NS)=DETI*( AA11*BB22-AA21*BB11 )
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=2
C
           DO L=2,LA
            SEDF(L,1,NS)=DELTI*DZC(2)*HP(L)*(SED(L,2,NS)-SEDS(L,2,NS))
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.EQ.3) THEN
        DO NS=1,NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISEDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.EQ.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),0.0,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.EQ.2) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
             DO K=1,KS
             DO L=2,LA
              LN=LNC(L)
              SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.GE.3) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
             DO K=1,KS
             DO L=2,LA
              LN=LNC(L)
              STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
             END DO
            END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loops
C
           K=3
           DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
           END DO
C
           K=2
           DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)-SEDF(L,K,NS)*WVEL
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
           END DO
C
           DO L=2,LA
            PROBDEP=0.
            WESE=0.
            TAUBC=QQ(L,0)/CTURB2
            UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
            VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
            CURANG=ATAN2(VTMP,UTMP)
            TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $           +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
            TAUB2=MAX(TAUB2,0.)
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
            IF(TAUB.GT.TAURB(L,KBT(L))) THEN
c **  mass erosion
               WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
               WESE=MIN(WESE,WESEMX)
              ELSE
               IF(TAUB.GT.TAURS(L,KBT(L))) THEN
c **  surface erosion
                  WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
                  WESE=MIN(WESE,WESEMX)
                  TAUE=(TAUB-TAURS(L,KBT(L)))/TAURS(L,KBT(L))
                  TAUE=MAX(TAUE,0.0)
                  WESE=WESE*(TAUE**TEXP(NS))
                 ELSE
c **  no erosion
                  WESE=0.0
               END IF
            END IF
c **  set probability of deposition
            IF(TAUB.LT.TAUD(NS)) PROBDEP=(TAUD(NS)-TAUB)/TAUD(NS)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL
            SED(L,1,NS)=CRIGHT/CLEFT
            SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE
            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)
c            IF(SEDBTMP.LT.0.0) THEN
c              SEDF(L,0,NS)=0.
c              SEDBTMP=SEDB1(L,KBT(L),NS)
c              SED(L,1,NS)=SEDS(L,1,NS)-SEDF(L,1,NS)*WVEL
c            END IF
c            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
c     $                        +S2TL*SEDB1(L,KBT(L),NS)
c            SEDB(L,KBT(L),NS)=SEDBTMP
            IF(SEDBTMP.LT.0.0) THEN
              SEDF(L,0,NS)=DELTI*SEDB1(L,KBT(L),NS)
              SEDBTMP=0.0
              SED(L,1,NS)=SEDS(L,1,NS)+(SEDF(L,0,NS)-SEDF(L,1,NS))*WVEL
            END IF
            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
     $                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     $                  +SEDVDRD*MIN(SEDF(L,0,NS),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.EQ.3
C
           IF(ISTOPT(6).EQ.1) THEN
C
           DO K=1,2
            DO L=2,LA
             CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
             GRADSED=(SED(L,K+1,NS)-SED(L,K,NS))/(DZC(K+1)+DZC(K))
             SEDAVG=0.5*(SED(L,K+1,NS)-SED(L,K,NS)+1.E-16)
             WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            END DO
           END DO
C
C     TVAR1S=LOWER DIAGONAL
           DO L=2,LA
           TVAR1S(L,1)=0
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
            END DO
           END DO
C     TVAR1N=UPPER DIAGONAL
           DO L=2,LA
           TVAR1N(L,KC)=0
           END DO
           DO K=1,KS
            DO L=2,LA
             TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1W=MAIN DIAGONAL
           DO L=2,LA
           TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
           TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
           END DO
           DO K=2,KS
            DO L=2,LA
             TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     $                -MIN(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1E=RIGHT HAND SIDE
           DO K=1,KC
            DO L=2,LA
             TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SED(L,K,NS)
            END DO
           END DO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
           DO L=2,LA
            TVAR3S(L)=TVAR1W(L,1)
           END DO
           DO L=2,LA
            TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
             TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
             TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     $                   TVAR3S(L)
            END DO
           END DO
           DO K=KS,1,-1
            DO L=2,LA
             TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            END DO
           END DO
           DO K=1,KC
            DO L=2,LA
             SED(L,K,NS)=TVAR2N(L,K)
            END DO
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=3
C
           DO L=2,LA
            SEDF(L,2,NS)=DELTI*DZC(3)*HP(L)*(SED(L,3,NS)-SEDS(L,3,NS))
           END DO
           DO L=2,LA
            SEDF(L,1,NS)=DELTI*DZC(2)*HP(L)*(SED(L,2,NS)-SEDS(L,2,NS))
     $                 +SEDF(L,2,NS)
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  COHESIVE SEDIMENT, KC.GT.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(ISTRAN(6).GE.1.AND.KC.GT.3) THEN
        DO NS=1,NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISEDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
c             WSETA(L,K,NS)=WSEDO(NS)
              WSETA(L,K,NS)=xLwsedo(L)  ! Ji, 3/30/01
             END DO
             END DO
           END IF
C
           IF(ISEDVW.EQ.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),0.0,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.EQ.2) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*DZIC(K+1)*HPI(L)*SQRT(STRESS)/VKC
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
             DO K=1,KS
             DO L=2,LA
              LN=LNC(L)
              SHEAR=HPI(L)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),SHEAR,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.GE.3) THEN
             K=0
             DO L=2,LA
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
             DO K=1,KS
             DO L=2,LA
              LN=LNC(L)
              STRESS=AV(L,K)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
              WSETA(L,K,NS)=CSEDSET(SED(L,K+1,NS),STRESS,ISEDVW)
             END DO
             END DO
            END IF
C
           IF(ISEDVW.GE.3) THEN  !Ji, hardwired by John for station c (39,29), 1/21/00
             K=0
              L=857
               TAUBC=QQ(L,0)/CTURB2
               UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
               VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
               CURANG=ATAN2(VTMP,UTMP)
               TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $              +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
               TAUB2=MAX(TAUB2,0.)
               STRESSS(K)=SQRT(TAUB2)
             DO K=1,KS
              L=857
              LN=LNC(L)
              STRESSS(K)=AV(L,K)*SQRT( DZIGSD4(K) )
     $              *SQRT( (U(L+1,K+1)-U(L+1,K)+U(L,K+1)-U(L,K))**2
     $                    +(V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
             END DO
            END IF
c
c      WRITE(11,6111)TIME,(WSETA(857,K,1),K=0,KS)
c      WRITE(41,6111)TIME,(STRESSS(K),K=0,KS)
 6111 FORMAT(F10.2,10E12.4)
C
C----------------------------------------------------------------------C
C
c **  horizontal loops
c
           K=KC
           DO L=2,LA
            SEDF(L,K,NS)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
           END DO
C
           DO K=KS,2,-1
           DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)-SEDF(L,K,NS)*WVEL
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
           END DO
           END DO
C
           DO L=2,LA                         ! # 11
            PROBDEP=0.
            WESE=0.
            TAUBC=QQ(L,0)/CTURB2
            UTMP=0.5*STCUV(L)*(U(L+1,1)+U(L,1))+1.E-12
            VTMP=0.5*STCUV(L)*(V(LNC(L),1)+V(L,1))
            CURANG=ATAN2(VTMP,UTMP)
            TAUB2=TAUBC*TAUBC+0.5*(QQWV2(L)*QQWV2(L))
     $           +FOURDPI*TAUBC*QQWV2(L)*COS(CURANG-WACCWE(L))
            TAUB2=MAX(TAUB2,0.)
            TAUB=SQRT(TAUB2)
c
            if(iswan.ge.1) then
c           if(L.eq.135) write(309,351) N,time,taub,tau(L)  ! Sta A
351   format(i5,999f14.5)
            taub=tau(L)            ! use wavetau's shear stress, Ji, 3/7/00
            endif
c
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SEDB(L,KBT(L),NS)  ! ensure not too much mass erosion at this time step
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX  !Ji, why, hardwired?, 2/26/00
            WESEMX=MAX(WESEMX,0.)             ! Ji, mass erosion is turned off by letting TAURB=1.0E6
            IF(TAUB.GT.TAURB(L,KBT(L))) THEN  ! Taurb = bulk erossion critical shear stress
c **  mass erosion
               WESE=CTMPDRY(L)*WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)   ! not really used
               WESE=MIN(WESE,WESEMX)
              ELSE
               TAUE=0.
               IF(TAUB.GT.TAURS(L,KBT(L))) THEN  ! Ji, 3/31/01
c               IF(TAUB.GT.XLTAUR(L,KBT(L))) THEN
c **  surface erosion
                  WESE=CTMPDRY(L)*WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)
                  WESE=MIN(WESE,WESEMX)
                  TAUE=(TAUB-TAURS(L,KBT(L)))/TAURS(L,KBT(L))
                  TAUE=MAX(TAUE,0.0)
                  WESE=WESE*(TAUE**TEXP(NS))
                 ELSE
c **  no erosion
                  WESE=0.0
               END IF
            END IF
c
c         wese=0.0 ! hardwired  for testing, Ji, 10/25/00
c
c **  set probability of deposition
c           IF(TAUB.LT.TAUD(NS)) PROBDEP=(TAUD(NS)-TAUB)/TAUD(NS) !@, Ji, sediment property, 2/7/00
            IF(TAUB.LT.XLTAUD(L)) PROBDEP=(XLTAUD(L)-TAUB)/XLTAUD(L)
            IF(SED(L,1,NS).GT.SEDMDGM) PROBDEP=1.
c
c         probdep=0.0 ! hardwired  for testing, Ji, 10/25/00
c
            WSETMP=PROBDEP*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
c           CRIGHT=MAX(SED(L,1,NS),0.)+(WESE-SEDF(L,1,NS))*WVEL  !it might cause conservation problem when there
c                                                                      ! is not enough sed on the bed
            xx1=wese*wvel                          ! max possible erosion at this time step, in g/m**3
c           xx2=sedb1(L,KBT(L),ns)*HPI(L)*DZIC(1)  ! available sed on bed for erosion, in g/m**3
            xx2=0
            do kx=1,KBT(L)
            xx2=xx2+sedb1(L,Kx,ns)
            enddo
            xx2=xx2*HPI(L)*DZIC(1)  ! available sed on bed for erosion, in g/m**3
            xx3=min(xx1,xx2)                       ! ensure sed conversation, no over erosion
            CRIGHT=MAX(SED(L,1,NS),0.)+xx3-SEDF(L,1,NS)*WVEL
c           SED(L,1,NS)=CRIGHT/CLEFT                ! John's way
c           SEDF(L,0,NS)=-WSETMP*SED(L,1,NS)+WESE   ! John's way
            xxSED1=CRIGHT/CLEFT                     ! = SED(L,1,NS) at (N+1) step, my way, see notes, 10/20/00
            SEDF(L,0,NS)=SEDF(L,1,NS)+(xxsed1-SED(L,1,NS))/WVEL          ! my way
            SED(L,1,NS)=xxSED1                                           ! my way, sed at (N+1)
c
            SEDBTMP=SEDB1(L,KBT(L),NS)-DELT*SEDF(L,0,NS)  ! checked by analytical solution
            IF(kbt(L).eq.1)  sedbtmp=max(sedbtmp,0.0) ! round off error could cause it, Ji
c
            SEDB1(L,KBT(L),NS)=S3TL*SEDB(L,KBT(L),NS)
     $                        +S2TL*SEDB1(L,KBT(L),NS)
            SEDB(L,KBT(L),NS)=SEDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SEDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.)
     $                  +SEDVDRD*MIN(SEDF(L,0,NS),0.) )
C
c     IF(L.EQ.1169.AND.NS.EQ.1) THEN   ! (20,35), Station E = 15
c       WRITE(21,6111)TIME,TAUBC,QQWV2(L),TAUB,WESE,WSETMP,SED(L,1,1) ! Ji, 10/19/00
c       WRITE(21,6119) KBT(L),TIME,sed(L,1,NS),
c    *         sedb(L,KBT(L), NS),sedf(L,0,NS),TAUB,TAURS(L,KBT(L))   ! ssedtox.dia2
c6119  format(i4,f14.4,999e12.4)
c     END IF
C
           END DO                        ! #11
C
c      WRITE(31,6111)TIME,(SEDF(857,K,1),K=0,KS)
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC.GT.3
C
           IF(ISTOPT(6).EQ.1) THEN
C
           DO K=1,KS
            DO L=2,LA
             CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
             GRADSED=(SED(L,K+1,NS)-SED(L,K,NS))/(DZC(K+1)+DZC(K))
             SEDAVG=0.5*(SED(L,K+1,NS)+SED(L,K,NS)+1.E-16)
             WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            END DO
           END DO
C
C     TVAR1S=LOWER DIAGONAL
           DO L=2,LA
           TVAR1S(L,1)=0
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
            END DO
           END DO
C     TVAR1N=UPPER DIAGONAL
           DO L=2,LA
           TVAR1N(L,KC)=0
           END DO
           DO K=1,KS
            DO L=2,LA
             TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1W=MAIN DIAGONAL
           DO L=2,LA
           TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
           TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
           END DO
           DO K=2,KS
            DO L=2,LA
             TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     $                -MIN(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1E=RIGHT HAND SIDE
           DO K=1,KC
            DO L=2,LA
             TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SED(L,K,NS)
            END DO
           END DO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
           DO L=2,LA
            TVAR3S(L)=TVAR1W(L,1)
           END DO
           DO L=2,LA
            TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
             TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
             TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     $                   TVAR3S(L)
            END DO
           END DO
           DO K=KS,1,-1
            DO L=2,LA
             TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            END DO
           END DO
           DO K=1,KC
            DO L=2,LA
             SED(L,K,NS)=TVAR2N(L,K)
            END DO
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC.GT.3
C
           DO L=2,LA
            SEDF(L,KS,NS)=DELTI*DZC(KC)*HP(L)*
     $                   (SED(L,KC,NS)-SEDS(L,KC,NS))
           END DO
C
           DO K=KS-1,1,-1
            DO L=2,LA
             SEDF(L,K,NS)=DELTI*DZC(K+1)*HP(L)*
     $                   (SED(L,K+1,NS)-SEDS(L,K+1,NS))+SEDF(L,K+1,NS)
            END DO
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  UNCOMMENT THE FOLLOWING FOR COHESIVE SEDIMENT DIAGNOSTICS
cdiag           SEDMX=-1.E+12
cdiag           SEDBMX=-1.E+12
cdiag           SEDFMX=-1.E+12
cdiag           SEDMN=1.E+12
cdiag           SEDBMN=1.E+12
cdiag           SEDFMN=1.E+12
cdiag           DO K=1,KC
cdiag            DO L=2,LA
cdiag             IF(SED(L,K,NS).GT.SEDMX) THEN
cdiag               LMX=L
cdiag               KMX=K
cdiag               SEDMX=SED(L,K,NS)
cdiag             END IF
cdiag             IF(SED(L,K,NS).LT.SEDMN) THEN
cdiag               LMN=L
cdiag               KMN=K
cdiag               SEDMN=SED(L,K,NS)
cdiag             END IF
cdiag            END DO
cdiag           END DO
cdiag           DO L=2,LA
cdiag            IF(SEDB(L,NS).GT.SEDBMX) THEN
cdiag              LBMX=L
cdiag              SEDBMX=SEDB(L,NS)
cdiag            END IF
cdiag            IF(SEDB(L,NS).LT.SEDBMN) THEN
cdiag              LBMN=L
cdiag              SEDBMN=SEDB(L,NS)
cdiag            END IF
cdiag            IF(SEDF(L,0,NS).GT.SEDFMX) THEN
cdiag              LFMX=L
cdiag              SEDFMX=SEDF(L,0,NS)
cdiag            END IF
cdiag            IF(SEDF(L,0,NS).LT.SEDFMN) THEN
cdiag              LFMN=L
cdiag              SEDFMN=SEDF(L,0,NS)
cdiag            END IF
cdiag           END DO
C
cdiag      L=LMX
cdiag      K=KMX
cdiag      WRITE(1,101)N,NS,IL(L),JL(L),K,SED(L,K,NS),SEDS(L,K,NS)
cdiag      L=LMN
cdiag      K=KMN
cdiag      WRITE(1,102)N,NS,IL(L),JL(L),K,SED(L,K,NS),SEDS(L,K,NS)
cdiag      L=LBMX
cdiag      WRITE(1,103)N,NS,IL(L),JL(L),SEDB(L,NS),SEDBS(L,NS),SEDF(L,O,NS)
cdiag      L=LBMN
cdiag      WRITE(1,104)N,NS,IL(L),JL(L),SEDB(L,NS),SEDBS(L,NS),SEDF(L,O,NS)
cdiag      L=LFMX
cdiag      WRITE(1,105)N,NS,IL(L),JL(L),SEDF(L,0,NS)
cdiag      L=LFMN
cdiag      WRITE(1,106)N,NS,IL(L),JL(L),SEDF(L,0,NS)
C
cdiag  101 FORMAT(' N,NS,I,J,K,SEDMX,SEDSMX = ',5I5,4E13.4)
cdiag  102 FORMAT(' N,NS,I,J,K,SEDMN,SEDSMN = ',5I5,4E13.4)
cdiag  103 FORMAT(' N,NS,I,J,SEDBMX,SEDBSMX = ',4I5,4E13.4)
cdiag  104 FORMAT(' N,NS,I,J,SEDBMN,SEDBSMN = ',4I5,4E13.4)
cdiag  105 FORMAT(' N,NS,I,J,SEDFMX,SEDFSMX = ',4I5,4E13.4)
cdiag  106 FORMAT(' N,NS,I,J,SEDFMN,SEDFSMN = ',4I5,4E13.4)
C
      DO NS=1,NSED
      DO K=1,KC
       DO L=2,LA
        IF(SED(L,K,NS).LT.0.) THEN
!          WRITE(1,107)N,NS,IL(L),JL(L),K,SED(L,K,NS)
        END IF
       END DO
      END DO
      END DO
C
      DO NS=1,NSED
      DO L=2,LA
       IF(SEDB(L,KBT(L),NS).LT.0.) THEN
!         WRITE(1,108)N,NS,IL(L),JL(L),SEDB(L,KBT(L),NS),SEDF(L,0,NS)
       END IF
      END DO
      END DO
C
  107 FORMAT(' N,NS,I,J,K,NEGSED = ',5I8,4E13.4)
  108 FORMAT(' N,NS,I,J,NEGSEDB,SEDF = ',4I8,4E13.4)
C
C**********************************************************************C
C
C **  SET ARMORING PARAMETERS
C
      IF(ISTRAN(7).GE.1) THEN
C
      IF(ISNDAL.EQ.0) THEN
        DO K=1,KB
        DO L=2,LA
          SIGPHI(L,K)=0.
        END DO
        END DO
      END IF
C
      IF(ISNDAL.EQ.1) THEN
C
C **  SET MAXIMUM DIAMETER AND PHI SIZE
C
      SNDDMX=0.
      DO NX=1,NSND
       NS=NSED+NX
       SNDDMX=MAX(SNDDMX,SEDDIA(NS))
       SEDPHI(NS)=-LOG(1000.*SEDDIA(NS))/LOG(2.)
      END DO
C
C **  SET MEAN PHI (STORE IN SNDBT AND STORE NONCOHESIVE
C **  VOLUME FRACTION IS SEDDIA50)
C
      DO K=1,KB
      DO L=2,LA
        SNDBT(L,K)=0.
        SEDDIA50(L,K)=0.
        SIGPHI(L,K)=0.
      END DO
      END DO
C
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        SNDBT(L,K)=SNDBT(L,K)+SEDPHI(NS)*VFRBED(L,K,NS)
        SEDDIA50(L,K)=SEDDIA50(L,K)+VFRBED(L,K,NS)
       END DO
       END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       IF(SEDDIA50(L,K).LE.0.) SEDDIA50(L,K)=1
       SNDBT(L,K)=SNDBT(L,K)/SEDDIA50(L,K)
      END DO
      END DO
C
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        SIGPHI(L,K)=SIGPHI(L,K)+( ( SEDPHI(NS)-SNDBT(L,K) )**2 )*
     &                          VFRBED(L,K,NS)/SEDDIA50(L,K)
       END DO
       END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
       SIGPHI(L,K)=SQRT(SIGPHI(L,K))
      END DO
      END DO
C
C **  SET MEAN D50
C
      DO K=1,KB
      DO L=2,LA
        SEDDIA50(L,K)=0.
        SNDBT(L,K)=0.
      END DO
      END DO
C
      DO NX=1,NSND
       NS=NSED+NX
       DO K=1,KB
       DO L=2,LA
        SEDDIA50(L,K)=SEDDIA50(L,K)+SEDDIA(NS)*SNDB(L,K,NX)
        SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
       END DO
       END DO
      END DO
C
      DO K=1,KB
      DO L=2,LA
        SEDDIA50(L,K)=0.
        IF(SNDBT(L,K).GT.0.) SEDDIA50(L,K)=SEDDIA50(L,K)/SNDBT(L,K)
      END DO
      END DO
C
      END IF
C
      END IF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.1) THEN
        DO NX=1,NSND
         NS=NX+NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
           DIASED=SEDDIA(NS)
           DIASED3=3.*DIASED
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           K=0
C
           IF(ISNDVW.EQ.0) THEN
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
           END IF
C
           IF(ISNDVW.GE.1) THEN
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)*
     $                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
             END DO
           END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loop
c
           DO L=2,LA
            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
            TAUB=MAX(QQ(L,0),QQMIN)/CTURB2
            TAUB2=TAUB*TAUB+0.5*(QQWV2(L)*QQWV2(L))
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),Nx)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
c **  set rouse parameter and equilibruim concentration
            USTAR=SQRT(TAUB)
            ROUSE=WSETA(L,0,NS)/(VKC*USTAR)
            ZEQ=DIASED3*HPI(L)
            ZEQMIN=0.5*DZC(1)
            ZEQ=MIN(ZEQ,ZEQMIN)
            ZEQI=1./ZEQ
            WSFAC=2.*(1.+ROUSE)/(2.+ROUSE*(1.-ZEQ))
            D50TMP=SEDDIA50(L,KBT(L))
            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
     $                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            IF(ROUSE.LT.0.999.OR.ROUSE.GT.1.001) THEN
               TOP=(ZEQ**(ROUSE-1.))-1.
               BOT=(1.-ROUSE)*(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
              ELSE
               TOP=LOG(ZEQI)
               BOT=(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
            END IF
c **  set resuspension flux
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,K,NS)*SNDEQ
c **  set deposition velocity
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
c            IF(SNDBTMP.LT.0.0) THEN
c              SNDF(L,0,NX)=0.
c              SNDBTMP=SNDB1(L,KBT(L),NX)
c              SND(L,1,NX)=SNDS(L,1,NX)-SNDF(L,1,NX)*WVEL
c            END IF
c            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
c     $                        +S2TL*SNDB1(L,KBT(L),NX)
c            SNDB(L,KBT(L),NX)=SNDBTMP
            IF(SNDBTMP.LT.0.0) THEN
              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            END IF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     $                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SNDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     $                  +SNDVDRD*MIN(SNDF(L,0,NX),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.2) THEN
        DO NX=1,NSND
         NS=NX+NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
           DIASED=SEDDIA(NS)
           DIASED3=3.*DIASED
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISNDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
             END DO
           END IF
C
           IF(ISNDVW.GE.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)*
     $                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
             END DO
             END DO
           END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loops
c
           K=2
           DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
           END DO
C
           DO L=2,LA
            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
            TAUB=MAX(QQ(L,0),QQMIN)/CTURB2
            TAUB2=TAUB*TAUB+0.5*(QQWV2(L)*QQWV2(L))
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),Nx)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
c **  set rouse parameter and equilibruim concentration
            USTAR=SQRT(TAUB)
            ROUSE=WSETA(L,0,NS)/(VKC*USTAR)
            ZEQ=DIASED3*HPI(L)*DZIC(1)
            ZEQMIN=0.5*DZC(1)
            ZEQ=MIN(ZEQ,ZEQMIN)
            ZEQI=1./ZEQ
            WSFAC=1
            D50TMP=SEDDIA50(L,KBT(L))
            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
     $                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            IF(ROUSE.LT.0.999.OR.ROUSE.GT.1.001) THEN
               TOP=(ZEQ**(ROUSE-1.))-1.
               BOT=(1.-ROUSE)*(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
              ELSE
               TOP=LOG(ZEQI)
               BOT=(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
            END IF
c **  set resuspension flux
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,K,NS)*SNDEQ
c **  set deposition velocity
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
c            IF(SNDBTMP.LT.0.0) THEN
c              SNDF(L,0,NX)=0.
c              SNDBTMP=SNDB1(L,KBT(L),NX)
c              SND(L,1,NX)=SNDS(L,1,NX)-SNDF(L,1,NX)*WVEL
c            END IF
c            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
c     $                        +S2TL*SNDB1(L,KBT(L),NX)
c            SNDB(L,KBT(L),NX)=SNDBTMP
            IF(SNDBTMP.LT.0.0) THEN
              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            END IF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     $                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SNDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     $                  +SNDVDRD*MIN(SNDF(L,0,NX),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.EQ.2
C
           IF(ISTOPT(7).EQ.1) THEN
C
           DO L=2,LA
            CRNUM=1.+DELT*WSETA(L,1,NS)*HPI(L)*DZIC(2)
            GRADSED=(SND(L,2,NX)-SND(L,1,NX))/(DZC(2)+DZC(1))
            SEDAVG=0.5*(SND(L,2,NX)+SND(L,1,NX)+1.E-16)
            WSETA(L,1,NS)=-CRNUM*DZC(2)*WSETA(L,1,NS)*GRADSED/SEDAVG
           END DO
C
           DO L=2,LA
            AA11=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
            AA12=-MAX(WSETA(L,1,NS),0.)
            AA21=MIN(WSETA(L,1,NS),0.)
            AA22=DELTI*DZC(2)*HP(L)+MAX(WSETA(L,1,NS),0.)
            BB11=DELTI*DZC(1)*HP(L)*SND(L,1,NX)
            BB22=DELTI*DZC(2)*HP(L)*SND(L,2,NX)
            DETI=1./(AA11*AA22-AA12*AA21)
            SND(L,1,NX)=DETI*( BB11*AA22-BB22*AA12 )
            SND(L,2,NX)=DETI*( AA11*BB22-AA21*BB11 )
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=2
C
           DO L=2,LA
            SNDF(L,1,NX)=DELTI*DZC(2)*HP(L)*(SND(L,2,NX)-SNDS(L,2,NX))
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.EQ.3) THEN
        DO NX=1,NSND
         NS=NX+NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
           DIASED=SEDDIA(NS)
           DIASED3=3.*DIASED
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISNDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
             END DO
           END IF
C
           IF(ISNDVW.GE.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)*
     $                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
             END DO
             END DO
           END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loops
c
           K=3
           DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
           END DO
C
           K=2
           DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
           END DO
c
           DO L=2,LA
            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
            TAUB=MAX(QQ(L,0),QQMIN)/CTURB2
            TAUB2=TAUB*TAUB+0.5*(QQWV2(L)*QQWV2(L))
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),Nx)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
c **  set rouse parameter and equilibruim concentration
            USTAR=SQRT(TAUB)
            ROUSE=WSETA(L,0,NS)/(VKC*USTAR)
            ZEQ=DIASED3*HPI(L)*DZIC(1)
            ZEQMIN=0.5*DZC(1)
            ZEQ=MIN(ZEQ,ZEQMIN)
            ZEQI=1./ZEQ
            WSFAC=1.
            D50TMP=SEDDIA50(L,KBT(L))
            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
     $                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            IF(ROUSE.LT.0.999.OR.ROUSE.GT.1.001) THEN
               TOP=(ZEQ**(ROUSE-1.))-1.
               BOT=(1.-ROUSE)*(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
              ELSE
               TOP=LOG(ZEQI)
               BOT=(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
            END IF
c **  set resuspension flux
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,K,NS)*SNDEQ
c **  set deposition velocity
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
c            IF(SNDBTMP.LT.0.0) THEN
c              SNDF(L,0,NX)=0.
c              SNDBTMP=SNDB1(L,KBT(L),NX)
c              SND(L,1,NX)=SNDS(L,1,NX)-SNDF(L,1,NX)*WVEL
c            END IF
c            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
c     $                        +S2TL*SNDB1(L,KBT(L),NX)
c            SNDB(L,KBT(L),NX)=SNDBTMP
            IF(SNDBTMP.LT.0.0) THEN
              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            END IF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     $                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SNDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     $                  +SNDVDRD*MIN(SNDF(L,0,NX),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.EQ.3
C
           IF(ISTOPT(7).EQ.1) THEN
C
           DO K=1,2
            DO L=2,LA
             CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
             GRADSED=(SND(L,K+1,NX)-SND(L,K,NX))/(DZC(K+1)+DZC(K))
             SEDAVG=0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)
             WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            END DO
           END DO
C
C     TVAR1S=LOWER DIAGONAL
           DO L=2,LA
           TVAR1S(L,1)=0
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
            END DO
           END DO
C     TVAR1N=UPPER DIAGONAL
           DO L=2,LA
           TVAR1N(L,KC)=0
           END DO
           DO K=1,KS
            DO L=2,LA
             TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1W=MAIN DIAGONAL
           DO L=2,LA
           TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
           TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
           END DO
           DO K=2,KS
            DO L=2,LA
             TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     $                -MIN(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1E=RIGHT HAND SIDE
           DO K=1,KC
            DO L=2,LA
             TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SND(L,K,NX)
            END DO
           END DO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
           DO L=2,LA
            TVAR3S(L)=TVAR1W(L,1)
           END DO
           DO L=2,LA
            TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
             TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
             TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     $                   TVAR3S(L)
            END DO
           END DO
           DO K=KS,1,-1
            DO L=2,LA
             TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            END DO
           END DO
           DO K=1,KC
            DO L=2,LA
             SND(L,K,NX)=TVAR2N(L,K)
            END DO
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC=3
C
           DO L=2,LA
            SNDF(L,2,NX)=DELTI*DZC(3)*HP(L)*(SND(L,3,NX)-SNDS(L,3,NX))
           END DO
           DO L=2,LA
            SNDF(L,1,NX)=DELTI*DZC(2)*HP(L)*(SND(L,2,NX)-SNDS(L,2,NX))
     $                 +SNDF(L,2,NX)
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  NONCOHESIVE SEDIMENT, KC.GT.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(ISTRAN(7).GE.1.AND.KC.GT.3) THEN
        DO NX=1,NSND
         NS=NX+NSED
         DSEDGMM=1./(1.E6*SSG(NS))
         IF(WSEDO(NS).GE.1.E-15) THEN
           DIASED=SEDDIA(NS)
           DIASED3=3.*DIASED
C
C----------------------------------------------------------------------C
C
C **  set settling velocities
C
           IF(ISNDVW.EQ.0) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)
             END DO
             END DO
           END IF
C
           IF(ISNDVW.GE.1) THEN
             DO K=0,KS
             DO L=2,LA
              WSETA(L,K,NS)=WSEDO(NS)*
     $                      CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)
             END DO
             END DO
           END IF
C
C----------------------------------------------------------------------C
C
c **  horizontal loop
c
           K=KC
           DO L=2,LA
            SNDF(L,K,NX)=0.
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
           END DO
C
           DO K=KS,2,-1
           DO L=2,LA
            WVEL=DELT*HPI(L)*DZIC(K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL
            SND(L,K,NX)=CRIGHT/CLEFT
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)
           END DO
           END DO
c
           DO L=2,LA
            SNDF(L,1,NX)=0.
            PROBDEP=0.
            WESE=0.
            TAUB=MAX(QQ(L,0),QQMIN)/CTURB2
            TAUB2=TAUB*TAUB+0.5*(QQWV2(L)*QQWV2(L))
            TAUB=SQRT(TAUB2)
c **  set maximum erosion rate
            WESEMX=0.5*DELTI*CTMPDRY(L)*SNDB(L,KBT(L),Nx)
            IF(KBT(L).NE.1) WESEMX=3.*WESEMX
            WESEMX=MAX(WESEMX,0.)
c **  set rouse parameter and equilibruim concentration
            USTAR=SQRT(TAUB)
            ROUSE=WSETA(L,0,NS)/(VKC*USTAR)
            ZEQ=DIASED3*HPI(L)*DZIC(1)
            ZEQMIN=0.5*DZC(1)
            ZEQ=MIN(ZEQ,ZEQMIN)
            ZEQI=1./ZEQ
            WSFAC=1.
            D50TMP=SEDDIA50(L,KBT(L))
            IF(ISNDAL.EQ.1) D50TMP=DIASED
            SIGP=SIGPHI(L,KBT(L))
            SNDEQB=CSNDEQC(DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUB,
     $                     D50TMP,SIGP,SNDDMX,ISNDEQ(NS))
            IF(ROUSE.LT.0.999.OR.ROUSE.GT.1.001) THEN
               TOP=(ZEQ**(ROUSE-1.))-1.
               BOT=(1.-ROUSE)*(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
              ELSE
               TOP=LOG(ZEQI)
               BOT=(ZEQI-1.)
               SNDEQ=SNDEQB*TOP/BOT
               SNDEQ=VFRBED(L,KBT(L),NS)*MAX(SNDEQ,0.)
            END IF
c **  set resuspension flux
            WESE=WSFAC*CTMPDRY(L)*WSETA(L,K,NS)*SNDEQ
c **  set deposition velocity
            WSETMP=WSFAC*WSETA(L,0,NS)
            WVEL=DELT*HPI(L)*DZIC(1)
            CLEFT=1.+WSETMP*WVEL
            CRIGHT=MAX(SND(L,1,NX),0.)+(WESE-SNDF(L,1,NX))*WVEL
            SND(L,1,NX)=CRIGHT/CLEFT
            SNDF(L,0,NX)=-WSETMP*SND(L,1,NX)+WESE
            SNDBTMP=SNDB1(L,KBT(L),NX)-DELT*SNDF(L,0,NX)
c            IF(SNDBTMP.LT.0.0) THEN
c              SNDF(L,0,NX)=0.
c              SNDBTMP=SNDB1(L,KBT(L),NX)
c              SND(L,1,NX)=SNDS(L,1,NX)-SNDF(L,1,NX)*WVEL
c            END IF
c            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
c     $                        +S2TL*SNDB1(L,KBT(L),NX)
c            SNDB(L,KBT(L),NX)=SNDBTMP
            IF(SNDBTMP.LT.0.0) THEN
              SNDF(L,0,NX)=DELTI*SNDB1(L,KBT(L),NX)
              SNDBTMP=0.0
              SND(L,1,NX)=SNDS(L,1,NX)+(SNDF(L,0,NX)-SNDF(L,1,NX))*WVEL
            END IF
            SNDB1(L,KBT(L),NX)=S3TL*SNDB(L,KBT(L),NX)
     $                        +S2TL*SNDB1(L,KBT(L),NX)
            SNDB(L,KBT(L),NX)=SNDBTMP
            QSBDTOP(L)=QSBDTOP(L)-DSEDGMM*SNDF(L,0,NS)
            QWBDTOP(L)=QWBDTOP(L)-DSEDGMM*
     $                 ( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.)
     $                  +SNDVDRD*MIN(SNDF(L,0,NX),0.) )
           END DO
C
C----------------------------------------------------------------------C
C
C **  ANTI-DIFFUSION OF NONCOHESIVE SEDIMENT  KC.GT.3
C
           IF(ISTOPT(7).EQ.1) THEN
C
           DO K=1,KS
            DO L=2,LA
             CRNUM=1.+DELT*WSETA(L,K,NS)*HPI(L)*DZIC(K+1)
             GRADSED=(SND(L,K+1,NX)-SND(L,K,NX))/(DZC(K+1)+DZC(K))
             SEDAVG=0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)
             WSETA(L,K,NS)=-CRNUM*DZC(K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            END DO
           END DO
C
C     TVAR1S=LOWER DIAGONAL
           DO L=2,LA
           TVAR1S(L,1)=0
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
            END DO
           END DO
C     TVAR1N=UPPER DIAGONAL
           DO L=2,LA
           TVAR1N(L,KC)=0
           END DO
           DO K=1,KS
            DO L=2,LA
             TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1W=MAIN DIAGONAL
           DO L=2,LA
           TVAR1W(L,1)=DELTI*DZC(1)*HP(L)-MIN(WSETA(L,1,NS),0.)
           TVAR1W(L,KC)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)
           END DO
           DO K=2,KS
            DO L=2,LA
             TVAR1W(L,K)=DELTI*DZC(KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)
     $                -MIN(WSETA(L,K,NS),0.)
            END DO
           END DO
C     TVAR1E=RIGHT HAND SIDE
           DO K=1,KC
            DO L=2,LA
             TVAR1E(L,K)=DELTI*DZC(KC)*HP(L)*SND(L,K,NX)
            END DO
           END DO
C
C     TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
           DO L=2,LA
            TVAR3S(L)=TVAR1W(L,1)
           END DO
           DO L=2,LA
            TVAR2N(L,1)=TVAR1E(L,1)/TVAR3S(L)
           END DO
           DO K=2,KC
            DO L=2,LA
             TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
             TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
             TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/
     $                   TVAR3S(L)
            END DO
           END DO
           DO K=KS,1,-1
            DO L=2,LA
             TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            END DO
           END DO
           DO K=1,KC
            DO L=2,LA
             SND(L,K,NX)=TVAR2N(L,K)
            END DO
           END DO
C
           END IF
C
C----------------------------------------------------------------------C
C
C **  FINAL FLUX KC.GE.3
C
           DO L=2,LA
            SNDF(L,KS,NX)=DELTI*DZC(KC)*HP(L)*
     $                   (SND(L,KC,NX)-SNDS(L,KC,NX))
           END DO
C
           DO K=KS-1,1,-1
            DO L=2,LA
             SNDF(L,K,NX)=DELTI*DZC(K+1)*HP(L)*
     $                   (SND(L,K+1,NX)-SNDS(L,K+1,NX))+SNDF(L,K+1,NX)
            END DO
           END DO
C
C----------------------------------------------------------------------C
C
         END IF
        END DO
      END IF
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN WATER COLUMN
C
C **  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
C **  TOXPARW(NS,NT) = PARTITION COEFFICIENT IN WATER COLUMN
C **  TOXPFTW(L,K,NT) = TOTAL PARTICULATE FRACTION IN WATER COLUMN
C                       USED AS TEMPORARY VARIBLE IN THIS AND
C                       FOLLOWING CODE BLOCK
C
      IF(TOXPARW(1,1).LT.1.E-6) RETURN
      IF(ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KC
         DO L=2,LA
          TOXPFW(L,K,NS,NT)=0.
         END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO K=1,KC
           DO L=2,LA
            TOXPFW(L,K,NS,NT)=SED(L,K,NS)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO K=1,KC
           DO L=2,LA
            TMPVAL=1.
            IF(SED(L,K,NS).GT.0.) TMPVAL=SED(L,K,NS)**TMPEXP
            TOXPFW(L,K,NS,NT)=TMPVAL*SED(L,K,NS)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO K=1,KC
           DO L=2,LA
            TOXPFW(L,K,NS,NT)=SND(L,K,NX)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO K=1,KC
           DO L=2,LA
            TMPVAL=1.
            IF(SND(L,K,NX).GT.0.) TMPVAL=SND(L,K,NX)**TMPEXP
            TOXPFW(L,K,NS,NT)=TMPVAL*SND(L,K,NX)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO K=1,KC
        DO L=2,LA
         TOXPFTW(L,K,NT)=0.
        END DO
       END DO
       DO NS=1,NSED+NSND
        DO K=1,KC
         DO L=2,LA
          TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)+TOXPFW(L,K,NS,NT)
         END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KC
         DO L=2,LA
          TOXPFW(L,K,NS,NT)=TOXPFW(L,K,NS,NT)/(1.+TOXPFTW(L,K,NT))
         END DO
        END DO
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS
C **  IN SEDIMENT BED
C
C **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED
C **  TOXPARB(NS,NT) = PARTITION COEFFICIENT IN SEDIMENT BED
C **  TOXPFTB(L,NT) = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
C                       USED AS TEMPORARY VARIBLE IN THIS AND
C                       FOLLOWING CODE BLOCK
C
      IF(ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        DO K=1,KB
        DO L=2,LA
          TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        TOXPFTB(L,K,NT)=0.
       END DO
       END DO
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         IF(SEDBALL(L,K).GT.0.0) THEN
           HBEDTMP=HBED(L,K)+1.E-9
           TOXPFB(L,K,NS,NT)=TOXPFB(L,K,NS,NT)/
     $                   (PORBED(L,K)*HBEDTMP+TOXPFTB(L,K,NT))
          ELSE
           TOXPFB(L,K,NS,NT)=0.
         END IF
        END DO
        END DO
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  CALCULATE PARTICULATE TOXIC CONTAMINANT SETTLING
C **  AND BED EXCHANGE FLUXES
C
C **  TOXF(L,KC,NT) = TOXIC CONTAMINANT SETTLING AND BED EXCHANGE
C                        FLUX.  USED AS TEMPORARY VARIABLE IN THIS
C                        AND FOLLOWING CODE BLOCK
C
      IF(ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX
       DO K=0,KC
        DO L=2,LA
         TOXF(L,K,NT)=0.
        END DO
       END DO
      END DO
C
      IF(KC.GE.2) THEN
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO K=1,KS
           DO L=2,LA
            TOXF(L,K,NT)=TOXF(L,K,NT)+SEDF(L,K,NS)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO K=1,KS
           DO L=2,LA
            TMPVAL=1.
            IF(SED(L,K+1,NS).GT.0.) TMPVAL=SED(L,K+1,NS)**TMPEXP
            TOXF(L,K,NT)=TOXF(L,K,NT)+TMPVAL*SEDF(L,K,NS)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO K=1,KS
           DO L=2,LA
            TOXF(L,K,NT)=TOXF(L,K,NT)+SNDF(L,K,NX)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO K=1,KS
           DO L=2,LA
            TMPVAL=1.
            IF(SND(L,K+1,NX).GT.0.) TMPVAL=SND(L,K+1,NX)**TMPEXP
            TOXF(L,K,NT)=TOXF(L,K,NT)+TMPVAL*SNDF(L,K,NX)*TOXPARW(NS,NT)
           END DO
          END DO
        END IF
       END DO
       END IF
      END DO
      END IF
C
      IF(KC.GE.2) THEN
      DO NT=1,NTOX
       DO K=1,KS
        DO L=2,LA
         TOXF(L,K,NT)=TOXF(L,K,NT)/(1.+TOXPFTW(L,K+1,NT))
        END DO
       END DO
      END DO
      END IF
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO L=2,LA
           TOXF(L,0,NT)=TOXF(L,0,NT)+
     $                   MIN(SEDF(L,0,NS),0.)*TOXPARW(NS,NT)
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO L=2,LA
           TMPVAL=1.
           IF(SED(L,1,NS).GT.0.) TMPVAL=SED(L,1,NS)**TMPEXP
           TOXF(L,0,NT)=TOXF(L,0,NT)+
     $         TMPVAL*MIN(SEDF(L,0,NS),0.)*TOXPARW(NS,NT)
          END DO
        END IF
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        IF(ITXPARW(NS,NT).EQ.0) THEN
          DO L=2,LA
           TOXF(L,0,NT)=TOXF(L,0,NT)+
     $          MIN(SNDF(L,0,NX),0.)*TOXPARW(NS,NT)
          END DO
        END IF
        IF(ITXPARW(NS,NT).EQ.1) THEN
          TMPEXP=CONPARW(NS,NT)
          DO L=2,LA
           TMPVAL=1.
           IF(SND(L,1,NX).GT.0.) TMPVAL=SND(L,1,NX)**TMPEXP
           TOXF(L,0,NT)=TOXF(L,0,NT)+
     $          TMPVAL*MIN(SNDF(L,0,NX),0.)*TOXPARW(NS,NT)
          END DO
        END IF
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO L=2,LA
        TOXF(L,0,NT)=TOXF(L,0,NT)/(1.+TOXPFTW(L,1,NT))
       END DO
      END DO
C
      DO NT=1,NTOX
       DO L=2,LA
        TOXFB(L,KB,NT)=0.
       END DO
      END DO
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        DO L=2,LA
         TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)
     &                 +MAX(SEDF(L,0,NS),0.)*TOXPARB(NS,NT)
        END DO
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        DO L=2,LA
         TOXFB(L,KBT(L),NT)=TOXFB(L,KBT(L),NT)
     &                +MAX(SNDF(L,0,NX),0.)*TOXPARB(NS,NT)
        END DO
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO L=2,LA
        IF(SEDBALL(L,KBT(L)).GT.0.0) THEN
          HBEDTMP=HBED(L,KBT(L))+1.E-9
          TOXFB(L,KBT(L),NT)
     $    =TOXFB(L,KBT(L),NT)/(PORBED(L,KBT(L))*HBEDTMP+
     $     TOXPFTB(L,KBT(L),NT)+1.0e-8)
         ELSE
          TOXFB(L,KBT(L),NT)=0.
        END IF
       END DO
      END DO
C
C ** DIAGNOSTICS OF FLUX
C
       IF(N.EQ.1)THEN
       OPEN(2,file='toxflx.dia')
       CLOSE(2,STATUS='DELETE')
       OPEN(2,file='toxflx.dia')
       DO L=2,LA
        WRITE(2,2222)IL(L),JL(L),HBED(L,KBT(L)),TOXPFTB(L,KBT(L),1),
     $         TOXFB(L,KBT(L),1),TOXF(L,0,1)
       END DO
       CLOSE(2)
       END IF
c
      END IF
C
C**********************************************************************C
C
C **  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN AND BED
C
      IF(ISTRAN(5).GE.1) THEN
C
      DO NT=1,NTOX
       DO K=1,KC
        DO L=2,LA
          TOXPFTW(L,K,NT)=TOXPFTW(L,K,NT)/(1.+TOXPFTW(L,K,NT))
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        IF(SEDBALL(L,K).GT.0.0) THEN
          HBEDTMP=HBED(L,K)+1.E-9
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     $                   /(PORBED(L,K)*HBEDTMP+TOXPFTB(L,K,NT))
         ELSE
          TOXPFTB(L,K,NT)=1.
        END IF
       END DO
       END DO
      END DO
C
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=1 (SINGLE LAYER IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.1) THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBT(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBT(L),NT)
          BB11=WVEL*TOX(L,1,NT)
          BB22=DELTI*TOXB1(L,KBT(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     $                      +S2TL*TOXB1(L,KBT(L),NT)
          TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
  676 FORMAT('N,I,J,TW,TB,TB1,TF,TFB=',3I5,5E13.4)
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=2 (TWO LAYERS IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.2) THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
         K=2
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBT(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBT(L),NT)
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBT(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     $                      +S2TL*TOXB1(L,KBT(L),NT)
          TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC=3 (THREE LAYERS IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.EQ.3) THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
         K=3
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         K=2
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBT(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBT(L),NT)
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBT(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     $                      +S2TL*TOXB1(L,KBT(L),NT)
          TOXB(L,KBT(L),NT)=DETI*(AA11*BB22-AA21*BB11)
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT, KC.GE.3 (THREE OR MORE LAYERS IN VERTICAL)
C
      IF(ISTRAN(5).GE.1.AND.KC.GE.3) THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
         K=KC
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
         END DO
C
         DO K=KS,2,-1
          DO L=2,LA
           WVEL=DELTI*HP(L)*DZC(K)
           CLEFT=WVEL-TOXF(L,K-1,NT)
           CRIGHT=WVEL*TOX(L,K,NT)-TOXF(L,K,NT)*TOX(L,K+1,NT)
           TOX(L,K,NT)=CRIGHT/CLEFT
          END DO
         END DO
C
         DO L=2,LA
          WVEL=DELTI*HP(L)*DZC(1)
          AA11=WVEL-TOXF(L,0,NT)
          AA12=-TOXFB(L,KBT(L),NT)
          AA21=TOXF(L,0,NT)
          AA22=DELTI+TOXFB(L,KBT(L),NT)
          BB11=WVEL*TOX(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)
          BB22=DELTI*TOXB1(L,KBT(L),NT)
          DETI=1./(AA11*AA22-AA12*AA21)
          TOX(L,1,NT)=DETI*(BB11*AA22-BB22*AA12)
          TOXB1(L,KBT(L),NT)=S3TL*TOXB(L,KBT(L),NT)
     $                      +S2TL*TOXB1(L,KBT(L),NT)
          TOXBTMP=DETI*(AA11*BB22-AA21*BB11)
          IF(TOXBTMP.LT.0.0) THEN
             TOXBTMP=TOXB1(L,KBT(L),NT)
             TOX(L,1,NT)=TOXS(L,1,NT)-TOXF(L,1,NT)*TOX(L,2,NT)/WVEL
          END IF
          TOXB(L,KBT(L),NT)=TOXBTMP
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
C **  UPDATE SEDIMENT BED PHYSICAL PROPERTIES
C
      IF(ISTRAN(6).GE.1.OR.ISTRAN(7).GE.1) THEN
C
C----------------------------------------------------------------------C
C
C ** CONSTANT POROSITY BED
c
      if(N.eq.1) then    ! save initial hbed, Ji, 10/21/00
      do L=2,LA
      do K=1,kb
      hbed00(L,K)=HBED(L,K)
      do NS=1,NSED
      sedb00(L,K,NS)=sedb1(L,K,NS)
      enddo
      enddo
      do k=kb+1,kbm
      hbed00(L,K)=Hbed00(L,kb)   ! the rest = the inital top layer's hbed
      do NS=1,NSED
      sedb00(L,K,NS)=sedb1(L,Kb,NS)
      enddo
      enddo
      enddo
      endif
c


      IF(IBMECH.EQ.0) THEN
c
c ** update top layer thickness
c  VOIDCON1=water volume/solid volume in the bed = porosity/(1-porosity), John & Ji, 4/3/00
c
	voidcon1=1.0  ! =1.0 -> porosity=0.5
c
      DO L=2,LA
       HBEDTMP=HBED1(L,KBT(L))+DELT*VOIDCON1*QSBDTOP(L) ! Bug, VOIDCON1 not defined, Ji
       QWBDTOP(L)=-BEDPORC*DELTI*(HBEDTMP-HBED1(L,KBT(L)))+QWBDTOP(L)
       HBED1(L,KBT(L))=S3TL*HBED(L,KBT(L))+S2TL*HBED1(L,KBT(L))
       HBED(L,KBT(L))=HBEDTMP
c      write(601,6011) N,L,KBT(L),qsbdtop(L),qwbdtop(L),hbedtmp,
c    +  (hbed1(L,Kx),Hbed(L,kx),kx=1,kbt(L))
6011  format(3i6,999e12.4)
      if(kbt(L).eq.1) HBED(L,1)=max(HBED(L,1),0.0)  ! Ji, 10/29/00
      END DO
c
      IF(KB.GT.1) THEN
c
c Case 1:  ** add new top layer
C
      DO L=2,LA
      sedbmax=sedb00(L,kbt(L),1)       ! hardwired for one class of sediment, Ji, 10/29/00
      IF(sedb(L,KBT(L), 1).gt.sedbMAX.and.KBT(L).LT.KBM) THEN
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
         sedb(L,KBT(L)+1,ns)=sedb(L,KBT(L),ns)-sedbMAX
         sedb(L,KBT(L),ns)=sedbMAX
         END DO
         END IF
C
c        IF(ISTRAN(7).GE.1) THEN
C
c        IF(ISTRAN(5).GE.1) THEN
C
         KBT(L)=KBT(L)+1
c
       END IF
      END DO
c
c Case 2:  ** rezone with new top layer added next time step
c
c     DO L=2,LA
c     sedbmax=sedb00(L,kbt(L),1)
c     IF(sedb(L,KBT(L), 1).gt.sedbMAX.         ! hardwired for one class of sediment, Ji, 10/29/00
c    +    and.KBT(L).eq.KBM.and.kbm.ge.2) THEN
C
c        HBED(L,1)=HBED(L,1)+HBED(L,2)         ! Ji, 10/29/00, john's way is completely wrong!
c        DO K=2,KBT(L)-1
c        HBED(L,K)=HBED(L,K+1)
c        END DO
c        HBED(L,KBT(L))=0.0
C
c        IF(ISTRAN(6).GE.1) THEN
c        DO NS=1,NSED
c        SEDB(L,1,NS)=SEDB(L,1,NS)+SEDB(L,2,NS)
c        DO K=2,KBT(L)-1
c        SEDB(L,K,NS)=SEDB(L,K+1,NS)
c        END DO
c        SEDB(L,KBT(L),NS)=0
c        END DO
c        END IF
c
c Note, sedb & hbed must be consistent, check it! Ji, 10/29/00
c        do k=1,kbt(L)
c        diff= hbed(L,k)-sden(1)*sedb(L,K,1)
c        write(801,8011) N,L,K,diff,hbed(L,k),sden(1)*sedb(L,K,1)
c8011  format(3i6,999e12.4)
c        enddo
C
c        IF(ISTRAN(7).GE.1) THEN            ! bug, but not corrected, Ji, 10/29/00
C
c        IF(ISTRAN(5).GE.1) THEN            ! bug, but not corrected, Ji, 10/29/00
C
c     KBT(L)=KBT(L)-1
C
c      END IF
c     END DO
c
c Case 3:  ** remove top layer
c
      DO L=2,LA
      sedbmax=sedb00(L,kbt(L),1)      ! hardwired for one class of sediment, Ji, 10/29/00
      IF(sedb(L,KBT(L), 1).lt.0.0.and.KBT(L).ge.2) THEN
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
!         SEDB(L,KBT(L),  NS)=sedbmax  ! Move next layer to the top layer and reset critical erosion stress
          SEDB(L,KBT(L)-1,NS)=SEDB(L,KBT(L)-1,NS)+SEDB(L,KBT(L),NS)
          SEDB(L,KBT(L),  NS)=0.0
         END DO
         END IF
C
         KBT(L)=KBT(L)-1
C
       END IF
      END DO
c
      END IF
      END IF
C
C----------------------------------------------------------------------C
C
C ** SIMPLE CONSOLIDATING BED
C
      IF(IBMECH.EQ.1) THEN
c
c ** determine time difference and update void ratio
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
        VDRBED2(L,K)=VDRBED1(L,K)
        TMP=(VDRBED1(L,K)-SEDVDRM)/(SEDVDRD-SEDVDRM)
        DELTVDR=-LOG(TMP)/SEDVRDT
        DELTVDR=DELTVDR+DELT
        FACTOR=EXP(-SEDVRDT*DELTVDR)
c        VDRTMP=SEDVDRM+(SEDVDRD-SEDVDRM)*FAC    ! Bug, fac not defined, Ji
        VDRTMP=SEDVDRM+(SEDVDRD-SEDVDRM)*FACtor  !?, JI
        VDRBED1(L,K)=S3TL*VDRBED(L,K)+S2TL*VDRBED1(L,K)
        VDRBED(L,K)=VDRTMP
       END IF
      END DO
      END DO
c
c ** update porosity
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
        PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
        PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
       END IF
      END DO
      END DO
c
c ** update lower layer thickness
c
      DO L=2,LA
       QWTRBED(L,0)=0.
      END DO
C
      DO K=1,KB
      DO L=2,LA
       KBTM1=KBT(L)-1
       IF(K.LE.KBTM1)THEN
        HBEDTMP=(1.+VDRBED(L,K))*HBED1(L,K)/(1.+VDRBED1(L,K))
        QWTRBED(L,K)=QWTRBED(L,K-1)-DELTI*(PORBED(L,K)*HBEDTMP)
     $             +DELTI*(VDRBED2(L,K)*HBED1(L,K)/(1.+VDRBED2(L,K)))
        HBED1(L,K)=S3TL*HBED(L,K)+S2TL*HBED1(L,K)
        HBED(L,K)=HBEDTMP
       END IF
      END DO
      END DO
c
c ** update top layer thickness
c
      DO L=2,LA
       HBEDTMP=(HBED1(L,KBT(L))/(1.+VDRBED1(L,KBT(L))))
     $        +DELT*QSBDTOP(L)
       HBEDTMP=HBEDTMP*(1.+VDRBED(L,KBT(L)))
       QWBDTOP(L)=QWTRBED(L,KBT(L)-1)-DELTI*(PORBED(L,KBT(L))*HBEDTMP)
     $   +DELTI*(VDRBED2(L,KBT(L))*HBED1(L,KBT(L))
     $   /(1.+VDRBED2(L,KBT(L))))+QWBDTOP(L)
       HBED1(L,KBT(L))=S3TL*HBED(L,KBT(L))+S2TL*HBED1(L,KBT(L))
       HBED(L,KBT(L))=HBEDTMP
      END DO
c
c ** add new top layer
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).GT.HBEDMAX) THEN
       IF(KBT(L).LT.KBM) THEN
C
         HOLDTOP=HBED(L,KBT(L))
         KBT(L)=KBT(L)+1
         HBED(L,KBT(L))=HOLDTOP-HBEDMAX
         HBED(L,KBT(L)-1)=HBEDMAX
         FKBT=HBED(L,KBT(L))/HOLDTOP
         FKBTM=HBED(L,KBT(L)-1)/HOLDTOP
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDBOLD=SEDB(L,KBT(L),NS)
          SEDB(L,KBT(L),NS)=FKBT*SEDBOLD
          SEDB(L,KBT(L)-1,NS)=FKBTM*SEDBOLD
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDBOLD=SNDB(L,KBT(L),NS)
          SNDB(L,KBT(L),NS)=FKBT*SNDBOLD
          SNDB(L,KBT(L)-1,NS)=FKBTM*SNDBOLD
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
         DO NT=1,NTOX
          TOXBOLD=TOXB(L,KBT(L),NT)
          TOXB(L,KBT(L),NT)=FKBT*TOXBOLD
          TOXB(L,KBT(L)-1,NT)=FKBTM*TOXBOLD
         END DO
         END IF
C
       END IF
       END IF
      END DO
c
c ** rezone with new top layer added next time step
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).GT.HBEDMAX) THEN
       IF(KBT(L).EQ.KBM.AND.KBM.GT.1) THEN
C
         HBED(L,1)=HBED(L,1)+HBED(L,2)
         IF(KBM.EQ.2) THEN
            HBED(L,2)=0
C            KBT(L)=1
         END IF
         IF(KBM.GT.2) THEN
           DO K=2,KBT(L)-1
            HBED(L,K)=HBED(L,K-1)
           END DO
           HBED(L,KBT(L))=0
C           KBT(L)=KBT(L)-1
         END IF
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDB(L,1,NS)=SEDB(L,1,NS)+SEDB(L,2,NS)
          IF(KBM.EQ.2) THEN
             SEDB(L,2,NS)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             SEDB(L,K,NS)=SEDB(L,K-1,NS)
            END DO
            SEDB(L,KBT(L),NS)=0
          END IF
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDB(L,1,NS)=SNDB(L,1,NS)+SNDB(L,2,NS)
          IF(KBM.EQ.2) THEN
             SNDB(L,2,NS)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             SNDB(L,K,NS)=SEDB(L,K-1,NS)
            END DO
            SNDB(L,KBT(L),NS)=0
          END IF
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
          DO NT=1,NTOX
          TOXB(L,1,NT)=TOXB(L,1,NT)+TOXB(L,2,NT)
          IF(KBM.EQ.2) THEN
             TOXB(L,2,NT)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             TOXB(L,K,NS)=TOXB(L,K-1,NT)
            END DO
            TOXB(L,KBT(L),NT)=0
          END IF
         END DO
         END IF
C
         IF(KBM.EQ.2) THEN
            KBT(L)=1
         END IF
         IF(KBM.GT.2) THEN
           KBT(L)=KBT(L)-1
         END IF
C
       END IF
       END IF
      END DO
c
c ** remove top layer
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).LT.0.) THEN
       IF(KBT(L).GT.1) THEN
C
         HBED(L,KBT(L)-1)=HBED(L,KBT(L)-1)+HBED(L,KBT(L))
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDB(L,KBT(L)-1,NS)=SEDB(L,KBT(L)-1,NS)+SEDB(L,KBT(L),NS)
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDB(L,KBT(L)-1,NS)=SNDB(L,KBT(L)-1,NS)+SNDB(L,KBT(L),NS)
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
         DO NT=1,NTOX
           TOXB(L,KBT(L)-1,NT)=TOXB(L,KBT(L)-1,NT)+TOXB(L,KBT(L),NT)
         END DO
         END IF
C
         KBT(L)=KBT(L)-1
C
       END IF
       END IF
      END DO
c
      END IF
C
C----------------------------------------------------------------------C
C
C ** FINITE STRAIN CONSOLIDATING BED
C
      IF(IBMECH.EQ.2) THEN
C
      WDENKGM3=1.E3
      WDENGMM3=1.E6
c
c ++  calculate transformed thickness of bed layers
c     dzbtr = transformed thickness of bed layer
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
       END IF
      END DO
      END DO
c
c ++  calculate effective stress
c     and hydraulic conductivity divied by (1+void)
c     hydcn = hydraulic conductivity divided by
c             (1+void ratio)    hydcn = K/(1+void)
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         STRSE(L,K)=FSTRSE(VDRBED(L,K))
         HYDCN(L,K)=FHYDCN(VDRBED(L,K))/(G*WDENKGM3*DZBTR(L,K))
       END IF
      END DO
      END DO
c
      DO K=1,KB
      DO L=2,LA
       KBTM1=KBT(L)-1
       IF(K.LE.KBTM1)THEN
         TMPVAL=( 1./HYDCN(L,K) )
     $       +( 1./HYDCN(L,K+1) )
         CCOEF(L,K)=2./TMPVAL
       END IF
      END DO
      END DO
C
c ++  calculate pressure components
c
      TMPVAL=G*WDENKGM3
c
      ZETATOP(L)=0.
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         ZETATOP(L)=ZETATOP(L)+DZBTR(L,K)*SGSM1(L,K)
       END IF
      END DO
      END DO
c
      STRSEM(L,1)=DZBTR(L,1)*SGSM1(L,1)
c
      DO K=2,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         STRSEM(L,K)=STRSEM(L,K-1)+DZBTR(L,K)*SGSM1(L,K)
       END IF
      END DO
      END DO
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         STRSEM(L,K)=(STRSEM(L,K)-0.5*DZBTR(L,K)*SGSM1(L,K))
       END IF
      END DO
      END DO
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         STRSEM(L,K)=TMPVAL*(ZETATOP(L)-STRSEM(L,K))
       END IF
      END DO
      END DO
c
c     prese = excess pore pressure in layer
c     presh = hydrostatic pressure in layer
c     prest = total pressure in layer
c     strst = total stress in layer
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         PRESE(L,K)=-STRSE(L,K)+STRSEM(L,K)
         PRESH(L,K)=G*WDENKGM3*(ZBEDGT(L)-ZBEDC(L,K))
         PREST(L,K)=PRESE(L,K)+PRESH(L,K)
         STRST(L,K)=PREST(L,K)+STRSE(L,K)
       END IF
      END DO
      END DO
c
c ++  calculate water volume fluxes
c     qwatrbed = specific water volume flux at top of layer k
c
      DO L=2,LA
       QWTRBED(L,0)=0.
      END DO
C
      DO K=1,KB
      DO L=2,LA
       KBTM1=KBT(L)-1
       IF(K.LE.KBTM1)THEN
         QWTRBED(L,K)=-CCOEF(L,K)*(PRESE(L,K+1)-PRESE(L,K))
       END IF
      END DO
      END DO
c
      DO L=2,LA
       QWTRBED(L,KBT(L))=CCOEF(L,KBT(L))*PRESE(L,KBT(L))
      END DO
c
c ++  calculate void ratios
c     vdrbed =  void ratio of bed layer
c
      VDRBED(L,1)=VDRBED1(L,1)-DELT*QWTRBED(L,1)/DZBTR1(L,1)
c
      DO K=2,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
       VDRBED(L,K)=VDRBED1(L,K)
     $         -DELT*(QWTRBED(L,K)-QWTRBED(L,K-1))/DZBTR1(L,K)
       END IF
      END DO
      END DO
c
	voidd=0.0 !bug, Ji, 10/10/00
	wdepsdd=0.0 !bug, Ji, 10/10/00
      DO L=2,LA
      VDRBED(L,KBT(L))=VDRBED(L,KBT(L))
     $         -DELT*(VDRBED1(L,KBT(L))-voidd)*wdepsdd/DZBTR1(L,KBT(L)) !bug, Ji, voidd, wdepsdd, is not defined, 3/5/00
      END DO
      write(6,*) "Voiddd= ", voidd
	stop
c
c ++  update layer thickness
c     hbed = bed layer thickness
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
        HBED(L,K)=HBED1(L,K)*(1.+VDRBED(L,K))/(1.+VDRBED1(L,K))
       END IF
      END DO
      END DO
c
      DO L=2,LA
       HBED(L,KBT(L))=HBED(L,KBT(L))
     $      +DELT*(1.+VDRBED(L,KBT(L)))*(wdepsdd-werosdd) !JI, bug, wdepsdd, werosdd not defined, 3/5/00
      write(6,*) "wdepsdd, werosdd= ", wdepsdd, werosdd
      END DO
c
c     zbedc = vertical coordinate of the center of bed layer
c     zbedg = vertical coordinate at top of bed layer
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         ZBEDG(L,K)=ZBEDG(L,K-1)+HBED(L,K)
       END IF
      END DO
      END DO
c
      DO L=2,LA
       ZBEDGT(L)=ZBEDG(L,KBT(L))
      END DO
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         ZBEDC(L,K)=0.5*(ZBEDG(L,K)+ZBEDG(L,K-1))
       END IF
      END DO
      END DO
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
         DZBTR(L,K)=HBED(L,K)/(1.+VDRBED(L,K))
       END IF
      END DO
      END DO
c
c ** update porosity
c
      DO K=1,KB
      DO L=2,LA
       IF(K.LE.KBT(L))THEN
        PORBED(L,K)=VDRBED(L,K)/(1.+VDRBED(L,K))
        PORBED1(L,K)=VDRBED1(L,K)/(1.+VDRBED1(L,K))
       END IF
      END DO
      END DO
c
c ** add new top layer
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).GT.HBEDMAX) THEN
       IF(KBT(L).LT.KBM) THEN
C
         HOLDTOP=HBED(L,KBT(L))
         KBT(L)=KBT(L)+1
         HBED(L,KBT(L))=HOLDTOP-HBEDMAX
         HBED(L,KBT(L)-1)=HBEDMAX
         FKBT=HBED(L,KBT(L))/HOLDTOP
         FKBTM=HBED(L,KBT(L)-1)/HOLDTOP
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDBOLD=SEDB(L,KBT(L),NS)
          SEDB(L,KBT(L),NS)=FKBT*SEDBOLD
          SEDB(L,KBT(L)-1,NS)=FKBTM*SEDBOLD
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDBOLD=SNDB(L,KBT(L),NS)
          SNDB(L,KBT(L),NS)=FKBT*SNDBOLD
          SNDB(L,KBT(L)-1,NS)=FKBTM*SNDBOLD
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
         DO NT=1,NTOX
          TOXBOLD=TOXB(L,KBT(L),NT)
          TOXB(L,KBT(L),NT)=FKBT*TOXBOLD
          TOXB(L,KBT(L)-1,NT)=FKBTM*TOXBOLD
         END DO
         END IF
C
       END IF
       END IF
      END DO
c
c ** rezone with new top layer added next time step
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).GT.HBEDMAX) THEN
       IF(KBT(L).EQ.KBM.AND.KBM.GT.1) THEN
C
         HBED(L,1)=HBED(L,1)+HBED(L,2)
         IF(KBM.EQ.2) THEN
            HBED(L,2)=0
C            KBT(L)=1
         END IF
         IF(KBM.GT.2) THEN
           DO K=2,KBT(L)-1
            HBED(L,K)=HBED(L,K-1)
           END DO
           HBED(L,KBT(L))=0
C           KBT(L)=KBT(L)-1
         END IF
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDB(L,1,NS)=SEDB(L,1,NS)+SEDB(L,2,NS)
          IF(KBM.EQ.2) THEN
             SEDB(L,2,NS)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             SEDB(L,K,NS)=SEDB(L,K-1,NS)
            END DO
            SEDB(L,KBT(L),NS)=0
          END IF
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDB(L,1,NS)=SNDB(L,1,NS)+SNDB(L,2,NS)
          IF(KBM.EQ.2) THEN
             SNDB(L,2,NS)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             SNDB(L,K,NS)=SEDB(L,K-1,NS)
            END DO
            SNDB(L,KBT(L),NS)=0
          END IF
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
          DO NT=1,NTOX
          TOXB(L,1,NT)=TOXB(L,1,NT)+TOXB(L,2,NT)
          IF(KBM.EQ.2) THEN
             TOXB(L,2,NT)=0
          END IF
          IF(KBM.GT.2) THEN
            DO K=2,KBT(L)-1
             TOXB(L,K,NS)=TOXB(L,K-1,NT)
            END DO
            TOXB(L,KBT(L),NT)=0
          END IF
         END DO
         END IF
C
         IF(KBM.EQ.2) THEN
            KBT(L)=1
         END IF
         IF(KBM.GT.2) THEN
           KBT(L)=KBT(L)-1
         END IF
C
       END IF
       END IF
      END DO
c
c ** remove top layer
c
      DO L=2,LA
       IF(HBED(L,KBT(L)).LT.0.) THEN
       IF(KBT(L).GT.1) THEN
C
         HBED(L,KBT(L)-1)=HBED(L,KBT(L)-1)+HBED(L,KBT(L))
C
         IF(ISTRAN(6).GE.1) THEN
         DO NS=1,NSED
          SEDB(L,KBT(L)-1,NS)=SEDB(L,KBT(L)-1,NS)+SEDB(L,KBT(L),NS)
         END DO
         END IF
C
         IF(ISTRAN(7).GE.1) THEN
         DO NS=1,NSND
          SNDB(L,KBT(L)-1,NS)=SNDB(L,KBT(L)-1,NS)+SNDB(L,KBT(L),NS)
         END DO
         END IF
C
         IF(ISTRAN(5).GE.1) THEN
         DO NT=1,NTOX
           TOXB(L,KBT(L)-1,NT)=TOXB(L,KBT(L)-1,NT)+TOXB(L,KBT(L),NT)
         END DO
         END IF
C
         KBT(L)=KBT(L)-1
C
       END IF
       END IF
      END DO
c
      END IF
C
C----------------------------------------------------------------------C
C
      END IF
C
C**********************************************************************C
C
C **  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
C
      DO NT=1,NTOX
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=0.
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       IF(ISTRAN(6).GE.1) THEN
       DO NS=1,NSED
        DO K=1,KB
        DO L=2,LA
         TOXPFB(L,K,NS,NT)=SEDB(L,K,NS)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
       IF(ISTRAN(7).GE.1) THEN
       DO NX=1,NSND
        NS=NX+NSED
        DO K=1,KB
        DO L=2,LA
          TOXPFB(L,K,NS,NT)=SNDB(L,K,NX)*TOXPARB(NS,NT)
        END DO
        END DO
       END DO
       END IF
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        TOXPFTB(L,K,NT)=0.
       END DO
       END DO
       DO NS=1,NSED+NSND
        DO K=1,KB
        DO L=2,LA
         TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
        END DO
        END DO
       END DO
      END DO
C
      DO NT=1,NTOX
       DO K=1,KB
       DO L=2,LA
        IF(HBED(L,K).GT.0.0) THEN
          TOXPFTB(L,K,NT)=TOXPFTB(L,K,NT)
     &                   /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
         ELSE
          TOXPFTB(L,K,NT)=1.
        END IF
       END DO
       END DO
      END DO
C
C**********************************************************************C
C
C **  DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER COLUMN LAYER
C
      IF(ISTRAN(5).GE.1) THEN
        DO NT=1,NTOX
C
        DO L=2,LA
         KBTP1=KBT(L)+1
         ALOW(L,1)=0.
         CUPP(L,KBTP1)=0.
         DO K=1,KBT(L)-1
          CUPP(L,K)=-DIFTOX(NT)*(PORBED(L,K)+PORBED(L,K+1))/
     $                     (HBED(L,K)+HBED(L,K+1))
         END DO
         CUPP(L,KBT(L))=-2.*DIFTOX(NT)*PORBED(L,KBT(L))/HBED(L,KBT(L))
         DO K=2,KBTP1
          ALOW(L,K)=CUPP(L,K-1)
         END DO
         DO K=1,KBT(L)
          BMNN(L,K)=DELTI*HBED(L,K)/(1.-TOXPFTB(L,K,NT))
         END DO
         BMNN(L,KBTP1)=DELTI*DZC(1)*HP(L)/(1.-TOXPFTW(L,1,NT))
         DO K=1,KBTP1
          BMNN(L,K)=BMNN(L,K)-ALOW(L,K)-CUPP(L,K)
         END DO
         DO K=1,KBT(L)
          RRHS(L,K)=DELTI*TOXB(L,K,NT)
         END DO
         RRHS(L,KBTP1)=DELTI*DZC(1)*HP(L)*TOX(L,1,NT)
        END DO
C
        DO L=2,LA
         KBTP1=KBT(L)+1
         BETTMP=BMNN(L,1)
         TOXTMP(L,1)=RRHS(L,1)/(BETTMP+1.0e-8)  ! J.S.
         DO KK=2,KBTP1
          GAMTMP(L,KK)=CUPP(L,KK-1)/(BETTMP+1.0e-8)
          BETTMP=BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
          TOXTMP(L,KK)=(RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))
     &		/(BETTMP+1.0e-8)
         END DO
         DO KK=KBT(L),1,-1
          TOXTMP(L,KK)=TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)
         END DO
        END DO
C
        DO L=2,LA
         KBTP1=KBT(L)+1
         DO K=1,KBT(L)
          TOXB(L,K,NT)=HBED(L,K)*TOXTMP(L,K)/(1.-TOXPFTB(L,K,NT))
         END DO
         TOX(L,1,NT)=TOXTMP(L,KBTP1)/(1.-TOXPFTW(L,1,NT))
        END DO
C
        END DO
       END IF
C
C**********************************************************************C
C
C **  TOXIC CONTAMINANT REACTIONS
C
      IF(ISTRAN(5).GE.1) THEN
        DO NT=1,NTOX
C
C----------------------------------------------------------------------C
C
C **     NOTES:
C
C        BULK DECAY COEFFICIENT
C
C          RKTOXWT=RKTOXW(NT)    !*(TEM(L,K)-TKTOXW(NT))**ALPTOX
C          RKTOXBT=RKTOXB(NT)    !*(TEM(L,1)-TRTOXB(NT))**ALPTOX
C
C        VOLITIZATION
C
C          VOLTOX(NT)
C          RMOLTX(NT)=MOLECULAR WEIGHT
C
C        PHOTOLOSIS
C
C          RKTOXP(NT)=BASE RATE
C          SKTOXP(NT)=SOLAR RADIATION AT BASE RATE
C
         DO K=1,KC
         DO L=2,LA
          CDECAYW(L,K)=1./(1.+DELT*RKTOXW(NT))
         END DO
         END DO
C
         DO K=1,KC
         DO L=2,LA
          TOX(L,K,NT)=CDECAYW(L,K)*TOX(L,K,NT)
         END DO
         END DO
C
         DO K=1,KB
         DO L=2,LA
          CDECAYB(L,K)=1./(1.+DELT*RKTOXB(NT))
         END DO
         END DO
C
         DO K=1,KB
         DO L=2,LA
          TOXB(L,K,NT)=CDECAYB(L,K)*TOXB(L,K,NT)
         END DO
         END DO
C
C----------------------------------------------------------------------C
C
        END DO
      END IF
C
C**********************************************************************C
C
      CLOSE(1)
      CLOSE(11)
      CLOSE(21)
      CLOSE(31)
      CLOSE(41)
C
C**********************************************************************C
C
      RETURN
      END

