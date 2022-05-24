C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WASP4
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
C
C **  SUBROUTINE WASPOUT WRITES OUTPUT FILES PROVIDING ADVECTIVE AND
C **  DIFFUSIVE TRANSPORT FIELDS FOR THE WASP4  WATER QUALITY MODEL
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      DIMENSION LUTMP(KCM*LCM),LDTMP(KCM*LCM),
     $          QTMP(KCM*LCM)
C
C**********************************************************************C
C
C **  WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C **  THE VALUE OF X IN THE F10.X FORMATS MAY NEED TO BE CHANGED
C **  FROM PROBLEM TO PROBLEM.  A PRELIMINARY RUN USING E10.3
C **  CAN BE USED TO SPEED THE ADJUSTMENT  
C   
C**********************************************************************C
C
C **  READ CONTROL DATA FOR WRITING TO WASP COMPATIABLE FILES     
C
C----------------------------------------------------------------------C
C
      SVPT=1. 
      IF(NTSMMT.LT.NTSPTC)SVPT=0.
C
      IF(JSWASP.EQ.1) THEN
      OPEN(1,FILE='efdc.wsp',STATUS='UNKNOWN')
C
C1**  READ CELL VOLUME PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
C
C2**  READ DIFFUSION PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) NRFLD,SCALR,CONVR
C
C3**  READ ADVECTION PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ
C
C4**  READ SEDIMENT VOLUME DEPTH
C
      READ(1,1)
      READ(1,1)
      READ(1,*) DEPSED
C
      CLOSE(1)
      END IF
C
    1 FORMAT (80X)
C
C**********************************************************************C
C
C **  WRITE HORIZONTAL POSITION AND LAYER FILE waspp.out 
C **  WRITE INITIAL VOLUME FILE waspc.out    
C
C **  file waspc.out is consistent with data group C specifications 
C **  on page 172 of the wasp4 manual PB88-185095, Jan 1988
C   
C **  file waspp.out defines the layer (1 is surface water layer, with
C **  layer numbering increasing with depth in water column) and
C **  horizontal positions in lon,lat or utme, utmn of the water
C **  quality (long term transport) cells or segements
C
C----------------------------------------------------------------------C
C
      IF (JSWASP.EQ.1) THEN
        OPEN(90,FILE='waspp.out',STATUS='UNKNOWN')
        OPEN(93,FILE='waspc.out',STATUS='UNKNOWN')
        CLOSE(90,STATUS='DELETE')
        CLOSE(93,STATUS='DELETE')
        OPEN(90,FILE='waspp.out',STATUS='UNKNOWN')
        OPEN(93,FILE='waspc.out',STATUS='UNKNOWN')
c       IVOPT=2
c       IBEDV=0
        WRITE(93,1031)IVOPT,IBEDV
c       SCALV=1.
c       CONVV=1.
        WRITE(93,1032)SCALV,CONVV
c       VMULT=0.
c       VEXP=0.
c       DMULT=0.
c       DEXP=0.
        LCLTM2=LCLT-2
        LWASP=0
        IF (KC.GT.1) THEN
          LTYPE=1
          KWASP=1
           DO LT=2,LALT
           LWASP=LWASP+1
           LBELOW=LWASP+LCLTM2
           I=ILLT(LT)
           J=JLLT(LT)
           L=LIJ(I,J)
           VOLUME=DXYP(L)*HLPF(L)*DZC(KC)
           WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
           WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
           END DO
          LTYPE=2
           DO K=KS,2,-1
           KWASP=KC-K+1
            DO LT=2,LALT
            LWASP=LWASP+1
            LBELOW=LWASP+LCLTM2
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            VOLUME=DXYP(L)*HLPF(L)*DZC(K)
            WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
            WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
            END DO
           END DO
        END IF 
        LTYPE=2
        IF (KC.EQ.1) LTYPE=1
        KWASP=KC
         DO LT=2,LALT
         LWASP=LWASP+1
C        LBELOW=0
         LBELOW=LWASP+LCLTM2
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         VOLUME=DXYP(L)*HLPF(L)*DZC(1)
         WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
         END DO
        LTYPE=3
        KWASP=KC+1
         DO LT=2,LALT
         LWASP=LWASP+1
         LBELOW=0
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         VOLUME=DXYP(L)*DEPSED
         WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
         END DO
        CLOSE(90)
        CLOSE(93)
      END IF
C
 1001 FORMAT(4I5,2F10.4)
 1031 FORMAT(2I5)
 1032 FORMAT(2F10.4)
c1033 FORMAT(3I10,5F10.2)
 1033 FORMAT(3I10,5E10.3)
C
C**********************************************************************C
C
C **  WRITE DIFFUSIVE AND DISPERSIVE TRANSPORT FILE waspb.out 
C
C **  file waspb.out is consistent with data group B specifications 
C **  on page 170 of the wasp4 manual PB88-185095, Jan 1988
C   
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1) THEN
        OPEN(91,FILE='waspb.out',STATUS='UNKNOWN')
        CLOSE(91,STATUS='DELETE')
        OPEN(91,FILE='waspb.out',STATUS='UNKNOWN')
c       NRFLD=1
        WRITE(91,1011)NRFLD
        NTEX=NTS/NTSMMT
c       SCALR=1.
c       CONVR=1.
        WRITE(91,1012)NTEX,SCALR,CONVR
        CLOSE(91)
      END IF
C
      OPEN(91,FILE='waspb.out',ACCESS='APPEND',STATUS='UNKNOWN')
C
      LCLTM2=LCLT-2
      NORSH=0
      NORSV=0 
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       NORSH=NORSH+INT(SUB(L))+INT(SVB(L))
       NORSV=NORSV+INT(SPB(L))
       END DO
      NORS=KC*NORSH+KS*NORSV
      WRITE(91,1013)NORS
C
      UNITY=1.
      DO K=KC,1,-1
      KMUL=KC-K
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SUB(L).EQ.1.) THEN
         LWASP=LT-1+KMUL*LCLTM2
         LWASPW=LWASP-1       
         LW=L-1
         ADDLW=DYU(L)*AHULPF(L,K)*DZC(K)*0.5*(HLPF(L)
     $       +HLPF(LW))*DXIU(L)
         WRITE(91,1014)ADDLW,UNITY,LWASPW,LWASP
       END IF
       END DO
      END DO
C
      UNITY=1.
      DO K=KC,1,-1
      KMUL=KC-K
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SVB(L).EQ.1.) THEN
         LWASP=LT-1+KMUL*LCLTM2
         LSLT=LSCLT(LT)
         LWASPS=LSLT-1+KMUL*LCLTM2
         LS=LSC(L)
         ADDLS=DXV(L)*AHVLPF(L,K)*DZC(K)*0.5*(HLPF(L)
     $       +HLPF(LS))*DYIV(L)
         WRITE(91,1014)ADDLS,UNITY,LWASPS,LWASP
       END IF
       END DO
      END DO
C
      IF (KC.GT.1) THEN
        UNITY=1.
        DO K=KS,1,-1
        KMUL1=KS-K
        KMUL2=KMUL1+1
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         IF (SPB(L).EQ.1.) THEN
           LWASP=LT-1+KMUL1*LCLTM2
           LBELOW=LT-1+KMUL2*LCLTM2
           ADDL=DXYP(L)*ABLPF(L,K)*DZIG(K)
           WRITE(91,1014)ADDL,UNITY,LWASP,LBELOW
         END IF 
         END DO
        END DO
      END IF 
C
      NBRK=6
      WRITE(91,1015)NBRK
C
      TSTOP=(DT*FLOAT(N)+TBEGIN*TCON)
      TSTART=TSTOP-DT*FLOAT(NTSMMT)
      TSTOP=TSTOP/86400.
      TSTART=TSTART/86400.
      TSMALL=1.E-5
      D1=0.
      T1=0.-2*TSMALL
      D2=0.
      T2=TSTART-TSMALL
      D3=1.
      T3=TSTART+TSMALL
      D4=1.
      T4=TSTOP-TSMALL
      D5=0.
      T5=TSTOP+TSMALL
      D6=0.
      T6=2*TSMALL+(DT*FLOAT(NTS)+TBEGIN*TCON)/86400.
      WRITE(91,1016)D1,T1,D2,T2,D3,T3,D4,T4
      WRITE(91,1016)D5,T5,D6,T6
C
      CLOSE(91)
C
 1011 FORMAT(I5)
 1012 FORMAT(I5,2F10.4)
 1013 FORMAT(I5)
c1014 FORMAT(2F10.0,2I5)
 1014 FORMAT(2E10.3,2I5)
 1015 FORMAT(I5)
 1016 FORMAT(4(2F10.5))
 1017 FORMAT(16I5)
C
C**********************************************************************C
C
C **  WRITE ADVECTIVE TRANSPORT FILE waspd.out    
C
C **  file waspd.out is consistent with data group D.1 specifications 
C **  on page 174 of the wasp4 manual PB88-185095, Jan 1988
C   
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1) THEN
        OPEN(92,FILE='waspd.out',STATUS='UNKNOWN')
        CLOSE(92,STATUS='DELETE')
        OPEN(92,FILE='waspd.out',STATUS='UNKNOWN')
c       IQOPT=1
c       NFIELD=1
        WRITE(92,1021)IQOPT,NFIELD
        NINQ=NTS/NTSMMT
c       SCALQ=1
c       CONVQ=1
        WRITE(92,1022)NINQ,SCALQ,CONVQ
        CLOSE(92)
      END IF
C
      OPEN(92,FILE='waspd.out',ACCESS='APPEND',STATUS='UNKNOWN')
      LCLTM2=LCLT-2
      NOQSH=0
      NOQSV=0
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       NOQSH=NOQSH+INT(SUB(L))+INT(SVB(L))
       NOQSV=NOQSV+INT(SWB(L))
       END DO
      NOQS=KC*NOQSH+KS*NOQSV
      WRITE(92,1023)NOQS
C
      LL=0
      DO K=KC,1,-1
      KMUL=KC-K
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SUB(L).EQ.1.) THEN
         LL=LL+1
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LDTMP(LL)-1
         QTMP(LL)=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
       END IF
       END DO
      END DO
C
      DO L=1,LL,4
      WRITE(92,1024)QTMP(L),  LUTMP(L),  LDTMP(L),
     $               QTMP(L+1),LUTMP(L+1),LDTMP(L+1),
     $               QTMP(L+2),LUTMP(L+2),LDTMP(L+2),
     $               QTMP(L+3),LUTMP(L+3),LDTMP(L+3)
      END DO
C
      LL=0
      DO K=KC,1,-1
      KMUL=KC-K
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SVB(L).EQ.1.) THEN
         LL=LL+1
         LSLT=LSCLT(LT)
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LSLT-1+KMUL*LCLTM2
         QTMP(LL)=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
       END IF
       END DO
      END DO
C
      DO L=1,LL,4
      WRITE(92,1024)QTMP(L),  LUTMP(L),  LDTMP(L),
     $               QTMP(L+1),LUTMP(L+1),LDTMP(L+1),
     $               QTMP(L+2),LUTMP(L+2),LDTMP(L+2),
     $               QTMP(L+3),LUTMP(L+3),LDTMP(L+3)
      END DO
C
      IF (KC.GT.1) THEN
        LL=0
        DO K=KS,1,-1
        KMUL1=KS-K
        KMUL2=KMUL1+1
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         IF (SWB(L).EQ.1.) THEN
           LL=LL+1
           LUTMP(LL)=LT-1+KMUL1*LCLTM2
           LDTMP(LL)=LT-1+KMUL2*LCLTM2
           QTMP(LL)=-DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
         END IF 
         END DO
        END DO
      END IF
C
      IF (KC.GT.1) THEN
        DO L=1,LL,4
        WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),
     $               QTMP(L+1),LUTMP(L+1),LDTMP(L+1),
     $               QTMP(L+2),LUTMP(L+2),LDTMP(L+2),
     $               QTMP(L+3),LUTMP(L+3),LDTMP(L+3)
        END DO
      END IF
C
      NBRKQ=6
      WRITE(92,1025)NBRKQ
      WRITE(92,1026)D1,T1,D2,T2,D3,T3,D4,T4
      WRITE(92,1026)D5,T5,D6,T6
C
      CLOSE(92)
C
 1021 FORMAT(2I5)
 1022 FORMAT(I5,2F10.4)
 1023 FORMAT(I5)
c1024 FORMAT(4(F10.0,2I5))
 1024 FORMAT(4(E10.3,2I5))
 1025 FORMAT(I5)
 1026 FORMAT(4(2F10.5))
C
C**********************************************************************C
C
C **  WRITE TO DYNHYD.HYD EMULATION FILES waspdh.out AND waspdhu.out
C
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1)THEN
        OPEN(90,FILE='waspdhd.out',STATUS='UNKNOWN')
        OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')
        OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',
     $                      FORM='UNFORMATTED')
        CLOSE(90,STATUS='DELETE')
        CLOSE(94,STATUS='DELETE')
        CLOSE(95,STATUS='DELETE')
        OPEN(90,FILE='waspdhd.out',STATUS='UNKNOWN')
        OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')
        OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',
     $                      FORM='UNFORMATTED')
        KCLC=KC*LCLT
        LCLTM2=LCLT-2
         DO KL=1,KCLC 
         NCHNC(KL)=0
         END DO
         DO M=1,10
          DO KL=1,KCLC 
          LCHNC(KL,M)=0
          END DO
         END DO
        NJUN=KC*(LCLT-2)
        NCHNH=0
        NCHNV=0
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         NCHNH=NCHNH+INT(SUB(L))+INT(SVB(L))
         NCHNV=NCHNV+INT(SWB(L))
         END DO
        NCHN=KC*NCHNH+KS*NCHNV
        ISTMP=0
        NODYN=NFLTMT
        TZERO=TBEGIN*TCON/86400.
        WRITE(90,901)NJUN,NCHN
        WRITE(94,941)NJUN,NCHN,DT,ISTMP,NTS,ISTMP,NODYN,TZERO
        WRITE(95)NJUN,NCHN,DT,ISTMP,NTS,ISTMP,NODYN,TZERO
C
C **  CHANNEL DATA
C
        RMNDUM=0.
        LCHN=0
         DO K=KC,1,-1
         KMUL=KC-K
          DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          IF (SUB(L).EQ.1.) THEN
            LDTM=LT-1+KMUL*LCLTM2
            LUTM=LDTM-1
            RLENTH=DXU(L)
            WIDTH=DYU(L)
            LCHN=LCHN+1
            NCHNC(LDTM)=NCHNC(LDTM)+1
            NCHNC(LUTM)=NCHNC(LUTM)+1
            LCHNC(LDTM,NCHNC(LDTM))=LCHN
            LCHNC(LUTM,NCHNC(LUTM))=LCHN
            WRITE(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            WRITE(95)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
          END IF
          END DO
         END DO
         DO K=KC,1,-1
         KMUL=KC-K
          DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          IF (SVB(L).EQ.1.) THEN
            LSLT=LSCLT(LT)
            LDTM=LT-1+KMUL*LCLTM2
            LUTM=LSLT-1+KMUL*LCLTM2
            RLENTH=DYV(L)
            WIDTH=DXV(L)
            LCHN=LCHN+1
            NCHNC(LDTM)=NCHNC(LDTM)+1
            NCHNC(LUTM)=NCHNC(LUTM)+1
            LCHNC(LDTM,NCHNC(LDTM))=LCHN
            LCHNC(LUTM,NCHNC(LUTM))=LCHN
            WRITE(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            WRITE(95)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
          END IF
          END DO
         END DO
         IF(KC.GT.1)THEN
           DO K=KS,1,-1
           KMUL1=KS-K
           KMUL2=KMUL1+1
            DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SWB(L).EQ.1.) THEN
              LUTM=LT-1+KMUL1*LCLTM2
              LDTM=LT-1+KMUL2*LCLTM2
              RLENTH=HLPF(L)*DZG(K)
              WIDTH=SQRT(DXYP(L))
              LCHN=LCHN+1
              NCHNC(LDTM)=NCHNC(LDTM)+1
              NCHNC(LUTM)=NCHNC(LUTM)+1
              LCHNC(LDTM,NCHNC(LDTM))=LCHN
              LCHNC(LUTM,NCHNC(LUTM))=LCHN
              WRITE(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              WRITE(95)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
            END IF
           END DO
          END DO
        END IF
C
C **  JUNCTION DATA
C
        DO K=KC,1,-1
        KMUL=KC-K
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         LCELL=LT-1+KMUL*LCLTM2
         WRITE(90,904)LCELL,DXYP(L),(LCHNC(LCELL,M),M=1,10)
         WRITE(94,944)DXYP(L),(LCHNC(LCELL,M),M=1,10)
         WRITE(95)DXYP(L),(LCHNC(LCELL,M),M=1,10)
         END DO
        END DO
C
        CLOSE(90)
        CLOSE(94)
        CLOSE(95)
      END IF
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP, VOLUME AND FLOW DATA
C
      OPEN(94,FILE='waspdh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      OPEN(95,FILE='waspdhu.out',ACCESS='APPEND',STATUS='UNKNOWN',
     $                      FORM='UNFORMATTED')
C
      LCLTM2=LCLT-2
      IZERO=0
      RZERO=0
C
      NSTEP=N-NTSMMT
      WRITE(94,945)NSTEP
C
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       LN=LNC(L)
       VOLUM=DXYP(L)*HLPF(L)*DZC(K)
       QIN=QSUMELPF(L)*DZC(K)
       FLOWXI=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
       FLOWYI=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
       FLOWZI=DXYP(L)*(WLPF(L,K-1)+SVPT*WVPT(L,K-1))
       FLOWXO=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
       FLOWYO=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
       FLOWZO=DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
       QQSUM=QIN+FLOWXI+FLOWYI+FLOWZI-FLOWXO-FLOWYO-FLOWZO
       DEPTH=HLPF(L)*DZC(K)
       VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K)
     $         +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
       VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K)
     $         +VHLPF(LN,K)+SVPT*VVPT(LN,K))/HLPF(L)
       VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1)
     $         +WLPF(L,K)+SVPT*WVPT(L,K))
       VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
       WRITE(94,946)VOLUM,QIN,QSUM,DEPTH,VELMAG
       WRITE(95)VOLUM,QIN,QQSUM,DEPTH,VELMAG 
       END DO
      END DO
C
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SUB(L).EQ.1.) THEN
         FLOWX=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
         WRITE(94,946)FLOWX
         WRITE(95)FLOWX
       END IF  
       END DO
      END DO
C
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SVB(L).EQ.1.) THEN
         FLOWY=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
         WRITE(94,946)FLOWY
         WRITE(95)FLOWY
       END IF 
       END DO
      END DO
C
      IF (KC.GT.1) THEN
      DO K=KS,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       IF (SWB(L).EQ.1) THEN
         FLOWZ=-DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
         WRITE(94,946)FLOWZ
         WRITE(95)FLOWZ
       END IF  
       END DO
      END DO
      END IF
C
      CLOSE(94)
      CLOSE(95)
C----------------------------------------------------------------------C
C
  901 FORMAT(2I5,E12.4,4I5,E12.4)
  902 FORMAT(I5,2X,3E12.4,2I5)
  903 FORMAT(3E12.4,2I5)
  904 FORMAT(I5,2X,E12.4,10I5)
  905 FORMAT(I5)
  906 FORMAT(5E12.4)
  941 FORMAT(2I5,E12.4,4I5,E12.4)
  942 FORMAT(3E12.4,2I5)
  943 FORMAT(3E12.4,2I5)
  944 FORMAT(E12.4,10I5)
  945 FORMAT(I5)
  946 FORMAT(5E12.4)
C
C**********************************************************************C
C
      JSWASP=0
C
      RETURN
      END
