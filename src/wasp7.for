C
C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WASP7
C
c     mike mortons wasp5mrm	  
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 3 June 94
c ==========
c Revisions:
c ==========
c  M. Morton 06/06/94: this version writes dispersion to the WASPDH.OUT
c     hydrodynamic file instead of WASPB.OUT.  A modified version of
c     WASP5 is necessary to utilize these dispersions.
c  M. Morton 06/07/94: writes hydrodynamic information and dispersion to
c     an unformatted binary file, WASPDHU.OUT.  The correct files
c     for WASP data groups B, C, and D are:
c        Data Group B use WASPB.MRM  (do not use WASPB.OUT)
c        Data Group C use WASPC.OUT
c        Data Group D use WASPD.MRM  (do not use WASPD.OUT)
c===========
C
C **  SUBROUTINE WASP5 WRITES OUTPUT FILES PROVIDING ADVECTIVE AND
C **  DIFFUSIVE TRANSPORT FIELDS FOR THE WASP5 WATER QUALITY MODEL
C
C**********************************************************************C
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
C
      DIMENSION LUTMP((KCM+1)*LCM),LDTMP((KCM+1)*LCM),
     $          QTMP((KCM+1)*LCM)
      CHARACTER*50 TITLEB,TITLEC
      CHARACTER*12 HYDFIL
C
      TITLEB='DATA GROUP B: EXCHANGE COEFFICIENTS'
      TITLEC='DATA GROUP C: VOLUMES'
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
C **  READ CONTROL DATA FOR WRITING TO WASP COMPATIBLE FILES     
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
      READ(1,*) NRFLD,SCALR,CONVR,ISNKH
C
C3**  READ ADVECTION PARAMETERS
C
      READ(1,1)
      READ(1,1)
      READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,isdhd
C
C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
C
      READ(1,1)
      READ(1,1)
      READ(1,*) DEPSED,TDINTS,SEDIFF, wss1, wss2, wss3                      !mrm
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
C **  on page 11 of the wasp5.1 manual part B, Sept 1993
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
        WRITE(93,1031)IVOPT,IBEDV,TDINTS,TITLEC
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
           dmult=hlpf(l)*dzc(kc)
           VOLUME=DXYP(L)*HLPF(L)*DZC(KC)
           IF(NTSMMT.LT.NTSPTC) THEN
             dmult=HMP(l)*dzc(kc)
             VOLUME=DXYP(L)*HMP(l)*DZC(KC)
           END IF
           WRITE(90,1001)LWASP,KWASP,I,J,L,KC
           WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP,I,J,L,KC
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
            dmult=hlpf(l)*dzc(k)                                             !mrm
            VOLUME=DXYP(L)*HLPF(L)*DZC(K)
            IF(NTSMMT.LT.NTSPTC) THEN
              dmult=HMP(l)*dzc(kc)
              VOLUME=DXYP(L)*HMP(l)*DZC(KC)
            END IF
            WRITE(90,1001)LWASP,KWASP,I,J,L,K
            WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP,I,J,L,KC
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
         dmult=hlpf(l)*dzc(1)                                               !mrm
         VOLUME=DXYP(L)*HLPF(L)*DZC(1)
         IF(NTSMMT.LT.NTSPTC) THEN
              dmult=HMP(l)*dzc(kc)
              VOLUME=DXYP(L)*HMP(l)*DZC(KC)
         END IF
         IONE=1
         WRITE(90,1001)LWASP,KWASP,I,J,L,IONE
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP,I,J,L,ione
        END DO
        LTYPE=3
        KWASP=KC+1
         DXYSUM=0.
         LWSPTMP=LWASP+1
         DO LT=2,LALT
         LWSPTMP=LWSPTMP+1
         END DO
c The following the lower benthic layer.  All upper benthic layer segments   !mrm
c have this layer immediately below them:                                    !mrm
         DO LT=2,LALT
         LWASP=LWASP+1
         LBELOW=LWSPTMP
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         DXYSUM=DXYSUM+DXYP(L)
         VOLUME=DXYP(L)*DEPSED
         IZERO=0
         WRITE(90,1001)LWASP,KWASP,I,J,L,IZERO
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                depsed,DEXP,I,J,L,izero                                            !mrm
c     $                DMULT,DEXP                                            !mrm
         END DO
c Next do the lower benthic layer:                                           !mrm
        LTYPE=4
        KWASP=KC+2
         LWASP=LWASP+1
         LBELOW=0
         dmult=depsed                                                       !mrm
         VOLUME=DXYSUM*DEPSED
         IM1=-1
         WRITE(90,1001)LWASP,KWASP,I,J,L,IM1
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP,I,J,L,IM1
        CLOSE(90)
        CLOSE(93)
      END IF
C
 1001 FORMAT(6I5,2F10.4)
 1031 FORMAT(2I5,F10.4,10X,A50)
 1032 FORMAT(2F10.4)
c FORMAT 1033 as commented out is troublesome ... better change shown
c1033 FORMAT(3I10,5F10.2)
 1033 FORMAT(3I10,F10.1,4F10.3,'   !',4i5)
C
C**********************************************************************C
C
C **  WRITE DIFFUSIVE AND DISPERSIVE TRANSPORT FILE waspb.out 
C
C **  file waspb.out is consistent with data group B specifications          !mrm
C **  on page 8 of the wasp5.1 manual part B, Sept 1993
C   
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1) THEN
        OPEN(91,FILE='waspb.out',STATUS='UNKNOWN')
        CLOSE(91,STATUS='DELETE')
        OPEN(91,FILE='waspb.out',STATUS='UNKNOWN')
c       NRFLD=1
        WRITE(91,1011)NRFLD,TITLEB
        NTEX=NTS/NTSMMT
c       SCALR=1.
c       CONVR=1.
        WRITE(91,1012)NTEX,SCALR,CONVR
        CLOSE(91)
c      END IF                                                                !mrm
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
         NORSH=NORSH+INT(SUBO(L))+INT(SVBO(L))
         NORSV=NORSV+INT(SPB(L))
       END DO
      NORS=ISNKH*KC*NORSH+KS*NORSV
      WRITE(91,1013)NORS
C
      IF(ISNKH.EQ.1) THEN
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
     $         +HLPF(LW))*DXIU(L)
              WRITE(91,1014) ADDLW,UNITY,LWASPW,LWASP
            END IF
          END DO
        END DO
      END IF
C
      IF(ISNKH.EQ.1) THEN
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
     $         +HLPF(LS))*DYIV(L)
              WRITE(91,1014) ADDLS,UNITY,LWASPS,LWASP
            END IF
          END DO
        END DO
      END IF
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
              WRITE(91,1014) ADDL,UNITY,LWASP,LBELOW
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
C **  ADD PORE WATER EXCHANGE FIELD ON LAST CALL
C
c      IF (N.GE.NTS) THEN                                                    !mrm
        OPEN(91,FILE='waspb.out',ACCESS='APPEND',STATUS='UNKNOWN')
C
        NTEX=1
        SCALR=1.
        CONVR=1.
        WRITE(91,1012)NTEX,SCALR,CONVR
        NORSV=0  
        DO LT=2,LALT
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        NORSV=NORSV+INT(SPB(L))
        END DO
        WRITE(91,1013)NORSV
        IF (KC.GE.1) THEN
          KMUL2=KC+1
          UNITY=1.
          DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          IF (SPB(L).EQ.1.) THEN
            LWASP=LT-1+KC*LCLTM2
            LBELOW=LT-1+KMUL2*LCLTM2
            ADDL=2.*DXYP(L)*SEDIFF/DEPSED
            WRITE(91,1014) ADDL,UNITY,LWASP,LBELOW
          END IF 
          END DO
        END IF
C
        NBRK=6
        WRITE(91,1015)NBRK
C
        TSTART=TBEGIN*TCON
        TSTOP=DT*FLOAT(NTS)+TBEGIN*TCON
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
        IBPTMP=0
        WRITE(91,1017)IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     $                IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     $                IBPTMP,IBPTMP,IBPTMP,IBPTMP,
     $                IBPTMP,IBPTMP,IBPTMP,IBPTMP
C        
        CLOSE(91)
c      END IF                                                                !mrm
      end if                                                                 !mrm
C
 1011 FORMAT(I5,10X,A50)
 1012 FORMAT(I5,2F10.4)
 1013 FORMAT(I5)
c1014 FORMAT(2F10.0,2I5,'   !',3i5,3x,a3)
 1014 FORMAT(2E10.3,2I5,F10.3,'   !',3i5,3x,a3)
 1015 FORMAT(I5)
 1016 FORMAT(4(e10.3,F10.5))
 1017 FORMAT(16I5)
C
C**********************************************************************C
C
C **  WRITE ADVECTIVE TRANSPORT FILE waspd.out    
C
C **  file waspd.out is consistent with data group D.1 specifications 
C **  on page 13 of the wasp5.1 manual part B, Sept 1993
C **  this file is written only if ISWASPD=1
C   
C----------------------------------------------------------------------C
C
C!!!!!!!!!!CHANGES ON NEXT 2 LINES
      IF(ISWASPD.EQ.1) THEN
C
       IF(JSWASP.EQ.1) THEN
        OPEN(92,FILE='waspd.out',STATUS='UNKNOWN')
        CLOSE(92,STATUS='DELETE')
        OPEN(92,FILE='waspd.out',STATUS='UNKNOWN')
c       IQOPT=1
c       NFIELD=1
        WRITE(92,1021)IQOPT,NFIELD,HYDFIL
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
C!!!!!!!!!!!!!!CHANGES ON NEXT 3 LINES
       NOQSH=NOQSH+INT(SUBO(L))+INT(SVBO(L))
cjh       IF (IJCTLT(I+1,J).EQ.6) NOQSH=NOQSH+1
cjh       IF (IJCTLT(I,J+1).EQ.6) NOQSH=NOQSH+1
       IF (IJCTLT(I+1,J).EQ.8) NOQSH=NOQSH+1
       IF (IJCTLT(I,J+1).EQ.8) NOQSH=NOQSH+1
       NOQSV=NOQSV+INT(SWB(L))
       END DO
      NOQS=KC*NOQSH+KS*NOQSV
      WRITE(92,1023)NOQS
C
      LL=0
C
      DO K=KC,1,-1
      KMUL=KC-K
c
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C!!!!!!!!!!!!!CHANGES ON NEXT 15 LINES
       IF (SUBO(L).EQ.1.) THEN
         LL=LL+1
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LDTMP(LL)-1
cjh         IF (IJCTLT(I-1,J).EQ.6) LUTMP(LL)=0
         IF (IJCTLT(I-1,J).EQ.8) LUTMP(LL)=0
         QTMP(LL)=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
       END IF
cjh       IF (IJCTLT(I+1,J).EQ.6) THEN
       IF (IJCTLT(I+1,J).EQ.8) THEN
         IF (SUBO(L+1).EQ.1.) THEN
           LL=LL+1
           LDTMP(LL)=0
           LUTMP(LL)=LT-1+KMUL*LCLTM2
           QTMP(LL)=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
         END IF
       END IF
       END DO
C
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
C!!!!!!!!!!!!!CHANGES ON NEXT 16 LINES
       IF (SVBO(L).EQ.1.) THEN
         LL=LL+1
         LSLT=LSCLT(LT)
         LDTMP(LL)=LT-1+KMUL*LCLTM2
         LUTMP(LL)=LSLT-1+KMUL*LCLTM2
cjh         IF (IJCTLT(I,J-1).EQ.6) LUTMP(LL)=0
         IF (IJCTLT(I,J-1).EQ.8) LUTMP(LL)=0
         QTMP(LL)=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
       END IF
cjh       IF (IJCTLT(I,J+1).EQ.6) THEN
       IF (IJCTLT(I,J+1).EQ.8) THEN
         LN=LNC(L)
         IF (SVBO(LN).EQ.1) THEN
           LL=LL+1
           LDTMP(LL)=0
           LUTMP(LL)=LT-1+KMUL*LCLTM2
           QTMP(LL)=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
         END IF
       END IF
       END DO
C
      END DO
C
      IF (KC.GT.1) THEN
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
      DO L=1,LL,4
      WRITE(92,1024) QTMP(L),  LUTMP(L),  LDTMP(L),
     $               QTMP(L+1),LUTMP(L+1),LDTMP(L+1),
     $               QTMP(L+2),LUTMP(L+2),LDTMP(L+2),
     $               QTMP(L+3),LUTMP(L+3),LDTMP(L+3)
      END DO
C
      NBRKQ=6
      WRITE(92,1025)NBRKQ
      WRITE(92,1026)D1,T1,D2,T2,D3,T3,D4,T4
      WRITE(92,1026)D5,T5,D6,T6
C
      CLOSE(92)
C!!!!!!!!!!CHANGES ON NEXT 2 LINES
C
      END IF
C
 1021 FORMAT(2I5,A12)
 1022 FORMAT(I5,2F10.4)
 1023 FORMAT(I5)
c1024 FORMAT(4(F10.0,2I5))
 1024 FORMAT(1p,4(E10.3,2I5))
 1025 FORMAT(I5)
 1026 FORMAT(4(2F10.5))

C**********************************************************************C     !mrm
c M.R. Morton's version of WASP Data Group D                                 !mrm
C **  WRITE ADVECTIVE TRANSPORT FILE waspd.mrm                               !mrm
C----------------------------------------------------------------------C     !mrm
      if (jswasp .eq. 1) then                                                !mrm
        open(92,file='waspd.mrm',status='unknown')                           !mrm
        write(92,2020) iqopt,nfield,hydfil                                   !mrm
        ll=0                                                                 !mrm
        ninq=0                                                               !mrm
        scalq=1.0                                                            !mrm
        convq=1.0/86400.0                                                    !mrm
c Data Block D.1 (Advective Flows) is not needed since HYD file is used:     !mrm
c        write(92,2021) ninq,scalq,convq                                     !mrm
c Data Block D.2 (Pore Water Flows) not needed:                              !mrm
        write(92,2022) ninq,scalq,convq                                      !mrm
c Data Block D.3 (Sediment #1 Transport Field):                              !mrm
        ninq=1                                                               !mrm
        write(92,2023) ninq,scalq,convq                                      !mrm
        if (kc.gt.1) then                                                    !mrm
          do k=ks,0,-1                                                       !mrm
            kmul1=ks-k                                                       !mrm
            kmul2=kmul1+1                                                    !mrm
            do lt=2,lalt                                                     !mrm
c              write(6,6999)k,kmul1,kmul2,lt
c              call F_FlushNow(6) 
              i=illt(lt)                                                     !mrm
c              write(6,6999)lt,i 
c              call F_FlushNow(6) 
              j=jllt(lt)                                                     !mrm
c              write(6,6999)lt,j
c              call F_FlushNow(6) 
              l=lij(i,j)                                                     !mrm
c              write(6,6999)k,kmul1,kmul2,lt,i,j,l,IL(L),JL(L),swb(L) 
c              call F_FlushNow(6) 
              if (swb(l).eq.1.) then                                         !mrm
                ll=ll+1                                                      !mrm
                lutmp(ll)=lt-1+kmul1*lcltm2                                  !mrm
                ldtmp(ll)=lt-1+kmul2*lcltm2                                  !mrm
c qtmp array will hold the plan view area of each cell:                      !mrm
                qtmp(ll)= dxyp(l)                                            !mrm
              end if                                                         !mrm
            end do                                                           !mrm
          end do                                                             !mrm
        end if                                                               !mrm
c                                                                            !mrm                                               
 6999 format(9i5,f5.1)                                                 
 6996 format(9i5,f5.1)                                                 
c                                                               
        write(92,2030) ll                                                    !mrm
        do l=1,ll,4                                                          !mrm
          write(92,1024) qtmp(l),  lutmp(l),  ldtmp(l),                      !mrm
     $                 qtmp(l+1),lutmp(l+1),ldtmp(l+1),                      !mrm
     $                 qtmp(l+2),lutmp(l+2),ldtmp(l+2),                      !mrm
     $                 qtmp(l+3),lutmp(l+3),ldtmp(l+3)                       !mrm
        end do                                                               !mrm
        nbrkq=2                                                              !mrm
        t1=1.0                                                               !mrm
        t2=366.0                                                             !mrm
        write(92,2030) nbrkq                                                 !mrm
        write(92,2031) wss1,t1,wss1,t2                                       !mrm
c Data Block D.4 (Sediment #2 Transport Field):                              !mrm
        ninq=1                                                               !mrm
        write(92,2024) ninq,scalq,convq                                      !mrm
        write(92,2030) ll                                                    !mrm
        do l=1,ll,4                                                          !mrm
          write(92,1024) qtmp(l),  lutmp(l),  ldtmp(l),                      !mrm
     $                 qtmp(l+1),lutmp(l+1),ldtmp(l+1),                      !mrm
     $                 qtmp(l+2),lutmp(l+2),ldtmp(l+2),                      !mrm
     $                 qtmp(l+3),lutmp(l+3),ldtmp(l+3)                       !mrm
        end do                                                               !mrm
        nbrkq=2                                                              !mrm
        t1=1.0                                                               !mrm
        t2=366.0                                                             !mrm
        write(92,2030) nbrkq                                                 !mrm
        write(92,2031) wss2,t1,wss2,t2                                       !mrm
c Data Block D.5 (Sediment #3 Transport Field):                              !mrm
        ninq=1                                                               !mrm
        write(92,2025) ninq,scalq,convq                                      !mrm
        write(92,2030) ll                                                    !mrm
        do l=1,ll,4                                                          !mrm
          write(92,1024) qtmp(l),  lutmp(l),  ldtmp(l),                      !mrm
     $                 qtmp(l+1),lutmp(l+1),ldtmp(l+1),                      !mrm
     $                 qtmp(l+2),lutmp(l+2),ldtmp(l+2),                      !mrm
     $                 qtmp(l+3),lutmp(l+3),ldtmp(l+3)                       !mrm
        end do                                                               !mrm
        nbrkq=2                                                              !mrm
        t1=1.0                                                               !mrm
        t2=366.0                                                             !mrm
        write(92,2030) nbrkq                                                 !mrm
        write(92,2031) wss3,t1,wss3,t2                                       !mrm
c add system bypass array to bottom of data group D:                         !mrm
        write(92,1017)ibptmp,ibptmp,ibptmp,ibptmp,                           !mrm
     +                ibptmp,ibptmp,ibptmp,ibptmp,                           !mrm
     +                ibptmp,ibptmp,ibptmp,ibptmp,                           !mrm
     +                ibptmp,ibptmp,ibptmp,ibptmp                            !mrm
        close(92)                                                            !mrm
      end if                                                                 !mrm
 2020 format(2i5,a12,'    Data Group D: Flows')                              !mrm
 2021 FORMAT(1p,I5,2e10.3,'    Data Block D.1 Advective Flows')              !mrm
 2022 FORMAT(1p,I5,2e10.3,'    Data Block D.2 Pore Water Flows')             !mrm
 2023 FORMAT(1p,I5,2e10.3,'    Data Block D.3 Sed. #1 Transport Field')      !mrm
 2024 FORMAT(1p,I5,2e10.3,'    Data Block D.4 Sed. #2 Transport Field')      !mrm
 2025 FORMAT(1p,I5,2e10.3,'    Data Block D.5 Sed. #3 Transport Field')      !mrm
 2030 format(i5)                                                             !mrm
 2031 format(2(e10.3,f10.5))                                                 !mrm

C
C**********************************************************************C
C
C **  WRITE TO EXTERNAL HYDRO FILE waspdh.out AND DIAGNOSTIC VERSION
C **  OF SAME FILE waspdhd.out 
C
C----------------------------------------------------------------------C
C
      IF(JSWASP.EQ.1)THEN
        OPEN(90,FILE='waspdhd.out',STATUS='UNKNOWN')
        if(iqopt.eq.3) OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')           !mrm
        if(iqopt.eq.4) OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',          !mrm
     $                      FORM='UNFORMATTED')
        open(96,file='waspb.mrm',status='unknown')                           !mrm
        CLOSE(90,STATUS='DELETE')
        if(iqopt.eq.3) CLOSE(94,STATUS='DELETE')                             !mrm
        if(iqopt.eq.4) CLOSE(95,STATUS='DELETE')                             !mrm
        close(96,status='delete')                                            !mrm
        OPEN(90,FILE='waspdhd.out',STATUS='UNKNOWN')
        if(iqopt.eq.3) OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')           !mrm
        if(iqopt.eq.4) OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',          !mrm
     $                      FORM='UNFORMATTED')
        open(96,file='waspb.mrm',status='unknown')                           !mrm
        write(96,1011) nrfld,titleb                                          !mrm
        ntexx=1                                                              !mrm
        write(96,1012) ntexx,scalr,convr                                     !mrm

c Write WASP5 Hydrodynamic File Data Record 1, Data Options:                 !mrm
c  NJUN = number of segments connected by flows from the hyd. file           !mrm
c  NCHN = number of interfacial flow pairs from the hyd. file                !mrm
c  DTWASP = WASP5 time step (seconds)                                        !mrm
c  TZERO = begin time step for hyd. file (seconds)                           !mrm
c  TENDHYD = end time step for hyd. file (seconds)                           !mrm
c  ISTMP = control switch, 0=time variable segment depths and velocities     !mrm
c          are read; 1=time variable segment depths and velocities are not   !mrm
c          read.                                                             !mrm
c
C       KCLC=KC*LCLT
C       LCLTM2=LCLT-2
C        DO KL=1,KCLC 
C        NCHNC(KL)=0
C        END DO
C        DO M=1,10
C         DO KL=1,KCLC 
C         LCHNC(KL,M)=0
C         END DO
C        END DO
        NJUN=KC*(LCLT-2)
        NCHNH=0
        NCHNV=0
C!!!!!!!!CHANGES NEXT 13 LINES
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          NCHNH=NCHNH+INT(SUBO(L))
cjh          IF (IJCTLT(I+1,J).EQ.6) THEN
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (SUBO(L+1).EQ.1.) NCHNH=NCHNH+1
          END IF 
          NCHNH=NCHNH+INT(SVBO(L))
cjh          IF (IJCTLT(I,J+1).EQ.6) THEN
          IF (IJCTLT(I,J+1).EQ.8) THEN
            IF (SVBO(LNC(L)).EQ.1.) NCHNH=NCHNH+1
          END IF 
          NCHNV=NCHNV+INT(SWB(L))
        END DO
        NCHN=KC*NCHNH+KS*NCHNV
        ISTMP=0
        NODYN=NFLTMT
        nodyn=nodyn                                                          !mrm
        dtwasp = dt * float(ntsmmt)                                          !mrm
        TZERO=TBEGIN*TCON
        TENDHYD=TZERO+NTS*DT
        WRITE(90,901)NJUN,NCHN
        if(iqopt.eq.3) then                                                  !mrm
          WRITE(94,941) NJUN,NCHN, dtwasp, TZERO,TENDHYD,ISTMP               !mrm
c        WRITE(94,941) NJUN,NCHN,DT,TZERO,TENDHYD,ISTMP
        end if                                                               !mrm
        if(iqopt.eq.4) then                                                  !mrm
          WRITE(95) njun,nchn, dtwasp, tzero,tendhyd,istmp                   !mrm
        end if
        write(96,1013) nchn                                                  !mrm
C
C **  CHANNEL DATA
C
c Write WASP5 Hydrodynamic File Data Record 2, Segment Interface Pairs:      !mrm
c   WASP expects to see boundary segments designated as "0".                 !mrm
c
        RMNDUM=0.
        LCHN=0
        DO K=KC,1,-1
          KMUL=KC-K
C!!!!!!!!!!!!!!CHANGES ON NEXT 38 LINES
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SUBO(L).EQ.1.) THEN
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LDTM-1
cjh              IF (IJCTLT(I-1,J).EQ.6) LUTM=0
              IF (IJCTLT(I-1,J).EQ.8) LUTM=0
              RLENTH=DXU(L)
              WIDTH=DYU(L)
              LCHN=LCHN+1
C             NCHNC(LDTM)=NCHNC(LDTM)+1
C             NCHNC(LUTM)=NCHNC(LUTM)+1
C             LCHNC(LDTM,NCHNC(LDTM))=LCHN
C             LCHNC(LUTM,NCHNC(LUTM))=LCHN
              if (isdhd. eq. 1) WRITE(90,902)LCHN,RLENTH,WIDTH,
     +          RMNDUM,LUTM,LDTM
              if (isdhd .eq. 2) write(90,'(2i5)') lutm,ldtm
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              if(iqopt.eq.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'u 0'                         !mrm
              if(iqopt.eq.4) WRITE(95) LUTM,LDTM                             !mrm
              write(96,1014) unity,unity,lutm,ldtm,unity,I,J,K,'u 0'                     !mrm
            END IF
cjh            IF (IJCTLT(I+1,J).EQ.6) THEN
            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (SUBO(L+1).EQ.1.) THEN
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                RLENTH=DXU(L+1)
                WIDTH=DYU(L+1)
                LCHN=LCHN+1
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                if (isdhd .eq. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +            RMNDUM,LUTM,LDTM
              if (isdhd .eq. 2) write(90,'(2i5)') lutm,ldtm
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                if(iqopt.eq.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'u+1'                       !mrm
                if(iqopt.eq.4) WRITE(95) LUTM,LDTM                           !mrm
                unity=1.0
                write(96,1014) unity,unity,lutm,ldtm,unity,I,J,K,'u+1'                   !mrm
              END IF
            END IF
          END DO
c        END DO
c        DO K=KC,1,-1
c        KMUL=KC-K
C!!!!!!!!CHANGES NEXT 41 LINES
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SVBO(L).EQ.1.) THEN
              LSLT=LSCLT(LT)
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LSLT-1+KMUL*LCLTM2
cjh              IF(IJCTLT(I,J-1).EQ.6) LUTM=0
              IF(IJCTLT(I,J-1).EQ.8) LUTM=0
              RLENTH=DYV(L)
              WIDTH=DXV(L)
              LCHN=LCHN+1
C             NCHNC(LDTM)=NCHNC(LDTM)+1
C             NCHNC(LUTM)=NCHNC(LUTM)+1
C             LCHNC(LDTM,NCHNC(LDTM))=LCHN
C             LCHNC(LUTM,NCHNC(LUTM))=LCHN
              if (isdhd .eq. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +          RMNDUM,LUTM,LDTM
              if (isdhd .eq. 2) write(90,'(2i5)') lutm,ldtm
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              if(iqopt.eq.3) WRITE(94,9941) LUTM,LDTM,I,J,K,'v 0'                         !mrm
              if(iqopt.eq.4) WRITE(95) LUTM,LDTM                             !mrm
              write(96,1014) unity,unity,lutm,ldtm,unity,I,J,K,'v 0'                     !mrm
            END IF
c            IF (IJCTLT(I,J+1).EQ.6) THEN
            IF (IJCTLT(I,J+1).EQ.8) THEN
              LN=LNC(L)
              IF (SVBO(LN).EQ.1.) THEN
                LSLT=LSCLT(LT)
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
                RLENTH=DYV(LN)
                WIDTH=DXV(LN)
                LCHN=LCHN+1
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                if (isdhd .eq. 1) WRITE(90,902) LCHN,RLENTH,WIDTH,
     +            RMNDUM,LUTM,LDTM
                if (isdhd .eq. 2) write(90,'(2i5)') lutm,ldtm
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                if(iqopt.eq.3) WRITE(94,9941) LUTM,LDTM,unity,I,
     +                                       J,K,'v+1'                       !mrm
                if(iqopt.eq.4) WRITE(95) LUTM,LDTM                           !mrm
                write(96,1014) unity,unity,lutm,ldtm,unity,I,J,K,'v+1'                   !mrm
              END IF
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
C               NCHNC(LDTM)=NCHNC(LDTM)+1
C               NCHNC(LUTM)=NCHNC(LUTM)+1
C               LCHNC(LDTM,NCHNC(LDTM))=LCHN
C               LCHNC(LUTM,NCHNC(LUTM))=LCHN
                WRITE(90,902)LCHN,RLENTH,WIDTH,RMNDUM,LUTM,LDTM
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                if(iqopt.eq.3) WRITE(94,9941)LUTM,LDTM,I,J,K,'w 0'                        !mrm
                if(iqopt.eq.4) WRITE(95) LUTM,LDTM                           !mrm
                write(96,1014) unity,unity,lutm,ldtm,unity,I,J,K,'w 0'                   !mrm
              END IF
            END DO
          END DO
c
c write out time series of zero dispersion coefficients:
c
          d1=0.0                                                             !mrm
          t1=tzero/tcon                                                      !mrm
          d2=0.0                                                             !mrm
          t2=tendhyd/tcon                                                    !mrm
          nbrkq=2                                                            !mrm
          write(96,905) nbrkq                                                !mrm
          write(96,1016) d1,t1, d2,t2                                        !mrm
c
c For exchange between the lower water surface layer and the upper           !mrm
c benthic layer, do the following:                                           !mrm
c
          write(96,1012) ntexx,scalr,convr                                   !mrm
          ntexx=0
          DO K=1,1                                                           !mrm
            DO LT=2,LALT                                                     !mrm
              I=ILLT(LT)                                                     !mrm
              J=JLLT(LT)                                                     !mrm
              L=LIJ(I,J)                                                     !mrm
              IF (SWB(L).EQ.1.) THEN                                         !mrm
                ntexx=ntexx+1                                                !mrm
              END IF
            END DO
          END DO
          write(96,1013) ntexx                                               !mrm
          DO K=1,1                                                           !mrm
            KMUL1=KS-K                                                       !mrm
            KMUL2=KMUL1+1                                                    !mrm
            kmul3=kmul2+1                                                    !mrm
            DO LT=2,LALT                                                     !mrm
              I=ILLT(LT)                                                     !mrm
              J=JLLT(LT)                                                     !mrm
              L=LIJ(I,J)                                                     !mrm
              IF (SWB(L).EQ.1.) THEN                                         !mrm
                LUTM=LT-1+KMUL2*LCLTM2                                       !mrm
                LDTM=LT-1+kmul3*LCLTM2                                       !mrm
                write(96,1014) dxyp(l),depsed,lutm,ldtm                      !mrm
              END IF
            END DO
          END DO
c
c write out time series of water-benthic exchange dispersion coefficients:   !mrm
c
          d1=sediff                                                          !mrm
          t1=tzero/tcon                                                      !mrm
          d2=sediff                                                          !mrm
          t2=tendhyd/tcon                                                    !mrm
          nbrkq=2                                                            !mrm
          write(96,905) nbrkq                                                !mrm
          write(96,1016) d1,t1, d2,t2                                        !mrm
          ibptmp=0                                                           !mrm
          write(96,1017)ibptmp,ibptmp,ibptmp,ibptmp,                         !mrm
     +                  ibptmp,ibptmp,ibptmp,ibptmp,                         !mrm
     +                  ibptmp,ibptmp,ibptmp,ibptmp,                         !mrm
     +                  ibptmp,ibptmp,ibptmp,ibptmp                          !mrm
        END IF
C
C **  JUNCTION DATA WITH INITIAL CONDITIONS
C
c Write WASP5 Hydrodynamic File Data Record 3, Initial Segment Properties:   !mrm
c
        veltmp=0.
        dumvol=0.
        DO K=KC,1,-1
        KMUL=KC-K
         DO LT=2,LALT
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         LCELL=LT-1+KMUL*LCLTM2
         DEPTMP=HLPF(L)*DZC(K)
         VOLTMP=DEPTMP*DXYP(L)
         IF(NTSMMT.LT.NTSPTC) THEN
           DEPTMP=HMP(L)*DZC(K)
           VOLTMP=DEPTMP*DXYP(L)
         END IF
         if (isdhd .eq. 1) WRITE(90,904) LCELL,VOLTMP,I,J
         if(iqopt.eq.3) WRITE(94,9440) VOLTMP,DEPTMP,VELTMP           !mrm
         if(iqopt.eq.4) WRITE(95) voltmp,deptmp,veltmp                       !mrm
         END DO
        END DO
C
        CLOSE(90)
        if(iqopt.eq.3) CLOSE(94)                                             !mrm
        if(iqopt.eq.4) CLOSE(95)                                             !mrm
        close(96)                                                            !mrm
      END IF
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP, VOLUME AND FLOW DATA
C
      OPEN(90,FILE='waspdhd.out',ACCESS='APPEND',STATUS='UNKNOWN')
      if(iqopt.eq.3) then                                                    !mrm
        OPEN(94,FILE='waspdh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      end if                                                                 !mrm
      if(iqopt.eq.4) then                                                    !mrm
        OPEN(95,FILE='waspdhu.out',ACCESS='APPEND',STATUS='UNKNOWN',
     $                      FORM='UNFORMATTED')
      end if                                                                 !mrm
      LCLTM2=LCLT-2
      IZERO=0
      RZERO=0
      izero=izero                                                            !mrm
      rzero=rzero                                                            !mrm
C
C     NSTEP=N-NTSMMT
C     WRITE(94,945)NSTEP
C
C Write WASP5 Hydrodynamic File Data Record 4, BQ(J) flow in interface       !mrm
C pair "J":                                                                  !mrm
c
c Advection and dispersion in the X-direction:
c
      lchnum=0
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
c +++++ following lines by M. Morton to input dispersion to HYD file:        !mrm
       addlw=0.0                                                             !mrm
       if (sub(l).eq.1.) then                                                !mrm
         lw=l-1                                                              !mrm
         addlw=dyu(l)*ahulpf(l,k)*dzc(k)*0.5*(hlpf(l)                        !mrm
     $         +hlpf(lw))*dxiu(l)                                            !mrm
       end if                                                                !mrm
c +++++ above added by M. Morton                                             !mrm
C!!!!!!!!CHANGES NEXT 12 LINES
       IF (SUBO(L).EQ.1.) THEN
         TMPVAL=UHLPF(L,K)+SVPT*UVPT(L,K)
         FLOWX=DYU(L)*TMPVAL*DZC(K)
         UDDXTMP=2.*TMPVAL*DXIU(L)/(HLPF(L)+HLPF(L-1))
         IMTMP=I-1
         lchnum=lchnum+1
         IDRTMP=1
         if (isdhd .eq. 1) WRITE(90,944) FLOWX,IMTMP,I,J,K
         if(iqopt.eq.3) WRITE(94,9946) FLOWX,UDDXTMP,addlw,IDRTMP                    !mrm/JMH
c         WRITE(94,946) FLOWX                                                !mrm
         if(iqopt.eq.4) WRITE(95) FLOWX,UDDXTMP,addlw,IDRTMP                        !mrm/JMH
       END IF
cjh       IF (IJCTLT(I+1,J).EQ.6) THEN
       IF (IJCTLT(I+1,J).EQ.8) THEN
         IF (SUBO(L+1).EQ.1.) THEN
           TMPVAL=UHLPF(L+1,K)+SVPT*UVPT(L+1,K)
           FLOWX=DYU(L+1)*TMPVAL*DZC(K)
           UDDXTMP=2.*TMPVAL*DXIU(L+1)/(HLPF(L+1)+HLPF(L))
           IPTMP=I+1
           lchnum=lchnum+1
           IDRTMP=1
           if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWX,I,IPTMP,J,K
           if(iqopt.eq.3) WRITE(94,9946) FLOWX,UDDXTMP,addlw,IDRTMP                  !mrm/JMH
c           WRITE(94,946) FLOWX                                              !mrm
           if(iqopt.eq.4) WRITE(95) FLOWX,UDDXTMP,addlw,IDRTMP                      !mrm/JMH
         END IF
       END IF  
       END DO
c     END DO
C
c Advection and dispersion in the Y-direction:
c
c     DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
c +++++ following lines by M. Morton to input dispersion to HYD file:        !mrm
       addls=0.0                                                             !mrm
       if (svb(l).eq.1.) then                                                !mrm
         ls=lsc(l)                                                           !mrm
         addls=dxv(l)*ahvlpf(l,k)*dzc(k)*0.5*(hlpf(l)                        !mrm
     $         +hlpf(ls))*dyiv(l)                                            !mrm
       end if                                                                !mrm
c +++++ above added by M. Morton                                             !mrm
C!!!!!!!CHANGES NEXT 13 LINES
       IF (SVBO(L).EQ.1.) THEN
         TMPVAL=VHLPF(L,K)+SVPT*VVPT(L,K)
         FLOWY=DXV(L)*TMPVAL*DZC(K)
         VDDYTMP=2.*TMPVAL*DYIV(L)/(HLPF(L)+HLPF(LSC(L)))
         JMTMP=J-1
         lchnum=lchnum+1
         IDRTMP=2
         if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWY,I,JMTMP,J,K
         if(iqopt.eq.3) WRITE(94,9946) FLOWY,VDDYTMP,addls,IDRTMP                          !mrm
c         WRITE(94,946) FLOWY                                               !mrm
         if(iqopt.eq.4) WRITE(95) FLOWY,VDDYTMP,addls,IDRTMP                              !mrm
       END IF
cjh       IF (IJCTLT(I,J+1).EQ.6) THEN
       IF (IJCTLT(I,J+1).EQ.8) THEN
         LN=LNC(L)
         IF (SVBO(LN).EQ.1.) THEN
           TMPVAL=VHLPF(LN,K)+SVPT*VVPT(LN,K)
           FLOWY=DXV(LN)*TMPVAL*DZC(K)
           VDDYTMP=2.*TMPVAL*DYIV(LN)/(HLPF(LN)+HLPF(L))
           JPTMP=J+1
           lchnum=lchnum+1
           IDRTMP=2
           if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWY,I,J,JPTMP,K
           if(iqopt.eq.3) WRITE(94,9946) FLOWY,VDDYTMP,addls,IDRTMP                          !mrm
c           WRITE(94,946) FLOWY                                              !mrm
           if(iqopt.eq.4) WRITE(95) FLOWY,VDDYTMP,addls,IDRTMP                              !mrm
         END IF
       END IF 
       END DO
      END DO
C
c Advection and dispersion in the Z-direction:
c
      IF (KC.GT.1) THEN
      DO K=KS,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
c +++++ following lines by M. Morton to input dispersion to HYD file:        !mrm
       addl=0.0                                                              !mrm
       if (spb(l).eq.1.) then                                                !mrm
         addl=dxyp(l)*ablpf(l,k)*dzig(k)                                     !mrm
         if (isdhd .eq. 2) write(90, '(4i5,e13.4)') i,j,k,l,ablpf(l,k)
       end if                                                                !mrm
c +++++ above added by M. Morton                                             !mrm
       IF (SWB(L).EQ.1) THEN
         TMPVAL=WLPF(L,K)+SVPT*WVPT(L,K)
         FLOWZ=-DXYP(L)*TMPVAL
         WDDZTMP=TMPVAL*DZIG(K)/HLPF(L)
         KPTMP=K+1
         IDRTMP=3
         lchnum=lchnum+1
         if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWZ,I,J,K,KPTMP
         if(iqopt.eq.3) WRITE(94,9946) FLOWZ,WDDZTMP,addl,IDRTMP
c         WRITE(94,946) FLOWZ
         if(iqopt.eq.4) WRITE(95) FLOWZ,WDDZTMP,addl,IDRTMP                                !mrm
       END IF  
       END DO
      END DO
      END IF
C
C Write WASP5 Hydrodynamic File Data Record 5, Segment Properties:           !mrm
C
      qqsum=0.
      LCELTMP=0
      DO K=KC,1,-1
       DO LT=2,LALT
       LCELTMP=LCELTMP+1
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
       LN=LNC(L)
       VOLUM=DXYP(L)*HLPF(L)*DZC(K)
       IF(NTSMMT.LT.NTSPTC) VOLUM=DXYP(L)*HP(L)*DZC(K)
C      QIN=QSUMELPF(L)*DZC(K)
C      FLOWXI=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
C      FLOWYI=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
C      FLOWZI=DXYP(L)*(WLPF(L,K-1)+SVPT*WVPT(L,K-1))
C      FLOWXO=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
C      FLOWYO=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
C      FLOWZO=DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
C      QQSUM=QIN+FLOWXI+FLOWYI+FLOWZI-FLOWXO-FLOWYO-FLOWZO
       DEPTH=HLPF(L)*DZC(K)
C       IF(NTSMMT.LT.NTSPTC) DEPTH=HP(L)*DZC(K)
       VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K)
     $         +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
       VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K)
     $         +VHLPF(LN,K)+SVPT*VVPT(LN,K))/HLPF(L)
       VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1)
     $         +WLPF(L,K)+SVPT*WVPT(L,K))
       VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
       if (isdhd .eq. 1) WRITE(90,902) LCELTMP,VOLUM,I,J,K
C       if(iqopt.eq.3) WRITE(94,946) VOLUM,qqsum,DEPTH,VELMAG                 !mrm
       if(iqopt.eq.3) WRITE(94,946) VOLUM,DEPTH,VELMAG                        !mrm
       if(iqopt.eq.4) WRITE(95) VOLUM, DEPTH, VELMAG                         !mrm
       END DO
      END DO
C
      CLOSE(90)
      if(iqopt.eq.3) CLOSE(94)                                               !mrm
      if(iqopt.eq.4) CLOSE(95)                                               !mrm
C
C----------------------------------------------------------------------C
C
  901 FORMAT(2I5,E12.5,4I5,E12.5)
  902 FORMAT(I5,2X,3F20.8,3I5)
  903 FORMAT(3E12.5,2I5)
  904 FORMAT(I5,2X,F20.8,10I5)
  905 FORMAT(I5)
  906 FORMAT(5E12.5)
  941 FORMAT(2I5,3F20.8,I5)
  942 FORMAT(3E12.5,2I5)
  943 FORMAT(3E12.5,2I5)
  944 FORMAT(I5,2x,F20.8,10I5)
 9440 FORMAT(4F20.8)
  945 FORMAT(I5)
  946 FORMAT(4e17.9)
 9946 FORMAT(3e17.9,I5)
 9941 FORMAT(2I5,'    !',3i5,3x,a3)
C
C**********************************************************************C
C
      JSWASP=0
C
      RETURN
      END
