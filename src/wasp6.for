C
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WASP6
C
C **  LAST MODIFIED BY JOHN HAMRICK ON 10 april 1999
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
      READ(1,*) DEPSED,TDINTS,SEDIFF, wss1, wss2, wss3
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
           dmult=HLPF(l)*DZC(kc)                             
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
            dmult=HLPF(l)*DZC(k)                              
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
         dmult=HLPF(l)*DZC(1)                                
         VOLUME=DXYP(L)*HLPF(L)*DZC(1)
         WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
        END DO
        LTYPE=3
        KWASP=KC+1
         DXYSUM=0.
         LWSPTMP=LWASP+1
         DO LT=2,LALT
         LWSPTMP=LWSPTMP+1
         END DO
c The following the lower benthic layer.  All upper benthic layer segments 
c have this layer immediately below them:                     
         DO LT=2,LALT
         LWASP=LWASP+1
         LBELOW=LWSPTMP
         I=ILLT(LT)
         J=JLLT(LT)
         L=LIJ(I,J)
         DXYSUM=DXYSUM+DXYP(L)
         VOLUME=DXYP(L)*DEPSED
         WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                depsed,DEXP                             
c     $                DMULT,DEXP                             
         END DO
c Next do the lower benthic layer:                            
        LTYPE=4
        KWASP=KC+2
         LWASP=LWASP+1
         LBELOW=0
         dmult=depsed                                        
         VOLUME=DXYSUM*DEPSED
         WRITE(90,1001)LWASP,KWASP,I,J,DLON(L),DLAT(L)
         WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP,
     $                DMULT,DEXP
        CLOSE(90)
        CLOSE(93)
      END IF
C
 1001 FORMAT(4I5,2F10.4)
 1031 FORMAT(2I5,F10.4,10X,A50)
 1032 FORMAT(2F10.4)
c FORMAT 1033 as commented out is troublesome ... better change shown
c1033 FORMAT(3I10,5F10.2)
 1033 FORMAT(3I10,F10.1,4F10.3)
C
C**********************************************************************C
C
C **  WRITE DIFFUSIVE AND DISPERSIVE TRANSPORT FILE waspb.out 
C
C **  file waSPB.out is consistent with data group B specifications  
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
c      END IF                                                 
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
c      IF (N.GE.NTS) THEN                                     
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
c      END IF                                                 
      end if                                                  
C
 1011 FORMAT(I5,10X,A50)
 1012 FORMAT(I5,2F10.4)
 1013 FORMAT(I5)
c1014 FORMAT(2F10.0,2I5)
 1014 FORMAT(2E10.3,2I5)
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
 1024 FORMAT(4(E10.3,2I5))
 1025 FORMAT(I5)
 1026 FORMAT(4(2F10.5))

C**********************************************************************C
c M.R. Morton's version of WASP Data Group D                  
C **  WRITE ADVECTIVE TRANSPORT FILE waspd.mrm                           
C----------------------------------------------------------------------C
      if (JSWASP .eq. 1) then                                   
        open(92,file='waspd.mrm',status='unknown')              
        write(92,2020) IQOPT,NFIELD,HYDFIL                      
        NINQ=0                                                  
        SCALQ=1.0                                               
        CONVQ=1.0/86400.0                                       
c Data Block D.1 (Advective Flows) is not needed since HYD file is used: 
c        write(92,2021) NINQ,SCALQ,CONVQ                        
c Data Block D.2 (Pore Water Flows) not needed:                 
        write(92,2022) NINQ,SCALQ,CONVQ                         
c Data Block D.3 (Sediment #1 Transport Field):                 
        NINQ=1                                                  
        write(92,2023) NINQ,SCALQ,CONVQ                         
        if (KC.gt.1) then                                       
          do K=KS,0,-1                                          
            KMUL1=KS-K                                          
            KMUL2=KMUL1+1                                       
            do lt=2,LALT                                        
              i=ILLT(lt)                                        
              j=JLLT(lt)                                        
              l=LIJ(i,j) 
              if (SWB(l).eq.1.) then                            
                ll=ll+1                                         
                LUTMP(ll)=lt-1+KMUL1*LCLTM2                     
                LDTMP(ll)=lt-1+KMUL2*LCLTM2                     
c QTMP array will hold the plan view area of each cell:         
                QTMP(ll)= DXYP(l)                               
              end if                                            
            end do                                              
          end do                                                
        end if
c 
        write(92,2030) ll                                       
        do l=1,ll,4                                             
          write(92,1024) QTMP(l),  LUTMP(l),  LDTMP(l),         
     $                 QTMP(l+1),LUTMP(l+1),LDTMP(l+1),         
     $                 QTMP(l+2),LUTMP(l+2),LDTMP(l+2),         
     $                 QTMP(l+3),LUTMP(l+3),LDTMP(l+3)          
        end do                                                  
        nbrkq=2                                                 
        t1=1.0                                                  
        t2=366.0                                                
        write(92,2030) nbrkq                                    
        write(92,2031) wss1,t1,wss1,t2                          
c Data Block D.4 (Sediment #2 Transport Field):                 
        NINQ=1                                                  
        write(92,2024) NINQ,SCALQ,CONVQ                         
        write(92,2030) ll                                       
        do l=1,ll,4                                             
          write(92,1024) QTMP(l),  LUTMP(l),  LDTMP(l),         
     $                 QTMP(l+1),LUTMP(l+1),LDTMP(l+1),         
     $                 QTMP(l+2),LUTMP(l+2),LDTMP(l+2),         
     $                 QTMP(l+3),LUTMP(l+3),LDTMP(l+3)          
        end do                                                  
        nbrkq=2                                                 
        t1=1.0                                                  
        t2=366.0                                                
        write(92,2030) nbrkq                                    
        write(92,2031) wss2,t1,wss2,t2                          
c Data Block D.5 (Sediment #3 Transport Field):                 
        NINQ=1                                                  
        write(92,2025) NINQ,SCALQ,CONVQ                         
        write(92,2030) ll                                       
        do l=1,ll,4                                             
          write(92,1024) QTMP(l),  LUTMP(l),  LDTMP(l),         
     $                 QTMP(l+1),LUTMP(l+1),LDTMP(l+1),         
     $                 QTMP(l+2),LUTMP(l+2),LDTMP(l+2),         
     $                 QTMP(l+3),LUTMP(l+3),LDTMP(l+3)          
        end do                                                  
        nbrkq=2                                                 
        t1=1.0                                                  
        t2=366.0                                                
        write(92,2030) nbrkq                                    
        write(92,2031) wss3,t1,wss3,t2                          
c add system bypass array to bottom of data group D:            
        write(92,1017)ibptmp,ibptmp,ibptmp,ibptmp,              
     +                ibptmp,ibptmp,ibptmp,ibptmp,              
     +                ibptmp,ibptmp,ibptmp,ibptmp,              
     +                ibptmp,ibptmp,ibptmp,ibptmp               
        close(92)                                               
      end if                                                    
 2020 format(2i5,a12,'  Data Group D: Flows')                 
 2021 FORMAT(I5,2e10.3,'  Data Block D.1 Advective Flows')    
 2022 FORMAT(I5,2e10.3,'  Data Block D.2 Pore Water Flows')   
 2023 FORMAT(I5,2e10.3,'  Data Block D.3 Sediment #1 Transport Field')
 2024 FORMAT(I5,2e10.3,'  Data Block D.4 Sediment #2 Transport Field')
 2025 FORMAT(I5,2e10.3,'  Data Block D.5 Sediment #3 Transport Field')
 2030 format(i5)                                                      
 2031 format(2(e10.3,f10.5))                                          

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
        if(IQOPT.eq.3) OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')   
        if(IQOPT.eq.4) OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',  
     $                      FORM='UNFORMATTED')
        open(96,file='waSPB.mrm',status='unknown')              
        CLOSE(90,STATUS='DELETE')
        if(IQOPT.eq.3) CLOSE(94,STATUS='DELETE')                
        if(IQOPT.eq.4) CLOSE(95,STATUS='DELETE')                
        close(96,status='delete')                               
        OPEN(90,FILE='waspdhd.out',STATUS='UNKNOWN')
        if(IQOPT.eq.3) OPEN(94,FILE='waspdh.out',STATUS='UNKNOWN')   
        if(IQOPT.eq.4) OPEN(95,FILE='waspdhu.out',STATUS='UNKNOWN',  
     $                      FORM='UNFORMATTED')
        open(96,file='waSPB.mrm',status='unknown')                   
        write(96,1011) nrfld,titleb                                  
        ntexx=1                                                      
        write(96,1012) ntexx,scalr,convr                             

c Write WASP5 Hydrodynamic File Data Record 1, Data Options:         
c  NJUN = number of segments connected by flows from the hyd. file   
c  NCHN = number of interfacial flow pairs from the hyd. file        
c  DTWASP = WASP5 time step (seconds)                                
c  TZERO = begin time step for hyd. file (seconds)                   
c  TENDHYD = end time step for hyd. file (seconds)                   
c  ISTMP = control switch, 0=time variable segment depths and velocities   
c          are read; 1=time variable segment depths and velocities are not 
c          read.                                                     
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
        nodyn=nodyn                                                  
        dtwasp = dt * float(ntsmmt)                                  
        TZERO=TBEGIN*TCON
        TENDHYD=TZERO+NTS*DT
        WRITE(90,901)NJUN,NCHN
        if(IQOPT.eq.3) then                                          
          WRITE(94,941) NJUN,NCHN, dtwasp, TZERO,TENDHYD,ISTMP       
c        WRITE(94,941) NJUN,NCHN,DT,TZERO,TENDHYD,ISTMP
        end if                                                       
        if(IQOPT.eq.4) then                                          
          WRITE(95) njun,nchn, dtwasp, tzero,tendhyd,istmp           
        end if
        write(96,1013) nchn                                          
C
C **  CHANNEL DATA
C
c Write WASP5 Hydrodynamic File Data Record 2, Segment Interface Pairs:
c   WASP expects to see boundary segments designated as "0".         
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
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              if(IQOPT.eq.3) WRITE(94,941) LUTM,LDTM                 
              if(IQOPT.eq.4) WRITE(95) LUTM,LDTM                     
              write(96,1014) unity,unity,lutm,ldtm                   
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
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                if(IQOPT.eq.3) WRITE(94,941) LUTM,LDTM               
                if(IQOPT.eq.4) WRITE(95) LUTM,LDTM                   
                unity=1.0
                write(96,1014) unity,unity,lutm,ldtm                 
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
cjh             IF(IJCTLT(I,J-1).EQ.6) LUTM=0
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
C             WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
              if(IQOPT.eq.3) WRITE(94,941) LUTM,LDTM                 
              if(IQOPT.eq.4) WRITE(95) LUTM,LDTM                     
              write(96,1014) unity,unity,lutm,ldtm                   
            END IF
cjh            IF (IJCTLT(I,J+1).EQ.6) THEN
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
C               WRITE(94,942)RLENTH,WIDTH,RMNDUM,LUTM,LDTM
                if(IQOPT.eq.3) WRITE(94,941) LUTM,LDTM               
                if(IQOPT.eq.4) WRITE(95) LUTM,LDTM                   
                write(96,1014) unity,unity,lutm,ldtm                 
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
                if(IQOPT.eq.3) WRITE(94,941)LUTM,LDTM                
                if(IQOPT.eq.4) WRITE(95) LUTM,LDTM                   
                write(96,1014) unity,unity,lutm,ldtm                 
              END IF
            END DO
          END DO
c
c write out time series of zero dispersion coefficients:
c
          d1=0.0                                                     
          t1=tzero/tcon                                              
          d2=0.0                                                     
          t2=tendhyd/tcon                                            
          nbrkq=2                                                    
          write(96,905) nbrkq                                        
          write(96,1016) d1,t1, d2,t2                                
c
c For exchange between the lower water surface layer and the upper   
c benthic layer, do the following:                                   
c
          write(96,1012) ntexx,scalr,convr                           
          ntexx=0
          DO K=1,1                                                   
            DO LT=2,LALT                                             
              I=ILLT(LT)                                             
              J=JLLT(LT)                                             
              L=LIJ(I,J)                                             
              IF (SWB(L).EQ.1.) THEN                                 
                ntexx=ntexx+1                                        
              END IF
            END DO
          END DO
          write(96,1013) ntexx                                       
          DO K=1,1                                                   
            KMUL1=KS-K                                               
            KMUL2=KMUL1+1                                            
            kmul3=KMUL2+1                                            
            DO LT=2,LALT                                             
              I=ILLT(LT)                                             
              J=JLLT(LT)                                             
              L=LIJ(I,J)                                             
              IF (SWB(L).EQ.1.) THEN                                 
                LUTM=LT-1+KMUL2*LCLTM2                               
                LDTM=LT-1+kmul3*LCLTM2                               
                write(96,1014) DXYP(l),depsed,lutm,ldtm              
              END IF
            END DO
          END DO
c
c write out time series of water-benthic exchange dispersion coefficients: 
c
          d1=sediff                                                  
          t1=tzero/tcon                                              
          d2=sediff                                                  
          t2=tendhyd/tcon                                            
          nbrkq=2                                                    
          write(96,905) nbrkq                                        
          write(96,1016) d1,t1, d2,t2                                
          ibptmp=0                                                   
          write(96,1017)ibptmp,ibptmp,ibptmp,ibptmp,                 
     +                  ibptmp,ibptmp,ibptmp,ibptmp,                 
     +                  ibptmp,ibptmp,ibptmp,ibptmp,                 
     +                  ibptmp,ibptmp,ibptmp,ibptmp                  
        END IF
C
C **  JUNCTION DATA WITH INITIAL CONDITIONS
C
c Write WASP5 Hydrodynamic File Data Record 3, Initial Segment Properties: 
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
         if (isdhd .eq. 1) WRITE(90,904) LCELL,VOLTMP,I,J
         if(IQOPT.eq.3) WRITE(94,9440) VOLTMP,DUMVOL,DEPTMP,VELTMP   
         if(IQOPT.eq.4) WRITE(95) voltmp,deptmp,veltmp               
         END DO
        END DO
C
        CLOSE(90)
        if(IQOPT.eq.3) CLOSE(94)                                     
        if(IQOPT.eq.4) CLOSE(95)                                     
        close(96)                                                    
      END IF
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP, VOLUME AND FLOW DATA
C
      OPEN(90,FILE='waspdhd.out',ACCESS='APPEND',STATUS='UNKNOWN')
      if(IQOPT.eq.3) then                                            
        OPEN(94,FILE='waspdh.out',ACCESS='APPEND',STATUS='UNKNOWN')
      end if                                                         
      if(IQOPT.eq.4) then                                            
        OPEN(95,FILE='waspdhu.out',ACCESS='APPEND',STATUS='UNKNOWN',
     $                      FORM='UNFORMATTED')
      end if                                                         
      LCLTM2=LCLT-2
      IZERO=0
      RZERO=0
      izero=izero                                                    
      rzero=rzero                                                    
C
C     NSTEP=N-NTSMMT
C     WRITE(94,945)NSTEP
C
C Write WASP5 Hydrodynamic File Data Record 4, BQ(J) flow in interface  
C pair "J":                                                          
      lchnum=0
      DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
c +++++ following lines by M. Morton to input dispersion to HYD file:
       addlw=0.0                                                     
       if (SUB(l).eq.1.) then                                        
         lw=l-1                                                      
         addlw=DYU(l)*AHULPF(l,k)*DZC(k)*0.5*(HLPF(l)                
     $         +HLPF(lw))*DXIU(l)                                    
       end if                                                        
c +++++ above added by M. Morton                                     
C!!!!!!!!CHANGES NEXT 12 LINES
       IF (SUBO(L).EQ.1.) THEN
         FLOWX=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
         IMTMP=I-1
         lchnum=lchnum+1
         if (isdhd .eq. 1) WRITE(90,944) FLOWX,IMTMP,I,J,K
         if(IQOPT.eq.3) WRITE(94,946) FLOWX, addlw                   
c         WRITE(94,946) FLOWX                                        
         if(IQOPT.eq.4) WRITE(95) FLOWX, addlw                       
       END IF
cjh       IF (IJCTLT(I+1,J).EQ.6) THEN
       IF (IJCTLT(I+1,J).EQ.8) THEN
         IF (SUBO(L+1).EQ.1.) THEN
           FLOWX=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
           IPTMP=I+1
           lchnum=lchnum+1
           if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWX,I,IPTMP,J,K
           if(IQOPT.eq.3) WRITE(94,946) FLOWX, addlw                 
c           WRITE(94,946) FLOWX                                      
           if(IQOPT.eq.4) WRITE(95) FLOWX, addlw                     
         END IF
       END IF  
       END DO
c     END DO
C
c     DO K=KC,1,-1
       DO LT=2,LALT
       I=ILLT(LT)
       J=JLLT(LT)
       L=LIJ(I,J)
c +++++ following lines by M. Morton to input dispersion to HYD file:      
       addls=0.0                                                     
       if (SVB(l).eq.1.) then                                        
         ls=LSC(l)                                                   
         addls=DXV(l)*AHVLPF(l,k)*DZC(k)*0.5*(HLPF(l)                
     $         +HLPF(ls))*DYIV(l)                                    
       end if                                                        
c +++++ above added by M. Morton                                     
C!!!!!!!CHANGES NEXT 13 LINES
       IF (SVBO(L).EQ.1.) THEN
         FLOWY=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
         JMTMP=J-1
         lchnum=lchnum+1
         if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWY,I,JMTMP,J,K
         if(IQOPT.eq.3) WRITE(94,946) FLOWY, addls                  
c         WRITE(94,946) FLOWY                                       
         if(IQOPT.eq.4) WRITE(95) FLOWY, addls                      
       END IF
cjh       IF (IJCTLT(I,J+1).EQ.6) THEN
       IF (IJCTLT(I,J+1).EQ.8) THEN
         LN=LNC(L)
         IF (SVBO(LN).EQ.1.) THEN
           FLOWY=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
           JPTMP=J+1
           lchnum=lchnum+1
           if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWY,I,J,JPTMP,K
           if(IQOPT.eq.3) WRITE(94,946) FLOWY, addls                 
c           WRITE(94,946) FLOWY                                      
           if(IQOPT.eq.4) WRITE(95) FLOWY, addls                     
         END IF
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
c +++++ following lines by M. Morton to input dispersion to HYD file:
       addl=0.0                                                      
       if (SPB(l).eq.1.) then                                        
         addl=DXYP(l)*ABLPF(l,k)*DZIG(k)                             
       end if                                                        
c +++++ above added by M. Morton                                     
       IF (SWB(L).EQ.1) THEN
         FLOWZ=-DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
         KPTMP=K+1
         lchnum=lchnum+1
         if (isdhd .eq. 1) WRITE(90,944) lchnum,FLOWZ,I,J,K,KPTMP
         if(IQOPT.eq.3) WRITE(94,946) FLOWZ, addl                    
c         WRITE(94,946) FLOWZ
         if(IQOPT.eq.4) WRITE(95) FLOWZ, addl                        
       END IF  
       END DO
      END DO
      END IF
C
C Write WASP5 Hydrodynamic File Data Record 5, Segment Properties:   
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
C      QIN=QSUMELPF(L)*DZC(K)
C      FLOWXI=DYU(L)*(UHLPF(L,K)+SVPT*UVPT(L,K))*DZC(K)
C      FLOWYI=DXV(L)*(VHLPF(L,K)+SVPT*VVPT(L,K))*DZC(K)
C      FLOWZI=DXYP(L)*(WLPF(L,K-1)+SVPT*WVPT(L,K-1))
C      FLOWXO=DYU(L+1)*(UHLPF(L+1,K)+SVPT*UVPT(L+1,K))*DZC(K)
C      FLOWYO=DXV(LN)*(VHLPF(LN,K)+SVPT*VVPT(LN,K))*DZC(K)
C      FLOWZO=DXYP(L)*(WLPF(L,K)+SVPT*WVPT(L,K))
C      QQSUM=QIN+FLOWXI+FLOWYI+FLOWZI-FLOWXO-FLOWYO-FLOWZO
       DEPTH=HLPF(L)*DZC(K)
       VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K)
     $         +UHLPF(L+1,K)+SVPT*UVPT(L+1,K))/HLPF(L)
       VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K)
     $         +VHLPF(LN,K)+SVPT*VVPT(LN,K))/HLPF(L)
       VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1)
     $         +WLPF(L,K)+SVPT*WVPT(L,K))
       VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
       if (isdhd .eq. 1) WRITE(90,902) LCELTMP,VOLUM,I,J,K
       if(IQOPT.eq.3) WRITE(94,946) VOLUM,qqsum,DEPTH,VELMAG         
       if(IQOPT.eq.4) WRITE(95) VOLUM, DEPTH, VELMAG                 
       END DO
      END DO
C
      CLOSE(90)
      if(IQOPT.eq.3) CLOSE(94)                                       
      if(IQOPT.eq.4) CLOSE(95)                                       
C
C----------------------------------------------------------------------C
C
  901 FORMAT(2I5,E12.4,4I5,E12.4)
  902 FORMAT(I5,2X,3F20.8,3I5)
  903 FORMAT(3E12.4,2I5)
  904 FORMAT(I5,2X,F20.8,10I5)
  905 FORMAT(I5)
  906 FORMAT(5E12.4)
  941 FORMAT(2I5,3F20.8,I5)
  942 FORMAT(3E12.4,2I5)
  943 FORMAT(3E12.4,2I5)
  944 FORMAT(I5,2x,F20.8,10I5)
 9440 FORMAT(4F20.8)
  945 FORMAT(I5)
  946 FORMAT(4e17.9)
C
C**********************************************************************C
C
      JSWASP=0
C
      RETURN
      END
