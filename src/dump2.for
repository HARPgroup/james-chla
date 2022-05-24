      SUBROUTINE DUMP2
c 
C  save EFDC output for GrADS, jeff ji, 10/20/98
C
      INCLUDE 'efdc.par'
      INCLUDE 'efdc.cmn'
	character fname*10
c
cji      write(6,*) ITIMES,iaverage,tcon
C ------------------------------------------------------------------------
      IF (IDUMP2.GE.1.AND.IDUMP2.LT.100) THEN       ! # 99
      IDUMP2=IDUMP2+100
c     write(6,*) "? dump2 ", idump2
      IUDUMP=281               ! hyts.bin
      IUDUM2=282               ! hy3d.bin
      open(IUDUMP,file='hyts.bin',form='binary')       ! time series
      open(IUDUM2,file='hy3d.bin',form='binary')       ! 4D arrays for GrADS
C hyts.bin
      DO J=1,MLTMSR
      ESAVE(J)=0.0
      TAUSAVE(J)=0.0
      SEBTSAVE(J)=0.0
      SNBTSAVE(J)=0.0
      VOLBSAVE(J)=0.0
      tbxxSAVE(J)=0.0
      tbyySAVE(J)=0.0
      DO NT=1,NTXM
      TBXSAVE(J,NT)=0.0
      TBXPSAVE(J,NT)=0.0
      ENDDO
      ENDDO
      DO J=1,MLTMSR
      DO K=1,KC
      UZSAVE(J,K)=0.0
      VZSAVE(J,K)=0.0
      WZSAVE(J,K)=0.0
      SZSAVE(J,K)=0.0
      TZSAVE(J,K)=0.0
      DYESAVE(J,K)=0.0
      SEDTSAVE(J,K)=0.0
      SNDTSAVE(J,K)=0.0
        fxSAVE(J,K)=0.0
        fySAVE(J,K)=0.0
      DO NT=1,NTXM
      TOXSAVE(J,K,NT)=0.0
      TOXPSAVE(J,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C for hy3d.bin
      DO L=1,LC
      DUMPP (L)=0.0
      DUMPtau(L)=0.0
      DUMPSEBT(L)=0.0
      DUMPSNBT(L)=0.0
      DUMPVOLB(L)=0.0
      DUMPtbxx(L)=0.0
      DUMPtbyy(L)=0.0
      DO NT=1,NTXM
      DUMPTBX(L,NT)=0.0
      DUMPTBXP(L,NT)=0.0
      ENDDO
      ENDDO
      DO K=1,KC
      DO L=1,LC
      DUMPU (L,K)=0.0
      DUMPV (L,K)=0.0
      DUMPW (L,K)=0.0
      DUMPS (L,K)=0.0
      DUMPT (L,K)=0.0
      DUMPDYE(L,K)=0.0
      DUMPSEDT(L,K)=0.0
      DUMPSNDT(L,K)=0.0
      DUMPfx  (L,K)=0.0
      DUMPfy  (L,K)=0.0
      ENDDO
      ENDDO
      DO NT=1,NTXM
      DO K=1,KC
      DO L=1,LC
      DUMPTOX(L,K,NT)=0.0
      DUMPTOXP(L,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C
c save location and time information to both binary files
      DO II=1,2
c     DO II=1,1
      IF(II.EQ.1) IUDUMx=IUDUMP        ! hydstx1.bin
      IF(II.EQ.2) IUDUMx=IUDUM2        !       2
      WRITE(IUDUMx) IC,JC,KC,LC,NTOX,NSED,NSND,LIJ,L2I,L2J,LNC,IGrADS
      WRITE(IUDUMx) IDUMPe   ,IDUMPu   ,IDUMPv
     1,IDUMPw   ,IDUMPs   ,IDUMPt   ,
     1 IDUMPd   ,IDUMPc   ,IDUMPn   ,IDUMPx
     1,IDUMPb   ,IDUMPfx   ,IDUMPfy
      WRITE(IUDUMx) MLTMSR,ITIMES,DT,
     &         (ILTMSR(I),JLTMSR(I),I=1,MLTMSR)                     ! time series
      WRITE(IUDUMx) JHM,(IHIST(I,2),I=1,JHM),IAVERAGE,TBEGIN           ! 4D
      WRITE(IUDUMx) DLON,DLAT,CUE,CVE,CUN,CVN                       ! lxly.inp
      WRITE(IUDUMx) DXP,DYP,BELV,DZC                                ! dxdy.inp
      ENDDO
C
      ENDIF     ! # 99
C
C-----------------------------------------------------------------------
C--  hyts.bin, time series *save(NT,*) arraeys ---------------------
c----------------------------------------------------------------------
c
c
c   !!! Time series fluxes are not calculated !!, jeff ji, 10/26/98
c
c
c
       IF(N.LT.NBTMSR.OR.N.GT.NSTMSR) GO TO 1000
c
       DO NT=1,MLTMSR
       I=ILTMSR(NT)
       J=JLTMSR(NT)
       L=LIJ(I,J)
       LN=LNC(L)
       LS=LSC(L)
c      IF(IDUMPe.GE.1) ESAVE(NT)=ESAVE(NT)+GI*P(L)         ! surface tides, mean=0
       IF(IDUMPe.GE.1) ESAVE(NT)=ESAVE(NT)+GI*P(L)-BELV(L) ! total water depth
       IF(IDUMPc.GE.1.or.IDUMPn.GE.1) TAUSAVE(NT)=TAUSAVE(NT)+TAU(L)
       DO K=1,KC
       IF(IDUMPu.GE.1) UZSAVE(NT,K)=UZSAVE(NT,K)+0.5*(U(L,K)+U(L+1,K))
c      IF(IDUMPu.GE.1) UZSAVE(NT,K)=UZSAVE(NT,K)+     U(L,K)            ! No shift
       IF(IDUMPv.GE.1) VZSAVE(NT,K)=VZSAVE(NT,K)+0.5*(V(L,K)+V(LN,K))
c      IF(IDUMPv.GE.1) VZSAVE(NT,K)=VZSAVE(NT,K)+     V(L,K)            ! No shift
       IF(IDUMPw.GE.1) WZSAVE(NT,K)=WZSAVE(NT,K)+0.5*(W(L,K)+W(L,K-1))
       IF(IDUMPs.GE.1) SZSAVE(NT,K)=SZSAVE(NT,K)+SAL(L,K)
       IF(IDUMPt.GE.1) TZSAVE(NT,K)=TZSAVE(NT,K)+TEM(L,K)
       IF(IDUMPd.GE.1) DYESAVE(NT,K)=DYESAVE(NT,K)+DYE(L,K)
       IF(IDUMPc.GE.1) SEDTSAVE(NT,K)=SEDTSAVE(NT,K)+SEDT(L,K)
       IF(IDUMPn.GE.1) SNDTSAVE(NT,K)=SNDTSAVE(NT,K)+SNDT(L,K)
       IF(IDUMPx.GE.1) THEN
       DO NT2=1,NTOX
       TOXSAVE(NT,K,NT2)=TOXSAVE(NT,K,NT2)+TOX(L,K,NT2)
       TOXPSAVE(NT,K,NT2)=TOXPSAVE(NT,K,NT2)+TOXPFTW(L,K,NT2)
       ENDDO
       ENDIF
       IF(IDUMPfx.GE.1)
     & fxSAVE(NT,K)=fxSAVE(NT,K)+U(L,K)*DZC(K)*0.25*                  ! No shift
     & (GI*P(L)-BELV(L)+GI*P(L-1)-BELV(L-1))*(DYP(L)+DYP(L-1))
       if(IDUMPfy.gt.0)
     & fySAVE(NT,K)=fySAVE(NT,K)+V(L,K)*DZC(K)*0.25*                  ! No shift
     & (GI*P(L)-BELV(L)+GI*P(LS)-BELV(LS))*(DXP(L)+DXP(LS))

       ENDDO   ! K
       ENDDO   ! NT
c
       IF(IDUMPb.GE.1) THEN
       DO NT=1,MLTMSR
       I=ILTMSR(NT)
       J=JLTMSR(NT)
       L=LIJ(I,J)
       do K=1,kb
       SEBTSAVE(NT)=SEBTSAVE(NT)+SEDBT(L,k)  ! for test case, ji, 12/11/00
       enddo
c      SNBTSAVE(NT)=SNBTSAVE(NT)+SNDBT(L)
c      VOLBSAVE(NT)=VOLBSAVE(NT)+VOLBW2(L)
c      tbxxSAVE(NT)=tbxxSAVE(NT)+tbx   (L)
c      tbyySAVE(NT)=tbyySAVE(NT)+tby   (L)
       ENDDO
c      DO NT=1,MLTMSR
c      I=ILTMSR(NT)
c      J=JLTMSR(NT)
c      L=LIJ(I,J)
c      DO NT2=1,NTOX
c      TBXSAVE(NT,NT2)=TBXSAVE(NT,NT2)+TOXB(L,NT2)
c      TBXPSAVE(NT,NT2)=TBXPSAVE(NT,NT2)+TOXPFTB(L,NT2)
c      ENDDO
c      ENDDO
       ENDIF
C
C-------------------------------------------------
C
C
      IF(MOD(N,ITIMES).NE.0) GO TO 1000
C
C - Average, since it's not taking too much CPU, IF statement not used here,
      DO NT=1,MLTMSR
       ESAVE(NT)=ESAVE(NT)/ITIMES
       TAUSAVE(NT)=TAUSAVE(NT)/ITIMES
       SEBTSAVE(NT)=SEBTSAVE(NT)/ITIMES
       SNBTSAVE(NT)=SNBTSAVE(NT)/ITIMES
       VOLBSAVE(NT)=VOLBSAVE(NT)/ITIMES
       tbxxSAVE(NT)=tbxxSAVE(NT)/ITIMES
       tbyySAVE(NT)=tbyySAVE(NT)/ITIMES
       ENDDO
c
       DO NT=1,MLTMSR
       DO NT2=1,NTOX
       TBXSAVE(NT,NT2)=TBXSAVE(NT,NT2)/ITIMES
       TBXPSAVE(NT,NT2)=TBXPSAVE(NT,NT2)/ITIMES
       ENDDO
       ENDDO
c
       DO NT=1,MLTMSR
       DO K=1,KC
       UZSAVE(NT,K)=UZSAVE(NT,K)/ITIMES
       VZSAVE(NT,K)=VZSAVE(NT,K)/ITIMES
       WZSAVE(NT,K)=WZSAVE(NT,K)/ITIMES
       SZSAVE(NT,K)=SZSAVE(NT,K)/ITIMES
       TZSAVE(NT,K)=TZSAVE(NT,K)/ITIMES
       DYESAVE(NT,K)=DYESAVE(NT,K)/ITIMES
       SEDTSAVE(NT,K)=SEDTSAVE(NT,K)/ITIMES
       SNDTSAVE(NT,K)=SNDTSAVE(NT,K)/ITIMES
         fxSAVE(NT,K)=  fxSAVE(NT,K)/ITIMES
         fySAVE(NT,K)=  fySAVE(NT,K)/ITIMES
       ENDDO
       ENDDO
c
       DO NT=1,MLTMSR
       DO NT2=1,NTOX
       DO K=1,KC
       TOXSAVE(NT,K,NT2)=TOXSAVE(NT,K,NT2)/ITIMES
       TOXPSAVE(NT,K,NT2)=TOXPSAVE(NT,K,NT2)/ITIMES
       ENDDO
       ENDDO
       ENDDO
C
C-------- WRITE ------
C
      TMIDDLE=(DT*(N-0.5*ITIMES))/TCON+TBEGIN
      WRITE(IUDUMP) TMIDDLE
c
c      write(301,3001) tmiddle, (sedtsave(1,k),k=1,kc)
3001  format(999f9.3)
C
      DO J=1,MLTMSR
      IF(IDUMPe.GE.1) WRITE(IUDUMP) ESAVE(J)
      IF(IDUMPc.GE.1.or.IDUMPn.GE.1) WRITE(IUDUMP) TAUSAVE(J)
      IF(IDUMPb.GE.1) THEN
      WRITE(IUDUMP) SEBTSAVE(J),SNBTSAVE(J),VOLBSAVE(J)
     *,tbxxSAVE(J),tbyySAVE(J)
      DO NT=1,NTOX
      WRITE(IUDUMP) TBXSAVE(J,NT),TBXPSAVE(J,NT)
      ENDDO
      ENDIF
      DO K=1,KC
      IF(IDUMPu.GE.1) WRITE(IUDUMP) UZSAVE(J,K)
      IF(IDUMPv.GE.1) WRITE(IUDUMP) VZSAVE(J,K)
      IF(IDUMPw.GE.1) WRITE(IUDUMP) WZSAVE(J,K)
      IF(IDUMPs.GE.1) WRITE(IUDUMP) SZSAVE(J,K)
      IF(IDUMPt.GE.1) WRITE(IUDUMP) TZSAVE(J,K)
      IF(IDUMPd.GE.1) WRITE(IUDUMP) DYESAVE(J,K)
      IF(IDUMPc.GE.1) WRITE(IUDUMP) SEDTSAVE(J,K)
      IF(IDUMPn.GE.1) WRITE(IUDUMP) SNDTSAVE(J,K)
      IF(IDUMPfx.GE.1) WRITE(IUDUMP) fxSAVE(J,K)
      IF(IDUMPfy.GE.1) WRITE(IUDUMP) fySAVE(J,K)
      IF(IDUMPx.GE.1) THEN
      DO NT=1,NTOX
      WRITE(IUDUMP) TOXSAVE(J,K,NT),TOXPSAVE(J,K,NT)
      ENDDO
      ENDIF
      ENDDO      ! K
      ENDDO      ! J
C
C*********************************************************
C
C      VSTOR=ESUM
C      EM=(EM-EMI)*SKILLI
C      APEC=(APE-APEI)*SKILLI
C      TSUM=TSUM/VOLUME
C      SSUM=SSUM/VOLUME
C
C      WRITE(IUDUMP) VSTOR,EM,APEC,TSUM,SSUM
C     ENDIF
C
      DO J=1,MLTMSR
      ESAVE(J)=0.0
      TAUSAVE(J)=0.0
      SEBTSAVE(J)=0.0
      SNBTSAVE(J)=0.0
      VOLBSAVE(J)=0.0
      tbxxSAVE(J)=0.0
      tbyySAVE(J)=0.0
      DO NT=1,NTOX
      TBXSAVE(J,NT)=0.0
      TBXPSAVE(J,NT)=0.0
      ENDDO
      ENDDO
C
      DO J=1,MLTMSR
      DO K=1,KC
      UZSAVE(J,K)=0.0
      VZSAVE(J,K)=0.0
      WZSAVE(J,K)=0.0
      SZSAVE(J,K)=0.0
      TZSAVE(J,K)=0.0
      DYESAVE(J,K)=0.0
      SEDTSAVE(J,K)=0.0
      SNDTSAVE(J,K)=0.0
        fxSAVE(J,K)=0.0
        fySAVE(J,K)=0.0
      DO NT=1,NTOX
      TOXSAVE(J,K,NT)=0.0
      TOXPSAVE(J,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C
 1000 CONTINUE
C-----------------------------------------------------------------------
C-  hy3d.bin, save 3D data DUMP*(LCM,*) arreya ----------------------
C-----------------------------------------------------------------------
      DO 2000 JHIST=1,JHM
      IF(N.LT.IHIST(JHIST,1).OR.N.GT.IHIST(JHIST,2)) GO TO 2000
c Water Depth
       DO L=2,LA
       IF(IDUMPe.GE.1) DUMPP (L)=DUMPP (L)+ P(L) ! dumpp will be finalized later by adjusting and averaging, must see mean for definition!
       IF(IDUMPc.GE.1.or.IDUMPn.GE.1) DUMPtau(L)=DUMPtau(L)+tau(L)
       ENDDO
C U,V,W,S,T,Dye,SEDT,SNDT in water column
c
c This kind of sigma level time average might have problem, especially in tidal waters,
c since the spacial locations change with time!!  jeff ji, 3/6/98
c
       DO K=1,KC
       DO L=2,LA
       IF(IDUMPu.GE.1) DUMPU (L,K)=DUMPU (L,K)+U (L,K)
       IF(IDUMPv.GE.1) DUMPV (L,K)=DUMPV (L,K)+V (L,K)
       IF(IDUMPw.GE.1) DUMPW (L,K)=DUMPW (L,K)+W (L,K)
       IF(IDUMPs.GE.1) DUMPS (L,K)=DUMPS (L,K)+SAL(L,K)
       IF(IDUMPt.GE.1) DUMPT (L,K)=DUMPT (L,K)+TEM(L,K)
       IF(IDUMPd.GE.1) DUMPDYE(L,K)=DUMPDYE(L,K)+DYE(L,K)
       IF(IDUMPc.GE.1) DUMPsedt(L,K)=DUMPsedt(L,K)+sedt(L,K)
       IF(IDUMPn.GE.1) DUMPsndt(L,K)=DUMPsndt(L,K)+sndt(L,K)
       if(IDUMPfx.gt.0)
     & DUMPfx(L,K)=DUMPfx(L,K)+U(L,K)*DZC(K)*0.25*
     & (GI*P(L)-BELV(L)+GI*P(L-1)-BELV(L-1))*(DYP(L)+DYP(L-1))
       if(IDUMPfy.gt.0) then
       LS=LSC(L)
       DUMPfy(L,K)=DUMPfy(L,K)+V(L,K)*DZC(K)*0.25*
     & (GI*P(L)-BELV(L)+GI*P(LS)-BELV(LS))*(DXP(L)+DXP(LS))
       endif
       ENDDO
       ENDDO
c TOX in water column
       IF(IDUMPx.EQ.0) go to 901
       DO NT=1,NTOX
       DO K=1,KC
       DO L=2,LA
       DUMPtox(L,K,NT)=DUMPtox(L,K,NT)+tox(L,K,NT)
       DUMPtoxp(L,K,NT)=DUMPtoxp(L,K,NT)+toxpftw(L,K,NT)
       ENDDO
       ENDDO
       ENDDO
 901   continue
c BED
       if(IDUMPb.EQ.0) GO TO 902
       DO L=2,LA
        do k=1,kb
        DUMPsebt(L)=DUMPsebt(L)+sedbt(L,K) ! for test case, Ji, 12/11/00
        enddo
c       DUMPsnbt(L)=DUMPsnbt(L)+sndbt(L)
c       DUMPvolb(L)=DUMPvolb(L)+volbw2(L)
c       DUMPtbxx(L)=DUMPtbxx(L)+tbx   (L)
c       DUMPtbyy(L)=DUMPtbyy(L)+tby   (L)
       ENDDO
c      DO NT=1,NTOX
c      DO L=2,LA
c       DUMPtbx(L,NT)=DUMPtbx(L,NT)+toxb(L,NT)
c       DUMPtbxp(L,NT)=DUMPtbxp(L,NT)+toxpftb(L,NT)
c      ENDDO
c      ENDDO
 902   CONTINUE
C
      IF(N.NE.IHIST(JHIST,2)) GO TO 2000
C
C ---- Average and Adjust output for graphics ----------------
c Water Depth
       DO L=2,LA
       IF(IDUMPe.GE.1) DUMPP (L)=GI*DUMPP (L)/IAVERAGE-BELV(L)
       IF(IDUMPc.GE.1.or.IDUMPn.GE.1) DUMPtau(L)=DUMPtau(L)/IAVERAGE
       ENDDO
c write & check water depth
c     TMIDDLE=(DT*(N-0.5*IAVERAGE))/TCON+TBEGIN
c     write(191,*) " Tmiddle=  ", tmiddle
c     write(191,199) (I,I=1,IC)
c199   format(4x,2x,999i4)
c      DO J=1,JC
c      write(191,198) J,(dumpp(LIJ(I,J)),I=1,IC)
c198   format(i4,2x,999f4.1)
c      END DO
c      write(191,199) (I,I=1,IC)
C U,V,W,S,T,Dye,SEDT,SNDT in water column
       DO K=1,KC
       DO L=2,LA
       LN=LNC(L)
c      IF(IDUMPu.GE.1) DUMPU (L,K)=0.5*(DUMPU (L,K)+DUMPU(L+1,K))/IAVERAGE
c      IF(IDUMPv.GE.1) DUMPV (L,K)=0.5*(DUMPV (L,K)+DUMPV(LN ,K))/IAVERAGE
       IF(IDUMPu.GE.1) DUMPU (L,K)=     DUMPU (L,K)              /
     *  IAVERAGE
       IF(IDUMPv.GE.1) DUMPV (L,K)=     DUMPV (L,K)              /
     *  IAVERAGE
       IF(IDUMPw.GE.1) DUMPW (L,K)=0.5*(DUMPW (L,K)+DUMPW(L,K-1))/
     *  IAVERAGE
       IF(IDUMPs.GE.1) DUMPS (L,K)=DUMPS (L,K)/IAVERAGE
       IF(IDUMPt.GE.1) DUMPT (L,K)=DUMPT (L,K)/IAVERAGE
       IF(IDUMPd.GE.1) DUMPDYE(L,K)=DUMPDYE(L,K)/IAVERAGE
       IF(IDUMPc.GE.1) DUMPsedt(L,K)=DUMPsedt(L,K)/IAVERAGE
       IF(IDUMPn.GE.1) DUMPsndt(L,K)=DUMPsndt(L,K)/IAVERAGE
       IF(IDUMPfx.GE.1) DUMPfx  (L,K)=DUMPfx  (L,K)/IAVERAGE
       IF(IDUMPfy.GE.1) DUMPfy  (L,K)=DUMPfy  (L,K)/IAVERAGE
       ENDDO
       ENDDO
c TOX in water column
       IF(IDUMPx.EQ.0) go to 911
       DO NT=1,NTOX
       DO K=1,KC
       DO L=2,LA
       DUMPtox(L,K,NT)=DUMPtox(L,K,NT)/IAVERAGE
       DUMPtoxp(L,K,NT)=DUMPtoxp(L,K,NT)/IAVERAGE
       ENDDO
       ENDDO
       ENDDO
 911   continue
c BED
       if(IDUMPb.EQ.0) GO TO 912
       DO L=2,LA
        DUMPsebt(L)=DUMPsebt(L)/IAVERAGE
        DUMPsnbt(L)=DUMPsnbt(L)/IAVERAGE
        DUMPvolb(L)=DUMPvolb(L)/IAVERAGE
        DUMPtbxx(L)=DUMPtbxx(L)/IAVERAGE
        DUMPtbyy(L)=DUMPtbyy(L)/IAVERAGE
       ENDDO
       DO NT=1,NTOX
       DO L=2,LA
        DUMPtbx(L,NT)=DUMPtbx(L,NT)/IAVERAGE
        DUMPtbxp(L,NT)=DUMPtbxp(L,NT)/IAVERAGE
       ENDDO
       ENDDO
 912   CONTINUE
c ---
c     write(229,*) tmiddle
c     do L=1,LC
c     write(229,487) L,dumps(L,kc)
c487  format(i5,f12.3)
c     enddo
c write --
      TMIDDLE=(DT*(N-0.5*IAVERAGE))/TCON+TBEGIN
      WRITE(IUDUM2) TMIDDLE
      IF(IDUMPe.GE.1) WRITE(IUDUM2) DUMPP
      IF(IDUMPc.GE.1.or.IDUMPn.GE.1) WRITE(IUDUM2) DUMPtau
      IF(IDUMPu.GE.1) WRITE(IUDUM2) DUMPU
      IF(IDUMPv.GE.1) WRITE(IUDUM2) DUMPV
      IF(IDUMPw.GE.1) WRITE(IUDUM2) DUMPW
      IF(IDUMPs.GE.1) WRITE(IUDUM2) DUMPS
      IF(IDUMPt.GE.1) WRITE(IUDUM2) DUMPT
      IF(IDUMPd.GE.1) WRITE(IUDUM2) DUMPDYE
      IF(IDUMPc.GE.1) WRITE(IUDUM2) DUMPSEDT
      IF(IDUMPn.GE.1) WRITE(IUDUM2) DUMPSNDT
      IF(IDUMPfx.GE.1) WRITE(IUDUM2) DUMPfx
      IF(IDUMPfy.GE.1) WRITE(IUDUM2) DUMPfy
      IF(IDUMPx.GE.1) WRITE(IUDUM2) DUMPTOX,DUMPTOXP
      IF(IDUMPb.GE.1) WRITE(IUDUM2) DUMPSEBT  !,DUMPSNBT,DUMPVOLB,DUMPTBX
c    1                             ,DUMPTBXP,DUMPtbxx,DUMPtbyy                   ! Test, 12/11/00
C
      DO L=1,LC
      DUMPP (L)=0.0
      DUMPtau(L)=0.0
      DUMPSEBT(L)=0.0
      DUMPSNBT(L)=0.0
      DUMPVOLB(L)=0.0
      DUMPtbxx(L)=0.0
      DUMPtbyy(L)=0.0
      DO NT=1,NTOX
      DUMPTBX(L,NT)=0.0
      DUMPTBXP(L,NT)=0.0
      ENDDO
      ENDDO
      DO K=1,KC
      DO L=1,LC
      DUMPU (L,K)=0.0
      DUMPV (L,K)=0.0
      DUMPW (L,K)=0.0
      DUMPS (L,K)=0.0
      DUMPT (L,K)=0.0
      DUMPDYE(L,K)=0.0
      DUMPSEDT(L,K)=0.0
      DUMPSNDT(L,K)=0.0
      DUMPfx  (L,K)=0.0
      DUMPfy  (L,K)=0.0
      ENDDO
      ENDDO
      DO NT=1,NTOX
      DO K=1,KC
      DO L=1,LC
      DUMPTOX(L,K,NT)=0.0
      DUMPTOXP(L,K,NT)=0.0
      ENDDO
      ENDDO
      ENDDO
C added for Test hot start, hardwired for 2001 test run, Ji, 12/6/01
c      if(abs(tmiddle-(120.0-0.5)).le.0.1) then
c	fname='maskd.inp'
c     write(6,*) tmiddle, n
c      call cellmask(fname)    ! hardwired for 2001 test run, Ji, 12/6/01
c      endif
c      if(abs(tmiddle-(305.0-0.5)).le.0.1) then
c      fname='mask.inp'
c     write(6,*) tmiddle, n
c      call cellmask(fname)    ! hardwired for 2001 test run, Ji, 12/6/01
c      endif
 2000 CONTINUE
C
      RETURN
      END
