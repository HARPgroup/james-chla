C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C


      SUBROUTINE SMRIN1(IWQDT,DT,IC,JC,KC,LA,ITIMES,NTSPTC)
C
C**********************************************************************C
C
C **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
C
C**********************************************************************C
C
!      INCLUDE 'efdc.par'
!      INCLUDE 'efdc.cmn'
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'     
C
c     PARAMETER (SMCW2=2.739726E-5) ! moved to efdc.cmn, Ji, 9/20/99
c      DIMENSION SMKPON(NSMGM),SMKPOP(NSMGM),SMKPOC(NSMGM) ! Moved to efdc.cmn, Ji, 9/20/99
c      DIMENSION SMTHKN(NSMGM),SMTHKP(NSMGM),SMTHKC(NSMGM)
      INTEGER isSKIP
      CHARACTER TITLE(3)*79, ccmrm*1
C
      IF(isTBIN.EQ.0) THEN
      OPEN(1,FILE='wqsdts1.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE='wqsdts2.out',STATUS='UNKNOWN')
      CLOSE(1,STATUS='DELETE')
C
      OPEN(1,FILE='wqsdts1.out',STATUS='UNKNOWN')
	WRITe(1,*)IWQTS,KC
      WRITE(1,1996)
      CLOSE(1)
C
      OPEN(1,FILE='wqsdts2.out',STATUS='UNKNOWN')
	WRITe(1,*)IWQTS,KC
      WRITE(1,1997)
      CLOSE(1)
      
      ELSE
      OPEN(1,FILE='wqsdts1.bin',STATUS='UNKNOWN') 
      CLOSE(1,STATUS='DELETE')   
      OPEN(550,FILE='wqsdts1.bin',acccess='stream')
      ENDIF
c
      SMCW2=2.739726E-5
C
 1996 FORMAT('c   I    J      TIME       nh41      nh42',
     $       '        fnh4       no31       no32      fno3',
     $       '        po41       po42      fpo4d      h2s1',
     $       '        h2s2      bfo2       bfcod       si1',
     $       '        si2       fsad        smt        bst',
     $       '        pon        pop        poc')
C
 1997 FORMAT('c   I    J      TIME       csod      nsod',
     $       '       d1po4       d1si        ss       jnit',
     $       '       jden      jaqh2s     jgch4       dgfn',
     $       '       dgfp       dgfc       dfn1       dfn2',
     $       '       dfn3       dfp1       dfp2       dfp3',
     $       '       dfc1       dfc2       dfc3')
C
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
C
      OPEN(1,FILE='wq3dsd.inp',STATUS='UNKNOWN')
C
C Read first line in WQ3DSD.INP file.  If first character is '#', then
C this is the new version with annotated comments added (i.e., uses the
C SKIPCOMM subroutine to skip comment lines.  Comment lines begin with
C a "C", "c", or "#" character in column 1.  If "#" is not found as the
C first character in the file, then the old method of reading the
C WQ3DSD.INP file is used to preserve backward compatability.
C
      isSKIP = 0
      READ(1,'(A1)') CCMRM
      BACKSPACE(1)
      IF (CCMRM .EQ. '#') isSKIP = 1
      CCMRM = '#'
C
C01 read main title cards:
C
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) (TITLE(M), M=1,3)
      WRITE(2,999)
      WRITE(2,5100) (TITLE(M), M=1,3)
C
C02 I/O control variables and temperature related variables
C
c     READ(1,999)
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) ISMZ,ISMICI,ISMRST,ISMHYST,ISMZB,isSDBIN
c      if (sodmult .le. 0.0) sodmult = 1.0
C
      WRITE(2,53)'* # of zones for spat. vary. parameters IN SPM =',ISMZ
      IF (ISMZ.GT.NSMZ) STOP 'ERROR!! ISMZ should be <= NSMZ'
      IF (ISMICI.EQ.1) THEN
        WRITE(2,50)'* Spatially/temporally-varying ICs from wqsdici.inp'
cxh        OPEN(INSMICI,FILE='SM-ICI.IN',STATUS='OLD')
CXH wqsdici.inp REPLACES SM-ICI.IN
       ELSE IF (ISMICI.EQ.2) THEN
        WRITE(2,50)'* Spatially/temporally-varying ICs from wqsdrst.inp'
cxh        OPEN(INSMRST,FILE='SM-RST.IN',STATUS='OLD')
CXH wqsdrst.inp REPLACES SM-RST.IN
       ELSE
        WRITE(2,50)'* Spatially/temporally constant initial conditions'
      END IF
      IF (ISMRST.EQ.1) THEN
        WRITE(2,50)'* Write spatial distributions to wqsdrst.out'
CXH        OPEN(ISMORST,FILE='SM-RST.OUT',STATUS='UNKNOWN')
CXH wqsdrst.out REPLACES SM-RST.OUT
cxh        WRITE(ISMORST,888)
C
       ELSE
        WRITE(2,50)'* No writing to ISMORST                           '
      END IF
      IF (ISMHYST.EQ.1) THEN
        WRITE(2,50)'* Hysteresis in benthic mixing is activated       '
       ELSE
        WRITE(2,50)'* Hysteresis in benthic mixing is not activated   '
      END IF
      IF (ISMZB.EQ.1) THEN
        WRITE(2,50)'* Diagnostic output for FUNC ZBRENT (zbrent.log)  '
        OPEN(99,FILE='zbrent.log',STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE='zbrent.log',STATUS='UNKNOWN')
        WRITE(99,53)'   ITNWQ    L    I    J         SOD          '
        CLOSE(99)
       ELSE
        WRITE(2,50)'* No diagnostic output for FUNC ZBRENT            '
      END IF
!
! isSDBIN > 0 turns on binary file output for benthic flux rates
!
      IF (isSDBIN .GT. 0) THEN  ! =0, Ji, 9/20/99
        CALL WQzero4(LA)
!         CALL INITbin4
      END IF
!
! transfer efdcwin.inp's parameters to override wq3dsd.inp's
! 
!      ISMICI=0
!      if(isresti.eq.1) IsmICI=2
!
      ISMTS=IWQTS
      SMTSDT=24.0*ITIMES/NTSPTC
      ISMTSB = IWQTSB !NINT(TSMTSB/DTD)
      ISMTSE = IWQTSE ! NINT(TSMTSE/DTD)
   !   ISMTSDT = NINT(SMTSDT*3600.0/DT)

  
C
C04
      DO M=1,IWQTS 
      LsmTS(M)=LWQTS(M)
!      DO NW=1,NTSsmV
!      ICsmTS(NW,M)=1
!      ENDDO
      ENDDO
c

!
  999 FORMAT(1X)
 5100 FORMAT(A79)
 5101 FORMAT(10I8)
 5103 FORMAT(10F8.4)
 5104 FORMAT(I8, 3F8.4)
   50 FORMAT(A50)
   51 FORMAT(A27, 3(F8.4,2X))
   52 FORMAT((A45, E10.4))
   53 FORMAT((A48, I10))
   55 FORMAT(A31, 2I5)
   84 FORMAT(3(A26,F10.4,A5,/), 2(A26,I8,A10,/))
C05
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      WRITE(2,999)
      READ(1,*) SMDIFT
      WRITE(2,52)'* Diff coeff (m^2/s) for sed temperature   = ',SMDIFT
      SMDIFT = SMDIFT*8.64E4
C
C06 Spatially constant parameters for spliting depositional fluxes of algae
C
c     READ(1,999)
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
C
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      WRITE(2,999)
      READ(1,*) SMFNBC(1),SMFNBC(2),SMFNBC(3),SMFNBD(1),SMFNBD(2),
     *  SMFNBD(3),SMFNBG(1),SMFNBG(2),SMFNBG(3)
      WRITE(2,50)'* Cyanobacteria-N split into G1, G2 & G3 classes  '
      WRITE(2,51)' : (FNBc1, FNBc2, FNBc3) = ', (SMFNBC(M),M=1,3)
      WRITE(2,50)'* Diatoms-N split into G1, G2 & G3 classes        '
      WRITE(2,51)' : (FNBd1, FNBd2, FNBd3) = ', (SMFNBD(M),M=1,3)
      WRITE(2,50)'* Blue-green algae-N split into G1, G2, G3 classes'
      WRITE(2,51)' : (FNBg1, FNBg2, FNBg3) = ', (SMFNBG(M),M=1,3)
      SUMNBC=SMFNBC(1)+SMFNBC(2)+SMFNBC(3)
      SUMNBD=SMFNBD(1)+SMFNBD(2)+SMFNBD(3)
      SUMNBG=SMFNBG(1)+SMFNBG(2)+SMFNBG(3)
      IF (SUMNBC.LT.0.9999.OR.SUMNBC.GT.1.0001)
     *  STOP 'ERROR!! SMFNBC(1)+SMFNBC(2)+SMFNBC(3) should be 1'
      IF (SUMNBD.LT.0.9999.OR.SUMNBD.GT.1.0001)
     *  STOP 'ERROR!! SMFNBD(1)+SMFNBD(2)+SMFNBD(3) should be 1'
      IF (SUMNBG.LT.0.9999.OR.SUMNBG.GT.1.0001)
     *  STOP 'ERROR!! SMFNBG(1)+SMFNBG(2)+SMFNBG(3) should be 1'
C07
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      WRITE(2,999)
      READ(1,*) SMFPBC(1),SMFPBC(2),SMFPBC(3),SMFPBD(1),SMFPBD(2),
     *  SMFPBD(3),SMFPBG(1),SMFPBG(2),SMFPBG(3)
      WRITE(2,50)'* Cyanobacteria-P split into G1, G2 & G3 classes  '
      WRITE(2,51)' : (FPBc1, FPBc2, FPBc3) = ', (SMFPBC(M),M=1,3)
      WRITE(2,50)'* Diatoms-P split into G1, G2 & G3 classes        '
      WRITE(2,51)' : (FPBd1, FPBd2, FPBd3) = ', (SMFPBD(M),M=1,3)
      WRITE(2,50)'* Blue-green algae-P split into G1, G2, G3 classes'
      WRITE(2,51)' : (FPBg1, FPBg2, FPBg3) = ', (SMFPBG(M),M=1,3)
      SUMPBC=SMFPBC(1)+SMFPBC(2)+SMFPBC(3)
      SUMPBD=SMFPBD(1)+SMFPBD(2)+SMFPBD(3)
      SUMPBG=SMFPBG(1)+SMFPBG(2)+SMFPBG(3)
      IF (SUMPBC.LT.0.9999.OR.SUMPBC.GT.1.0001)
     *  STOP 'ERROR!! SMFPBC(1)+SMFPBC(2)+SMFPBC(3) should be 1'
      IF (SUMPBD.LT.0.9999.OR.SUMPBD.GT.1.0001)
     *  STOP 'ERROR!! SMFPBD(1)+SMFPBD(2)+SMFPBD(3) should be 1'
      IF (SUMPBG.LT.0.9999.OR.SUMNBG.GT.1.0001)
     *  STOP 'ERROR!! SMFPBG(1)+SMFPBG(2)+SMFPBG(3) should be 1'
C08
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      WRITE(2,999)
      READ(1,*) SMFCBC(1),SMFCBC(2),SMFCBC(3),SMFCBD(1),SMFCBD(2),
     *  SMFCBD(3),SMFCBG(1),SMFCBG(2),SMFCBG(3)
      WRITE(2,50)'* Cyanobacteria-C split into G1, G2 & G3 classes  '
      WRITE(2,51)' : (FCBc1, FCBc2, FCBc3) = ', (SMFCBC(M),M=1,3)
      WRITE(2,50)'* Diatoms-C split into G1, G2 & G3 classes        '
      WRITE(2,51)' : (FCBd1, FCBd2, FCBd3) = ', (SMFCBD(M),M=1,3)
      WRITE(2,50)'* Blue-green algae-C split into G1, G2, G3 classes'
      WRITE(2,51)' : (FCBg1, FCBg2, FCBg3) = ', (SMFCBG(M),M=1,3)
      SUMCBC=SMFCBC(1)+SMFCBC(2)+SMFCBC(3)
      SUMCBD=SMFCBD(1)+SMFCBD(2)+SMFCBD(3)
      SUMCBG=SMFCBG(1)+SMFCBG(2)+SMFCBG(3)
      IF (SUMCBC.LT.0.9999.OR.SUMCBC.GT.1.0001)
     *  STOP 'ERROR!! SMFCBC(1)+SMFCBC(2)+SMFCBC(3) should be 1'
      IF (SUMCBD.LT.0.9999.OR.SUMCBD.GT.1.0001)
     *  STOP 'ERROR!! SMFPBD(1)+SMFCBD(2)+SMFCBD(3) should be 1'
      IF (SUMCBG.LT.0.9999.OR.SUMCBG.GT.1.0001)
     *  STOP 'ERROR!! SMFCBG(1)+SMFCBG(2)+SMFCBG(3) should be 1'
C
C09 Spatially constant parameters for DIAGENESIS
C
c     READ(1,999)
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMKPON(1),SMKPON(2),SMKPON(3),SMKPOP(1),SMKPOP(2),
     *  SMKPOP(3),SMKPOC(1),SMKPOC(2),SMKPOC(3)
      WRITE(2,50)'* Diagenesis rate at 20oC in Layer 2 (/day)       '
      WRITE(2,51)' : (KPON1,KPON2,KPON3)   = ', (SMKPON(M),M=1,3)
      WRITE(2,51)' : (KPOP1,KPOP2,KPOP3)   = ', (SMKPOP(M),M=1,3)
      WRITE(2,51)' : (KPOC1,KPOC2,KPOC3)   = ', (SMKPOC(M),M=1,3)
C10
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMTHKN(1),SMTHKN(2),SMTHKN(3),SMTHKP(1),SMTHKP(2),
     *  SMTHKP(3),SMTHKC(1),SMTHKC(2),SMTHKC(3)
      WRITE(2,50)'* Temperature effect on diagenesis rate           '
      WRITE(2,51)' : (ThKN1,ThKN2,ThKN3)   = ', (SMTHKN(M),M=1,3)
      WRITE(2,51)' : (ThKP1,ThKP2,ThKP3)   = ', (SMTHKP(M),M=1,3)
      WRITE(2,51)' : (ThKC1,ThKC2,ThKC3)   = ', (SMTHKC(M),M=1,3)
C
C11 Spatially constant parameters common for SEDIMENT FLUX
C
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMM1,SMM2,SMTHDD,SMTHDP,SMPOCR,SMKMDP,SMKBST,
     *  XSMDPMIN,SMRBIBT
C12
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMO2BS,SMTDMBS,SMTCMBS
      WRITE(2,50)'* Solid concentrations (kg/L) in Layers 1 and 2   '
      WRITE(2,51)' : (rM1, rM2)            = ', SMM1,SMM2
      WRITE(2,50)'* Temp effect on mixing in dissolved & particulate'
      WRITE(2,51)' : (ThDd, ThDp)          = ', SMTHDD,SMTHDP
      WRITE(2,52)'* Half-sat. const of O2 for particle mixing= ',SMKMDP
     *          ,': First-order decay rate for stress (/day) = ',SMKBST
     *          ,'* Ratio of bio-irrigation to bioturbation  = ',SMRBIBT
     *          ,'* Reference conc (gC/m^3) for GPOC(1)      = ',SMPOCR
     *         ,'* Minimum diffusion coeff (m^2/day)        = ',XSMDPMIN
     *          ,'* Critical O2 (g/m^3) for benth. hysteresis= ',SMO2BS
     *          ,': Time lag (days) for max stress to be kept= ',SMTDMBS
     *          ,': Time duration (d) above which hysteresis = ',SMTCMBS
      ISMTDMBS = NINT(SMTDMBS/DTWQ)
      ISMTCMBS = NINT(SMTCMBS/DTWQ)
      SM1OKMDP = 1.0/SMKMDP
      SMBST1 = 1.0 / (1.0 + SMKBST*DTWQ)
C
C13 Spatially constant parameters for NH4, NO3 & PO4 FLUX
C
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMP1NH4,SMP2NH4,SMKMNH4,SMKMO2N,SMTHNH4,SMTHNO3,
     *  SMP2PO4,SMCO2PO4
      WRITE(2,50)'* Partition coeff bet/ dissolved and sorbed NH4   '
      WRITE(2,51)' : (P1NH4, P2NH4)        = ', SMP1NH4,SMP2NH4
      WRITE(2,50)'* Half-sat. const for nitri. (gN/m^3, gO2/m^3)    '
      WRITE(2,51)' : (KMNH4, KMNH4O2)      = ', SMKMNH4,SMKMO2N
      WRITE(2,50)'* Temp effect on KNH4 & KNO3                      '
      WRITE(2,51)' : (ThNH4, ThNO3)        = ', SMTHNH4,SMTHNO3
      WRITE(2,52)'* Partition coef for PO4 in Layer 2 (L/kg) = ',SMP2PO4
     *         ,': Critical DO (mg/L) for PO4 sorption      = ',SMCO2PO4
      SMFD1NH4 = 1.0 / (1.0 + SMM1*SMP1NH4)
      SMFP1NH4 = 1.0 - SMFD1NH4
      SMFD2NH4 = 1.0 / (1.0 + SMM2*SMP2NH4)
      SMFP2NH4 = 1.0 - SMFD2NH4
      SMKMO2N = SMKMO2N * 2.0
      SMFD2PO4 = 1.0 / (1.0 + SMM2*SMP2PO4)
      SMFP2PO4 = 1.0 - SMFD2PO4
C
C14 Spatially constant parameters for H2S/CH4 FLUX and SOD
C
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMP1H2S,SMP2H2S,SMKD1HS,SMKP1HS,SMTHH2S,SMKMH2S,
     *  SMKCH4,SMTHCH4,SMCSHSCH
      WRITE(2,50)'* Partition coeff for H2S in Layer 1 (L/kg)       '
      WRITE(2,51)' : (P1H2S, P2H2S)        = ', SMP1H2S,SMP2H2S
      WRITE(2,50)'* Reaction vel (m/d) for dissol & part. in Layer 1'
      WRITE(2,51)' : (KH2Sd1, KH2Sp1)      = ', SMKD1HS,SMKP1HS
      WRITE(2,52)'* Critical sal (ppt) for H2S/CH4 oxidation = ',
     *  SMCSHSCH
      WRITE(2,52)'* Temperature effect on H2S oxidation rate = ',SMTHH2S
     *          ,': Oxygen effect (mg/L) on H2S oxidation    = ',SMKMH2S
      WRITE(2,52)'* Methane oxidation reaction velocity (m/d)= ',SMKCH4
     *          ,': Temperature effect on CH4 oxidation rate = ',SMTHCH4
      SMFD1H2S = 1.0 / (1.0 + SMM1*SMP1H2S)
      SMFP1H2S = 1.0 - SMFD1H2S
      SMFD2H2S = 1.0 / (1.0 + SMM2*SMP2H2S)
      SMFP2H2S = 1.0 - SMFD2H2S
      XSMK1H2S = (SMKD1HS*SMKD1HS*SMFD1H2S + SMKP1HS*SMKP1HS*SMFP1H2S)
     *  / (2.0*SMKMH2S)
!
c     J.S. Check this part. The following block should removed. !!!
      SMFD1NH4 = 1.0 / (1.0 + SMM1*SMP1NH4)
      SMFP1NH4 = 1.0 - SMFD1NH4
      SMFD2NH4 = 1.0 / (1.0 + SMM2*SMP2NH4)
      SMFP2NH4 = 1.0 - SMFD2NH4
      SMKMO2N = SMKMO2N * 2.0
      SMFD2PO4 = 1.0 / (1.0 + SMM2*SMP2PO4)
      SMFP2PO4 = 1.0 - SMFD2PO4
c
      SMFD1H2S = 1.0 / (1.0 + SMM1*SMP1H2S)
      SMFP1H2S = 1.0 - SMFD1H2S
      SMFD2H2S = 1.0 / (1.0 + SMM2*SMP2H2S)
      SMFP2H2S = 1.0 - SMFD2H2S
      XSMK1H2S = (SMKD1HS*SMKD1HS*SMFD1H2S + SMKP1HS*SMKP1HS*SMFP1H2S)
     *  / (2.0*SMKMH2S)

!

C15
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMO2C,SMO2NO3,SMO2NH4
      WRITE(2,52)'* Stoichi coef for C used by H2S ox (gO2/gC)=',SMO2C
     *          ,': Stoichi coef for C used by denitr (gO2/gN)=',SMO2NO3
     *          ,': Stoichi coef for O2 used by nitri (gO2/gN)=',SMO2NH4
C
C16 Spatially constant parameters for SILICA
C
      WRITE(2,999)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      IF (isSKIP .EQ. 0) READ(1,999)
      READ(1,*) SMKSI,SMTHSI,SMKMPSI,SMSISAT,SMP2SI,SMDP1SI,SMCO2SI,
     *  SMJDSI
      WRITE(2,52)'* PSi dissol. rate at 20C in Layer 2 (/d)  = ',SMKSI
     *          ,': Temperature effect on PSi dissolution    = ',SMTHSI
     *          ,': Sat. conc. in pore water (g Si/m^3)      = ',SMSISAT
     *          ,'* Partition coef for Si in Layer 2 (L/kg)  = ',SMP2SI
     *          ,': Incremental in part. coef in Layer1, P1Si= ',SMDP1SI
     *          ,': Critical DO (mg/L) for Si sorption       = ',SMCO2SI
     *          ,'* Detrital flux (g/m^2/d) except diatoms   = ',SMJDSI
     *          ,'* Dissolution half-sat constant (g Si/m^3) = ',SMKMPSI
      SMFD2SI = 1.0 / (1.0 + SMM2*SMP2SI)
      SMFP2SI = 1.0 - SMFD2SI
C
C Set up look-up table for temperature dependency over -5oC to 35oC
C
      DO IT=1,NWQTD
        STEMP = REAL(IT-1)*0.1 - 4.95  !It's -5C -> 50C, since NWQTD=550, Ji, 8/29/02
c        STEMP = REAL(M-1)*0.1 - 4.95
        TT20 = STEMP-20.0
        DO M=1,3
          SMTDND(IT,M) = SMKPON(M) * SMTHKN(M)**TT20
          SMTDPD(IT,M) = SMKPOP(M) * SMTHKP(M)**TT20
          SMTDCD(IT,M) = SMKPOC(M) * SMTHKC(M)**TT20
        END DO
        SMTDDP(IT) = SMTHDP**TT20
        SMTDDD(IT) = SMTHDD**TT20
        SMTDNH4(IT) = SMTHNH4**TT20
        SMTDNO3(IT) = SMTHNO3**TT20
        SMK1H2S(IT)= XSMK1H2S * SMTHH2S**TT20
        SMTD1CH4(IT) = 0.97656**TT20 * 20.0
        SMTD2CH4(IT) = SMKCH4 * SMTHCH4**TT20
c	write(302,*) IT,SMTD2CH4(IT),SMKCH4,SMTHCH4,tt20
        SMTDSI(IT) = SMKSI * SMTHSI**TT20
      END DO
C
C17 Spatially constant ICs for sediment state variables (g/m^3)
C
      WRITE(2,998)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      READ(1,*) SMPON(1,1),SMPON(1,2),SMPON(1,3),SMPOP(1,1),
     *  SMPOP(1,2),SMPOP(1,3),SMPOC(1,1),SMPOC(1,2),SMPOC(1,3)
     *  ,SEDNF, SEDFF
      IF (ISMICI.NE.1 .AND. ISMICI.NE.2)
     *  WRITE(2,5105) SMPON(1,1),SMPON(1,2),SMPON(1,3),SMPOP(1,1),
     *  SMPOP(1,2),SMPOP(1,3),SMPOC(1,1),SMPOC(1,2),SMPOC(1,3)
C18
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      READ(1,*) SM1NH4(1),SM2NH4(1),SM2NO3(1),SM2PO4(1),SM2H2S(1),
     *  SMPSI(1),SM2SI(1),SMBST(1),SMT(1)
      IF (ISMICI.NE.1 .AND. ISMICI.NE.2) THEN
        WRITE(2,5105) SM1NH4(1),SM2NH4(1),SM2NO3(1),SM2PO4(1),SM2H2S(1),
     *    SMPSI(1),SM2SI(1),SMBST(1),SMT(1)
        DO L=2,LA
          DO M=1,NSMG
            SMPON(L,M)=SMPON(1,M)
            SMPOP(L,M)=SMPOP(1,M)
            SMPOC(L,M)=SMPOC(1,M)
          END DO
          SM1NH4(L)=SM1NH4(1)
          SM2NH4(L)=SM2NH4(1)
          SM2NO3(L)=SM2NO3(1)
          SM2PO4(L)=SM2PO4(1)
          SM2H2S(L)=SM2H2S(1)
          SMPSI(L) =SMPSI(1)
          SM2SI(L) =SM2SI(1)
          SMBST(L) =SMBST(1)
          SMT(L)   =SMT(1)
        END DO
      END IF
C
C19 Spatially varying physical & rate velocity: convert W2 from cm/yr to m/d
C   SMDIFT in m^2/d
C
      WRITE(2,998)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      DO I=1,ISMZ
        READ(1,*) MM,SMHSED(I),SMW2(I),SMDD(I),SMDP(I),SMKNH4(I),
     *    SMK1NO3(I),SMK2NO3(I),SMDP1PO4(I) !     , sodmult(i),psload10
        WRITE(2,56) MM,SMHSED(I),SMW2(I),SMDD(I),SMDP(I),SMKNH4(I),
     *    SMK1NO3(I),SMK2NO3(I),SMDP1PO4(I), sodmult(i)
        SMW2(I) = SMW2(I)*SMCW2
        SMDTOH(I) = DTWQ/SMHSED(I)
        SMHODT(I) = SMHSED(I)/DTWQ
        SMDP(I) = SMDP(I) / (SMHSED(I)*SMPOCR+ 1.E-18)
        SMDD(I) = SMDD(I) / (SMHSED(I)+ 1.E-18)
        SMKNH4(I) = SMKNH4(I)*SMKNH4(I) * SMKMNH4
        SMK1NO3(I) = SMK1NO3(I)*SMK1NO3(I)
        SM1DIFT(I) = SMDIFT * SMDTOH(I)/(SMHSED(I)+ 1.E-18)
        SM2DIFT(I) = 1.0 / (1.0 + SM1DIFT(I))
        SMW2DTOH(I) = 1.0 + SMW2(I)*SMDTOH(I)
        SMW2PHODT(I) = SMW2(I) + SMHODT(I)
        SMDPMIN(I) = XSMDPMIN / (SMHSED(I)+ 1.E-18)
      END DO
C
C20 Spatially varying parameters: distribution coefficients for RPOM
C
      WRITE(2,998)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      IF (isSKIP .GT. 0) CALL SKIPCOMM(1,CCMRM)
      READ(1,5100) TITLE(1)
      WRITE(2,5100) TITLE(1)
      DO I=1,ISMZ
        READ(1,*) MM,SMFNR(I,1),SMFNR(I,2),SMFNR(I,3),SMFPR(I,1),
     *    SMFPR(I,2),SMFPR(I,3),SMFCR(I,1),SMFCR(I,2),SMFCR(I,3)
        WRITE(2,54) MM,SMFNR(I,1),SMFNR(I,2),SMFNR(I,3),SMFPR(I,1),
     *    SMFPR(I,2),SMFPR(I,3),SMFCR(I,1),SMFCR(I,2),SMFCR(I,3)
        SUMNBC=SMFNR(I,1)+SMFNR(I,2)+SMFNR(I,3)
        SUMNBD=SMFPR(I,1)+SMFPR(I,2)+SMFPR(I,3)
        SUMNBG=SMFCR(I,1)+SMFCR(I,2)+SMFCR(I,3)
        IF (SUMNBC.LT.0.9999.OR.SUMNBC.GT.1.0001)
     *    STOP 'ERROR!! SMFNR(I,1)+SMFNR(I,2)+SMFNR(I,3) should be 1'
        IF (SUMNBD.LT.0.9999.OR.SUMNBD.GT.1.0001)
     *    STOP 'ERROR!! SMFPR(I,1)+SMFPR(I,2)+SMFPR(I,3) should be 1'
        IF (SUMNBG.LT.0.9999.OR.SUMNBG.GT.1.0001)
     *    STOP 'ERROR!! SMFCR(I,1)+SMFCR(I,2)+SMFCR(I,3) should be 1'
      END DO
C
      CLOSE(1)
 6666 format(a30)
C
  998 FORMAT(80X)
 5105 FORMAT(10F8.2)
   54 FORMAT(I8, 10F8.3)
   56 FORMAT(I8, 3F8.3, E8.1, 6F8.3)
C
C Read in mapping infor. for spatially-varying SED parameters (unit #7).
C
      DO L=2,LA
       ISMZMAP(L)=1
      END DO
C
c      IF(ISWQSMAP.EQ.1) THEN
      if (ismz .gt. 1) then
C
      OPEN(1,FILE='wqsdmap.inp',STATUS='UNKNOWN')
C
        WRITE(2,999)
        READ(1,90) (TITLE(M), M=1,3)
        WRITE(2,90) (TITLE(M), M=1,3)
        READ(1,999)
        WRITE(2,999)
        WRITE(2,92)
        IN=0
        IJC=IC*JC
        DO M=1,IJC
          READ(1,*,END=1111) I,J,ISMZX
          IN=IN+1
          IF (LIJW(I,J).LT.1) THEN
            PRINT*, 'i, j, IJCT(i,j) = ', I,J,LIJW(I,J)
            STOP 'ERROR!! invalid (i,j) in FILE wqsdmap.inp'
          END IF
          L = LIJW(I,J)
          ISMZMAP(L)=ISMZX
          WRITE(2,91) L,I,J,ISMZMAP(L)
        END DO
 1111   CONTINUE
        IF (IN.NE.(LA-1)) THEN
          PRINT*, 'all active sed. cells should be mapped for SED par.'
          STOP 'ERROR!! number of lines in FILE wqsdmap.inp =\ (LA-1)'
        END IF

        CLOSE(1)

      END IF
C
      CLOSE(2)
C
   90 FORMAT(A79)
   91 FORMAT(15I5)
   92 FORMAT('    L    I    J    ISMZMAP')
C
      RETURN
      END
