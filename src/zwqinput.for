!**********************************************************************C
! 
      SUBROUTINE RWQINP(IWQDT,DT,IC,JC,KC,LA,NTSPTC)
! 
!**********************************************************************C
!
! **  LAST MODIFIED BY JOHN HAMRICK ON  7 APRIL 1997
!
! **  MODIFIED BY J.S. ON 5/5/98 TO ADD A MACALGAL.
!**********************************************************************C
!
! Read in from the unit #8
!: I/O control variables
!: spatially and temporally constant REAL parameters
!
!    NTSPTC = number timestep per referecen cycle
!    KC     = total layers
!    LA     = total cell in 1D arrary
!    DT     = time step in second
!**********************************************************************C
!
      INCLUDE 'wq.par'
      INCLUDE 'wqcom.cmn'
!
      PARAMETER (CONV1=1.0E3,CONV2=8.64E4)
      DIMENSION XPSL(NWQVM),XDSL(NWQVM)
      DIMENSION WQTMC1(NWQZM),WQTMC2(NWQZM),WQTMD1(NWQZM),WQTMD2(NWQZM),
     &  WQTMG1(NWQZM),WQTMG2(NWQZM),WQTMM1(NWQZM),WQTMM2(NWQZM),
     $  WQTMp1(NWQZM), WQTMp2(NWQZM),
     &  WQKG1C(NWQZM),WQKG2C(NWQZM),WQKG1D(NWQZM),WQKG2D(NWQZM),
     &  WQKG1G(NWQZM),WQKG2G(NWQZM),WQKG1M(NWQZM),WQKG2M(NWQZM),
     $         WQKG1p(NWQZM), WQKG2p(NWQZM)
      INTEGER isSKIP
      CHARACTER TITLE(5)*79, ccmrm*1, ccmr*1
!
      OPEN(2,FILE='wq3d.out',STATUS='UNKNOWN',ACCESS='APPEND')
!
      OPEN(1,FILE='wq3dwc.inp',STATUS='UNKNOWN')
!
      DO I=1,21
	IWTRC(I)=1
	ENDDO
	IWTRC(20)=0
	IWTRC(21)=0
!C01
      CCMR='C'
	isSKIP=0
      CALL SKIPCOMM(1,CCMRM)
      READ(1,90) TITLE(1)
      WRITE(2,90) TITLE(1)

!C02 I/O conrol variables

      WRITE(2,90)'C02 CONTROL VARIABLE CARD'
      CALL SKIPCOMM(1,CCMRM)
!               1    2    3      4    5      6      7    8    9
      READ(1,*) NWQV,NWQZ,NWQPS,NWQTD,NWQTS,NTSWQV,NSMG,NSMZ,NTSSMV,
     $          NSMTS,NWQKDPT,ITRYCHLA 

C use efdc.par values, Ji, 9/17/99
c     NWQV   =  NWQVM           ! 1
      NWQZ   =  NWQZM           ! 2
      NWQPS  =  NWQPSM          ! 3
c     NWQTD  =  NWQTDM          ! 4 use 550
      NWQTS  =  NWQTSM          ! 5
c     NTSWQV =  NTSWQVM         ! 6
      NSMG   =  NSMGM           ! 7
      NSMZ   =  NSMZM           ! 8
      NTSSMV =  NTSSMVM         ! 9
      NSMTS  =  NSMTSM          ! 10
c check array dimensions, jeff Ji, 8/5/99
      if(NWQV  .gt.NWQVM  ) then   ! 1
      write(6,*) " Error:  NWQV  >  NWQVM  ", NWQV  , NWQVM
      stop
      endif
      if(NWQZ  .gt.NWQZM  ) then   ! 2
      write(6,*) " Error:  NWQZ  >  NWQZM  ", NWQZ  , NWQZM
      stop
      endif  
      if(NWQPS .gt.NWQPSM ) then   ! 3
      write(6,*) " Error:  NWQPS >  NWQPSM ", NWQPS , NWQPSM
      stop
      endif
      if(NWQTD .gt.NWQTDM ) then   ! 4
      write(6,*) " Error:  NWQTD >  NWQTDM ", NWQTD , NWQTDM
      stop
      endif
      if(NWQTS .gt.NWQTSM ) then   ! 5
      write(6,*) " Error:  NWQTS >  NWQTSM ", NWQTS , NWQTSM
      stop
      endif
      if(NTSWQV.gt.NWQVM) then   ! 6
      write(6,*) " Error:  NTSWQV>  NWQVM", NTSWQV, NWQVM
      stop
      endif
      if(NSMG  .gt.NSMGM  ) then   ! 7
      write(6,*) " Error:  NSMG  >  NSMGM  ", NSMG  , NSMGM
      stop
      endif
      if(NSMZ  .gt.NSMZM  ) then   ! 8 
      write(6,*) " Error:  NSMZ  >  NSMZM  ", NSMZ  , NSMZM
      stop
      endif
      if(NTSSMV.gt.NTSSMVM) then   ! 9
      write(6,*) " Error:  NTSSMV>  NTSSMVM", NTSSMV, NTSSMVM
      stop
      endif
      if(NSMTS .gt.NSMTSM ) then   ! 10
      write(6,*) " Error:  NSMTS >  NSMTSM ", NSMTS , NSMTSM
      stop
      endif
!
! Check array dimensions
!
!C03
      WRITE(2,*)'C03 '
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQDT,IWQM,IWQBEN,IWQSI,IWQFCB,IWQSRP,IWQSTOX, IWQKA,
     $   WQKRO,iSHADE
      IWQDT0=IWQDT
	NWQKDPT=IWQDT /2   ! Double check J.S. if not 3-level time steping
!
!C04
!
      WRITE(2,*)'C04'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQZ,IWQNC,IWQRST,NDMWQ,LDMWQ,NDDOAVG,NDLTAVG,IDNOTRVA
c
!      if(IC9.gt.0) then   ! Ji for grid.inp, 9/10/99
!      LDMWQ=LDM9
!      endif
c check array dimensions, jeff Ji, 8/5/99
!      if(LDMWQ .gt.(LCMWQ-2)) then   ! 1
!      write(6,*) " Error:  LDMWQ >(LCMWQ-2)", LDMWQ , LCMWQ
!      stop
!      endif
c5
      WRITE(2,*)'C05'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQICI,IWQAGR,IWQSTL,IWQSUN,IWQPSL,IWQNPL, isDIURDO,
     $          WQDIUDT
c     $          ISDIURDO,  ISDIUDIA, ICDIA, JCDIA
c     isDIURDO > 0 turns on diurnal DO output to binary file
c     WQDIUDT = time interval for writing to diurnal DO file (hours)
c
      IWQDIUDT = NINT(WQDIUDT*3600.0/DT)
      WRITE(2,83)': Frequency of Diurnal DO output (in DT unit) =',
     $  IWQDIUDT
C
      WRITE(2,83)'* IWQDT (DTWQ(d) = DT(s)*IWQDT/86400)        = ',IWQDT
      DTD = DT/86400.0
      DTWQ = DTD*REAL(IWQDT) !*REAL(NWQKDPT)             ! WQ model DT, usually = 2*DTD
      DTWQO2 = DTWQ*0.5
C
      IF (IWQM.EQ.1) THEN
        WRITE(2,80)'* Full version with 21 variables is activated     '
       ELSE IF (IWQM.EQ.2) THEN
        WRITE(2,80)'* Small version with 9 variables is activated     '
       ELSE
        STOP '** ERROR!!! Invalid IWQM value **'
      END IF
C
      IF (IWQBEN.EQ.1) THEN
        WRITE(2,80)'* Sediment process model is activated             '
       ELSE IF (IWQBEN.EQ.0) THEN
        WRITE(2,80)'* Spatially/temporally constant BF is specified   '
       ELSE IF (IWQBEN.EQ.2) THEN
        WRITE(2,80)'* Spatially and/or temporally-varying BF specified'
       ELSE
        STOP '** ERROR!!! Invalid IWQBEN value **'
      END IF
C
      IF (IWQSI.EQ.1) THEN
        WRITE(2,80)'* Silica state variables (SU & SA) are modeled    '
        IF (IWQM.EQ.2) STOP '** ERROR!!! Incompatible IWQM & IWQSI'
       ELSE
        WRITE(2,80)'* No Silica (SU & SU) limitation                  '
      END IF
C
      IF (IWQFCB.EQ.1) THEN
        WRITE(2,80)'* FCB (fecal coliform bacteria) is modeled        '
       ELSE
        WRITE(2,80)'* FCB (fecal coliform bacteria) is not modeled    '
      END IF
C
      IF (IWQSRP.EQ.1) THEN
        WRITE(2,80)'* TAM is used for sorption of PO4t/SA: model TAM  '
       ELSE IF (IWQSRP.EQ.2) THEN
        WRITE(2,80)'* TSS is used for sorption of PO4t/SA: model TSS  '
       ELSE
        WRITE(2,80)'* No sorption of PO4t/SA: may model TSS & no TAM  '
      END IF
C
      IF (IWQSTOX.EQ.1) THEN
        WRITE(2,80)'* Salinity toxicity is applied to cyanobacteria   '
       ELSE
        WRITE(2,80)'* No salinity toxicity: saltwater cyanobacteria   '
      END IF
C
      IF (IWQKA.EQ.0) THEN
        WRITE(2,80)'* User-specified constant reaeration set to WQKRO '
        WRITE(2,80)'*   reaeration due to wind set to zero            '
      END IF
      IF (IWQKA.EQ.1) THEN
        WRITE(2,80)'* User-specified constant reaeration set to WQKRO '
        WRITE(2,80)'*   reaeration due to wind added to WQKRO         '
      END IF
      IF (IWQKA.EQ.2) THEN
        WRITE(2,80)'* OConnor-Dobbins reaeration formula is used      '
      END IF
      IF (IWQKA.EQ.3) THEN
        WRITE(2,80)'* Owens & Gibbs (1964) reaeration formula is used '
      END IF
      IF (iSHADE .EQ. 0) THEN
        WRITE(2,80)'* iSHADE=0, shade factor set to 1.0 for all cells '
      END IF
      DO LL = 2, LA
        pSHADE(LL) = 1.0
      END DO
      IF (iSHADE .GT. 0) THEN
        WRITE(2,80)'* iSHADE>0, shade factors read from SHADEMAP.INP  '
        OPEN(3, FILE='shademap.inp', STATUS='UNKNOWN')
        call skipcomm(3, ccmrm)
9001    READ(3,*,END=9002) II, JJ, xmrm
        IF (II .LE. 0) GO TO 9002
        IF (LIJW(II,JJ).LT.1) THEN
          PRINT*, 'i, j, IJCT(i,j) = ', II,JJ,LIJW(II,JJ)
          STOP 'ERROR!! invalid (i,j) in file SHADEMAP.INP'
        END IF
        LL=LIJW(II,JJ)
        pSHADE(LL)=xmrm
9002    CLOSE(3)
      END IF
C
      WRITE(2,83)'* # of zones for spatially varying parameters =',IWQZ
      IF (IWQZ.GT.NWQZ) STOP 'ERROR!! IWQZ should be <= NWQZ'
C
      IF (IWQNC.EQ.1) THEN
        WRITE(2,80)'* Write negative conc. information to NEG-CONC.LOG'
       ELSE
        WRITE(2,80)'* No wrting of negative concentration information '
      END IF
C
      IF (IWQRST.GE.1) THEN
        WRITE(2,80)'* Write spatial distributions to IWQORST          '
       ELSE
        WRITE(2,80)'* No writing to IWQORST                           '
      END IF
C
      WRITE(2,999)
C
      IF (IWQICI.EQ.1) THEN
        WRITE(2,80)'* Spatially/temporally-varying ICs from INWQICI   '
       ELSE IF (IWQICI.GE.2) THEN
        WRITE(2,80)'* Spatially/temporally-varying ICs from INWQRST   '
       ELSE
        WRITE(2,80)'* Spatially/temporally constant initial conditions'
      END IF
C
      IF (IWQAGR.EQ.1) THEN
        WRITE(2,80)'* Spatially a/o temporally-varying algal kinetics '
       ELSE
        WRITE(2,80)'* Spatially/temporally constant algal kinetics    '
      END IF
C
      IF (IWQSTL.EQ.1) THEN
        WRITE(2,80)'* Spatially and/or temporally-varying settling vel'
       ELSE
        WRITE(2,80)'* Spatially/temporally constant settling velocity '
      END IF
C
      IF (IWQSUN.GE.1) THEN
        WRITE(2,80)'* Temporally-varying Io & FD                      '
       ELSE
        WRITE(2,80)'* Temporally constant Io & FD                     '
      END IF
C
c          IF (IWQPSL.EQ.1) THEN
c             WRITE(2,80)'* Temporally-varying point source input           '
c            ELSE
c             WRITE(2,80)'* Temporally constant point source input          '
c           END IF
C
      IF (IWQNPL.EQ.1) THEN
        WRITE(2,80)'* Spatially and/or temporally-varying NPS input   '
       ELSE
        WRITE(2,80)'* Spatially/temporally constant NPS input         '
      END IF
C6
      WRITE(2,*)'C6-------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQTS,TWQTSB,TWQTSE,WQTSDT, isWQAVG, isWQMIN, isWQMAX,
     *          isCOMP,isBIN, isTBIN,isFLUX,IAVGBIN,IRESTHR
      if(isBIN.eq.1) write(*,*)'Out put 3D resutls'
c isWQAVG > 0 turns on binary file output for WQ daily averages
c isWQMIN > 0 turns on binary file output for WQ daily minimums
c isWQMAX > 0 turns on binary file output for WQ daily maximums
c isCOMP  > 0 turns on binary file output for DO component analysis
!      
      IAVGBIN0=IAVGBIN	
	ISMTSDT=NTSPTC/24*IAVGBIN          ! in hours
	IAVGBIN=NTSPTC/24*IAVGBIN
      IRESTHR=NTSPTC/24*IRESTHR

      IF (isCOMP .GT. 0) THEN
        CALL WQzero3(LA,KC)
        CALL INITbin3(TBEGIN, DT,LA,KC)
      END IF
C
      IF (IWQTS.GT.NWQTS) THEN
        WRITE(2,80)'** IWQTS should be <= NWQTS **                    '
        IWQTS=NWQTS
      END IF
C
C
      WRITE(2,84)
     * '* Time-series output from ', TWQTSB, ' day ',
     * '                       to ', TWQTSE, ' day ',
     * '                    every ', WQTSDT, ' hour',
     * '                       at ', IWQTS,  ' locations',
     * ' Bin File switch isWQAVG =', isWQAVG,' (0=off)  ',
     * ' Bin File switch isWQMIN =', isWQMIN,' (0=off)  ',
     * ' Bin File switch isWQMAX =', isWQMAX,' (0=off)  ',
     * ' Bin File switch isCOMP  =', isCOMP, ' (0=off)  '
      WRITE(2,999)
C
C7
      WRITE(2,*)'C7 Print Controal --------------'
      CALL SKIPCOMM(1,CCMRM)
C
      IF (IWQTS.GE.1) THEN
        WRITE(2,80)': ICWQTS(I)=1, time-series output for variable I  '
        WRITE(2,80)': ICWQTS(I)\=1, no time-series output for var. I  '
        WRITE(2,999)
C
        DO M=1,IWQTS
          READ(1,*) II,JJ
		DO NW=1,22 
		 ICWQTS(NW,M)=1
		ENDDO 
          IF (LIJW(II,JJ).LT.1 ) THEN
            WRITE(2,86)  II,JJ,M
            WRITE(2,80)'ERROR!! invalid (i,j): time-series location'
            STOP
          END IF
          LWQTS(M)=LIJW(II,JJ)
          WRITE(2,94) II,JJ
        END DO
C
      END IF
C
      IWQTSB = NINT(TWQTSB/DTD)
      IWQTSE = NINT(TWQTSE/DTD)
      IWQTSDT = NINT(WQTSDT*3600.0/DT)
      WRITE(2,999)
      WRITE(2,83)': Time-series starting time step (in DT unit) =',
     *  IWQTSB
      WRITE(2,83)': Time-series ending time step (in DT unit)   =',
     *  IWQTSE
      WRITE(2,83)': Frequency of TS output  (in DT unit)        =',
     *  IWQTSDT
      IF (MOD(IWQTSDT,IWQDT).NE.0)
     *  STOP 'ERROR!! IWQTSDT should be multiple of IWQDT'
C
  999 FORMAT(1X)
   90 FORMAT(A79)
   91 FORMAT(10I8)
   92 FORMAT(10F8.4)
   93 FORMAT(I8,3F8.4)
   94 FORMAT(2I5, 13I5, /, 10X, 9I5)
   80 FORMAT(A50)
   81 FORMAT(A27, 4(F8.4,2X))
   82 FORMAT((A45, F8.4))
   83 FORMAT(A47, I10)
   84 FORMAT(3(A26,F10.4,A5,/), 5(A26,I8,A10,/))
   86 FORMAT(' I,J,M = ',3I10)
!
!C08: Constant parameters for ALGAE (see Table 3-1)
!
      WRITE(2,*)'C08-------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQKHNC,WQKHND,WQKHNG,WQKHNM,WQKHPC,WQKHPD,WQKHPG,WQKHPM,
     1 WQKHS,WQSTOX
      WRITE(2,80)'* Half-sat. cosntant (g/m^3) for nutrient uptake  '
      WRITE(2,81)' : (KHNc, KHPc)          = ', WQKHNC,WQKHPC
      WRITE(2,81)' : (KHNd, KHPd, KHS)     = ', WQKHND,WQKHPD,WQKHS
      WRITE(2,81)' : (KHNg, KHPg)          = ', WQKHND,WQKHPG
      WRITE(2,81)' : (KHNm, KHPm)          = ', WQKHNM,WQKHPM
      WRITE(2,82)'* Sal. where Microsystis growth is halved  = ', WQSTOX
      WQSTOX = WQSTOX*WQSTOX
!
!C09:constant parameters for ALGAE (see Table 3-1)
!    
      Write(2,*)'C09------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      WRITE(2,80)'* Carbon-to-Chl ratio (g C per mg Chl)            '
  
      DO I=1,IWQZ
      READ(1,*) WQKETSS(I),WQKECHL(I),WQKESAT(I),WQKEB(I),WQCHLC(I),
     * WQCHLD(I),WQCHLG(I),WQCHLM,WQDOPC,WQDOPD,WQDOPG,WQDOPM
      WRITE(2,80)'* Light extinc. coeff. due to TSS & Chl           '
      WRITE(2,81)' : KeTSS (/m per g/m^3)  = ', WQKETSS
      WRITE(2,81)' : KeChl (/m per mg/m^3) = ', WQKECHL
      IF (WQKECHL(I) .LT. 0.0) THEN
        WRITE(2,80) '* Use Riley (1956) equation for WQKECHL          '
        WRITE(2,80) ' : KeChl = 0.054*CHL**0.667 + 0.0088*CHL         '
      END IF
      WQCHLC(I)=1.0/(WQCHLC(I)+ 1.E-12)    ! be careful! inversed!, 9/15/02
      WQCHLD(I)=1.0/(WQCHLD(I)+ 1.E-12)
      WQCHLG(I)=1.0/(WQCHLG(I)+ 1.E-12)
      WQCHLM=1.0/(WQCHLM+ 1.E-12)  
      WRITE(2,81)' : (CChlc, CChld, CChlg) = ',    
     & WQCHLC(I),WQCHLD(I),WQCHLG(I)         
      ENDDO
      WRITE(2,80)'* Depth (m) of maximum algal growth               '
      WRITE(2,81)' : (DOPTc, DOPTd, DOPTg) = ', WQDOPC,WQDOPD,WQDOPG

!C10:constant parameters for ALGAE (see Table 3-1)
!   

	  
      Write(2,*)'C10------------------------------'   
	CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQI0,WQISMINC,WQISMIN,WQISMING,WQFD,WQCIA,WQCIB,WQCIC,
     &          WQCIM,PARADJ
      WRITE(2,82)'* Initial Io (ly/d) at water surface       = ',WQI0
     *          ,'  Minimum optimum solar radiation (ly/d)   = ',WQISMIN
     *          ,'  Fractional daylength                     = ',WQFD
     *          ,'  Weighting factor for rad. at current day = ',WQCIA
     *          ,'  Weighting factor for rad. at (-1) day    = ',WQCIB
     *          ,'  Weighting factor for rad. at (-2) days   = ',WQCIC
     *          ,'  Fraction of Solar Radiation that is PAR  = ',PARADJ
      WQI1 = WQI0
      WQI2 = WQI0
      WQI3 = WQI0
      WQI0opt = 0.0
 !     IF (IWQSUN .EQ. 2) THEN   ! Need revise J.S. !!!
 !       WQISMIN = 0.0
 !     END IF

!      DO I=1,IWQZ
!      REAC(I)=REAC(1)
!	ENDDO
C
c      READ(1,999)
c      READ(1,*) WQTMC,WQTMD,WQTMG,WQKG1C,WQKG2C,WQKG1D,WQKG2D,WQKG1G,
c     *  WQKG2G
c      WRITE(2,80)'* Optimum temperature for algal growth (degC)     '
c      WRITE(2,81)' : (TMc, TMd, TMg      ) = ', WQTMC,WQTMD,WQTMG
c      WRITE(2,80)'* Suboptimal temperature effect for algal growth  '
c      WRITE(2,81)' : (KTG1c, KTG1d, KTG1g) = ', WQKG1C,WQKG1D,WQKG1G
c      WRITE(2,80)'* Superoptimal temperature effect for algal growth'
c      WRITE(2,81)' : (KTG2c, KTG2d, KTG2g) = ', WQKG2C,WQKG2D,WQKG2G
c 
C11 
C
      Write(2,*)'C11------------------------------'  
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ
      READ(1,*)WQTMC1(I),WQTMC2(I),WQTMD1(I),WQTMD2(I),WQTMG1(I)
     & ,WQTMG2(I),WQTMM1(I),WQTMM2(I),WQTMp1(I), WQTMp2(I)
      ENDDO
C12
      Write(2,*)'C12------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ      
      READ(1,*)WQKG1C(I),WQKG2C(I),WQKG1D(I),WQKG2D(I),WQKG1G(I)
     &  ,WQKG2G(I),WQKG1M(I),WQKG2M(I),WQKG1p(I), WQKG2p(I)

      WRITE(2,80)'* Lower Optimum temp for algal growth (degC)     '
      WRITE(2,81)' : (TMc1, TMd1, TMg1   ) = ',
     &  WQTMC1(I),WQTMD1(I),WQTMG1(I)
      WRITE(2,80)'* Upper Optimum temp for algal growth (degC)     '
      WRITE(2,81)' : (TMc2, TMd2, TMg2   ) = ', 
     & WQTMC2(I),WQTMD2(I),WQTMG2(I)
      ENDDO
C13
C   
      Write(2,*)'C13------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQTRC,WQTRD,WQTRG,WQTRM,WQKTBC,WQKTBD,WQKTBG,WQKTBM
      WRITE(2,80)'* Reference temperature for algal metabolism (oC) '
      WRITE(2,81)' : (TRc, TRd, TRg)       = ', WQTRC,WQTRD,WQTRG
      WRITE(2,80)'* Temperature effect for algal metabolism         '
      WRITE(2,81)' : (KTBc, KTBd, KTBg)    = ', WQKTBC,WQKTBD,WQKTBG
C
C Constant parameters for CARBON (see Table 3-2)
C14
      Write(2,*)'C14------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFCRP,WQFCLP,WQFCDP,WQFCDC,WQFCDD,WQFCDG,
     * WQKHRC,WQKHRD,WQKHRG
      WRITE(2,80)'* Carbon distribution coeff for algal predation   '
      WRITE(2,81)' : (FCRP, FCLP, FCDP)    = ', WQFCRP,WQFCLP,WQFCDP
      WRITE(2,80)'* Carbon distribution coeff for algal metabolism  '
      WRITE(2,81)' : (FCDc, FCDd, FCDg)    = ', WQFCDC,WQFCDD,WQFCDG
      WRITE(2,80)'* Half-sat. constant (gO/m*3) for algal DOC excret'
      WRITE(2,81)' : (KHRc, KHRd, KHRg)    = ', WQKHRC,WQKHRD,WQKHRG
      CFCDCWQ = 1.0 - WQFCDC
      CFCDDWQ = 1.0 - WQFCDD
      CFCDGWQ = 1.0 - WQFCDG
      XC = ABS(1.0 - (WQFCRP+WQFCLP+WQFCDP))
      IF (XC .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FCRP+FCLP+FCDP not equal to 1.0'
        WRITE(2,*)
      END IF
C   
C
C15
      Write(2,*)'C15------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFCRPM,WQFCLPM,WQFCDPM,WQFCDM,WQKHRM
C16
      Write(2,*)'C16------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQKRC,WQKLC,WQKDC,WQKRCALG,WQKLCALG,WQKDCALG,
     *         WQKDCALM
      WRITE(2,80)'* Minimum dissolution rate (/day) of organic C    '
      WRITE(2,81)' : (KRC, KLC, KDC)       = ', WQKRC,WQKLC,WQKDC
      WRITE(2,80)'* Constant relating dissolution rate to algae     '
      WRITE(2,81)' : (KRCalg,KLCalg,KDCalg)= ', WQKRCALG,WQKLCALG,
     *  WQKDCALG
C17
      Write(2,*)'C17------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQTRHDR,WQTRMNL,WQKTHDR,WQKTMNL,WQKHORDO,WQKHDNN,
     *  WQAANOX
      WRITE(2,80)'* Reference temp for hydrolysis/mineralization(oC)'
      WRITE(2,81)' : (TRHDR, TRMNL)        = ', WQTRHDR,WQTRMNL
      WRITE(2,80)'* Temperature effect on hydrolysis/mineralization '
      WRITE(2,81)' : (KTHDR, KTMNL)        = ', WQKTHDR,WQKTMNL
      WRITE(2,80)'* Half-sat. constant for oxic resp/denitrification'
      WRITE(2,81)' : (KHORDO, KHDNN)       = ', WQKHORDO,WQKHDNN
      WRITE(2,80)'* Ration of denitrification to oxic DOC resp      '
      WRITE(2,81)' : (AANOX)               = ', WQAANOX
      WQAANOX = WQAANOX*WQKHORDO
C
C Constant parameters for PHOSPHORUS (Table 3-3)
C
C18
C  
C
      Write(2,*)'C18------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFPRP,WQFPLP,WQFPDP,WQFPIP,WQFPRC,WQFPRD,WQFPRG,
     * WQFPLC,WQFPLD,WQFPLG
      WRITE(2,80)'* Phosphorus distribution coef for algal predation'
      WRITE(2,81)' : (FPRP,FPLP,FPDP,FPIP) = ', WQFPRP,WQFPLP,WQFPDP,
     *  WQFPIP
      WRITE(2,80)'* Phosphorus dist coef of RPOP for algal metabolis'
      WRITE(2,81)' : (FPRc, FPRd, FPRg)    = ', WQFPRC,WQFPRD,WQFPRG
      WRITE(2,80)'* Phosphorus dist coef of LPOP for algal metabolis'
      WRITE(2,81)' : (FPLc, FPLd, FPLg)    = ', WQFPLC,WQFPLD,WQFPLG
      XP = ABS(1.0 - (WQFPRP+WQFPLP+WQFPDP+WQFPIP))
      IF (XP .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPRP+FPLP+FPDP+FPIP not equal to 1.0'
        WRITE(2,*)
      END IF
C19
      Write(2,*)'C19------------------------------' 
	CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFPRPM,WQFPLPM,WQFPDPM,WQFPIPM,WQFPRM,WQFPLM,WQAPCM

C20
C   
C
      Write(2,*)'C20------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFPDC,WQFPDD,WQFPDG,WQFPDM,WQFPIC,WQFPID,WQFPIG,
     * WQFPIM,WQKPO4P
      IF (IWQSRP.NE.1 .AND. IWQSRP.NE.2) THEN
        WQKPO4P = 0.0
        WRITE(2,80)': no sorption of PO4t/SA, so KPO4p is forced to 0 '
      END IF
      WRITE(2,80)'* Phosphorus dist coef of DOP for algal metabolism'
      WRITE(2,81)' : (FPDc, FPDd, FPDg)    = ', WQFPDC,WQFPDD,WQFPDG
      WRITE(2,80)'* Phosphorus dist coef of NH4 for algal metabolism'
      WRITE(2,81)' : (FPIc, FPId, FPIg)    = ', WQFPIC,WQFPID,WQFPIG
      WRITE(2,82)'* Paritition coeff for sorbed/dissolved PO4 =',WQKPO4P
      XPC = ABS(1.0 - (WQFPRC+WQFPLC+WQFPDC+WQFPIC))
      IF (XPC .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPRc+FPLc+FPDc+FPIc not equal to 1.0'
        WRITE(2,*)
      END IF
      XPD = ABS(1.0 - (WQFPRD+WQFPLD+WQFPDD+WQFPID))
      IF (XPD .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPRd+FPLd+FPDd+FPId not equal to 1.0'
        WRITE(2,*)
      END IF
      XPG = ABS(1.0 - (WQFPRG+WQFPLG+WQFPDG+WQFPIG))
      IF (XPG .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPRg+FPLg+FPDg+FPIg not equal to 1.0'
        WRITE(2,*)
      END IF
C21
      Write(2,*)'C21------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO L=1,IWQZ
      READ(1,*) WQKRP,WQKLP,WQKDP,WQKRPALG,WQKLPALG,WQKDPALG,
     * WQCP1PRM1(L), WQCP2PRM,WQCP3PRM
      WRITE(2,80)'* Minimum hydrolysis rate (/day) of organic P     '
      WRITE(2,81)' : (KRP, KLP, KDP)       = ', WQKRP,WQKLP,WQKDP
      WRITE(2,80)'* Constant relating hydrolysis rate to algae      '
      WRITE(2,81)' : (KRPalg,KLPalg,KDPalg)= ', WQKRPALG,WQKLPALG,
     *  WQKDPALG
      WRITE(2,80)'* Constant used in determining P-to-C ratio       '
      WRITE(2,81)' : (CPprm1,CPprm2,CPprm3)= ', WQCP1PRM,WQCP2PRM,
     *  WQCP3PRM
      ENDDO
C
C Constant parameters for NITROGEN (Table 3-4)
C
C22
C   
C
      Write(2,*)'C22------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFNRP,WQFNLP,WQFNDP,WQFNIP,WQFNRC,WQFNRD,WQFNRG,
     *  WQFNLC,WQFNLD,WQFNLG
      WRITE(2,80)'* Nitrogen distribution coeff for algal predation '
      WRITE(2,81)' : (FNRP,FNLP,FNDP,FNIP) = ', WQFNRP,WQFNLP,WQFNDP,
     *  WQFNIP
      WRITE(2,80)'* Nitrogen dist coef of RPON for algal metabolism '
      WRITE(2,81)' : (FNRc, FNRd, FNRg)    = ', WQFNRC,WQFNRD,WQFNRG
      WRITE(2,80)'* Nitrogen dist coef of LPON for algal metabolism '
      WRITE(2,81)' : (FNLc, FNLd, FNLg)    = ', WQFNLC,WQFNLD,WQFNLG
      XN = ABS(1.0 - (WQFNRP+WQFNLP+WQFNDP+WQFNIP))
      IF (XN .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FNRP+FNLP+FNDP+FNIP not equal to 1.0'
        WRITE(2,*)
      END IF
C23
      Write(2,*)'C23------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFNRPM,WQFNLPM,WQFNDPM,WQFNIPM,WQFNRM,WQFNLM
C24
C   
      Write(2,*)'C24------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO L=1,IWQZ
      READ(1,*) WQFNDC,WQFNDD,WQFNDG,WQFNDM,WQFNIC,WQFNID,WQFNIG,
     * WQFNIM,WQANCC1(L),WQANCD1(L),WQANCG1(L),WQANCM
      ENDDO
      WRITE(2,80)'* Nitrogen dist coef of DON for algal metabolism  '
      WRITE(2,81)' : (FNDc, FNDd, FNDg)    = ', WQFNDC,WQFNDD,WQFNDG
      WRITE(2,80)'* Nitrogen dist coef of NH4 for algal metabolism  '
      WRITE(2,81)' : (FNIc, FNId, FNIg)    = ', WQFNIC,WQFNID,WQFNIG
      WRITE(2,80)'* Nitrogen-to-carbon ratio in algae               '
      WRITE(2,81)' : (ANCc, ANCd, ANCg)    = ', WQANCC,WQANCD,WQANCG
      XNC = ABS(1.0 - (WQFNRC+WQFNLC+WQFNDC+WQFNIC))
      IF (XNC .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FNRc+FNLc+FNDc+FNIc not equal to 1.0'
        WRITE(2,*)
      END IF
      XND = ABS(1.0 - (WQFNRD+WQFNLD+WQFNDD+WQFNID))
      IF (XND .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FNRd+FNLd+FNDd+FNId not equal to 1.0'
        WRITE(2,*)
      END IF
      XNG = ABS(1.0 - (WQFNRG+WQFNLG+WQFNDG+WQFNIG))
      IF (XNG .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FNRg+FNLg+FNDg+FNIg not equal to 1.0'
        WRITE(2,*)
      END IF
C25
      Write(2,*)'C25------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQANDC,WQNITM,WQKHNDO,WQKHNN,WQTNIT,WQKN1,WQKN2
      WRITE(2,82)'* Mass NO3 reduced per DOC oxidized (gN/gC)= ',WQANDC
     *          ,'* Maximum nitrification rate (g N /m^3/d)  = ',WQNITM
     *          ,'  Reference temp for nitrification (degC)  = ',WQTNIT
      WRITE(2,80)'* Nitrification half-sat constant for DO & NH4    '
      WRITE(2,81)' : (KHNitDO, KHNitN)     = ', WQKHNDO,WQKHNN
      WRITE(2,80)'* Sub & super-optimum temp effect on nitrification'
      WRITE(2,81)' : (KNit1, KNit2)        = ', WQKN1,WQKN2
C26
      Write(2,*)'C26------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQKRN,WQKLN,WQKDN,WQKRNALG,WQKLNALG,WQKDNALG
      WRITE(2,80)'* Minimum hydrolysis rate (/day) of organic N     '
      WRITE(2,81)' : (KRN, KLN, KDN)       = ', WQKRN,WQKLN,WQKDN
      WRITE(2,80)'* Constant relating hydrolysis rate to algae      '
      WRITE(2,81)' : (KRNalg,KLNalg,KDNalg)= ', WQKRNALG,WQKLNALG,
     *  WQKDNALG
C27
C Constant parameters for SILICA (Table 3-5)
C
      Write(2,*)'C27------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFSPP,WQFSIP,WQFSPD,WQFSID,WQASCD,WQKSAP,WQKSU,
     *  WQTRSUA,WQKTSUA
      IF (IWQSRP.NE.1 .AND. IWQSRP.NE.2) THEN
        WQKSAP = 0.0
        WRITE(2,80)': no sorption of PO4t/SA, so KSAp is forced to 0  '
      END IF
      WRITE(2,80)'* Silica distribution coeff for diatom predation  '
      WRITE(2,81)' : (FSPP, FSIP)          = ', WQFSPP,WQFSIP
      WRITE(2,80)'* Silica distribution coeff for diatom metabolism '
      WRITE(2,81)' : (FSPd, FSId)          = ', WQFSPD,WQFSID
      WRITE(2,82)'* Silica-to-carbon ratio in diatoms        = ',WQASCD
     *          ,'* Paritition coeff for sorbed/dissolved SA = ',WQKSAP
     *          ,'* Dissolution rate (/d) of PSi             = ',WQKSU
     *          ,'  Reference temp for PSi dissolution (oC)  = ',WQTRSUA
     *          ,'  Temperature effect on PSi dissolution    = ',WQKTSUA
cxh   IF (WQFSPP+WQFSIP.NE.1.0) STOP 'ERROR!! invalid FSPP,FSIP values'
cxh   IF (WQFSPD+WQFSID.NE.1.0) STOP 'ERROR!! invalid FSPd,FSId values'
C
C28
C Constant parameters for COD & DO (Table 3-6)
C
      Write(2,*)'C28------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQAOCR,WQAONT,WQKTR,WQKHCOD,WQKCD,WQTRCOD,WQKTCOD,
     *   WQAOCRpm, WQAOCRrm
      WRITE(2,82)'* DO-to-carbon ratio in respiration        = ',WQAOCR
     *          ,': Mass DO consumed per NH4-N nitrified     = ',WQAONT
     *          ,': Proporn. constant for DO-reaeration (MKS)= ',WQKRO
     *          ,'  Temperature effect on DO-reaeration      = ',WQKTR
     *          ,'* Half-sat constant of DO for COD (gO2/m^3)= ',WQKHCOD
     *          ,': Oxidation rate of COD (/day)             = ',WQKCD
     *          ,'  Reference temp for COD oxidation (oC)    = ',WQTRCOD
     *          ,'  Temperature effect on COD oxidation      = ',WQKTCOD
     *        ,': DO-to-carbon ratio macroalgae photosynth = ',WQAOCRpm
     *        ,': DO-to-carbon ratio macroalgae respiration= ',WQAOCRrm
C
C29 Constant parameters for TAM & FCB (Table 3-7)
C
      Write(2,*)'C29------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQKHBMF,WQBFTAM,WQTTAM,WQKTAM,WQTAMDMX,WQKDOTAM,
     *  WQKFCB,WQTFCB
      WRITE(2,82)
     *  '* DO where TAM release is half anoxic one  = ',WQKHBMF
     * ,'  Anoxic release of TAM (mol/m^2/d)        = ',WQBFTAM
     * ,'  Reference temp for TAM release (oC)      = ',WQTTAM
     * ,'  Temperature effect on TAM release        = ',WQKTAM
     * ,': TAM solubility at anoxic cond. (mol/m^3) = ',WQTAMDMX
     * ,'  Constant relating TAM solubility to DO   = ',WQKDOTAM
     * ,'* First-order die-off rate at 20oC (/d)    = ',WQKFCB
     * ,'  Temperature effect on bacteria die-off   = ',WQTFCB
C
C Set up look-up table for temperature dependency over -15oC to 40oC
C

      DO M=1,NWQTD
        WTEMP =1.00*REAL(M-1)*0.1 - 4.95   ! -5C -> 50C, Ji, 8/29/02
c        WTEMP =1.00*REAL(M-1)*0.1 - 14.95
c       WTEMP =1.23*REAL(M-1)*0.1 - 14.95
c       WTEMP =1.4*REAL(M-1)*0.1 - 4.95
C       WTEMP =REAL(M-1)*0.1 - 4.95
C
c        IF (WTEMP.GT.WQTMC) THEN
c          WQTDGC(M) = EXP( WQKG1C*(WTEMP-WQTMC)*(WQTMC-WTEMP) )
c         ELSE
c          WQTDGC(M) = EXP( WQKG2C*(WTEMP-WQTMC)*(WQTMC-WTEMP) )
c        END IF
C
c        IF (WTEMP.GT.WQTMD) THEN
c          WQTDGD(M) = EXP( WQKG1D*(WTEMP-WQTMD)*(WQTMD-WTEMP) )
c         ELSE
c          WQTDGD(M) = EXP( WQKG2D*(WTEMP-WQTMD)*(WQTMD-WTEMP) )
c        END IF
C
c        IF (WTEMP.GT.WQTMG) THEN
c          WQTDGG(M) = EXP( WQKG1G*(WTEMP-WQTMG)*(WQTMG-WTEMP) )
c         ELSE
c          WQTDGG(M) = EXP( WQKG2G*(WTEMP-WQTMG)*(WQTMG-WTEMP) )
c        END IF
c Modified on 03/04/96 by JMH/MRM:
       DO I=1,IWQZ
        WQTDGC(M,I)=1.
        IF (WTEMP.LT.WQTMC1(I)) THEN
         WQTDGC(M,I)=EXP(-WQKG1C(I)*(WTEMP-WQTMC1(I))*(WTEMP-WQTMC1(I)))
        END IF
        IF (WTEMP.GT.WQTMC2(I)) THEN
         WQTDGC(M,I)=EXP(-WQKG2C(I)*(WTEMP-WQTMC2(I))*(WTEMP-WQTMC2(I)))
        END IF
C
        WQTDGD(M,I)=1.
        IF (WTEMP.LT.WQTMD1(I)) THEN
         WQTDGD(M,I)=EXP(-WQKG1D(I)*(WTEMP-WQTMD1(I))*(WTEMP-WQTMD1(I)))
        END IF
        IF (WTEMP.GT.WQTMD2(I)) THEN
         WQTDGD(M,I)=EXP(-WQKG2D(I)*(WTEMP-WQTMD2(I))*(WTEMP-WQTMD2(I)))
        END IF
C
        WQTDGG(M,I)=1.
        IF (WTEMP.LT.WQTMG1(I)) THEN
         WQTDGG(M,I)=EXP(-WQKG1G(I)*(WTEMP-WQTMG1(I))*(WTEMP-WQTMG1(I)))
        END IF
        IF (WTEMP.GT.WQTMG2(I)) THEN
         WQTDGG(M,I)=EXP(-WQKG2G(I)*(WTEMP-WQTMG2(I))*(WTEMP-WQTMG2(I)))
        END IF
        WQTDGG(M,I) = max(WQTDGG(M,I),T_Lim )
        WQTDGD(M,I) = max(WQTDGD(M,I),T_Lim)
        WQTDGC(M,I) = max(WQTDGC(M,I),T_Lim)
C J.S.  5/5/98
        WQTDGM(M,I)=1.
        IF(IDNOTRVA.GT.0)THEN
         IF (WTEMP.LT.WQTMM1(I)) THEN
         WQTDGM(M,I)=EXP(-WQKG1M(I)*(WTEMP-WQTMM1(I))*(WTEMP-WQTMM1(I)))
         END IF
         IF (WTEMP.GT.WQTMM2(I)) THEN
         WQTDGM(M,I)=EXP(-WQKG2M(I)*(WTEMP-WQTMM2(I))*(WTEMP-WQTMM2(I)))
         END IF
           WQTDRM(M)= EXP( WQKTBM*(WTEMP-WQTRM) )
        END IF

        

C J.S.
C MRM: 06/20/98
C  The following WQTDGP variable is a temperature related adjustment
C  to the predation and/or basal matabolism rate to allow diatoms
C  to bloom in winter (or other time of year).
        WQTDGp(M,I)=1.
        IF (WTEMP.LT.WQTMp1(I)) THEN
        WQTDGp(M,I)= EXP(-WQKG1p(I)*(WTEMP-WQTMp1(I))*(WTEMP-WQTMp1(I)))
        END IF
        IF (WTEMP.GT.WQTMD2(I)) THEN
        WQTDGp(M,I)= EXP(-WQKG2p(I)*(WTEMP-WQTMp2(I))*(WTEMP-WQTMp2(I)))
        END IF
        ENDDO
        ENDDO
c MRM
C
C        write(1,555)WTEMP,WQTDGC(M),WQTDGD(M),WQTDGG(M),WQTDGM(M) 
        DO M=1,NWQTD
        WTEMP =1.00*REAL(M-1)*0.1 - 4.95 
555     format(f7.2,4e12.4)
        WQTDRC(M) = EXP( WQKTBC*(WTEMP-WQTRC) )
        WQTDRD(M) = EXP( WQKTBD*(WTEMP-WQTRD) )
        WQTDRG(M) = EXP( WQKTBG*(WTEMP-WQTRG) )
        WQTDHDR(M) = EXP( WQKTHDR*(WTEMP-WQTRHDR) )
        WQTDMNL(M) = EXP( WQKTMNL*(WTEMP-WQTRMNL) )
C
        IF (WTEMP.GT.WQTNIT) THEN
          WQTDNIT(M) = WQNITM*EXP( WQKN1*(WTEMP-WQTNIT)*(WQTNIT-WTEMP) )
         ELSE
          WQTDNIT(M) = WQNITM*EXP( WQKN2*(WTEMP-WQTNIT)*(WQTNIT-WTEMP) )
        END IF
C
        WQKSUA(M) = WQKSU * EXP( WQKTSUA*(WTEMP-WQTRSUA) )
        WQKCOD(M) = WQKCD * EXP( WQKTCOD*(WTEMP-WQTRCOD) )
        TT20 = WTEMP-20.0
        WQTDKR(M) = WQKTR**TT20
        WQTDTAM(M) = WQKHBMF * WQBFTAM * EXP( WQKTAM*(WTEMP-WQTTAM) )
        WQTT = WQKFCB * WQTFCB**TT20 * DTWQO2                            ! (3-19)
        WQTD1FCB(M) = 1.0 - WQTT
        WQTD2FCB(M) = 1.0 / (1.0 + WQTT)
        ENDDO

!
! spatially/temporally constant ALGAL GROWTH, RESPIRATION & PREDATION RATES
!     
!30 45   
!
      Write(2,*)'C30------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ      
      READ(1,*) WQPMC(I),WQPMD(I),WQPMG(I),WQPMM(I),WQBMRC(I),
     *    WQBMRD(I),WQBMRG(I),WQBMRM(I),WQPRRC(I),WQPRRD(I),
     *    WQPRRG(I),WQPRRM(I) !,WQKEB(1)
       IF(ITRYCHLA.GE.2) THEN
       WQPMC(I)=WQPMC(I)*(1/WQCHLC(I))*1000    ! be careful! inversed!, 9/15/02
       WQPMD(I)=WQPMD(I)*(1/WQCHLD(I))*1000
       WQPMG(I)=WQPMG(I)*(1/WQCHLG(I))*1000
       ENDIF
       IF (IWQAGR.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,80)'* Algal growth rate (/day)                        '
        WRITE(2,21)' : (PMc, PMd, PMg)       = ', WQPMC(1),WQPMD(1),
     *    WQPMG(1)
        WRITE(2,80)'* Algal basal metabolism rate (/day)              '
        WRITE(2,21)' : (BMRc, BMRd, BMRg)    = ', WQBMRC(1),WQBMRD(1),
     *    WQBMRG(1)
        WRITE(2,80)'* Algal predation rate (/day)                     '
        WRITE(2,21)' : (PRRc, PRRd, PRRg)    = ', WQPRRC(1),WQPRRD(1),
     *    WQPRRG(1)
        WRITE(2,82)
     *    '* Base light extinction coefficient (/m)   = ',WQKEB(1)
C        write(2,82)
C     *    '* WQSDCOEF (Secchi depth = WQSDCOEF/WQKESS)= ',wqsdcoef(1)
        END IF
!          WQPMC(I)=WQPMC(1)
!          WQPMD(I)=WQPMD(1)
!          WQPMG(I)=WQPMG(1)
!          WQPMM(I)=WQPMM(1)
!          WQBMRC(I)=WQBMRC(1)
!          WQBMRD(I)=WQBMRD(1)
!          WQBMRG(I)=WQBMRG(1)
!          WQBMRM(I)=WQBMRM(1)
!          WQPRRC(I)=WQPRRC(1)
!          WQPRRD(I)=WQPRRD(1)
!          WQPRRG(I)=WQPRRG(1)
!          WQPRRM(I)=WQPRRM(1)
!          WQKEB(I)=WQKEB(1)
C          wqsdcoef(i)=wqsdcoef(1)
        END DO

!
!C31  (46)
! spatially/temporally constant SETTLING VELOCITIES and REAERATION FACTOR
!
      Write(2,*)'C31------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ
      READ(1,*)WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),WQWSLP(I),WQWSS(I),
     *      WQWSM,REAC(I)
      IF (IWQSTL.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,80)'* Algal settling rate (m/day)                     '
        WRITE(2,21)' : (WSc, WSd, WSg)       = ', WQWSC(1),WQWSD(1),
     *    WQWSG(1)
        WRITE(2,80)'* POM settling rate (m/day)                       '
        WRITE(2,21)' : (WSrp, WSlp)          = ', WQWSRP(1),WQWSLP(1)
        WRITE(2,80)'* Settling rate of particulate metal (m/day)      '
        WRITE(2,21)' : (WSs)                 = ', WQWSS(1)
      END IF
!        DO I=2,IWQZ
!          WQWSC(I)=WQWSC(1)
!          WQWSD(I)=WQWSD(1)
!          WQWSG(I)=WQWSG(1)
!          WQWSRP(I)=WQWSRP(1)
!          WQWSLP(I)=WQWSLP(1)
!          WQWSS(I)=WQWSS(1)
 !         REAC(I)=REAC(1)
        END DO


!
! spatially/temporally constant INITIAL CONDITIONS: WQCHLx=1/WQCHLx
! read data points & do internal interpolation?
!
!C32 (C44)
!
      Write(2,*)'C32------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) (WQV(1,1,NW), NW=1,6)
      READ(1,*) (WQV(1,1,NW), NW=7,13)
      READ(1,*) (WQV(1,1,NW), NW=14,NWQV),WQV(1,1,22),WQMCMIN
      Write(*,*)'Minimum algae = ', WQMCMIN
      IF (IWQICI.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,21)' : (Bc, Bd, Bg)          = ', (WQV(1,1,NW),NW=1,3)
        WRITE(2,21)' : (RPOC, LPOC, DOC)     = ', (WQV(1,1,NW),NW=4,6)
        WRITE(2,21)' : (RPOP,LPOP,DOP,PO4t)  = ', (WQV(1,1,NW),NW=7,10)
        WRITE(2,21)' : (RPON, LPON, DON)     = ', (WQV(1,1,NW),NW=11,13)
        WRITE(2,21)' : (NH4, NO3)            = ', (WQV(1,1,NW),NW=14,15)
        WRITE(2,21)' : (SU, SA, COD, DO)     = ', (WQV(1,1,NW),NW=16,19)
        WRITE(2,981)' : (TAM, FCB,MALG)      = ',
     *    (WQV(1,1,NW),NW=20,NWQV+1)
        WQCHL(1,1) = WQV(1,1,1)*WQCHLC(1) + WQV(1,1,2)*WQCHLD(1)
     *    + WQV(1,1,3)*WQCHLG(1)
C  IWQSRP  = switch for sediment sorption (1=TAM sorption; 2=sediment sorption)
C  wqkpo4p = KPO4p = partition coefficient for sorbed/dissolved PO4
        IF (IWQSRP.EQ.1) THEN
          O2WQ = MAX(WQV(1,1,19), 0.0)
          WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ), WQV(1,1,20) )
          WQTAMP(1,1) = WQV(1,1,20) - WQTAMD
          WQPO4D(1,1) = WQV(1,1,10) / (1.0 + WQKPO4P*WQTAMP(1,1))
          WQSAD(1,1)  = WQV(1,1,17) / (1.0 + WQKSAP*WQTAMP(1,1))
         ELSE IF (IWQSRP.EQ.2) THEN
          WQPO4D(1,1) = WQV(1,1,10) / (1.0 + WQKPO4P*SEDTWQ(1,1))   !=> Dissolved PO4t, Eq.(3-8c), Ji, 9/18/02
          WQSAD(1,1)  = WQV(1,1,17) / (1.0 + WQKSAP*SEDTWQ(1,1))    !=> Dissolved SA is also controlled by SEDT!
         ELSE
          WQPO4D(1,1) = WQV(1,1,10)
          WQSAD(1,1)  = WQV(1,1,17)
        END IF
        DO NW=1,NWQV
         TVARWQ=WQV(1,1,NW)
         DO K=1,KC
          WQV(LA,K,NW) = TVARWQ
          WQV(1 ,K,NW) = TVARWQ
          WQVO(LA,K,NW) = TVARWQ
          WQVO(1 ,K,NW) = TVARWQ
         END DO
        END DO
        DO NW=1,NWQV
         TVARWQ=WQV(1,1,NW)
           DO K=1,KC
           DO L=2,LA
            WQV(L,K,NW) = TVARWQ
           END DO
          END DO
        END DO
        XWQCHL = WQCHL(1,1)
        XWQTAMP = WQTAMP(1,1)
        XWQPO4D = WQPO4D(1,1)
        XWQSAD = WQSAD(1,1)
        DO K=1,KC
         WQCHL(LA,K) = XWQCHL
         WQTAMP(LA,K) = XWQTAMP
         WQPO4D(LA,K) = XWQPO4D
         WQSAD(LA,K) = XWQSAD
         WQCHL(1,K) = XWQCHL
         WQTAMP(1,K) = XWQTAMP
         WQPO4D(1,K) = XWQPO4D
         WQSAD(1,K) = XWQSAD
        END DO

         DO K=1,KC
          DO L=2,LA
           WQCHL(L,K) = XWQCHL
           WQTAMP(L,K) = XWQTAMP
           WQPO4D(L,K) = XWQPO4D
           WQSAD(L,K) = XWQSAD
          END DO
        END DO
      END IF

C   Modified by J.S.(5/5/98)
C
C   Read cell mapping file 'macalgmp.inp' and set initial condition.
C   All macalgae will only stay at the bottom layer.
C
c     IF(IDNOTRVA.GT.0) THEN
      IF(IDNOTRVA.GT.21) THEN  ! no macalgmp.inp when IDNOTRVA=1, Ji, 7/23/99
      DO L=1,LCMWQ 
        SMAC(L) =0
      END DO
      OPEN(3,FILE='macalgmp.inp',STATUS='UNKNOWN')
      DO NS=1,6
        READ(3,999)
      END DO
      READ(3,*)NMACELL
      READ(3,999)
      DO NS=1,NMACELL
        READ(3,*) II,JJ,TT_p
        IF (LIJW(II,JJ).LT.1) THEN
          PRINT*, 'i, j, LIJW(II,JJ) = ', II,JJ,LIJW(II,JJ)
           STOP 'ERROR!! invalid (i,j) in FILE wqcalgmp.inp'
        END IF
        LL=LIJW(II,JJ)
        SMAC(LL)=1.0
        WQV(LL,1,IDNOTRVA)=TT_p
        WQVO(LL,1,IDNOTRVA)=TT_p
      END DO
      WQV(1,1,IDNOTRVA)=0
      CLOSE(3)
c      write(6,6666)'read macalgmp.inp'
      END IF
 6666 format(a30)
!
!C33 (47)
! spatially/temporally constant BENTHIC FLUXES
!
      Write(2,*)'C33------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQBFPO4D(1),WQBFNH4(1),WQBFNO3(1),WQBFSAD(1),
     *  WQBFCOD(1),WQBFO2(1)
      IF (IWQBEN.EQ.0) THEN
        WRITE(2,21)' : (PO4d, NH4, NO3)      = ',WQBFPO4D(1),WQBFNH4(1),
     *    WQBFNO3(1)
        WRITE(2,21)' : (SAd, COD, DO)        = ',WQBFSAD(1),WQBFCOD(1),
     *    WQBFO2(1)
        DO L=2,LA
          WQBFPO4D(L)=WQBFPO4D(1)
          WQBFNH4(L)=WQBFNH4(1)
          WQBFNO3(L)=WQBFNO3(1)
          WQBFSAD(L)=WQBFSAD(1)
          WQBFCOD(L)=WQBFCOD(1)
          WQBFO2(L)=WQBFO2(1)
        END DO
      END IF
!
!C34 Read in total numer of time series open boundary condition
!
      WRITE(2,*)'C34------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,90) TITLE(1)
      WRITE(2,90) TITLE(1)
        DO M=1,5
          READ(1,90) TITLE(M)
          WRITE(2,90) TITLE(M)
        END DO
      READ(1,*) (NWQCSR(NW),NW=1,NWQV)
!
! check array dimensions (NWQCSRM: parameter set in WQ.PAR)
!
      do NW=1,NWQV
      if(NWQCSR(NW).gt.NWQCSRM) then   ! 1
      write(6,*) " Error:  NWQCSR(NW)>NWQCSRM", NW,NWQCSR(NW),NWQCSRM
      stop
      endif
      ENDDO
      WRITE(2,970)(NWQCSR(NW),NW=1,NWQV)
!
!C35
! 
! Read total number of Open boundary condition
!
      WRITE(2,*)'C35------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,90) TITLE(1)
      WRITE(2,90) TITLE(1)
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) NWQOBS,NWQOBW,NWQOBE,NWQOBN
      WRITE(2,23)'* # of open bdry cells on SOUTH bdry        = ',NWQOBS
      WRITE(2,23)'* # of open bdry cells on WEST bdry         = ',NWQOBW
      WRITE(2,23)'* # of open bdry cells on EAST bdry         = ',NWQOBE
      WRITE(2,23)'* # of open bdry cells on NORTH bdry        = ',NWQOBN
      IF (NWQOBS.GT.NBBSWM) STOP 'ERROR!! NWQOBS should <= NBBSWM'
      IF (NWQOBW.GT.NBBWWM) STOP 'ERROR!! NWQOBW should <= NBBWWM'
      IF (NWQOBE.GT.NBBEWM) STOP 'ERROR!! NWQOBE should <= NBBEWM'
      IF (NWQOBN.GT.NBBNWM) STOP 'ERROR!! NWQOBN should <= NBBNWM'
      WRITE(2,999)
!      WRITE(2,80)'* constant OBC at (ICBx(M),JCBx(M)) if IWQOBx(M)=0'
!      WRITE(2,80)': read time-series OBCs IWQOBx times if IWQOBx > 0'
!
!36  Read parameters for SOUTH OPEN BOUNDARY
!
      WRITE(2,*)'C36------------------------------'
	WRITE(2,*)'South open boundary condition'
      CALL SKIPCOMM(1,CCMRM)
      IF (NWQOBS.GT.0) THEN
        DO M=1,NWQOBS
          READ(1,*) IWQCBS(M),JWQCBS(M),(IWQOBS(M,NW),NW=1,NWQV)
          WRITE(2,969) IWQCBS(M),JWQCBS(M),(IWQOBS(M,NW),NW=1,NWQV)
        END DO
      END IF
!
!C37 WEST BOUNDARY
!
      WRITE(2,*)'C37------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      IF (NWQOBW.GT.0) THEN
        DO M=1,NWQOBW
          READ(1,*) IWQCBW(M),JWQCBW(M),(IWQOBW(M,NW),NW=1,NWQV)
          WRITE(2,969) IWQCBW(M),JWQCBW(M),(IWQOBW(M,NW),NW=1,NWQV)
	   do NW=1,NWQV
          NWQCSR(NW)=max(NWQCSR(NW),IWQOBw(M,NW))
!          WQOBCw(M,1,NW)=0.0
!          WQOBCw(M,2,NW)=0.0
         enddo
        END DO
      END IF

!
!C38 EAST BDRY
!
      WRITE(2,*)'C38------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      IF (NWQOBE.GT.0) THEN
        DO M=1,NWQOBE
          READ(1,*) IWQCBE(M),JWQCBE(M),(IWQOBE(M,NW),NW=1,NWQV)
          WRITE(2,969) IWQCBE(M),JWQCBE(M),(IWQOBE(M,NW),NW=1,NWQV)
        END DO
       END IF
!
!C39 NORTH BDRY
!
      WRITE(2,*)'C39------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      IF (NWQOBN.GT.0) THEN
        DO M=1,NWQOBN
          READ(1,*) IWQCBN(M),JWQCBN(M),(IWQOBN(M,NW),NW=1,NWQV)
          WRITE(2,969) IWQCBN(M),JWQCBN(M),(IWQOBN(M,NW),NW=1,NWQV)
        END DO
      END IF


!
! temporally-constant values for POINT SOURCE INPUT in (kg/d) except
! xPSQ (m^3/s), xPO2 (g/m^3), xPTAM (kmol/d), xPFCB (MPN/100mL).
!: In GEs, load is in (g/d) except TAM in (mol/d), FCB in (MPN/100mL)*(m^3/s).
!: To convert kg/d to g/d, (xW kg/d)*(10^3 g/kg) = (CONV1*xW g/d) with
!   CONV1=1.0E3.
!: For O2, (xPSQ m^3/s)*(xPO2 g/m^3)*(86400 s/d) = (CONV2*xPSQ*xPO2 g/d)
!   with CONV2=8.64E4.
!: For TAM, (xPTAM kmol/d)*(10^3 mol/kmol) = (CONV1*xPTAM mol/d).
!: For FCB, (xPFCB MPN/100mL)*(xPSQ m^3/s)*(86400 s/d) =
!   (CONV2*xPSQ*xPFCB (MPN/100mL)*m^3/d).
!
!
!40 Point source input ----- 
!      
      WRITE(2,*)'Card 40---------------------------------'
      CALL SKIPCOMM(1,'C')
      READ(1,*) IWQPS,NPSTMSR
!
! check array dimensions
!
      if(NPSTMSR.gt.NWQPSRM) then   ! 1
      write(6,*) " Error:  NPSTMSR >NWQPSRM ", NPSTMSR, NWQPSRM
      stop
      endif

      WRITE(2,23)'* Number of cells for point source input  = ',IWQPS
      WRITE(2,23)'* Number with variable point source input = ',NPSTMSR
!
! C41 location of point source
!
      WRITE(2,*)'Card 41---------------------------------'
      CALL SKIPCOMM(1,'C')
      DO M=1,IWQPS
        READ(1,*) I,J,K
	  write(2,*) 'IWPS  JWPS  WPSER '
        write(2,294) I,J,K
!c check array dimensions
          IF (LIJW(I,J).LT.1) then
            write(*,911) i,j
  911       format(/,' I,J = ', 2i5)
            STOP 'ERROR!! invalid (i,j) in FILE wq3dwc.inp for PSL'
          end if
!c
!cmrm M.R. Morton  02/20/1999
!cmrm Added a fix here so that multiple point source loads (PSL) can be added
!cmrm to the same I,J cell.  In original code, the PSL read in last
!cmrm overwrites any earlier PSLs in the same I,J cell.
!cmrm    M = point source number
!cmrm    K = Time series ID
!cmrm    L = cell number for location I,J
!cmrm    IWQPSC(L,K) = M, index pointer to constant PSL number
!cmrm    IWQPSV(L,K) = index pointer to time-variable PSL in WQPSL.INP
!cmrm    ICPSL(M) = saves I cell index of PSL M
!cmrm    JCPSL(M) = saves J cell index of PSL M
!cmrm    KCPSL(M) = saves K layer index of PSL M
!cmrm    MVPSL(M) = index pointer to time-variable PSL in WQPSL.INP
          L = LIJW(I,J)
          If(L.EQ.1) write(*,*)'Check PT ', M,I,J,K
!          IWQPSC(L,K)=M
!          IWQPSV(L,K)=ITMP
          ICPSL(M)=I
          JCPSL(M)=J
!         KCPSL(M)=K
          MVPSL(M)=K
 !         WQPSQC(M)=XPSQ
 !         DO NW=1,18
 !           WQWPSLC(M,NW) = XPSL(NW) * CONV1
 !         END DO
 !         WQTT = XPSQ*CONV2
 !         WQWPSLC(M,19) = XPSL(19) * WQTT
 !         WQWPSLC(M,20) = XPSL(20) * CONV1
 !         WQWPSLC(M,NWQV) = XPSL(NWQV) * WQTT
c        END IF
      END DO
C
c49 --
c
C spatially/temporally-constant values for NON-POINT SOURCE INPUT in (kg/d)
C except xDSQ (m^3/s), xNO2 (mg/L), xNTAM (kmol/d), xNFCB (MPN/100mL).
c mrm  Note: This group is now used for dry atmospheric deposition with
c      constituent units of g/m2/day except FCB which is MPN/m2/day.
C
!
!C42:Constant Dry Atmospheric Deposition
! 
      WRITE(2,*)'Card 42---------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) XDSQ,(XDSL(NW),NW=1,6)
      READ(1,*) (XDSL(NW),NW=7,13)
      READ(1,*) (XDSL(NW),NW=14,NWQV)
      IF (IWQNPL.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,21)' : (DSQ, Bc, Bd, Bg)     = ',XDSQ,(XDSL(NW),NW=1,3)
        WRITE(2,21)' : (RPOC, LPOC, DOC)     = ',(XDSL(NW),NW=4,6)
        WRITE(2,21)' : (RPOP,LPOP,DOP,PO4t)  = ',(XDSL(NW),NW=7,10)
        WRITE(2,21)' : (RPON, LPON, DON)     = ',(XDSL(NW),NW=11,13)
        WRITE(2,21)' : (NH4, NO3)            = ',(XDSL(NW),NW=14,15)
        WRITE(2,21)' : (SU, SA, COD, DO)     = ',(XDSL(NW),NW=16,19)
        WRITE(2,981)' : (TAM, FCB)            = ',(XDSL(NW),NW=20,NWQV)
        WQDSQ(1,1) = XDSQ
        DO NW=1,18
          WQWDSL(1,1,NW) = XDSL(NW) * CONV1
        END DO
        WQTT = XDSQ*CONV2
        WQWDSL(1,1,19) = XDSL(19) * WQTT
        WQWDSL(1,1,20) = XDSL(20) * CONV1
        WQWDSL(1,1,NWQV) = XDSL(NWQV) * WQTT
C
        DO L=2,LA
          WQDSQ(L,1)=WQDSQ(1,1)
          DO NW=1,NWQV
c M. Morton modified the line below so that constant atmospheric deposition
c can be added via this routine instead of constant NPS input which the
c original code called for and which was not particularly useful.
c Input data (xdsl) are in g/m2/day and are multiplied by the cell surface
c area (dxyp) to get g/day.  Atmospheric deposition only enters thru surface
c layer (kc):
c            WQWDSL(L,1,NW)=WQWDSL(1,1,NW)
            WQWDSL(L,kc,nw) = XDSL(NW) * DXYPWQ(L) !DXYPW(L)
          END DO
        END DO
      END IF
!
!43
!
! Wet Atmospheric Deposition concentrations (mg/L; MPN/L) 
! and will be multiplied by the rainfall flow rate into each grid cell
! to get a load in kg/day.
!
      WRITE(2,*)'Card 43---------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) (wqatm(NW),NW=1,6)
      READ(1,*) (wqatm(NW),NW=7,13)
      READ(1,*) (wqatm(NW),NW=14,NWQV)
      WRITE(2,999)
      WRITE(2,90) TITLE(1)
      WRITE(2, 21)' : (Bc, Bd, Bg)          = ',(wqatm(NW),NW=1,3)
      WRITE(2, 21)' : (RPOC, LPOC, DOC)     = ',(wqatm(NW),NW=4,6)
      WRITE(2, 21)' : (RPOP,LPOP,DOP,PO4t)  = ',(wqatm(NW),NW=7,10)
      WRITE(2, 21)' : (RPON, LPON, DON)     = ',(wqatm(NW),NW=11,13)
      WRITE(2, 21)' : (NH4, NO3)            = ',(wqatm(NW),NW=14,15)
      WRITE(2, 21)' : (SU, SA, COD, DO)     = ',(wqatm(NW),NW=16,19)
      WRITE(2,981)' : (TAM, FCB)            = ',(wqatm(NW),NW=20,NWQV)

C
c44
c
C Input/output file names for spatially and/or temporally varying parameters
C
C     READ(1,999)
C     WRITE(2,999)
      CALL SKIPCOMM(1,CCMRM)
      READ(1,90) TITLE(1)
      WRITE(2,999)
      WRITE(2,90) TITLE(1)
C
      READ(1,295) RSTOFN
      WRITE(2,85)'* Output file for restart writing         = ', RSTOFN
      IF (IWQRST.GE.1) THEN
        OPEN(99,FILE=RSTOFN,STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE=RSTOFN,STATUS='UNKNOWN')
CXH        OPEN(UNIT=IWQORST,FILE=RSTOFN,STATUS='UNKNOWN')
        WRITE(99,888)
        CLOSE(99)
       ELSE
        IF (RSTOFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQORST/RSTOFN'
      END IF
  888 FORMAT('    L    K',
     *'          Bc          Bd          Bg        RPOC        LPOC',
     *'         DOC        RPOP        LPOP         DOP        PO4t',
     *'        RPON        LPON         DON         NH4         NO3',
     *'          SU          SA         COD          O2         TAM',
     *'         FCB')
C
      READ(1,295) ICIFN
      WRITE(2,85)'* File for initial conditions             = ', ICIFN
c      IF (IWQICI.EQ.1) THEN
cxh        OPEN(UNIT=INWQICI,FILE=ICIFN,STATUS='OLD')
c       ELSE IF (IWQICI.EQ.2) THEN
cxh        OPEN(UNIT=INWQRST,FILE=ICIFN,STATUS='OLD')
c       ELSE
!        IF (ICIFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQICI/ICIFN'
c      END IF
C
C      READ(1,295) OBCFN
C      WRITE(2,85)'* File for downriver boundary conditions  = ', OBCFN
C      IF (IWQOBC.EQ.1) THEN
C        OPEN(UNIT=INWQOBC,FILE=OBCFN,STATUS='OLD')
C       ELSE
C        IF (OBCFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQOBC/OBCFN'
C      END IF
C
      READ(1,295) AGRFN
      WRITE(2,85)'* File for algal growth, resp., predatat. = ', AGRFN
      IF (IWQAGR.EQ.1) THEN
cxh        OPEN(UNIT=INWQAGR,FILE=AGRFN,STATUS='OLD')
       ELSE
        IF (AGRFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQAGR/AGRFN'
      END IF
C
      READ(1,295) STLFN
      WRITE(2,85)'* File for settling rates of algae, part. = ', STLFN
      IF (IWQSTL.EQ.1) THEN
cxh        OPEN(UNIT=INWQSTL,FILE=STLFN,STATUS='OLD')
       ELSE
        IF (STLFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQSTL/STLFN'
      END IF
C
      READ(1,295) SUNFN
      WRITE(2,85)'* File for Io, FD, Te, KT                 = ', SUNFN
      IF (IWQSUN.EQ.1) THEN
cxh        OPEN(UNIT=INWQSUN,FILE=SUNFN,STATUS='OLD')
       ELSE
        IF (SUNFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQSUN/SUNFN'
      END IF
C
      READ(1,295) BENFN
      WRITE(2,85)'* File for benthic flux                   = ', BENFN
      IF (IWQBEN.EQ.2) THEN
CXH        OPEN(UNIT=INWQBEN,FILE=BENFN,STATUS='OLD')
       ELSE
        IF (BENFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQBEN/BENFN'
      END IF
C
      READ(1,295) PSLFN
      WRITE(2,85)'* File for point source input             = ', PSLFN
c      IF (IWQPSL.EQ.1) THEN
CXH        OPEN(UNIT=INWQPSL,FILE=PSLFN,STATUS='OLD')
c       ELSE
c        IF (PSLFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQPSL/PSLFN'
c      END IF
C
      READ(1,295) NPLFN
      WRITE(2,85)'* File for NPS input including atm. input = ', NPLFN
      IF (IWQNPL.EQ.1) THEN
CXH        OPEN(UNIT=INWQNPL,FILE=NPLFN,STATUS='OLD')
       ELSE
        IF (NPLFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQNPL/NPLFN'
      END IF
C
      READ(1,295) NCOFN
      WRITE(2,85)'* Diagnostic file for negative concentrat = ', NCOFN
C
C45 Diagnose control
C
      SODMULT(1)=1
      REAC(1)=1
      CALL SKIPCOMM(1,CCMRM)
      Read(1,*)IDIA,SODMULT(1),BFHN,BFNO, BFP,
     &  DO_G,DO_R,AN_Lim,P_lim,T_Lim,REAC(1),F_AB,
     &  RN_R,RP_R,RC_R
     
	IF(IDIA.EQ.0) then
	SODMULT(1)=1.0
	BFHN =1.0
	BFNO =1.0
	BFP  =1.0
	DO_G =1.0
	DO_R =1.0
	AN_Lim =0.0
	P_lim  =0.0
	T_Lim  =0.0
!	REAC(1)=1.0
	F_AB =1.0
	RN_R =0.0
	RP_R =0.0
	RC_R =0.0
	ENDIF
	
      DO I=1,IWQZ
!      REAC(I)=REAC(1)
      SODMULT(I)=SODMULT(1)
	ENDDO 

c46 Forced Temperature
C
      CALL SKIPCOMM(1,CCMRM)	     
      read(1,*)IFOCETMP,TPHPSE,   TPMAG, TMPADD  
	    	     
      
      CLOSE(1)
C      write(6,6666)'read wq3dwc.inp'
C
      IF (IWQNC.EQ.1) THEN
        OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')
        WRITE(1,284)'* Negative concentration occurs:'
        CLOSE(1)
       ELSE
!        IF (NCOFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQNC/NCOFN'
      END IF
C
  294 FORMAT(2I4,2I3, 7F8.3, /, 14X, 7F8.3, /, 14X, 8F8.3)
  295 FORMAT(44X, A50)
   96 FORMAT(2I5, 13I5, /, 10X, 8I5)
  969 FORMAT(2I4,1X,21I3)
  970 FORMAT(1X,21I3)
   97 FORMAT(2I4, 6F8.3, /, 8X, 7F8.3, /, 8X, 8F8.3)
   98 FORMAT(6F8.4, /, 7F8.4, /, 8F8.4)
   99 FORMAT(7F8.4, /, 7F8.4, /, 8F8.4)
c mrm  21 FORMAT(A27, 4(F8.3,2X))
   21 FORMAT(A27, 1p, 4e11.3)
c mrm  981 FORMAT(A27, F8.3,2X,F8.2,2X,F8.2)
  981 FORMAT(A27, 1p, 3e11.3)
   23 FORMAT(A46, I5)
   85 FORMAT(A44, A50)
  284 FORMAT(A32, /, 'Name    ITNWQ    L    I    J    K       Conc')
c
c----------------------------------------------------------------------
      IWQWIN=0
!      call wqwin            ! Ji, 9/17/99
c--------------------------------------------------------------------
C
C
C Read in mapping information for spatially-varying parameters (unit #7).
C
      DO K=1,KC
       DO L=2,LA
        IWQZMAP(L,K)=1
       END DO
      END DO
C
      DO L=2,LA
       DO K=1,5
       FLXPND(L,K)=1.0
       ENDDO
      ENDDO
c mrm removed the following hardwire so wqwcmap.inp will be read whenever
c     IWQZ is greater than 1 which implies the map file must be read in:
c      IF(ISWQCMAP.EQ.1)THEN
      IF (IWQZ .GT. 1) THEN
C
       OPEN(1,FILE='wqwcmap.inp',STATUS='UNKNOWN')
C
       CALL SKIPCOMM(1,CCMRM)	
       IN=0
       DO M=2,LA
       READ(1,*,END=1111) I,J,IWQZX,P_TT,TN_T,DO_TT,w_tt,A_TT
!        write(*,'(A11,2I6,5F7.2)')'Par. Map = ',
!     *    M,IWQZX,P_TT,TN_T,DO_TT,w_tt,A_TT
        IF (LIJW(I,J).LT.1 ) THEN
        PRINT*, 'i, j, k, IJCT(i,j) = ', I,J,K,LIJW(I,J)
        STOP 'ERROR!! invalid (i,j) in FILE wqwcmap.inp'
        END IF
        IN=IN+1
        L = LIJW(I,J)
        FLXPND(L,1)=P_TT
        FLXPND(L,2)=TN_T
        FLXPND(L,3)=DO_TT
  !      FLXPND(L,4)=w_tt   ! adj vertical transport
  !      FLXPND(L,5)=A_TT   ! adj addy viscosity
        DO K=1,KC
        IWQZMAP(L,K)=IWQZX
        WRITE(2,31) L,I,J,K,IWQZMAP(L,K)
        ENDDO
       ENDDO 
 1111  CONTINUE
       CLOSE(1)
      ENDIF


c       OPEN(1,FILE='wqwcmap.inp',STATUS='UNKNOWN')

c       WRITE(2,999)
c       READ(1,30) (TITLE(M), M=1,3)
c       WRITE(2,30) (TITLE(M), M=1,3)
C      READ(1,999)
c       READ(1,999)
c       WRITE(2,999)
c       WRITE(2,32)
c       IN=0
c       IJKC=IC*JC*KC
c       DO M=1,IJKC
c        READ(1,*,END=1111) I,J,K,IWQZX
c        IN=IN+1
c        IF (LIJW(I,J).LT.1 ) THEN
c          PRINT*, 'i, j, k, IJCT(i,j) = ', I,J,K,LIJW(I,J)
c          STOP 'ERROR!! invalid (i,j) in FILE wqwcmap.inp'
c        END IF
c        L = LIJW(I,J)
c        IWQZMAP(L,K)=IWQZX
c        WRITE(2,31) L,I,J,K,IWQZMAP(L,K)
c       END DO
c 1111  CONTINUE
c       IF (IN.NE.(LA-1)*KC) THEN
c        PRINT*, 'all active water cells should be mapped for WQ par.'
c        STOP 'ERROR!! number of lines in FILE wqwcmap.inp =\ (LA-1)'
c       END IF
C
c       CLOSE(1)
C
c      END IF
C      write(6,6666)'read wqwcmap.inp'
c
c mrm Code below added by M. Morton 03/29/98:
C Read in mapping information for spatially-varying benthic fluxes.
c Formulated for Peconic Bay data which includes %mud for each cell as
c well as mapping to both mud and sand fluxes.  Subroutine RWQBEN2
c contains the code to interpolate the final flux for the cell based
c on the percent mud and the mud/sand fluxes.
C
      IF (iwqben .EQ. 2) THEN
c
       DO K=1,2
        DO L=2,LA
         ibenmap(L,K)=1
         xbenmud(L) = 50.0
        END DO
       END DO
C
       OPEN(1,FILE='wqbenmap.inp',STATUS='UNKNOWN')
C
       WRITE(2,999)
       DO M=1,4
        READ(1,30) TITLE(M)
        WRITE(2,30) TITLE(M)
       END DO
c       READ(1,999)
       WRITE(2,999)
       WRITE(2,33)
       IN=0
       IJKC=IC*JC
       DO M=1,IJKC
         READ(1,*,END=1112) I, J, xmud, izmud, izsand
         IN=IN+1
         IF (LIJW(I,J).LT.1 ) THEN
           PRINT*, 'i, j, k, IJCT(i,j) = ', I,J,LIJW(I,J)
           STOP 'ERROR!! invalid (i,j) in FILE wqbenmap.inp'
         END IF
         L = LIJW(I,J)
         ibenmap(L,1) = izmud
         ibenmap(L,2) = izsand
         xbenmud(L) = xmud / 100.0
         WRITE(2,34) L, I, J, xbenmud(L), ibenmap(L,1), ibenmap(L,2)
       END DO
 1112  CONTINUE
       IF (IN .NE. (LA-1)) THEN
         PRINT*, 'all active water cells should be mapped for WQ par.'
         STOP 'ERROR!! number of lines in FILE wqbenmap.inp <> (LA-1)'
       END IF
C
       CLOSE(1)
C
      END IF
C      write(6,6666)'read wqbenmap.inp'
C
cmrm Code above to previous cmrm added by M. Morton 03/29/98.
      CLOSE(2)
C
   30 FORMAT(A79)
   31 FORMAT(15I5)
   32 FORMAT('    L    I    J    L IWQZMAP')
   33 FORMAT('     L    I    J   MUD IZmud IZsand')
   34 format(' ',3i5, f6.2, 2i6)
C
      RETURN
      END

