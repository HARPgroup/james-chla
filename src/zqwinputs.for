!**********************************************************************C
!
      SUBROUTINE RWQINPS(IWQDT,DT,IC,JC,KC,LA,NTSPTC,IWQS)
!
!**********************************************************************C
!
! **  LAST MODIFIED BY JIAN SHEN ON  July 2008
!
! **  
!**********************************************************************C
! 
!  Water Quality input file for implified version
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
      PARAMETER (CONV1=1.0E3,CONV2=8.64E4,NACP=12)
      DIMENSION XPSL(NWQVM),XDSL(NWQVM),IACP(NACP),IACP1(NACP)
      DIMENSION WQTMC1(NWQZM),WQTMC2(NWQZM),WQTMD1(NWQZM),WQTMD2(NWQZM),
     &  WQTMG1(NWQZM),WQTMG2(NWQZM),WQTMM1(NWQZM),WQTMM2(NWQZM),
     $  WQTMp1(NWQZM), WQTMp2(NWQZM),
     &  WQKG1C(NWQZM),WQKG2C(NWQZM),WQKG1D(NWQZM),WQKG2D(NWQZM),
     &  WQKG1G(NWQZM),WQKG2G(NWQZM),WQKG1M(NWQZM),WQKG2M(NWQZM),
     $         WQKG1p(NWQZM), WQKG2p(NWQZM)
      INTEGER isSKIP
      CHARACTER TITLE(5)*79, ccmrm*1, ccmr*1
      DATA   IACP/3,5,6,8,9,10,12,13,14,15,18,19/
      DATA  IACP1/3,6,6,9,9,10,13,13,14,15,18,19/
!
      OPEN(1,FILE='wq3dwcs.inp',STATUS='UNKNOWN')
      OPEN(2,FILE='wq3dwcs.out',STATUS='UNKNOWN')
!
      DO I=1,21
	IWTRC(I)=0
	ENDDO
	DO I=1,NACP
	if(IWQS.eq.1) then
	IWTRC(IACP(I))=1
	else
	IWTRC(IACP1(I))=1	
	endif
	ENDDO
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
     $          NSMTS,NWQKDPT

      NWQZ   =  NWQZM           ! 2
      NWQPS  =  NWQPSM          ! 3
      NWQTS  =  NWQTSM          ! 5
      NSMG   =  NSMGM           ! 7
      NSMZ   =  NSMZM           ! 8
      NTSSMV =  NTSSMVM         ! 9
      NSMTS  =  NSMTSM          ! 10

! check array dimensions 

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
      READ(1,*) IWQDT,IWQBEN,IWQKA, WQKRO !,REAC(1)

      IWQM = 0  ! reduced version
      IWQSI =0  ! silica is not simulated
      IWQFCB =0 ! no FC is simulated
      IWQSRP =2 ! use sediment for PO4 adsorption
      IWQSTOX=0 ! no cyanobacteria salinity toxicity
	iSHADE =0 ! no shading
 

      DO LL = 2, LA      ! no shadal
        pSHADE(LL) = 1.0
      END DO
!
      if(mod(IWQDT,2).NE.0) IWQDT=IWQDT+1
	NWQKDPT=IWQDT /2       ! Double check J.S. if not 3-level time steping. WQ model is called every 2-timestep,need divid 2 here
      DTD = DT/86400.0
      DTWQ = DTD*REAL(IWQDT) ! WQ model DT, usually = 2*DTD
      DTWQO2 = DTWQ*0.5
!
!C04
!
      WRITE(2,*)'C04'
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQZ,IWQNC,IWQRST,LDMWQ,IDNOTRVA
      IF(IDNOTRVA.EQ.1) THEN
	IDNOTRVA =22
	ELSE
	IDNOTRVA =0
	ENDIF
      NDMWQ =1 
	NDDOAVG =0 
	NDLTAVG = 91
!
!C05
!      
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQICI,IWQPSL,IWQNPL
       
	IWQAGR = 0 ! constant kinetics for algae growth
      IWQSTL = 0
      IWQSUN = 2 !use hourly solar rad. from wqsun.inp file
	isDIURDO =0 ! not active diurnal D.O.

	WQDIUDT = 24 ! not used           ! Check J.S.
      IWQDIUDT = NINT(WQDIUDT*3600.0/DT)

C 
C06
C
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) IWQTS,TWQTSB,TWQTSE,WQTSDT,  isWQMIN, isWQMAX,
     *          isCOMP,isBIN, IAVGBIN,IRESTHR
!      
	ISMTSDT=NTSPTC/24*IAVGBIN          ! in hours
	IAVGBIN=NTSPTC/24*IAVGBIN
      IRESTHR=NTSPTC/24*IRESTHR
      isWQAVG =0   ! not used check J.S.

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
C7
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
      READ(1,*) WQKHNG,WQKHPG,WQKHNM,WQKHPM

      WQKHNC =WQKHNG
      WQKHND =WQKHNG
      WQKHPC =WQKHPG
	WQKHPD =WQKHPG
      WQKHS  =0.05
	WQSTOX =1.0

      WQSTOX = WQSTOX*WQSTOX
!
!C09:constant parameters for ALGAE (see Table 3-1)
!    
      Write(2,*)'C09------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ
      READ(1,*) WQKETSS(I),WQKECHL(I),WQKEB(I),WQCHLG(I),WQDOPG,
     + WQCHLM,WQDOPM, WQKEMA(I),REAC(I)
      WRITE(2,80)'* Light extinc. coeff. due to TSS & Chl           '
      WRITE(2,81)' : KeTSS (/m per g/m^3)  = ', WQKETSS
      WRITE(2,81)' : KeChl (/m per mg/m^3) = ', WQKECHL
      IF (WQKECHL(I) .LT. 0.0) THEN
        WRITE(2,80) '* Use Riley (1956) equation for WQKECHL          '
        WRITE(2,80) ' : KeChl = 0.054*CHL**0.667 + 0.0088*CHL         '
      END IF
      ENDDO
!
      WQCHLC(I) = 0.04 
     	WQCHLD(I) = 0.04
      DOPTc  = 1.0
      DOPTd  = 1.0
      WQDOPD =1
	WQDOPC =1
C

      WRITE(2,*)'* Carbon-to-Chl ratio (g C per mg Chl)',WQCHLG

      WRITE(2,*)'* Depth (m) of maximum algal growth   ',WQDOPG
!
      WQCHLC(I)=1.0/(WQCHLC(I)+ 1.E-12)    ! be careful! inversed!, 9/15/02
      WQCHLD(I)=1.0/(WQCHLD(I)+ 1.E-12)
      WQCHLG(I)=1.0/(WQCHLG(I)+ 1.E-12)
      WQCHLM=1.0/(WQCHLM+ 1.E-12)
!
!C10:constant parameters for ALGAE (see Table 3-1)
!     
      Write(2,*)'C10------------------------------'   
	CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQI0,WQISMIN,WQFD,WQCIA,WQCIB,WQCIC,WQCIM,PARADJ
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
!      IF (IWQSUN .EQ. 2) THEN   ! 2010 revised J.S. !!!
!        WQISMIN = 0.0
!      END IF
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
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQTMG1,WQTMG2, WQTMM1, WQTMM2,iTMgOP
      WRITE(2,82)'C11' 
      WQTMC1 = 10.0 
      WQTMC2 = 23.0  
      WQTMD1 = 15.0  
      WQTMD2 = 23.0
      WQTMp1 = 22.0  
	WQTMp2 = 28.0
C12
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQKG1G,WQKG2G,WQKG1M,WQKG2M
      WRITE(2,82)'C12'      
      WQKG1C = 0.003  
	WQKG2C = 0.005  
	WQKG1D = 0.001  
	WQKG2D = 0.005
      WQKG1p = 0.001  
      WQKG2p = 0.001
C13
C   
      Write(2,*)'C13------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQTRG,WQKTBG,WQTRM,WQKTBM
      WRITE(2,82)'C13'
      WQTRC = 20.0  
	WQTRD = 20.0  
      WQKTBC = 0.069 
	WQKTBD = 0.069 

C
C Constant parameters for CARBON (see Table 3-2)
C14
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFCLP,WQFCDP,WQFCDG, WQKHRG,ICBOD
      WRITE(2,82)'C14'
      WQFCRP =0.0
      WQFCDC =0.0
	WQFCDD =0.0
      WQKHRC =0.5
	WQKHRD =0.5
      CFCDCWQ = 1.0 - WQFCDC
      CFCDDWQ = 1.0 - WQFCDD
      CFCDGWQ = 1.0 - WQFCDG
      XC = ABS(1.0 - (WQFCLP+WQFCDP))
      IF (XC .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FCLP+FCDP not equal to 1.0'
        WRITE(2,*)
      END IF   
C
C15
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFCLPM,WQFCDPM,WQFCDM,WQKHRM
      WQFCRPM = 0.0
      WRITE(2,82)'C15'
C16
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQKLC,WQKDC
      WRITE(2,82)'C16'
     *         
      WQKRC =   0.0
      WQKRCALG =0.0
	WQKLCALG =0.0
      WQKDCALG =0.0
      WQKDCALM =0.0
C17
      Write(2,*)'C17------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQTRMNL,WQKTMNL,WQKHORDO,WQKHDNN, WQAANOX
      WRITE(2,82)'C17'
	WQTRHDR = 20
      WQKTHDR =0.069	
      WQAANOX = WQAANOX*WQKHORDO
C
C Constant parameters for PHOSPHORUS (Table 3-3)
C
C18
C  
C
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQFPLP,WQFPDP,WQFPIP, WQFPLG, WQFPDG,WQFPIG,WQKPO4P
      WRITE(2,82)'C18'
      WQFPRP =0.0
      WQFPRC =0.0
	WQFPRD =0.0
      WQFPRG =0.0
	WQFPLC =0.0
	WQFPLD =0.0

      XP = ABS(1.0 - (WQFPLP+WQFPDP+WQFPIP))
      IF (XP .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPLP+FPDP+FPIP not equal to 1.0'
        WRITE(2,*)
      END IF
      XPG = ABS(1.0 - (WQFPLG+WQFPDG+WQFPIG))
      IF (XPG .GT. 0.0001) THEN
        WRITE(2,*)
        WRITE(2,*) ' WARNING!  FPLg+FPDg+FPIg not equal to 1.0'
        WRITE(2,*)
      END IF

C19

	CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFPLPM,WQFPDPM,WQFPIPM,WQFPLM,WQFPDM,WQFPIM,WQAPCM
      WRITE(2,82)'C19'
      WQFPRPM =0.0
      WQFPRM  =0.0
C20
C   
c      CALL SKIPCOMM(1,CCMRM)
c
      WQFPDC = 0.0
	WQFPDD = 1.0
      WQFPIC = 0.0
	WQFPID = 1.0

      IF (IWQSRP.NE.1 .AND. IWQSRP.NE.2) THEN
        WQKPO4P = 0.0
        WRITE(2,80)': no sorption of PO4t/SA, so KPO4p is forced to 0 '
      END IF

C21
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQKLP,WQKDP,WQCP1PRM, WQCP2PRM,WQCP3PRM
      WRITE(2,82)'C21'
      WQKRP =0
      WQKRPALG =0.0
	WQKLPALG =0.0
	WQKDPALG =0.0
C
C Constant parameters for NITROGEN (Table 3-4)
C
C22
C   
C
      Write(2,*)'C22------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFNLP,WQFNDP,WQFNIP,WQFNLG,WQFNDG,WQFNIG,WQANCG
      WRITE(2,82)'C22'
	WQFNRP =0.0
      WQFNRC =0.0
	WQFNRD =0.0
      WQFNLC =0.0
	WQFNLD =0.0
      WQFNRG =0.0

      XN = ABS(1.0 - (WQFNLP+WQFNDP+WQFNIP))
      IF (XN .GT. 0.0001) THEN
        WRITE(*,*)
        WRITE(*,*) ' WARNING!  FNRP+FNLP+FNDP+FNIP not equal to 1.0'
        WRITE(*,*)
      END IF
      XNG = ABS(1.0 - (WQFNLG+WQFNDG+WQFNIG))
      IF (XNG .GT. 0.0001) THEN
        WRITE(*,*)
        WRITE(*,*) ' WARNING!  FNRg+FNLg+FNDg+FNIg not equal to 1.0'
        WRITE(*,*)
      END IF
C23
      Write(2,*)'C23------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*)WQFNLPM,WQFNDPM,WQFNIPM,WQFNLM,WQFNDM,WQFNIM,WQANCM,
     + WQFNCm
      WQFNRPM =0.0
      WQFNRM  =0.0
      WRITE(2,82)'C23'
C24
C   
c      Write(2,*)'C24------------------------------' 
c      CALL SKIPCOMM(1,CCMRM)
c      READ(1,*) WQFNDC,WQFNDD,WQFNDG,WQFNDM,WQFNIC,WQFNID,WQFNIG,
c     * WQFNIM,WQANCC,WQANCD,WQANCG,WQANCM
      WQFNDC =0.0
	WQFNDD =0.0
      WQFNIC =0.0
	WQFNID =0.0
      WQANCC =0.0
	WQANCD =0.0
	WQANCC =0.0
	WQANCD =0.0
C25
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQANDC,WQNITM,WQKHNDO,WQKHNN,WQTNIT,WQKN1,WQKN2
      WRITE(2,82)'C25'
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
      READ(1,*) WQKLN,WQKDN

      WQKRN =0.0
      WQKRNALG =0.0
	WQKLNALG =0.0
	WQKDNALG =0.0
C27
C Constant parameters for SILICA (Table 3-5)
C
c      Write(2,*)'C27------------------------------' 
c      CALL SKIPCOMM(1,CCMRM)
c      READ(1,*) WQFSPP,WQFSIP,WQFSPD,WQFSID,WQASCD,WQKSAP,WQKSU,
c     *  WQTRSUA,WQKTSUA
c      IF (IWQSRP.NE.1 .AND. IWQSRP.NE.2) THEN
c        WQKSAP = 0.0
c        WRITE(2,80)': no sorption of PO4t/SA, so KSAp is forced to 0  '
c      END IF
C
C28
C Constant parameters for COD & DO (Table 3-6)
C
      Write(2,*)'C28------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQAOCR,WQAONT,WQKTR,WQKHCOD,WQKCD,WQTRCOD,WQKTCOD,
     *   WQAOCRpm, WQAOCRrm
      WRITE(2,82)'C28'
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
c      Write(2,*)'C29------------------------------' 
c      CALL SKIPCOMM(1,CCMRM)
c      READ(1,*) WQKHBMF,WQBFTAM,WQTTAM,WQKTAM,WQTAMDMX,WQKDOTAM,
c     *  WQKFCB,WQTFCB
c      WRITE(2,82)
c     *  '* DO where TAM release is half anoxic one  = ',WQKHBMF
c     * ,'  Anoxic release of TAM (mol/m^2/d)        = ',WQBFTAM
c     * ,'  Reference temp for TAM release (oC)      = ',WQTTAM
c     * ,'  Temperature effect on TAM release        = ',WQKTAM
c     * ,': TAM solubility at anoxic cond. (mol/m^3) = ',WQTAMDMX
c     * ,'  Constant relating TAM solubility to DO   = ',WQKDOTAM
c     * ,'* First-order die-off rate at 20oC (/d)    = ',WQKFCB
c     * ,'  Temperature effect on bacteria die-off   = ',WQTFCB
C
C Set up look-up table for temperature dependency over -15oC to 40oC
C

      DO M=1,NWQTD
        WTEMP =1.00*REAL(M-1)*0.1 - 4.95   ! -5C -> 50C, Ji, 8/29/02
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
        WQTDGG(M,I)=max(WQTDGG(M,I),T_Lim)
C J.S.  5/5/98
        WQTDGM(M,I)=1.
        IF(IDNOTRVA.GT.0.or.iTMgop.EQ.1)THEN
         IF (WTEMP.LT.WQTMM1(I)) THEN
         WQTDGM(M,I)=EXP(-WQKG1M(I)*(WTEMP-WQTMM1(I))*(WTEMP-WQTMM1(I)))
         END IF
         IF (WTEMP.GT.WQTMm2(I)) THEN
         WQTDGM(M,I)=EXP(-WQKG2M(I)*(WTEMP-WQTMM2(I))*(WTEMP-WQTMM2(I)))
         END IF
         WQTDRM(M) = EXP( WQKTBM*(WTEMP-WQTRM) )
        END IF
C J.S.
C MRM: 06/20/98
C  The following WQTDGP variable is a temperature related adjustment
C  to the predation and/or basal matabolism rate to allow diatoms
C  to bloom in winter (or other time of year).
        WQTDGp(M,I)=1.
        IF (WTEMP.LT.WQTMp1(I)) THEN
        WQTDGp(M,I)=EXP(-WQKG1p(I)*(WTEMP-WQTMp1(I))*(WTEMP-WQTMp1(I)))
        END IF
        IF (WTEMP.GT.WQTMD2(I)) THEN
        WQTDGp(M,I)=EXP(-WQKG2p(I)*(WTEMP-WQTMp2(I))*(WTEMP-WQTMp2(I)))
        END IF
        
       ENDDO
       ENDDO
       
       DO M=1,NWQTD
        WTEMP =1.00*REAL(M-1)*0.1 - 4.95   ! -5C -> 50C, Ji, 8/29/02      
c MRM
C
C        write(1,555)WTEMP,WQTDGC(M),WQTDGD(M),WQTDGG(M),WQTDGM(M)
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
        END DO
!
! spatially/temporally constant ALGAL GROWTH, RESPIRATION & PREDATION RATES
!     
!30    
!
      Write(2,*)'C30------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ
      READ(1,*)WQPMG(I),WQBMRG(I),WQPRRG(I),WQPMM(I),WQBMRM(I),WQPRRM(I)
      WRITE(2,82)'C30'
      WQPMC(I) =0.0
	WQPMD(I) =0.0
      WQBMRC(I)=0.0
      WQBMRD(I)=0.0
      WQPRRC(I)=0.0
	WQPRRD(I)=0.0

      IF (IWQAGR.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,80)'* Algal growth rate (/day)                        '
        WRITE(2,21)' : (PMc, PMd, PMg)       = ', WQPMC(I),WQPMD(I),
     *    WQPMG(1)
        WRITE(2,80)'* Algal basal metabolism rate (/day)              '
        WRITE(2,21)' : (BMRc, BMRd, BMRg)    = ', WQBMRC(I),WQBMRD(I),
     *    WQBMRG(1)
        WRITE(2,80)'* Algal predation rate (/day)                     '
        WRITE(2,21)' : (PRRc, PRRd, PRRg)    = ', WQPRRC(I),WQPRRD(I),
     *    WQPRRG(1)
      END IF     
      ENDDO
        DO I=2,IWQZ
    !    WQPMC(I)=WQPMC(1)  
    !    WQPMD(I)=WQPMD(1)
    !    WQPMG(I)=WQPMG(1)
    !    WQPMM(I)=WQPMM(1)
    !    WQBMRC(I)=WQBMRC(1)
    !    WQBMRD(I)=WQBMRD(1)
    !    WQBMRG(I)=WQBMRG(1)
    !    WQBMRM(I)=WQBMRM(1)
    !    WQPRRC(I)=WQPRRC(1)
    !    WQPRRD(I)=WQPRRD(1)
    !    WQPRRG(I)=WQPRRG(1)
    !    WQPRRM(I)=WQPRRM(1)
         WQKEB(I)=WQKEB(1)
        END DO

!
!C31  (46)
! spatially/temporally constant SETTLING VELOCITIES and REAERATION FACTOR
!
      Write(2,*)'C31------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      DO I=1,IWQZ
      READ(1,*)WQWSG(I),WQWSLP(I),WQWSM,FSTOC,FSTOP,FSTON 
      WRITE(2,82)'31'
      WQWSC(I)=0.0
	WQWSD(I)=0.0
	WQWSRP(I)=0.0
      WQWSS(I) =0.0


      IF (IWQSTL.NE.1) THEN
        WRITE(2,999)
        WRITE(2,90) TITLE(1)
        WRITE(2,80)'* Algal settling rate (m/day)                     '
        WRITE(2,21)' : (WSg)       = ',  WQWSG(I)
        WRITE(2,80)'* POM settling rate (m/day)                       '
        WRITE(2,21)' : WSlp)          = ', WQWSLP(I)
      END IF
      ENDDO
      
 !       DO I=2,IWQZ
 !         WQWSC(I)=WQWSC(1)
 !         WQWSD(I)=WQWSD(1)
 !        WQWSG(I)=WQWSG(1)
 !         WQWSRP(I)=WQWSRP(1)
 !         WQWSLP(I)=WQWSLP(1)
 !         WQWSS(I)=WQWSS(1)
          ! REAC(I)=REAC(1)
 !       END DO
   

!
! spatially/temporally constant INITIAL CONDITIONS: WQCHLx=1/WQCHLx
! read data points & do internal interpolation?
!
!C32 (C44)
!
      Write(2,*)'C32------------------------------' 
      CALL SKIPCOMM(1,CCMRM)
      READ(1,*) WQV(1,1,3),WQV(1,1,5),WQV(1,1,6),WQV(1,1,8),WQV(1,1,9)
     $ ,WQV(1,1,12),WQV(1,1,13)
      READ(1,*) WQV(1,1,14),WQV(1,1,15),WQV(1,1,18),WQV(1,1,19)
     $ ,WQV(1,1,22),WQMCMIN

      WQV(1,1,1)=0.0
      WQV(1,1,2)=0.0
      WQV(1,1,4)=0.0
      WQV(1,1,7)=0.0
      WQV(1,1,11)=0.0
      WQV(1,1,16)=0.0
      WQV(1,1,17)=0.0
      WQV(1,1,20)=0.0
      WQV(1,1,21)=0.0

      IF (IWQICI.NE.1) THEN

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
        READ(3,*) II,JJ,tt_p
        IF (LIJW(II,JJ).LT.1) THEN
          PRINT*, 'i, j, LIJW(II,JJ) = ', II,JJ,LIJW(II,JJ)
           STOP 'ERROR!! invalid (i,j) in FILE wqcalgmp.inp'
        END IF
        LL=LIJW(II,JJ)
        SMAC(LL)=1.0
        WQV(LL,1,IDNOTRVA)=tt_p
        WQVO(LL,1,IDNOTRVA)=tt_p
!	  WQPMM(LL)=gm_t
!	  WQBMRM(LL)=rm_t
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
      READ(1,*) WQBFPO4D(1),WQBFNH4(1),WQBFNO3(1),
     *  WQBFCOD(1),WQBFO2(1)
      WQBFSAD(1)=0.0

      IF (IWQBEN.EQ.0) THEN
        WRITE(2,21)' : (PO4d, NH4, NO3)      = ',WQBFPO4D(1),WQBFNH4(1),
     *    WQBFNO3(1)
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
	DO NW=1,NWQV
	NWQCSR(NW)=0
	ENDDO
      READ(1,*) (NWQCSR(IACP(NW)), NW=1,NACP)
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
      DO M=1,NWQOBS
      DO NW=1,NWQV
      IWQOBS(M,NW)=0
	ENDDO
	ENDDO
      IF (NWQOBS.GT.0) THEN
       DO M=1,NWQOBS
        READ(1,*)IWQCBS(M),JWQCBS(M),(IWQOBS(M,IACP(NW)), NW=1,NACP)
        WRITE(2,969)IWQCBS(M),JWQCBS(M),(IWQOBS(M,IACP(NW)),NW=1,NACP)
       END DO
      END IF
!
!C37 WEST BOUNDARY
!
      WRITE(2,*)'C37------------------------------'
      CALL SKIPCOMM(1,CCMRM)
      DO M=1,NWQOBW
      DO NW=1,NWQV
      IWQOBW(M,NW)=0
	ENDDO
	ENDDO
      IF (NWQOBW.GT.0) THEN
       DO M=1,NWQOBW
        READ(1,*)IWQCBW(M),JWQCBW(M),(IWQOBW(M,IACP(NW)), NW=1,NACP)
        WRITE(2,969)IWQCBW(M),JWQCBW(M),(IWQOBW(M,IACP(NW)), NW=1,NACP)
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
      DO M=1,NWQOBE
      DO NW=1,NWQV
      IWQOBE(M,NW)=0
	ENDDO
      ENDDO
      IF (NWQOBE.GT.0) THEN
       DO M=1,NWQOBE
        READ(1,*)IWQCBE(M),JWQCBE(M),(IWQOBE(M,IACP(NW)), NW=1,NACP)
        WRITE(2,969)IWQCBE(M),JWQCBE(M),(IWQOBE(M,IACP(NW)), NW=1,NACP)
        END DO
       END IF
!
!C39 NORTH BDRY
!
      WRITE(2,*)'C39------------------------------'
      CALL SKIPCOMM(1,CCMRM)
	DO M=1,NWQOBN
      DO NW=1,NWQV
      IWQOBN(M,NW)=0
      ENDDO
	ENDDO

      IF (NWQOBN.GT.0) THEN
       DO M=1,NWQOBN
        READ(1,*)IWQCBN(M),JWQCBN(M),(IWQOBN(M,IACP(NW)), NW=1,NACP)
        WRITE(2,969)IWQCBN(M),JWQCBN(M),(IWQOBN(M,IACP(NW)), NW=1,NACP)
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
	  write(2,*) 'IWPS  JWPS  WPSER '      
      CALL SKIPCOMM(1,'C')
      DO M=1,IWQPS
        READ(1,*) I,J,K
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
C42
      SODMULT(1)=1
!     REAC(1)=1
      CALL SKIPCOMM(1,CCMRM)
      Read(1,*)IDIA,SODMULT(1),BFHN,BFNO, BFP,
     &  DO_G,DO_R,AN_Lim,P_lim,T_Lim,REAC_1,F_AB,
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
	
	WRITE(2,*)'Card 42'
	write(2,*)'SODMULT= ',SODMULT(1)
	write(2,*)'BFHN =   ',BFHN
	write(2,*)'BFNO = ',BFNO 
	write(2,*)'BFP  = ',BFP  
	write(2,*)'DO_G = ',DO_G 
	write(2,*)'DO_R = ',DO_R 
	write(2,*)'AN_Lim = ',AN_Lim 
	write(2,*)'P_lim  = ',P_lim  
	write(2,*)'T_Lim  = ',T_Lim  
	write(2,*)'F_AB = ',F_AB 
	write(2,*)'RN_R = ',RN_R 
	write(2,*)'RP_R = ',RP_R 
	write(2,*)'RC_R = ',RC_R 
	
      DO I=1,IWQZ
      IF(IDIA.GT.0)REAC(I)=REAC(I)*REAC_1
      SODMULT(I)=SODMULT(1)
	ENDDO 
      
c43 Forced Temperature
C
      CALL SKIPCOMM(1,CCMRM)	     
      read(1,*)IFOCETMP,TPHPSE,   TPMAG, TMPADD  
       WRITE(2,*)'Card 43'
       write(2,*)IFOCETMP,TPHPSE, TPMAG, TMPADD  

c49 --
c
C spatially/temporally-constant values for NON-POINT SOURCE INPUT in (kg/d)
C except xDSQ (m^3/s), xNO2 (mg/L), xNTAM (kmol/d), xNFCB (MPN/100mL).
c mrm  Note: This group is now used for dry atmospheric deposition with
c      constituent units of g/m2/day except FCB which is MPN/m2/day.
C
!
!C42:Constant Dry Atmospheric Deposition
!  Set to zero
      DO NW=1,NWQV
	XDSL(NW)=0.0
      WQWDSL(1,1,NW)=0.0
	ENDDO
      XDSQ=0

C      WRITE(2,*)'Card 42---------------------------------'
C      CALL SKIPCOMM(1,CCMRM)
C      READ(1,*) XDSQ,(XDSL(NW),NW=1,6)
C      READ(1,*) (XDSL(NW),NW=7,13)
C      READ(1,*) (XDSL(NW),NW=14,NWQV)
C      IF (IWQNPL.NE.1) THEN
C        WRITE(2,999)
C        WRITE(2,90) TITLE(1)
C        WRITE(2,21)' : (DSQ, Bc, Bd, Bg)     = ',XDSQ,(XDSL(NW),NW=1,3)
C        WRITE(2,21)' : (RPOC, LPOC, DOC)     = ',(XDSL(NW),NW=4,6)
C        WRITE(2,21)' : (RPOP,LPOP,DOP,PO4t)  = ',(XDSL(NW),NW=7,10)
C        WRITE(2,21)' : (RPON, LPON, DON)     = ',(XDSL(NW),NW=11,13)
C        WRITE(2,21)' : (NH4, NO3)            = ',(XDSL(NW),NW=14,15)
C        WRITE(2,21)' : (SU, SA, COD, DO)     = ',(XDSL(NW),NW=16,19)
C        WRITE(2,981)' : (TAM, FCB)            = ',(XDSL(NW),NW=20,NWQV)
C        DO NW=1,18
C          WQWDSL(1,1,NW) = XDSL(NW) * CONV1
C        END DO
C        WQTT = XDSQ*CONV2
C        WQWDSL(1,1,19) = XDSL(19) * WQTT
C        WQWDSL(1,1,20) = XDSL(20) * CONV1
C        WQWDSL(1,1,NWQV) = XDSL(NWQV) * WQTT
C
        WQDSQ(1,1) = XDSQ
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
C      END IF
!
!43
!
! Wet Atmospheric Deposition concentrations (mg/L; MPN/L) 
! and will be multiplied by the rainfall flow rate into each grid cell
! to get a load in kg/day.
!
      DO NW=1,NWQV
      wqatm(NW)=0.0
	ENDDO

c      WRITE(2,*)'Card 43---------------------------------'
c      CALL SKIPCOMM(1,CCMRM)
c      READ(1,*) (wqatm(NW),NW=1,6)
c      READ(1,*) (wqatm(NW),NW=7,13)
c      READ(1,*) (wqatm(NW),NW=14,NWQV)
c      WRITE(2,999)
c      WRITE(2,90) TITLE(1)
c      WRITE(2, 21)' : (Bc, Bd, Bg)          = ',(wqatm(NW),NW=1,3)
c      WRITE(2, 21)' : (RPOC, LPOC, DOC)     = ',(wqatm(NW),NW=4,6)
c      WRITE(2, 21)' : (RPOP,LPOP,DOP,PO4t)  = ',(wqatm(NW),NW=7,10)
c      WRITE(2, 21)' : (RPON, LPON, DON)     = ',(wqatm(NW),NW=11,13)
c      WRITE(2, 21)' : (NH4, NO3)            = ',(wqatm(NW),NW=14,15)
c      WRITE(2, 21)' : (SU, SA, COD, DO)     = ',(wqatm(NW),NW=16,19)
c      WRITE(2,981)' : (TAM, FCB)            = ',(wqatm(NW),NW=20,NWQV)
c
C
c51
c
C Input/output file names for spatially and/or temporally varying parameters
C
       RSTOFN='wqwcrst.out'

      IF (IWQRST.GE.1) THEN
        OPEN(99,FILE=RSTOFN,STATUS='UNKNOWN')
        CLOSE(99,STATUS='DELETE')
        OPEN(99,FILE=RSTOFN,STATUS='UNKNOWN')
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
      ICIFN='none'
      AGRFN='none'
	STLFN='none'
      SUNFN='none'
      BENFN='none'
	PSLFN='none'
      NPLFN='none'
	NCOFN='wq3dnc.log'
C
      IF (IWQNC.EQ.1) THEN
        OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=NCOFN,STATUS='UNKNOWN')
        WRITE(1,284)'* Negative concentration occurs:'
        CLOSE(1)
       ELSE
 !       IF (NCOFN(1:4).NE.'none') STOP 'ERROR!! invalid IWQNC/NCOFN'
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
c mrm removed the following hardwire so wqwcmap.inp will be read whenever
c     IWQZ is greater than 1 which implies the map file must be read in:
c      IF(ISWQCMAP.EQ.1)THEN
      DO L=2,LA
       DO K=1,5
       FLXPND(L,K)=1.0
       ENDDO
      ENDDO
      IF (IWQZ .GT. 1) THEN
C
       OPEN(1,FILE='wqwcmap.inp',STATUS='UNKNOWN')
C
       WRITE(2,999)
       READ(1,30) (TITLE(M), M=1,3)
       WRITE(2,30) (TITLE(M), M=1,3)
C      READ(1,999)
       READ(1,999)
       WRITE(2,999)
       WRITE(2,32)
       IN=0
       DO M=2,LA
       READ(1,*,END=1111) I,J,IWQZX,P_TT,TN_T,DO_TT,w_tt,A_TT
        write(*,'(A11,2I6,5F7.2)')'Par. Map = ',
     *    M,IWQZX,P_TT,TN_T,DO_TT,w_tt,A_TT
        IF (LIJW(I,J).LT.1 ) THEN
        PRINT*, 'i, j, k, IJCT(i,j) = ', I,J,K,LIJW(I,J)
        STOP 'ERROR!! invalid (i,j) in FILE wqwcmap.inp'
        END IF
        IN=IN+1
        L = LIJW(I,J)
        FLXPND(L,1)=P_TT
        FLXPND(L,2)=TN_T
        FLXPND(L,3)=DO_TT
        FLXPND(L,4)=w_tt
        FLXPND(L,5)=A_TT
        DO K=1,KC
        IWQZMAP(L,K)=IWQZX
        WRITE(2,31) L,I,J,K,IWQZMAP(L,K)
        ENDDO
       ENDDO 
 !      IN=0
 !      IJKC=IC*JC*KC
 !      DO M=1,IJKC
 !       READ(1,*,END=1111) I,J,K,IWQZX
 !       IN=IN+1
 !       IF (LIJW(I,J).LT.1 ) THEN
 !         PRINT*, 'i, j, k, IJCT(i,j) = ', I,J,K,LIJW(I,J)
 !         STOP 'ERROR!! invalid (i,j) in FILE wqwcmap.inp'
 !       END IF
 !       L = LIJW(I,J)
 !       IWQZMAP(L,K)=IWQZX
 !       WRITE(2,31) L,I,J,K,IWQZMAP(L,K)
 !      END DO
       
 1111  CONTINUE
!       IF (IN.NE.(LA-1)*KC) THEN
!        PRINT*, 'all active water cells should be mapped for WQ par.'
!        STOP 'ERROR!! number of lines in FILE wqwcmap.inp =\ (LA-1)'
!       END IF
C
       CLOSE(1)
C
      END IF
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

