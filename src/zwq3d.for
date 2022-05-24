C 
C**********************************************************************C
C**********************************************************************C
C**********************************************************************C
C
      SUBROUTINE WQ3D(ISTL,N,DT,TCON,TBEGIN,TIDALP,
     &  NTSPTC,IWQDT,LA,KC,IC,JC,IWQS,IOP_SAVE,iyear1,NCSTEP,SECDLAST)
C
C**********************************************************************C
C
C  CONTROL SUBROUTINE FOR WATER QUALITY MODEL
C
C  ORGINALLY CODED BY K.-Y. PARK
C  OPTIMIZED AND MODIFIED BY J. M. HAMRICK
C  LAST MODIFIED BY J. M. HAMRICK  7 APRIL 1997
C  Last modified by AEE 3/18/2007
C
!
!  WQ3D be called every 2 timestep
!
!  DT       = Time step
!  N        = N time step for the model hydrodynamc loop
!  TBEGIN   = Begining ofthe model run (cold start TBEGIN = 0)
!  TCON     = CONVERSION MULTIPLIER TO CHANGE TBEGIN TO SECONDS
!  KC       = Number of total vertical layers
!  LA       = Total active water cell (determined by Mapping 2D  to 1D)
!  IC       = Total horizontal cell (x-direction)
!  JC       = Total horizontal cell (y-direction)
!  NTSPTC   = NUMBER OF TIME STEPS PER REFERENCE TIME PERIOD (suggest to use 86400 for water quality model)
!  IWQDT    = number of water quality time steps per hydrodynamic time step (=2/or even number)
!             read in from input file 'WQ3DWC.INP'
!  TIDALP   = REFERENCE TIME PERIOD IN SEC (ie 44714.16s or 86400s) 
!  IWQS     = 1 simplified water quality version, 
!           = 0 full version
C**********************************************************************C
C
!      INCLUDE 'efdc.par'
      include 'wq.par'
!      INCLUDE 'efdc.cmn'
      INCLUDE 'wqcom.cmn'
      INTEGER, SAVE :: iniwq=1, IWQPT=0, ITHOUR,ITHRPT=0
      REAL*8 SECDLAST  
      DATA IWQTICI,IWQTAGR,IWQTSTL,IWQTSUN,IWQTBEN,IWQTPSL,IWQTNPL/7*0/
      DATA ISMTICI/0/
      iwqtsun=iwqtsun
      iwqtben=iwqtben
      iwqtpsl=iwqtpsl
      iwqtnpl=iwqtnpl
      
      ITHOUR=NTSPTC/24 
      
      IF(iniwq.eq.1) THEN
	CALL WQ3DINP(LA,KC,IC,JC,IWQDT,DT,NTSPTC,IWQS)
      CALL RWQSUN(DT,N,TBEGIN,TCON,NTSPTC,NCSTEP,SECDLAST)
	CALL WQMAPPING
	iniwq=0 
	ITHOUR=NTSPTC/24
	ENDIF

      TIMTMP=(DT*FLOAT(N))/86400
      IF(NCSTEP.GT.0) TIMTMP=SECDLAST/TCON 
!
! **  READ INITIAL CONDITIONS
!
      IF (IWQICI.EQ.1 .AND. ITNWQ.EQ.IWQTICI) then
	CALL RWQICI(IWQTICI,LA,KC)
	endif
!
! **  UPDATE TIME SERIES BOUNDARY CONDITIONS
!
      CALL RWQCSR(KC,DT,N,TCON,TBEGIN,NWQV,NTSPTC,NCSTEP,SECDLAST) 
!
! **  READ TIME/SPACE VARYING ALGAE PARAMETERS
!
      IF (IWQAGR.EQ.1 .AND. ITNWQ.EQ.IWQTAGR) CALL RWQAGR(IWQTAGR)
!
! **  READ TIME/SPACE VARYING SETTLING VELOCITIES
!
      IF (IWQSTL.EQ.1 .AND. ITNWQ.EQ.IWQTSTL) CALL RWQSTL(IWQTSTL)
!
! **  UPDATE SOLAR RADIATION INTENSITY
!   Update occurs only when the simulation day changes.
!   WQI1 = solar radiation on previous day
!   WQI2 = solar radiation two days ago
!   WQI3 = solar radiation three days ago
!
      ISTPDAY = INT((86400.0/TIDALP)*REAL(NTSPTC)) 
 !     IF (MOD(N, ISTPDAY) .EQ. 0) THEN
      IF (MOD(ITNWQ, ISTPDAY) .EQ. 0) THEN
        WQI3 = WQI2
        WQI2 = WQI1
        WQI1 = WQI0opt
        WQI0opt = 0.0
      END IF
!
! **  READ SOLAR RADIATION INTENSITY AND DAYLIGHT LENGTH
!
! Note:  IWQSUN=1 calls subroutine RWQSUN which reads the daily
!                 solar radiation data from file wqsunser.inp which
!                 are in units of Langleys/day. Note the file is same as 
!                 hourly input data. The time interval must be daily
!        IWQSUN=2 uses the hourly solar radiation data from ASER.INP
!                 and convert watts/m**2 to Langleys/day using 2.065.
!                 and adjust for photosynthetic active radiation by 0.43
!
      IF (IWQSUN.EQ.1)  THEN
        call RWQSUN(DT,N,TBEGIN,TCON,NTSPTC,NCSTEP,SECDLAST)
        WQI0=SUNDT   !SOLSRDT
        WQFD=0.5     !SOLFRDT   
      END IF
C
      IF (IWQSUN.EQ.2)  THEN
       call RWQSUN(DT,N,TBEGIN,TCON,NTSPTC,NCSTEP,SECDLAST)
!
!  PARadj = solar radiation multiplied by this factor to get the
!           photoactive available radiation (PAR) for algae growth = 0.43
!           0.43*2.065* SOLARAVG=0.88795*

        WQI0 = PARADJ*2.065*SUNDT  
        WQI0opt = MAX(WQI0opt, WQI0)
 !       write(*,*)'WQI0 =',WQI0,SUNDT,WQI0opt,PARADJ
        WQFD=1.
      END IF
! 
! **  READ BENTHIC FLUX IF REQUIRED
!
!  IF (IWQBEN.EQ.2 .AND. ITNWQ.EQ.IWQTBEN) CALL RWQBEN(IWQTBEN)          
!  Added by M. Morton 03/29/98
!  call spatially and time varying benthic flux here.  Only call RWQBEN2
!  if simulation time is >= the next time in the benthic file.
!
!   The subroutine rwqben2 should be revised for application in special case !
!
      if (iwqben .eq. 2) then                                               
        timtmp = (dt*float(n) + tcon*tbegin)/86400.0                        
        if (timtmp .ge. benday) then                                        
	  call rwqben2(timtmp,LA)                                              
        end if                                                               
      end if                                                                 
!
! **  UPDATE POINT SOURCE LOADINGS
!
      IF (IWQPSL.GE.1)
     & CALL RWQPSL(DT,N,TBEGIN,TCON,LA,KC,IWQS,NTSPTC,NCSTEP,SECDLAST)
!
! mrm ** Update atmospheric wet deposition loadings
!
      call RWQATM(LA,KC,DT,N,TBEGIN,TCON)
!
! **  READ SEDIMENT MODEL INITIAL CONDITION
!
      IF (IWQBEN.EQ.1) THEN
        IF (ISMICI.EQ.1 .AND. ITNWQ.EQ.ISMTICI) CALL RSMICI(ISMTICI,LA)
      END IF
!
! **  UPDATE TIME IN DAYS
!
!      ITNWQ = ITNWQ + IWQDT ! 2
      TINDAY = TINDAY + DTWQ
!
!   Keep track of D.O. transport component (i.e., the NEW D.O.
!   following the call to CALWQC minus OLD D.O. before the call).
!   First subtract the OLD D.O. here:
!
      DO K=1,KC
        DO L=2,LA
          xmrm = WQV(L,K,19)*DTWQ*DZCWQ(K)*HPWQ(L) 
          xDOtrn(L,K) = xDOtrn(L,K) - xmrm
          xDOall(L,K) = xDOall(L,K) - xmrm
        END DO
      END DO
!
! **  CALCULATE PHYSICAL TRANSPORT
! **  WQV(L,K,NW) SENT TO PHYSICAL TRANSPORT AND TRANSPORTED
! **  VALUE RETURNED IN WQV(L,K,NW)
!
	 CALL WQMAPPING
       CALL CALWQCN(2,WQV,LCMWQ,KCWM,NWQVM,IWQFCB,NWQV,IWTRC,F_AB)  

       CALL WWQNC2(LA,KC)  ! check negtice value due to transportation
!
!   Keep track of D.O. transport component (i.e., the NEW D.O.
!   following the call to CALWQC minus OLD D.O. before the call).
!   Now add the NEW D.O. here:
!
      DO K=1,KC
        DO L=2,LA
          xmrm = WQV(L,K,19)*DTWQ*DZCWQ(K)*HPWQ(L)
          xDOtrn(L,K) = xDOtrn(L,K) + xmrm
          xDOall(L,K) = xDOall(L,K) + xmrm
        END DO
      END DO
!
! **  LOAD WQV INTO WQVO FOR REACTION CALCULATION
!
       NMALG=0
       IF(IDNOTRVA.GT.0) NMALG=1
       DO NW=1,NWQV+NMALG
         DO K=1,KC
           DO L=2,LA
              WQVO(L,K,NW)=WQV(L,K,NW)
           END DO
         END DO
       END DO
!
! **  UPDATE WATER COLUMN KINETICS AND SEDIMENT MODEL
! **  OVER LONGER TIME INTERVALS THAN PHYSICAL TRANSPORT
! **  IF NWQKDPT .GT. 1
!
! **  WQ3D is called every even timestep. NEWKDPT is 1/2 of the IWQDT J.S. 11/28/2014
! **  ITNWQ = Count WQ timestep (note, IWQDT is even number of hynamic transport timestep
! **  IWQPT =  Avage interval count
! **  ITHRPT = Hourly average count.
! **  IWQS =0 full version, > 0 simplified version
!
      NWQKCNT=NWQKCNT+1

      IF (NWQKCNT.EQ.NWQKDPT) THEN
        NWQKCNT=0
        ITNWQ = ITNWQ + IWQDT   
        IWQPT = IWQPT + IWQDT    
        ITHRPT= ITHRPT + IWQDT 
!
! **    CALCULATE KINETIC SOURCES AND SINKS
!
        IF(IWQS.GE.0) THEN
         IF(IWQS.GE.1) THEN
	    CALL WQSKES(LA,KC,DT,N,TCON,TBEGIN,TCTMSR,NTSPTC,IWQS,SECDLAST)
         ELSE
          CALL WQSKE(LA,KC,DT,N,TCON,TBEGIN,TCTMSR,NTSPTC,SECDLAST)
	   ENDIF
         NMALG=0
         IF(IDNOTRVA.GT.0) NMALG=1
  !  After computing the reaction, WQVO=WQVO+WQV. Output is in the middle of N,N+IWQDT     J.S. 11/28/2014
         DO NW=1,NWQV+NMALG
          DO K=1,KC
           DO L=2,LA
              WQVO(L,K,NW)=WQVO(L,K,NW)*0.5
           END DO
          END DO
         END DO
        ENDIF
!
! **    DIAGNOSE NEGATIVE CONCENTRATIONS (due to kinetics)
!
        CALL WWQNC(LA,KC)  !
!
! **    WRITE TIME SERIES
!
!        IF (ITNWQ.GE.IWQTSB .AND. ITNWQ.LE.IWQTSE) THEN
!          IF (MOD(ITNWQ,IWQTSDT).EQ.0) then IWQPT
!         IF (MOD(ITNWQ,IAVGBIN).EQ.0) then
 !        write(*,*)'Time day=',IWQPT,IAVGBIN

         IF (mod(IWQPT,IAVGBIN).EQ.0) then       
 !         write(*,*)'Output day IWQPT,ITNWQ,IAVGBIN ='
 !         write(*,*)'        ', IWQPT,ITNWQ,IAVGBIN
          IWQPT=0
!		  CALL WWQTS(TINDAY,KC,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST)
	    if(isCOMP.gt. 0) then     ! output diagnostice restuls 
	    call WWQTSbin(LA,KC,DT,N,TCON,TBEGIN,TIDALP,NTSPTC,NCSTEP,SECDLAST)
	    endif
	   ENDIF
	   	   
         IF(mod(IWQPT,ITHOUR).EQ.0) THEN
          ITTTT=ITHOUR/IWQDT 
  !        write(*,*)'Check save max N,IWQPT,ITHOUR= ',N,IWQPT,ITHOUR
          if(isBIN.GE.1.or.isWQMIN.GE.1)then
            call WQDUMP3(N,KC,LA,DT,NTSPTC,TCON,TBEGIN,iyear1,
     &        NCSTEP,SECDLAST,IWQPT,ITTTT) 
          endif
         ENDIF	   
!        END IF
!
! **    CALL SEDIMENT MODEL
!  
         IF (IWQBEN.EQ.1.and.IWQS.GE.0) THEN  
!
          CALL SMMBE(LA,KC,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST)
!
          IF (ISMTS.GE.1) THEN
!
! **      WRITE SEDIMENT MODEL TIME SERIES
!
!            IF (ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE) THEN
       !       IF (MOD(ITNWQ,ISMTSDT).EQ.0) 
           IF (IWQPT.EQ.0) then   ! output using IAVGBIN interval
            CALL WSMTS(TINDAY,DT,N,TCON,TBEGIN,NTSPTC,NCSTEP,SECDLAST) 
           END IF
          END IF
!
! **      WRITE SEDIMENT MODEL FLUXES TO BINARY FILE:
!
          IF (ITNWQ.GE.ISMTSB .AND. ITNWQ.LE.ISMTSE) THEN
            CALL WSMTSbin(N,LA,iyear1)
          END IF
        END IF
!
!     Save binary file
!
!      if(IOP_SAVE.GT.0) 
!     &  Call zwqfield(N,DT,TCON,TBEGIN,IOP_SAVE,LA,KC)

      END IF

! **  END IF ON KINETIC AND SEDIMENT UPDATE

! **  WRITE RESTART FILES

       IF(IWQRST.GE.1.and.mod(N,NTSPTC*IWQRST).eq.0) then  ! write restart file reference cycle. AEE
         CALL WWQRST(LA,KC,DT,N,TCON,TBEGIN,NCSTEP,SECDLAST)
         CALL WSMRST(LA,DT,N,TCON,TBEGIN,NCSTEP,SECDLAST)
	 ENDIF


  100 FORMAT('  TIME = ',A11,' HH.MM.SS.hh')
  600 FORMAT(' ITNWQ,IWQTSDT,ITMP,IWQTSB,TWQTSE = ',5I5)
C
      RETURN
      END
